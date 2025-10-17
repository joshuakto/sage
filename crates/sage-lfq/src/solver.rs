use crate::MaxLfqConfig;
use nalgebra::DMatrix;
use sprs::CsMat;
use std::collections::{HashSet, VecDeque};

pub struct ProteinProfile {
    pub lfq_intensities: Vec<f32>,
    pub raw_intensities: Vec<f32>,
    pub n_peptides: usize,
    pub n_quantified_samples: usize,
}

pub struct ProteinSolver;

struct RatioMatrix {
    ratios: DMatrix<f32>,
    counts: DMatrix<usize>,
    active_samples: Vec<bool>,
}

impl ProteinSolver {
    pub fn quantify(
        peptide_submatrix: &CsMat<f32>,
        norm_factors: &[f64],
        config: &MaxLfqConfig,
    ) -> Option<ProteinProfile> {
        let ratio_matrix = build_ratio_matrix(peptide_submatrix, norm_factors)?;

        if !is_connected(&ratio_matrix) {
            return None;
        }

        let log_intensities = solve_least_squares(&ratio_matrix)?;
        if log_intensities
            .iter()
            .filter(|value| value.is_finite())
            .count()
            < config.min_samples_for_protein
        {
            return None;
        }

        let lfq_intensities = rescale_to_absolute(&log_intensities, peptide_submatrix);
        let raw_intensities = extract_raw_intensities(peptide_submatrix);

        Some(ProteinProfile {
            lfq_intensities,
            raw_intensities,
            n_peptides: peptide_submatrix.rows(),
            n_quantified_samples: log_intensities.len(),
        })
    }
}

fn build_ratio_matrix(submatrix: &CsMat<f32>, norm_factors: &[f64]) -> Option<RatioMatrix> {
    let n_samples = submatrix.cols();
    if n_samples == 0 {
        return Some(RatioMatrix {
            ratios: DMatrix::from_element(0, 0, 0.0),
            counts: DMatrix::from_element(0, 0, 0usize),
            active_samples: Vec::new(),
        });
    }

    let mut ratios = DMatrix::from_element(n_samples, n_samples, 0.0f32);
    let mut counts = DMatrix::from_element(n_samples, n_samples, 0usize);
    let mut active = vec![false; n_samples];

    for row in submatrix.outer_iterator() {
        let entries: Vec<(usize, f64)> = row
            .iter()
            .map(|(col, &value)| {
                (
                    col,
                    value as f64 - norm_factors.get(col).copied().unwrap_or(0.0),
                )
            })
            .collect();

        for (col, _) in &entries {
            active[*col] = true;
        }

        if entries.len() < 2 {
            continue;
        }

        for i in 0..entries.len() {
            for j in (i + 1)..entries.len() {
                let (col_i, value_i) = entries[i];
                let (col_j, value_j) = entries[j];
                let diff = (value_j - value_i) as f32;

                ratios[(col_i, col_j)] += diff;
                ratios[(col_j, col_i)] -= diff;
                counts[(col_i, col_j)] += 1;
                counts[(col_j, col_i)] += 1;
            }
        }
    }

    let mut has_ratio = false;
    for i in 0..n_samples {
        for j in 0..n_samples {
            let count = counts[(i, j)];
            if count > 0 {
                ratios[(i, j)] /= count as f32;
                has_ratio = true;
            }
        }
    }

    if !has_ratio {
        return None;
    }

    Some(RatioMatrix {
        ratios,
        counts,
        active_samples: active,
    })
}

fn is_connected(ratios: &RatioMatrix) -> bool {
    let n = ratios.ratios.nrows();
    if n == 0 {
        return true;
    }

    let active_indices: Vec<usize> = ratios
        .active_samples
        .iter()
        .enumerate()
        .filter_map(|(idx, active)| if *active { Some(idx) } else { None })
        .collect();

    if active_indices.is_empty() {
        return false;
    }

    if active_indices.len() == 1 {
        return true;
    }

    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();

    queue.push_back(active_indices[0]);
    visited.insert(active_indices[0]);

    while let Some(node) = queue.pop_front() {
        for neighbor in 0..n {
            if ratios.counts[(node, neighbor)] > 0 && !visited.contains(&neighbor) {
                visited.insert(neighbor);
                queue.push_back(neighbor);
            }
        }
    }

    active_indices.iter().all(|idx| visited.contains(idx))
}

fn solve_least_squares(ratios: &RatioMatrix) -> Option<Vec<f64>> {
    let n = ratios.ratios.nrows();
    if n == 0 {
        return Some(Vec::new());
    }

    let mut values = vec![None; n];
    let mut queue = VecDeque::new();

    let start = ratios
        .active_samples
        .iter()
        .enumerate()
        .find(|(_, active)| **active)
        .map(|(idx, _)| idx)?;

    values[start] = Some(0.0f64);
    queue.push_back(start);

    while let Some(node) = queue.pop_front() {
        let node_value = values[node].unwrap();
        for neighbor in 0..n {
            if ratios.counts[(node, neighbor)] == 0 {
                continue;
            }

            let diff = ratios.ratios[(node, neighbor)] as f64;
            let candidate = node_value + diff;

            match values[neighbor] {
                Some(current) => {
                    let avg = (current + candidate) / 2.0;
                    values[neighbor] = Some(avg);
                }
                None => {
                    values[neighbor] = Some(candidate);
                    queue.push_back(neighbor);
                }
            }
        }
    }

    if values.iter().any(|value| value.is_none()) {
        return None;
    }

    let mut result: Vec<f64> = values.into_iter().map(|value| value.unwrap()).collect();
    let mean = result.iter().copied().sum::<f64>() / result.len() as f64;
    for value in &mut result {
        *value -= mean;
    }

    Some(result)
}

fn rescale_to_absolute(log_intensities: &[f64], original_matrix: &CsMat<f32>) -> Vec<f32> {
    if log_intensities.is_empty() {
        return Vec::new();
    }

    let weights: Vec<f64> = log_intensities
        .iter()
        .map(|value| 2f64.powf(*value))
        .collect();
    let weight_sum: f64 = weights.iter().copied().sum();

    let mut total_linear = 0.0f64;
    for row in original_matrix.outer_iterator() {
        for (_, &value) in row.iter() {
            total_linear += 2f64.powf(value as f64);
        }
    }

    if weight_sum == 0.0 {
        return vec![0.0; log_intensities.len()];
    }

    weights
        .iter()
        .map(|weight| (weight / weight_sum * total_linear) as f32)
        .collect()
}

fn extract_raw_intensities(submatrix: &CsMat<f32>) -> Vec<f32> {
    let n_samples = submatrix.cols();
    let mut sums = vec![0.0f64; n_samples];
    let mut counts = vec![0usize; n_samples];

    for row in submatrix.outer_iterator() {
        for (col, &value) in row.iter() {
            sums[col] += 2f64.powf(value as f64);
            counts[col] += 1;
        }
    }

    sums.iter()
        .zip(counts.iter())
        .map(|(sum, count)| {
            if *count > 0 {
                (sum / *count as f64) as f32
            } else {
                0.0
            }
        })
        .collect()
}

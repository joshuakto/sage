use crate::MaxLfqConfig;
use nalgebra::DMatrix;
use sprs::CsMat;
use std::collections::{HashMap, HashSet, VecDeque};

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
        config: &MaxLfqConfig,
    ) -> Option<ProteinProfile> {
        let ratio_matrix = build_ratio_matrix(peptide_submatrix)?;

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
            n_quantified_samples: log_intensities
                .iter()
                .filter(|value| value.is_finite())
                .count(),
        })
    }
}

/// Build pairwise ratio matrix using median aggregation (MaxLFQ standard).
/// 
/// For each sample pair (i,j), computes the median of log-ratios across all
/// peptides observed in both samples. Median provides robustness against
/// outlier peptides compared to mean aggregation.
fn build_ratio_matrix(submatrix: &CsMat<f32>) -> Option<RatioMatrix> {
    let n_samples = submatrix.cols();
    if n_samples == 0 {
        return Some(RatioMatrix {
            ratios: DMatrix::from_element(0, 0, 0.0),
            counts: DMatrix::from_element(0, 0, 0usize),
            active_samples: Vec::new(),
        });
    }

    // Collect all pairwise ratios for median calculation
    let mut pairwise_ratios: HashMap<(usize, usize), Vec<f32>> = HashMap::new();
    let mut active = vec![false; n_samples];

    for row in submatrix.outer_iterator() {
        let entries: Vec<(usize, f64)> = row
            .iter()
            .map(|(col, &value)| (col, value as f64))
            .collect();

        for (col, _) in &entries {
            active[*col] = true;
        }

        if entries.len() < 2 {
            continue;
        }

        // Collect log-ratios for each sample pair
        for i in 0..entries.len() {
            for j in (i + 1)..entries.len() {
                let (col_i, value_i) = entries[i];
                let (col_j, value_j) = entries[j];
                let diff = (value_j - value_i) as f32;

                pairwise_ratios
                    .entry((col_i, col_j))
                    .or_insert_with(Vec::new)
                    .push(diff);
                pairwise_ratios
                    .entry((col_j, col_i))
                    .or_insert_with(Vec::new)
                    .push(-diff);
            }
        }
    }

    // Compute median for each pair (more robust than mean)
    let mut ratios = DMatrix::from_element(n_samples, n_samples, 0.0f32);
    let mut counts = DMatrix::from_element(n_samples, n_samples, 0usize);
    let mut has_ratio = false;

    for ((i, j), mut values) in pairwise_ratios {
        if !values.is_empty() {
            values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let median = compute_median(&values);
            
            ratios[(i, j)] = median;
            counts[(i, j)] = values.len();
            has_ratio = true;
        }
    }

    let active_count = active.iter().filter(|&&is_active| is_active).count();

    if !has_ratio {
        if active_count == 1 {
            return Some(RatioMatrix {
                ratios,
                counts,
                active_samples: active,
            });
        }

        return None;
    }

    Some(RatioMatrix {
        ratios,
        counts,
        active_samples: active,
    })
}

/// Compute median of sorted values.
fn compute_median(sorted_values: &[f32]) -> f32 {
    let n = sorted_values.len();
    if n == 0 {
        return 0.0;
    }
    if n % 2 == 0 {
        (sorted_values[n / 2 - 1] + sorted_values[n / 2]) / 2.0
    } else {
        sorted_values[n / 2]
    }
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

    if ratios
        .active_samples
        .iter()
        .enumerate()
        .any(|(idx, active)| *active && values[idx].is_none())
    {
        return None;
    }

    let mean = {
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for (idx, value) in values.iter().enumerate() {
            if ratios.active_samples[idx] {
                if let Some(v) = value {
                    sum += *v;
                    count += 1;
                }
            }
        }

        if count == 0 {
            0.0
        } else {
            sum / count as f64
        }
    };

    let mut result = Vec::with_capacity(values.len());
    for (idx, value) in values.into_iter().enumerate() {
        match value {
            Some(mut v) => {
                v -= mean;
                result.push(v);
            }
            None => {
                debug_assert!(!ratios.active_samples[idx]);
                result.push(f64::NAN);
            }
        }
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
    let weight_sum: f64 = weights
        .iter()
        .copied()
        .filter(|weight| weight.is_finite())
        .sum();

    let mut total_linear = 0.0f64;
    for row in original_matrix.outer_iterator() {
        for (_, &value) in row.iter() {
            total_linear += 2f64.powf(value as f64);
        }
    }

    if weight_sum == 0.0 {
        return vec![f32::NAN; log_intensities.len()];
    }

    weights
        .iter()
        .map(|weight| {
            if weight.is_finite() {
                (weight / weight_sum * total_linear) as f32
            } else {
                f32::NAN
            }
        })
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

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    #[test]
    fn preserves_consistent_fold_change_in_ratios() {
        let mut triplets = TriMat::new((2, 2));
        triplets.add_triplet(0, 0, 10.0);
        triplets.add_triplet(0, 1, 11.0);
        triplets.add_triplet(1, 0, 12.0);
        triplets.add_triplet(1, 1, 13.0);

        let matrix = triplets.to_csr();

        let ratio_matrix = build_ratio_matrix(&matrix).expect("ratio matrix");
        let log_intensities = solve_least_squares(&ratio_matrix).expect("lfq solution");

        let diff = log_intensities[1] - log_intensities[0];
        assert!((diff - 1.0).abs() < 1e-6, "expected fold-change of 1, got {diff}");
    }
    
    #[test]
    fn median_robust_to_outliers() {
        // Test that median aggregation is robust to outlier peptides
        let mut triplets = TriMat::new((5, 2)); // 5 peptides, 2 samples
        
        // 4 peptides show 1.0 log2 fold-change
        triplets.add_triplet(0, 0, 10.0);
        triplets.add_triplet(0, 1, 11.0);
        triplets.add_triplet(1, 0, 10.0);
        triplets.add_triplet(1, 1, 11.0);
        triplets.add_triplet(2, 0, 10.0);
        triplets.add_triplet(2, 1, 11.0);
        triplets.add_triplet(3, 0, 10.0);
        triplets.add_triplet(3, 1, 11.0);
        
        // 1 outlier peptide shows 5.0 log2 fold-change (contamination/mismatch)
        triplets.add_triplet(4, 0, 10.0);
        triplets.add_triplet(4, 1, 15.0);
        
        let matrix = triplets.to_csr();
        let ratio_matrix = build_ratio_matrix(&matrix).expect("ratio matrix");
        
        // Median should be ~1.0, not affected by outlier
        let median_ratio = ratio_matrix.ratios[(0, 1)];
        assert!((median_ratio - 1.0).abs() < 0.1, 
                "median should be ~1.0, got {median_ratio}");
    }
    
    #[test]
    fn single_sample_protein_quantification() {
        // Test that proteins quantified in only one sample are handled
        let config = MaxLfqConfig {
            min_peptides_per_ratio: 2,
            min_samples_for_protein: 1,
        };
        
        let mut triplets = TriMat::new((2, 4)); // 2 peptides, 4 samples
        triplets.add_triplet(0, 2, 10.0); // Only in sample 2
        triplets.add_triplet(1, 2, 12.0);
        
        let matrix = triplets.to_csr();
        let result = ProteinSolver::quantify(&matrix, &config);
        
        assert!(result.is_some(), "single-sample protein should be quantified");
        let profile = result.unwrap();
        assert_eq!(profile.n_quantified_samples, 1);
    }
    
    #[test]
    fn disconnected_samples_rejected() {
        // Test that proteins with no shared peptides between samples are rejected
        let config = MaxLfqConfig {
            min_peptides_per_ratio: 2,
            min_samples_for_protein: 1,
        };
        
        let mut triplets = TriMat::new((4, 4)); // 4 peptides, 4 samples
        // Peptides 0,1 only in samples 0,1
        triplets.add_triplet(0, 0, 10.0);
        triplets.add_triplet(0, 1, 11.0);
        triplets.add_triplet(1, 0, 10.0);
        triplets.add_triplet(1, 1, 11.0);
        
        // Peptides 2,3 only in samples 2,3 (disconnected!)
        triplets.add_triplet(2, 2, 10.0);
        triplets.add_triplet(2, 3, 11.0);
        triplets.add_triplet(3, 2, 10.0);
        triplets.add_triplet(3, 3, 11.0);
        
        let matrix = triplets.to_csr();
        let result = ProteinSolver::quantify(&matrix, &config);
        
        assert!(result.is_none(), "disconnected samples should be rejected");
    }
    
    #[test]
    fn median_calculation_correctness() {
        // Test median calculation for different scenarios
        assert_eq!(compute_median(&[1.0]), 1.0);
        assert_eq!(compute_median(&[1.0, 2.0, 3.0]), 2.0);
        assert_eq!(compute_median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_eq!(compute_median(&[]), 0.0);
    }
}

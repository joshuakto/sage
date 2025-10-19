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
        global_normalization: Option<&[f64]>,
        config: &MaxLfqConfig,
    ) -> Option<ProteinProfile> {
        let ratio_matrix = build_ratio_matrix(peptide_submatrix, global_normalization, config)?;

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
/// Applies global normalization offsets (if provided) and computes median
/// log-ratios for each sample pair. Enforces min_peptides_per_ratio threshold.
///
/// # Arguments
/// * `submatrix` - Peptide intensities (log-scale) for this protein
/// * `global_normalization` - Optional normalization offsets to subtract
/// * `config` - Configuration including min_peptides_per_ratio threshold
fn build_ratio_matrix(
    submatrix: &CsMat<f32>,
    global_normalization: Option<&[f64]>,
    config: &MaxLfqConfig,
) -> Option<RatioMatrix> {
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
            .map(|(col, &value)| {
                // Apply global normalization by SUBTRACTING offsets
                let normalized = global_normalization
                    .and_then(|norm| norm.get(col).copied())
                    .map(|offset| value as f64 - offset)
                    .unwrap_or(value as f64);
                (col, normalized)
            })
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

    // FIX REVIEWER ISSUE #1: Enforce min_peptides_per_ratio threshold
    // Invalidate edges with insufficient peptides
    let min_peptides = config.min_peptides_per_ratio;
    for i in 0..n_samples {
        for j in 0..n_samples {
            if counts[(i, j)] > 0 && counts[(i, j)] < min_peptides {
                counts[(i, j)] = 0;
                ratios[(i, j)] = 0.0;
            }
        }
    }

    // Recheck if any valid ratios remain after thresholding
    has_ratio = counts.iter().any(|&c| c > 0);

    let active_count = active.iter().filter(|&&is_active| is_active).count();

    // Special case: single-sample proteins are quantifiable
    if !has_ratio && active_count == 1 {
        return Some(RatioMatrix {
            ratios,
            counts,
            active_samples: active,
        });
    }

    // If no ratios but multiple samples, can't quantify
    if !has_ratio {
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

/// FIX REVIEWER ISSUE #2: Proper least-squares solver
///
/// Solves the MaxLFQ optimization problem by minimizing:
/// Σ(i<j) weight[i,j] * (x_i - x_j - ratio[i,j])²
///
/// This builds and solves the normal equations (A^T A)x = A^T b where:
/// - Each edge (i,j) contributes equation: x_i - x_j = ratio[i,j]
/// - Edges are weighted by peptide count for robustness
/// - Reference sample fixed at 0 to remove gauge freedom
fn solve_least_squares(ratios: &RatioMatrix) -> Option<Vec<f64>> {
    let n = ratios.ratios.nrows();
    if n == 0 {
        return Some(Vec::new());
    }

    // Find active samples and reference sample
    let active_indices: Vec<usize> = ratios
        .active_samples
        .iter()
        .enumerate()
        .filter(|(_, &active)| active)
        .map(|(idx, _)| idx)
        .collect();

    if active_indices.is_empty() {
        return None;
    }

    // Special case: single sample (no edges)
    // Just return 0 for that sample, NaN for others
    if active_indices.len() == 1 {
        let mut result = vec![f64::NAN; n];
        result[active_indices[0]] = 0.0;
        return Some(result);
    }

    // Build normal equations over active samples only to avoid singular matrix
    // Map: global index -> position in active_indices
    let m = active_indices.len();
    let mut index_map = vec![None; n];
    for (pos, &idx) in active_indices.iter().enumerate() {
        index_map[idx] = Some(pos);
    }

    let mut ata = DMatrix::zeros(m, m);
    let mut atb = nalgebra::DVector::zeros(m);

    // Add equation for each edge: x_j - x_i = ratio[i,j]
    // where ratio[i,j] = median(log(I_j) - log(I_i)) across peptides
    // Weighted by peptide count for robustness
    for i in 0..n {
        for j in (i + 1)..n {
            let count = ratios.counts[(i, j)];
            if count == 0 {
                continue;
            }

            // Map to active indices
            let Some(i_pos) = index_map[i] else { continue };
            let Some(j_pos) = index_map[j] else { continue };

            let weight = count as f64; // Weight by peptide count
            let ratio = ratios.ratios[(i, j)] as f64;

            // Equation: x_j - x_i = ratio[i,j]
            // In normal equations form:
            ata[(i_pos, i_pos)] += weight;
            ata[(j_pos, j_pos)] += weight;
            ata[(i_pos, j_pos)] -= weight;
            ata[(j_pos, i_pos)] -= weight;

            atb[i_pos] -= weight * ratio;
            atb[j_pos] += weight * ratio;
        }
    }

    // Fix gauge: set first active sample to 0
    // Replace first active sample's equation with x_ref = 0
    let reference_pos = 0; // First position in active array
    for j in 0..m {
        ata[(reference_pos, j)] = 0.0;
    }
    ata[(reference_pos, reference_pos)] = 1.0;
    atb[reference_pos] = 0.0;

    // Solve using LU decomposition (with fallback to QR)
    let solution = ata
        .clone()
        .lu()
        .solve(&atb)
        .or_else(|| ata.qr().solve(&atb))?;

    // Map solution back to full-sized vector (n samples)
    // Active samples get solved values, inactive samples get NaN
    let mut values = vec![f64::NAN; n];
    for (pos, &global_idx) in active_indices.iter().enumerate() {
        values[global_idx] = solution[pos];
    }
    
    // Center values (subtract mean of active samples)
    let mean = {
        let sum: f64 = active_indices.iter().map(|&idx| values[idx]).sum();
        sum / active_indices.len() as f64
    };

    for &idx in &active_indices {
        values[idx] -= mean;
    }

    Some(values)
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

        let config = MaxLfqConfig::default();
        let ratio_matrix = build_ratio_matrix(&matrix, None, &config).expect("ratio matrix");
        let log_intensities = solve_least_squares(&ratio_matrix).expect("lfq solution");

        let diff = log_intensities[1] - log_intensities[0];
        assert!(
            (diff - 1.0).abs() < 1e-6,
            "expected fold-change of 1, got {diff}"
        );
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
        let config = MaxLfqConfig::default();
        let ratio_matrix = build_ratio_matrix(&matrix, None, &config).expect("ratio matrix");

        // Median should be ~1.0, not affected by outlier
        let median_ratio = ratio_matrix.ratios[(0, 1)];
        assert!(
            (median_ratio - 1.0).abs() < 0.1,
            "median should be ~1.0, got {median_ratio}"
        );
    }
    
    #[test]
    fn single_sample_protein_quantification() {
        // Test that proteins quantified in only one sample are handled
        let config = MaxLfqConfig {
            min_peptides_per_ratio: 2,
            min_samples_for_protein: 1,
            ..Default::default()
        };

        let mut triplets = TriMat::new((2, 4)); // 2 peptides, 4 samples
        triplets.add_triplet(0, 2, 10.0); // Only in sample 2
        triplets.add_triplet(1, 2, 12.0);

        let matrix = triplets.to_csr();
        let result = ProteinSolver::quantify(&matrix, None, &config);

        assert!(
            result.is_some(),
            "single-sample protein should be quantified"
        );
        let profile = result.unwrap();
        assert_eq!(profile.n_quantified_samples, 1);
    }
    
    #[test]
    fn disconnected_samples_rejected() {
        // Test that proteins with no shared peptides between samples are rejected
        let config = MaxLfqConfig {
            min_peptides_per_ratio: 2,
            min_samples_for_protein: 1,
            ..Default::default()
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
        let result = ProteinSolver::quantify(&matrix, None, &config);

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

    #[test]
    fn applies_global_normalization_offsets() {
        // Test that global normalization offsets are correctly SUBTRACTED
        let mut triplets = TriMat::new((1, 2));
        triplets.add_triplet(0, 0, 10.0);
        triplets.add_triplet(0, 1, 12.0);

        let matrix = triplets.to_csr();
        let config = MaxLfqConfig {
            min_peptides_per_ratio: 1, // Allow single peptide
            ..Default::default()
        };
        let normalization = [0.0, 1.0]; // Sample 1 has 1.0 offset to subtract

        let ratio_matrix = build_ratio_matrix(&matrix, Some(&normalization), &config).expect("ratio");
        
        // After subtracting offset: sample 0 = 10.0, sample 1 = 12.0 - 1.0 = 11.0
        // Ratio should be 1.0 (not 2.0)
        assert!((ratio_matrix.ratios[(0, 1)] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn enforces_min_peptides_per_ratio_threshold() {
        // Test that sample pairs with insufficient peptides are excluded
        let config = MaxLfqConfig {
            min_peptides_per_ratio: 2,
            ..Default::default()
        };

        let mut triplets = TriMat::new((1, 3)); // 1 peptide, 3 samples
        triplets.add_triplet(0, 0, 10.0);
        triplets.add_triplet(0, 1, 11.0);
        triplets.add_triplet(0, 2, 12.0);

        let matrix = triplets.to_csr();
        let result = build_ratio_matrix(&matrix, None, &config);

        // With only 1 peptide per pair and threshold=2, all edges invalid
        // Should return None (cannot quantify with no valid edges)
        assert!(result.is_none(), "Should reject protein with insufficient peptides per ratio");
    }

    #[test]
    fn least_squares_handles_cyclic_inconsistencies() {
        // Test that proper least-squares solver minimizes error across cycles
        // Triangle with 2+ peptides per edge
        let config = MaxLfqConfig {
            min_peptides_per_ratio: 1, // Allow edges with 1+ peptides
            ..Default::default()
        };
        
        let mut triplets = TriMat::new((6, 3));
        // 2 peptides between samples 0,1 with diff ~1.0
        triplets.add_triplet(0, 0, 10.0);
        triplets.add_triplet(0, 1, 11.0);
        triplets.add_triplet(1, 0, 10.0);
        triplets.add_triplet(1, 1, 11.0);
        // 2 peptides between samples 1,2 with diff ~1.0
        triplets.add_triplet(2, 1, 11.0);
        triplets.add_triplet(2, 2, 12.0);
        triplets.add_triplet(3, 1, 11.0);
        triplets.add_triplet(3, 2, 12.0);
        // 2 peptides between samples 0,2 with diff ~1.5 (inconsistent)
        triplets.add_triplet(4, 0, 10.0);
        triplets.add_triplet(4, 2, 11.5);
        triplets.add_triplet(5, 0, 10.0);
        triplets.add_triplet(5, 2, 11.5);

        let matrix = triplets.to_csr();
        let ratio_matrix = build_ratio_matrix(&matrix, None, &config).expect("ratio");
        let log_intensities = solve_least_squares(&ratio_matrix).expect("solution");

        // Should find values that minimize total squared error
        // Not dependent on traversal order (unlike BFS)
        let diff_01 = log_intensities[1] - log_intensities[0];
        let diff_12 = log_intensities[2] - log_intensities[1];
        let diff_02 = log_intensities[2] - log_intensities[0];

        // Errors should be roughly balanced (not all error on one edge)
        let error_01 = (diff_01 - 1.0).abs();
        let error_12 = (diff_12 - 1.0).abs();
        let error_02 = (diff_02 - 1.5).abs();
        
        // Total squared error should be minimized
        assert!(error_01 < 0.3);
        assert!(error_12 < 0.3);
        assert!(error_02 < 0.3);
    }

    #[test]
    fn handles_protein_in_subset_of_samples() {
        // Bug fix test: protein appears in samples 0,1 of a 4-sample experiment
        // Should build system over active samples only, not singular matrix
        let config = MaxLfqConfig {
            min_peptides_per_ratio: 1,
            ..Default::default()
        };

        // 4 samples total, but protein only in samples 0 and 1
        let mut triplets = TriMat::new((2, 4));
        triplets.add_triplet(0, 0, 10.0); // Peptide 1 in samples 0,1
        triplets.add_triplet(0, 1, 11.0);
        triplets.add_triplet(1, 0, 10.0); // Peptide 2 in samples 0,1
        triplets.add_triplet(1, 1, 11.0);

        let matrix = triplets.to_csr();
        let ratio_matrix = build_ratio_matrix(&matrix, None, &config).expect("ratio");
        let log_intensities = solve_least_squares(&ratio_matrix).expect("solution");

        // Should succeed (not return None due to singular matrix)
        assert_eq!(log_intensities.len(), 4);
        
        // Samples 0,1 should have values
        assert!(!log_intensities[0].is_nan());
        assert!(!log_intensities[1].is_nan());
        
        // Samples 2,3 should be NaN (inactive)
        assert!(log_intensities[2].is_nan());
        assert!(log_intensities[3].is_nan());
        
        // Ratio between active samples should be preserved
        let diff = log_intensities[1] - log_intensities[0];
        assert!((diff - 1.0).abs() < 1e-6);
    }
}

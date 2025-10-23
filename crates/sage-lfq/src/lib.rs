mod delayed_normalization;
pub mod error;
pub mod matrix;
pub mod solver;

pub use error::MaxLfqError;
pub use matrix::IntensityMatrix;
pub use solver::ProteinProfile;

use solver::ProteinSolver;

/// Reason why a protein group was skipped during quantification
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SkipReason {
    /// Protein had disconnected peptide graphs (peptides don't connect all samples)
    DisconnectedGraph,
    /// Protein had fewer peptides than the minimum threshold
    InsufficientPeptides,
}

/// Information about a protein group that could not be quantified
#[derive(Debug, Clone)]
pub struct SkippedProtein {
    pub protein_ids: Vec<String>,
    pub reason: SkipReason,
}

/// Result of protein quantification including both successful and skipped proteins
#[derive(Debug, Clone)]
pub struct QuantificationResult {
    /// Successfully quantified protein groups
    pub quantified: Vec<ProteinQuantResult>,
    /// Protein groups that were skipped with reasons
    pub skipped: Vec<SkippedProtein>,
}

/// Trait representing the minimal interface required to build LFQ intensity matrices
/// from peptide-level quantification traces.
pub trait QuantTrace {
    fn intensities(&self) -> &[f64];
    fn peptide_index(&self) -> usize;
}

/// Main entry point for MaxLFQ quantification.
///
/// Quantifies protein groups using the MaxLFQ algorithm. Returns both successfully
/// quantified proteins and those that were skipped with explicit reasons.
///
/// # Returns
///
/// Returns `Ok(QuantificationResult)` containing:
/// - `quantified`: Successfully quantified protein groups
/// - `skipped`: Protein groups that could not be quantified, with reasons
///
/// # Errors
///
/// Returns `Err` only for unrecoverable failures:
/// - `InvalidDimensions`: Peptide indices out of bounds
/// - `NormalizationFailed`: Global normalization computation failed
/// - `InvalidReferenceSample`: Reference sample index invalid
/// - `OptimizationError`: Numerical solver failed
pub fn quantify_proteins<T: QuantTrace>(
    peptide_traces: &[T],
    protein_groups: &[(Vec<String>, Vec<usize>)], // (protein_ids, peptide_indices)
    config: MaxLfqConfig,
) -> Result<QuantificationResult, MaxLfqError> {
    let n_samples = peptide_traces
        .first()
        .map(|trace| trace.intensities().len())
        .unwrap_or(0);

    let intensity_matrix = matrix::IntensityMatrix::from_peptide_traces(peptide_traces, n_samples);

    // Compute global normalization offsets if enabled
    let normalization = if config.use_global_normalization {
        Some(delayed_normalization::compute_delayed_normalization(
            &intensity_matrix,
            config.reference_sample,
        )?)
    } else {
        None
    };

    let mut quantified = Vec::with_capacity(protein_groups.len());
    let mut skipped = Vec::new();

    for (protein_ids, peptide_indices) in protein_groups {
        if peptide_indices.is_empty() {
            continue;
        }

        if peptide_indices
            .iter()
            .any(|&idx| idx >= intensity_matrix.n_peptides)
        {
            return Err(MaxLfqError::InvalidDimensions);
        }

        let peptide_submatrix = intensity_matrix.get_protein_peptides(peptide_indices);

        if peptide_submatrix.rows() < config.min_peptides_per_ratio {
            // Skip proteins below threshold - should be filtered upstream,
            // but checked here as a safety measure for standalone library use
            skipped.push(SkippedProtein {
                protein_ids: protein_ids.clone(),
                reason: SkipReason::InsufficientPeptides,
            });
            continue;
        }

        match ProteinSolver::quantify(&peptide_submatrix, normalization.as_deref(), &config) {
            Some(profile) => {
                let mut sample_coverage = vec![false; intensity_matrix.n_samples];
                for row in peptide_submatrix.outer_iterator() {
                    for (col, _) in row.iter() {
                        sample_coverage[col] = true;
                    }
                }

                quantified.push(ProteinQuantResult {
                    protein_ids: protein_ids.clone(),
                    lfq_intensities: profile.lfq_intensities,
                    peptide_count: profile.n_peptides,
                    sample_coverage,
                });
            }
            None => {
                // Skip proteins with disconnected ratio graphs - common when
                // peptides only connect subsets of samples (e.g., separate batches)
                skipped.push(SkippedProtein {
                    protein_ids: protein_ids.clone(),
                    reason: SkipReason::DisconnectedGraph,
                });
            }
        }
    }

    Ok(QuantificationResult {
        quantified,
        skipped,
    })
}

#[derive(Debug, Clone)]
pub struct MaxLfqConfig {
    pub min_peptides_per_ratio: usize,   // Default: 2
    pub min_samples_for_protein: usize,  // Default: 1
    pub use_global_normalization: bool,  // Default: true
    pub reference_sample: Option<usize>, // Default: None (auto-select)
}

impl Default for MaxLfqConfig {
    fn default() -> Self {
        Self {
            min_peptides_per_ratio: 2,
            min_samples_for_protein: 1,
            use_global_normalization: true,
            reference_sample: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ProteinQuantResult {
    pub protein_ids: Vec<String>,
    pub lfq_intensities: Vec<f32>,
    pub peptide_count: usize,
    pub sample_coverage: Vec<bool>,
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MockTrace {
        intensities: Vec<f64>,
        peptide_idx: usize,
    }

    impl QuantTrace for MockTrace {
        fn intensities(&self) -> &[f64] {
            &self.intensities
        }

        fn peptide_index(&self) -> usize {
            self.peptide_idx
        }
    }

    #[test]
    fn disconnected_proteins_are_skipped_not_fatal() {
        // Create test data with 4 samples
        // Protein 1: Fully connected (peptides in all samples)
        // Protein 2: Disconnected (peptides only in samples 0-1 OR 2-3, no shared peptides)
        let traces = vec![
            // Protein 1 peptides (fully connected)
            MockTrace {
                intensities: vec![10.0, 11.0, 12.0, 13.0],
                peptide_idx: 0,
            },
            MockTrace {
                intensities: vec![20.0, 21.0, 22.0, 23.0],
                peptide_idx: 1,
            },
            // Protein 2 peptides (disconnected: pep2 only in samples 0-1, pep3 only in samples 2-3)
            MockTrace {
                intensities: vec![30.0, 31.0, f64::NAN, f64::NAN],
                peptide_idx: 2,
            },
            MockTrace {
                intensities: vec![f64::NAN, f64::NAN, 32.0, 33.0],
                peptide_idx: 3,
            },
        ];

        let protein_groups = vec![
            (vec!["Protein1".to_string()], vec![0, 1]), // Fully connected
            (vec!["Protein2".to_string()], vec![2, 3]), // Disconnected
        ];

        let config = MaxLfqConfig::default();

        // Should not panic or return error - disconnected protein should be skipped
        let result = quantify_proteins(&traces, &protein_groups, config);
        assert!(result.is_ok(), "Should not fail on disconnected proteins");

        let quant_result = result.unwrap();

        // Only the connected protein should be quantified
        assert_eq!(quant_result.quantified.len(), 1, "Should quantify 1 protein");
        assert_eq!(quant_result.quantified[0].protein_ids, vec!["Protein1"]);
        assert_eq!(quant_result.quantified[0].peptide_count, 2);

        // All samples should have coverage for the connected protein
        assert_eq!(quant_result.quantified[0].sample_coverage, vec![true, true, true, true]);

        // Verify the disconnected protein was skipped with correct reason
        assert_eq!(quant_result.skipped.len(), 1, "Should skip 1 protein");
        assert_eq!(quant_result.skipped[0].protein_ids, vec!["Protein2"]);
        assert_eq!(quant_result.skipped[0].reason, SkipReason::DisconnectedGraph);
    }
}

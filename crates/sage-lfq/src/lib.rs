mod delayed_normalization;
pub mod error;
pub mod matrix;
pub mod solver;

pub use error::MaxLfqError;
pub use matrix::IntensityMatrix;
pub use solver::ProteinProfile;

use solver::ProteinSolver;

/// Trait representing the minimal interface required to build LFQ intensity matrices
/// from peptide-level quantification traces.
pub trait QuantTrace {
    fn intensities(&self) -> &[f64];
    fn peptide_index(&self) -> usize;
}

/// Main entry point for MaxLFQ quantification
pub fn quantify_proteins<T: QuantTrace>(
    peptide_traces: &[T],
    protein_groups: &[(Vec<String>, Vec<usize>)], // (protein_ids, peptide_indices)
    config: MaxLfqConfig,
) -> Result<Vec<ProteinQuantResult>, MaxLfqError> {
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

    let mut results = Vec::with_capacity(protein_groups.len());

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
            return Err(MaxLfqError::InsufficientPeptides(protein_ids.join(",")));
        }

        match ProteinSolver::quantify(&peptide_submatrix, normalization.as_deref(), &config) {
            Some(profile) => {
                let mut sample_coverage = vec![false; intensity_matrix.n_samples];
                for row in peptide_submatrix.outer_iterator() {
                    for (col, _) in row.iter() {
                        sample_coverage[col] = true;
                    }
                }

                results.push(ProteinQuantResult {
                    protein_ids: protein_ids.clone(),
                    lfq_intensities: profile.lfq_intensities,
                    peptide_count: profile.n_peptides,
                    sample_coverage,
                });
            }
            None => return Err(MaxLfqError::DisconnectedProtein),
        }
    }

    Ok(results)
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

pub mod error;
pub mod matrix;
pub mod normalization;
pub mod solver;

pub use error::MaxLfqError;
pub use matrix::IntensityMatrix;
pub use normalization::NormalizationFactors;
pub use solver::ProteinProfile;

use sage_core::lfq::PeptideQuantTrace;
use solver::ProteinSolver;

/// Main entry point for MaxLFQ quantification
pub fn quantify_proteins(
    peptide_traces: &[PeptideQuantTrace],
    protein_groups: &[(Vec<String>, Vec<usize>)], // (protein_ids, peptide_indices)
    config: MaxLfqConfig,
) -> Result<Vec<ProteinQuantResult>, MaxLfqError> {
    let n_samples = peptide_traces
        .first()
        .map(|trace| trace.intensities.len())
        .unwrap_or(0);

    let intensity_matrix = matrix::IntensityMatrix::from_peptide_traces(peptide_traces, n_samples);
    
    // Note: Normalization is handled implicitly by MaxLFQ's least-squares optimization
    // Pre-normalization was removed as it inappropriately removed biological signal

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

        match ProteinSolver::quantify(&peptide_submatrix, &config) {
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
    pub min_peptides_per_ratio: usize,  // Default: 2
    pub min_samples_for_protein: usize, // Default: 1
}

#[derive(Debug, Clone)]
pub struct ProteinQuantResult {
    pub protein_ids: Vec<String>,
    pub lfq_intensities: Vec<f32>,
    pub peptide_count: usize,
    pub sample_coverage: Vec<bool>,
}

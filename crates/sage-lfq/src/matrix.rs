use sage_core::lfq::PeptideQuantTrace;
use sprs::{CsMat, TriMat};

pub struct IntensityMatrix {
    // CSR format for efficient row access
    pub matrix: CsMat<f32>,
    pub peptide_ids: Vec<String>,
    pub sample_names: Vec<String>,
    pub n_peptides: usize,
    pub n_samples: usize,
}

impl IntensityMatrix {
    /// Convert peptide traces to sparse log2 intensity matrix
    pub fn from_peptide_traces(traces: &[PeptideQuantTrace], n_samples: usize) -> Self {
        let mut triplets = TriMat::new((traces.len(), n_samples));

        for (row, trace) in traces.iter().enumerate() {
            for (col, &intensity) in trace.intensities.iter().take(n_samples).enumerate() {
                if intensity > 0.0 {
                    let value = (intensity as f32).log2();
                    triplets.add_triplet(row, col, value);
                }
            }
        }

        Self {
            matrix: triplets.to_csr(),
            peptide_ids: traces.iter().map(|t| format!("{}", t.peptide.0)).collect(),
            sample_names: (0..n_samples).map(|i| format!("S{}", i)).collect(),
            n_peptides: traces.len(),
            n_samples,
        }
    }

    /// Get peptides for a specific protein
    pub fn get_protein_peptides(&self, peptide_indices: &[usize]) -> CsMat<f32> {
        let mut triplets = TriMat::new((peptide_indices.len(), self.n_samples));

        for (row_pos, &src_row) in peptide_indices.iter().enumerate() {
            if let Some(row) = self.matrix.outer_view(src_row) {
                for (col, &value) in row.iter() {
                    triplets.add_triplet(row_pos, col, value);
                }
            }
        }

        triplets.to_csr()
    }

    /// Count shared peptides between two samples
    pub fn count_shared_peptides(&self, sample_i: usize, sample_j: usize) -> usize {
        if sample_i >= self.n_samples || sample_j >= self.n_samples {
            return 0;
        }

        let mut shared = 0usize;

        for row in self.matrix.outer_iterator() {
            let mut has_i = false;
            let mut has_j = false;

            for (col, _) in row.iter() {
                if col == sample_i {
                    has_i = true;
                } else if col == sample_j {
                    has_j = true;
                }

                if has_i && has_j {
                    shared += 1;
                    break;
                }
            }
        }

        shared
    }
}

use sage_core::database::PeptideIx;
use sage_core::lfq::{Peak, PeptideQuantTrace, PrecursorId};
use sage_core::ml::matrix::Matrix;
use sage_lfq::{quantify_proteins, MaxLfqConfig};

#[test]
fn test_crate_compiles() {
    let traces = vec![PeptideQuantTrace {
        precursor: PrecursorId::Combined(PeptideIx(0)),
        peptide: PeptideIx(0),
        decoy: false,
        peak: Peak {
            rt: 0,
            spectral_angle: 0.0,
            score: 0.0,
            q_value: 0.0,
        },
        intensities: vec![100.0, 200.0, 150.0],
        reference_file_id: 0,
        dot_product: Matrix::zeros(0, 0),
        spectral_angle: Matrix::zeros(0, 0),
        isotope_traces: Matrix::zeros(0, 0),
        raw_isotope_traces: Matrix::zeros(0, 0),
        isotopic_distribution: [0.0; 3],
        time_warps: vec![],
    }];

    let config = MaxLfqConfig {
        min_peptides_per_ratio: 1,
        min_samples_for_protein: 1,
    };

    let result = quantify_proteins(&traces, &[(vec!["P1".to_string()], vec![0])], config);
    assert!(result.is_ok());
}

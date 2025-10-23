use sage_lfq::{quantify_proteins, MaxLfqConfig, QuantTrace};

#[derive(Clone)]
struct TestTrace {
    intensities: Vec<f64>,
    peptide_index: usize,
}

impl QuantTrace for TestTrace {
    fn intensities(&self) -> &[f64] {
        &self.intensities
    }

    fn peptide_index(&self) -> usize {
        self.peptide_index
    }
}

#[test]
fn test_crate_compiles() {
    let traces = vec![TestTrace {
        intensities: vec![100.0, 200.0, 150.0],
        peptide_index: 0,
    }];

    let config = MaxLfqConfig {
        min_peptides_per_ratio: 1,
        min_samples_for_protein: 1,
        ..Default::default()
    };

    let result = quantify_proteins(&traces, &[(vec!["P1".to_string()], vec![0])], config);
    assert!(result.is_ok());
}

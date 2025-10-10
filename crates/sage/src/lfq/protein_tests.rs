#![cfg(all(test, feature = "protein-quant-prototype"))]

use std::sync::Arc;

use crate::{
    database::{IndexedDatabase, PeptideIx},
    enzyme::Position,
    lfq::{Peak, PeptideQuantTrace, PrecursorId, ProteinQuantTrace},
    ml::matrix::Matrix,
    peptide::Peptide,
};

trait IndexedDatabaseProteinsExt {
    fn proteins(&self, ix: PeptideIx) -> Vec<Arc<str>>;
}

impl IndexedDatabaseProteinsExt for IndexedDatabase {
    fn proteins(&self, ix: PeptideIx) -> Vec<Arc<str>> {
        self.peptides[ix.0 as usize].proteins.clone()
    }
}

fn fake_database(mut peptides: Vec<Peptide>) -> IndexedDatabase {
    peptides
        .iter_mut()
        .for_each(|pep| pep.proteins.sort_unstable());
    IndexedDatabase {
        peptides,
        fragments: Vec::new(),
        ion_kinds: Vec::new(),
        min_value: Vec::new(),
        potential_mods: Vec::new(),
        bucket_size: 1,
        generate_decoys: true,
        decoy_tag: "DECOY_".to_string(),
    }
}

fn fake_peptide(seq: &str, proteins: &[&str], decoy: bool) -> Peptide {
    Peptide {
        decoy,
        sequence: Arc::<[u8]>::from(seq.as_bytes()),
        modifications: vec![0.0; seq.len()],
        nterm: None,
        cterm: None,
        monoisotopic: 0.0,
        missed_cleavages: 0,
        semi_enzymatic: false,
        position: Position::Full,
        proteins: proteins.iter().map(|acc| Arc::<str>::from(*acc)).collect(),
    }
}

fn fake_trace(peptide: usize, decoy: bool, q_value: f32, intensities: &[f64]) -> PeptideQuantTrace {
    PeptideQuantTrace {
        precursor: PrecursorId::Combined(PeptideIx(peptide as u32)),
        peptide: PeptideIx(peptide as u32),
        decoy,
        peak: Peak {
            rt: 0,
            spectral_angle: 0.0,
            score: 0.0,
            q_value,
        },
        intensities: intensities.to_vec(),
        reference_file_id: 0,
        dot_product: Matrix::zeros(0, 0),
        spectral_angle: Matrix::zeros(0, 0),
        isotope_traces: Matrix::zeros(0, 0),
        raw_isotope_traces: Matrix::zeros(0, 0),
        isotopic_distribution: [0.0f32; 3],
        time_warps: vec![0isize; intensities.len()],
    }
}

#[test]
#[ignore = "ProteinQuantTrace helper not implemented yet"]
fn intensity_rollup_sums_targets() {
    let run_count = 3;
    let db = fake_database(vec![
        fake_peptide("PEPA", &["P12345"], false),
        fake_peptide("PEPB", &["P12345"], false),
    ]);
    let traces = vec![
        fake_trace(0, false, 0.004, &[100.0, 0.0, 30.0]),
        fake_trace(1, false, 0.006, &[50.0, 50.0, 0.0]),
    ];
    let proteins = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, 0.01);
    let expected_accessions = db.proteins(PeptideIx(0));
    assert_eq!(expected_accessions, db.proteins(PeptideIx(1)));

    let target = proteins
        .get(expected_accessions[0].as_ref())
        .expect("target protein missing");
    assert_eq!(target.accessions, expected_accessions);
    assert!(!target.decoy);
    assert_eq!(target.intensities, vec![150.0, 50.0, 30.0]);
    assert!((target.q_value - 0.004).abs() < f32::EPSILON);
    assert_eq!(target.total_peptide_count, 2);
    assert_eq!(target.passing_peptide_count, 2);
    assert_eq!(target.peptide_indices, vec![PeptideIx(0), PeptideIx(1)]);
    assert_eq!(target.run_coverage, vec![2, 1, 1]);
}

#[test]
#[ignore = "ProteinQuantTrace helper not implemented yet"]
fn decoy_peptides_are_tracked_separately() {
    let run_count = 3;
    let db = fake_database(vec![
        fake_peptide("PEPA", &["P12345"], false),
        fake_peptide("PEPB", &["DECOY_P12345"], true),
    ]);
    let traces = vec![
        fake_trace(0, false, 0.004, &[100.0, 10.0, 0.0]),
        fake_trace(1, true, 0.002, &[12.0, 8.0, 0.0]),
    ];
    let proteins = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, 0.01);

    let target_accession = db.proteins(PeptideIx(0));
    let target = proteins
        .get(target_accession[0].as_ref())
        .expect("target protein missing");
    assert!(!target.decoy);
    assert_eq!(target.total_peptide_count, 1);
    assert_eq!(target.passing_peptide_count, 1);
    assert_eq!(target.peptide_indices, vec![PeptideIx(0)]);
    assert_eq!(target.run_coverage, vec![1, 1, 0]);
    assert_eq!(target.intensities, vec![100.0, 10.0, 0.0]);

    let decoy_accession = db.proteins(PeptideIx(1));
    let decoy = proteins
        .get(decoy_accession[0].as_ref())
        .expect("decoy protein missing");
    assert!(decoy.decoy);
    assert_eq!(proteins.len(), 2);
    assert_eq!(decoy.accessions, decoy_accession);
    assert_eq!(decoy.total_peptide_count, 1);
    assert_eq!(decoy.passing_peptide_count, 1);
    assert_eq!(decoy.peptide_indices, vec![PeptideIx(1)]);
    assert_eq!(decoy.run_coverage, vec![1, 1, 0]);
    assert_eq!(decoy.intensities, vec![12.0, 8.0, 0.0]);
}

#[test]
#[ignore = "ProteinQuantTrace helper not implemented yet"]
fn q_value_filtering_excludes_high_q_peptides() {
    let run_count = 3;
    let db = fake_database(vec![
        fake_peptide("PEPA", &["P12345"], false),
        fake_peptide("PEPB", &["P12345"], false),
        fake_peptide("PEPC", &["P12345"], false),
    ]);
    let traces = vec![
        fake_trace(0, false, 0.004, &[100.0, 0.0, 30.0]),
        fake_trace(1, false, 0.006, &[50.0, 50.0, 0.0]),
        fake_trace(2, false, 0.2, &[400.0, 400.0, 400.0]),
    ];
    let proteins = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, 0.01);

    let accession = db.proteins(PeptideIx(0));
    let target = proteins
        .get(accession[0].as_ref())
        .expect("target protein missing");

    assert!(!target.decoy);
    assert_eq!(target.total_peptide_count, 3);
    assert_eq!(target.passing_peptide_count, 2);
    assert_eq!(target.peptide_indices, vec![PeptideIx(0), PeptideIx(1)]);
    assert_eq!(target.intensities, vec![150.0, 50.0, 30.0]);
    assert_eq!(target.run_coverage, vec![2, 1, 1]);
    assert!((target.q_value - 0.004).abs() < f32::EPSILON);
}

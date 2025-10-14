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

fn fake_trace_with_precursor(
    precursor: PrecursorId,
    peptide: usize,
    decoy: bool,
    q_value: f32,
    intensities: &[f64],
) -> PeptideQuantTrace {
    PeptideQuantTrace {
        precursor,
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

fn fake_trace(peptide: usize, decoy: bool, q_value: f32, intensities: &[f64]) -> PeptideQuantTrace {
    fake_trace_with_precursor(
        PrecursorId::Combined(PeptideIx(peptide as u32)),
        peptide,
        decoy,
        q_value,
        intensities,
    )
}

#[test]
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
    let mut expected_accessions = db.proteins(PeptideIx(0));
    expected_accessions.sort_unstable();
    expected_accessions.dedup();
    assert_eq!(expected_accessions, db.proteins(PeptideIx(1)));

    let target = proteins
        .get(&expected_accessions)
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

    let mut target_accession = db.proteins(PeptideIx(0));
    target_accession.sort_unstable();
    target_accession.dedup();
    let target = proteins
        .get(&target_accession)
        .expect("target protein missing");
    assert!(!target.decoy);
    assert_eq!(target.total_peptide_count, 1);
    assert_eq!(target.passing_peptide_count, 1);
    assert_eq!(target.peptide_indices, vec![PeptideIx(0)]);
    assert_eq!(target.run_coverage, vec![1, 1, 0]);
    assert_eq!(target.intensities, vec![100.0, 10.0, 0.0]);

    let mut decoy_accession = db.proteins(PeptideIx(1));
    decoy_accession.sort_unstable();
    decoy_accession.dedup();
    let decoy = proteins
        .get(&decoy_accession)
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
fn reversed_decoy_peptides_get_unique_protein_groups() {
    let run_count = 2;
    let db = fake_database(vec![
        fake_peptide("PEPTIDE", &["P12345"], false),
        fake_peptide("EDITPEP", &["P12345"], true),
    ]);
    let traces = vec![
        fake_trace(0, false, 0.004, &[25.0, 50.0]),
        fake_trace(1, true, 0.006, &[10.0, 5.0]),
    ];
    let proteins = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, 0.01);

    let mut target_accession = db.proteins(PeptideIx(0));
    target_accession.sort_unstable();
    target_accession.dedup();
    let target = proteins
        .get(&target_accession)
        .expect("target protein missing");
    assert_eq!(proteins.len(), 2);
    assert!(!target.decoy);
    assert_eq!(target.accessions, target_accession);
    assert_eq!(target.intensities, vec![25.0, 50.0]);

    let decoy_accession = vec![Arc::<str>::from(format!("{}{}", db.decoy_tag, "P12345"))];
    let decoy = proteins
        .get(&decoy_accession)
        .expect("decoy protein missing");
    assert!(decoy.decoy);
    assert_eq!(decoy.accessions, decoy_accession);
    assert_eq!(decoy.intensities, vec![10.0, 5.0]);
}

#[test]
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

    let mut accession = db.proteins(PeptideIx(0));
    accession.sort_unstable();
    accession.dedup();
    let target = proteins.get(&accession).expect("target protein missing");

    assert!(!target.decoy);
    assert_eq!(target.total_peptide_count, 3);
    assert_eq!(target.passing_peptide_count, 2);
    assert_eq!(target.peptide_indices, vec![PeptideIx(0), PeptideIx(1)]);
    assert_eq!(target.intensities, vec![150.0, 50.0, 30.0]);
    assert_eq!(target.run_coverage, vec![2, 1, 1]);
    assert!((target.q_value - 0.004).abs() < f32::EPSILON);
}

#[test]
fn run_coverage_counts_unique_peptides_per_run() {
    let run_count = 2;
    let db = fake_database(vec![fake_peptide("PEPU", &["P12345"], false)]);
    let traces = vec![
        fake_trace(0, false, 0.004, &[0.0, 10.0]),
        fake_trace(0, false, 0.003, &[5.0, 0.0]),
    ];

    let proteins = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, 0.01);

    let mut accession = db.proteins(PeptideIx(0));
    accession.sort_unstable();
    accession.dedup();
    let target = proteins.get(&accession).expect("target protein missing");

    assert_eq!(target.intensities, vec![5.0, 10.0]);
    assert_eq!(target.run_coverage, vec![1, 1]);
    assert_eq!(target.peptide_indices, vec![PeptideIx(0)]);
    assert_eq!(target.total_peptide_count, 2);
    // The passing peptide count mirrors the number of unique peptide indices.
    assert_eq!(target.passing_peptide_count, target.peptide_indices.len());
}

#[test]
fn passing_peptide_count_matches_unique_peptides_without_charge_combining() {
    let run_count = 1;
    let db = fake_database(vec![fake_peptide("PEPA", &["P12345"], false)]);
    let traces = vec![
        fake_trace_with_precursor(
            PrecursorId::Charged((PeptideIx(0), 2)),
            0,
            false,
            0.001,
            &[100.0],
        ),
        fake_trace_with_precursor(
            PrecursorId::Charged((PeptideIx(0), 3)),
            0,
            false,
            0.002,
            &[150.0],
        ),
        // The third trace represents a second quantified feature for the same
        // peptide/charge combination (e.g. an additional chromatographic peak
        // promoted during alignment). The protein rollup should still include
        // its intensity while ensuring the unique peptide count is not
        // inflated.
        fake_trace_with_precursor(
            PrecursorId::Charged((PeptideIx(0), 2)),
            0,
            false,
            0.003,
            &[200.0],
        ),
    ];

    let proteins = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, 0.01);

    let mut accession = db.proteins(PeptideIx(0));
    accession.sort_unstable();
    accession.dedup();
    let trace = proteins
        .get(&accession)
        .expect("expected protein entry to exist");

    assert_eq!(trace.total_peptide_count, 3);
    assert_eq!(trace.peptide_indices.len(), 1);
    assert_eq!(trace.passing_peptide_count, 1);
    assert_eq!(trace.intensities, vec![450.0]);
    assert_eq!(trace.run_coverage, vec![1]);
}

fn digest_protein_map(
    proteins: &std::collections::BTreeMap<Vec<Arc<str>>, ProteinQuantTrace>,
) -> String {
    use std::fmt::Write;

    let mut fingerprint = String::new();
    for (accessions, trace) in ProteinQuantTrace::ordered_groups(proteins) {
        if !fingerprint.is_empty() {
            fingerprint.push('\n');
        }

        let accession_list = accessions
            .iter()
            .map(|a| a.as_ref())
            .collect::<Vec<_>>()
            .join(",");
        let intensities = trace
            .intensities
            .iter()
            .map(|value| format!("{value:.6}"))
            .collect::<Vec<_>>()
            .join(",");
        let coverage = trace
            .run_coverage
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let peptides = trace
            .peptide_indices
            .iter()
            .map(|ix| ix.0.to_string())
            .collect::<Vec<_>>()
            .join(",");

        write!(
            &mut fingerprint,
            "{}|{}|{}|{}|{:.6}|{}|{}|{}",
            accession_list,
            trace.decoy,
            trace.total_peptide_count,
            trace.passing_peptide_count,
            trace.q_value,
            peptides,
            intensities,
            coverage,
        )
        .expect("failed to write digest");
    }

    fingerprint
}

#[test]
fn protein_grouping_is_deterministic_across_input_orders() {
    let run_count = 4;
    let db = fake_database(vec![
        fake_peptide("PEPA", &["P10000", "P20000"], false),
        fake_peptide("PEPB", &["P10000", "P20000"], false),
        fake_peptide("PEPC", &["P30000"], false),
        fake_peptide("PEPD", &["DECOY_P40000"], true),
    ]);

    let traces = vec![
        fake_trace(0, false, 0.004, &[10.0, 20.0, 0.0, 5.0]),
        fake_trace(1, false, 0.006, &[2.0, 4.0, 6.0, 8.0]),
        fake_trace(2, false, 0.008, &[100.0, 0.0, 50.0, 0.0]),
        fake_trace(3, true, 0.002, &[0.0, 25.0, 25.0, 25.0]),
    ];

    let canonical = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, 0.05);
    let canonical_digest = digest_protein_map(&canonical);

    let mut reversed = traces.clone();
    reversed.reverse();
    let reversed_digest = digest_protein_map(&ProteinQuantTrace::group_by_accession(
        &db, &reversed, run_count, 0.05,
    ));

    let mut rotated = traces.clone();
    rotated.rotate_left(2);
    let rotated_digest = digest_protein_map(&ProteinQuantTrace::group_by_accession(
        &db, &rotated, run_count, 0.05,
    ));

    assert_eq!(canonical_digest, reversed_digest);
    assert_eq!(canonical_digest, rotated_digest);
}

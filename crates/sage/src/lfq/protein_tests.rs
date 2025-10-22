use std::sync::Arc;

use crate::{
    database::{IndexedDatabase, PeptideIx},
    enzyme::Position,
    lfq::{
        quantify_protein_groups, MaxLfqConfig, Peak, PeptideQuantTrace, PrecursorId,
        ProteinQuantTrace,
    },
    ml::matrix::Matrix,
    peptide::Peptide,
};

trait IndexedDatabaseProteinsExt {
    fn proteins(&self, ix: PeptideIx) -> Vec<Arc<str>>;
}

#[test]
fn quantify_protein_groups_converts_peptides() {
    let run_count = 2;
    let db = fake_database(vec![
        fake_peptide("PEPA", &["P12345"], false),
        fake_peptide("PEPB", &["P12345"], false),
    ]);

    let traces = vec![
        fake_trace(0, false, 0.005, &[100.0, 50.0]),
        fake_trace(1, false, 0.007, &[10.0, 0.0]),
    ];

    let max_precursor_q = 0.01;
    let groups = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, max_precursor_q);

    let results = quantify_protein_groups(
        &traces,
        &groups,
        max_precursor_q,
        MaxLfqConfig {
            min_peptides_per_ratio: 1,
            min_samples_for_protein: 1,
            use_global_normalization: false,
            reference_sample: None,
        },
    )
    .expect("protein rollup");

    assert_eq!(results.len(), 1);
    let protein = &results[0];
    assert_eq!(protein.quant.protein_ids, vec!["P12345".to_string()]);
    assert_eq!(protein.quant.peptide_count, 2);
    assert_eq!(protein.quant.sample_coverage.len(), run_count);
    assert!(protein.quant.sample_coverage[0]);
}

#[test]
fn single_peptide_proteins_are_filtered_not_fatal() {
    let run_count = 2;
    let db = fake_database(vec![
        fake_peptide("SINGLE", &["P11111"], false),    // Single-peptide protein
        fake_peptide("MULTI_A", &["P22222"], false),   // Multi-peptide protein
        fake_peptide("MULTI_B", &["P22222"], false),
    ]);

    let traces = vec![
        fake_trace(0, false, 0.005, &[100.0, 50.0]),
        fake_trace(1, false, 0.006, &[200.0, 150.0]),
        fake_trace(2, false, 0.007, &[300.0, 250.0]),
    ];

    let max_precursor_q = 0.01;
    let groups = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, max_precursor_q);
    assert_eq!(groups.len(), 2); // Both proteins initially grouped

    // With min_peptides=2, single-peptide protein should be filtered, not cause error
    let results = quantify_protein_groups(
        &traces,
        &groups,
        max_precursor_q,
        MaxLfqConfig {
            min_peptides_per_ratio: 2,
            min_samples_for_protein: 1,
            use_global_normalization: false,
            reference_sample: None,
        },
    )
    .expect("should not fail on single-peptide proteins");

    // Only the multi-peptide protein should be quantified
    assert_eq!(results.len(), 1);
    let protein = &results[0];
    assert_eq!(protein.quant.protein_ids, vec!["P22222".to_string()]);
    assert_eq!(protein.quant.peptide_count, 2);
}

#[test]
fn peptide_counts_distinguish_total_vs_passing() {
    // Regression test for bug where pre-filtering traces caused
    // total_peptide_count to equal passing_peptide_count
    let run_count = 2;
    let db = fake_database(vec![
        fake_peptide("PEPA", &["P12345"], false),
        fake_peptide("PEPB", &["P12345"], false),
        fake_peptide("PEPC", &["P12345"], false),
    ]);

    let max_precursor_q = 0.01;
    let traces = vec![
        fake_trace(0, false, 0.001, &[100.0, 50.0]),   // Passing (q=0.001 <= 0.01)
        fake_trace(1, false, 0.005, &[200.0, 150.0]),  // Passing (q=0.005 <= 0.01)
        fake_trace(2, false, 0.15, &[300.0, 250.0]),   // Failing (q=0.15 > 0.01)
    ];

    // Pass ALL traces (including high q-value) to group_by_accession
    // This is how the runner SHOULD work after the fix
    let groups = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, max_precursor_q);

    assert_eq!(groups.len(), 1, "Should have 1 protein group");
    let protein = groups.values().next().expect("protein group");

    // Critical assertions: these will FAIL with the buggy pre-filtering
    assert_eq!(
        protein.total_peptide_count, 3,
        "Should count ALL peptides (including those above q-value threshold)"
    );
    assert_eq!(
        protein.passing_peptide_count, 2,
        "Should count only peptides passing q-value threshold"
    );
    assert!(
        protein.total_peptide_count > protein.passing_peptide_count,
        "Bug check: total MUST exceed passing when some peptides have high q-values. \
         If this fails, pre-filtering is removing high q-value traces before aggregation."
    );

    // Verify the passing peptides are correctly identified
    assert_eq!(protein.peptide_indices.len(), 2, "Should track 2 unique passing peptides");
}

impl IndexedDatabaseProteinsExt for IndexedDatabase {
    fn proteins(&self, ix: PeptideIx) -> Vec<Arc<str>> {
        self.peptides[ix.0 as usize].proteins.clone()
    }
}

fn build_fake_database(mut peptides: Vec<Peptide>, generate_decoys: bool) -> IndexedDatabase {
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
        generate_decoys,
        decoy_tag: "DECOY_".to_string(),
    }
}

fn fake_database(peptides: Vec<Peptide>) -> IndexedDatabase {
    build_fake_database(peptides, true)
}

fn fake_database_with_fasta_decoys(peptides: Vec<Peptide>) -> IndexedDatabase {
    build_fake_database(peptides, false)
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
fn fasta_provided_decoys_keep_traces_and_peptides_separate() {
    let run_count = 2;
    let db = fake_database_with_fasta_decoys(vec![
        fake_peptide("PEPTIDE", &["P12345"], false),
        fake_peptide("DECOY", &["P67890"], true),
    ]);
    let traces = vec![
        fake_trace(0, false, 0.004, &[75.0, 25.0]),
        fake_trace(0, true, 0.006, &[30.0, 0.0]),
        fake_trace(1, true, 0.005, &[12.0, 3.0]),
    ];

    let proteins = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, 0.01);

    let mut target_accession = db.proteins(PeptideIx(0));
    target_accession.sort_unstable();
    target_accession.dedup();
    let target = proteins
        .get(&target_accession)
        .expect("target protein missing");
    assert!(!target.decoy);
    assert_eq!(target.accessions, target_accession);
    assert_eq!(target.intensities, vec![75.0, 25.0]);
    assert_eq!(target.total_peptide_count, 1);
    assert_eq!(target.passing_peptide_count, 1);
    assert_eq!(target.peptide_indices, vec![PeptideIx(0)]);
    assert_eq!(target.run_coverage, vec![1, 1]);

    let trace_decoy_accession = vec![Arc::<str>::from(format!(
        "{}{}#TRACE",
        db.decoy_tag, "P12345"
    ))];
    let trace_decoy = proteins
        .get(&trace_decoy_accession)
        .expect("synthetic decoy protein missing");
    assert!(trace_decoy.decoy);
    assert_eq!(trace_decoy.accessions, trace_decoy_accession);
    assert_eq!(trace_decoy.intensities, vec![30.0, 0.0]);
    assert_eq!(trace_decoy.total_peptide_count, 1);
    assert_eq!(trace_decoy.passing_peptide_count, 1);
    assert_eq!(trace_decoy.peptide_indices, vec![PeptideIx(0)]);
    assert_eq!(trace_decoy.run_coverage, vec![1, 0]);

    let fasta_decoy_accession = vec![Arc::<str>::from(format!("{}{}", db.decoy_tag, "P67890"))];
    let fasta_decoy = proteins
        .get(&fasta_decoy_accession)
        .expect("FASTA decoy protein missing");
    assert!(fasta_decoy.decoy);
    assert_eq!(fasta_decoy.accessions, fasta_decoy_accession);
    assert_eq!(fasta_decoy.intensities, vec![12.0, 3.0]);
    assert_eq!(fasta_decoy.total_peptide_count, 1);
    assert_eq!(fasta_decoy.passing_peptide_count, 1);
    assert_eq!(fasta_decoy.peptide_indices, vec![PeptideIx(1)]);
    assert_eq!(fasta_decoy.run_coverage, vec![1, 1]);

    assert_eq!(proteins.len(), 3);
}

#[test]
fn synthetic_trace_decoys_do_not_merge_with_fasta_decoys() {
    let run_count = 2;
    let db = fake_database_with_fasta_decoys(vec![
        fake_peptide("TARGET", &["P12345"], false),
        fake_peptide("FASTA", &["P12345"], true),
    ]);
    let traces = vec![
        fake_trace(0, false, 0.004, &[80.0, 20.0]),
        fake_trace(0, true, 0.005, &[10.0, 0.0]),
        fake_trace(1, true, 0.006, &[5.0, 5.0]),
    ];

    let proteins = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, 0.01);

    let mut target_accession = db.proteins(PeptideIx(0));
    target_accession.sort_unstable();
    target_accession.dedup();
    let target = proteins
        .get(&target_accession)
        .expect("target protein missing");
    assert_eq!(target.intensities, vec![80.0, 20.0]);

    let trace_decoy_accession = vec![Arc::<str>::from(format!(
        "{}{}#TRACE",
        db.decoy_tag, "P12345"
    ))];
    let trace_decoy = proteins
        .get(&trace_decoy_accession)
        .expect("synthetic trace decoy missing");
    assert!(trace_decoy.decoy);
    assert_eq!(trace_decoy.intensities, vec![10.0, 0.0]);

    let fasta_decoy_accession = vec![Arc::<str>::from(format!("{}{}", db.decoy_tag, "P12345"))];
    let fasta_decoy = proteins
        .get(&fasta_decoy_accession)
        .expect("FASTA decoy missing");
    assert!(fasta_decoy.decoy);
    assert_eq!(fasta_decoy.intensities, vec![5.0, 5.0]);

    assert_eq!(proteins.len(), 3);
}

#[test]
fn decoy_traces_for_target_peptides_form_separate_groups() {
    let run_count = 2;
    let db = fake_database(vec![fake_peptide("PEPA", &["P12345"], false)]);
    let traces = vec![
        fake_trace(0, false, 0.004, &[100.0, 25.0]),
        fake_trace(0, true, 0.006, &[12.5, 7.5]),
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
    assert_eq!(target.intensities, vec![100.0, 25.0]);
    assert_eq!(target.total_peptide_count, 1);
    assert_eq!(target.passing_peptide_count, 1);

    let decoy_accession = vec![Arc::<str>::from(format!(
        "{}{}#TRACE",
        db.decoy_tag, "P12345"
    ))];
    let decoy = proteins
        .get(&decoy_accession)
        .expect("decoy protein missing");
    assert!(decoy.decoy);
    assert_eq!(decoy.accessions, decoy_accession);
    assert_eq!(decoy.intensities, vec![12.5, 7.5]);
    assert_eq!(decoy.total_peptide_count, 1);
    assert_eq!(decoy.passing_peptide_count, 1);
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
fn multiple_charge_states_are_quantified_by_maxlfq() {
    let run_count = 2;
    let db = fake_database(vec![
        fake_peptide("PEPA", &["P12345"], false),
        fake_peptide("PEPB", &["P12345"], false),
    ]);
    
    // Create traces with multiple charge states per peptide
    let traces = vec![
        // PEPA with charge 2
        fake_trace_with_precursor(
            PrecursorId::Charged((PeptideIx(0), 2)),
            0,
            false,
            0.001,
            &[100.0, 50.0],
        ),
        // PEPA with charge 3 (additional evidence for same peptide)
        fake_trace_with_precursor(
            PrecursorId::Charged((PeptideIx(0), 3)),
            0,
            false,
            0.002,
            &[80.0, 40.0],
        ),
        // PEPB with charge 2
        fake_trace_with_precursor(
            PrecursorId::Charged((PeptideIx(1), 2)),
            1,
            false,
            0.003,
            &[200.0, 100.0],
        ),
        // PEPB with charge 3
        fake_trace_with_precursor(
            PrecursorId::Charged((PeptideIx(1), 3)),
            1,
            false,
            0.004,
            &[150.0, 75.0],
        ),
    ];

    let max_precursor_q = 0.01;
    let groups = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, max_precursor_q);

    // Verify protein rollup aggregated all charge states
    let mut accession = db.proteins(PeptideIx(0));
    accession.sort_unstable();
    let protein_trace = groups.get(&accession).expect("protein missing");
    assert_eq!(protein_trace.intensities, vec![530.0, 265.0]); // Sum of all traces
    assert_eq!(protein_trace.total_peptide_count, 4); // 4 traces total
    assert_eq!(protein_trace.passing_peptide_count, 2); // 2 unique peptides

    // Now verify MaxLFQ receives all 4 traces (not just 2)
    let results = quantify_protein_groups(
        &traces,
        &groups,
        max_precursor_q,
        MaxLfqConfig {
            min_peptides_per_ratio: 2,
            min_samples_for_protein: 1,
            use_global_normalization: false,
            reference_sample: None,
        },
    )
    .expect("quantification should succeed");
    
    assert_eq!(results.len(), 1);
    let protein = &results[0];
    assert_eq!(protein.quant.protein_ids, vec!["P12345".to_string()]);
    
    // The peptide_count should be 4 (all charge state traces), not 2
    // This verifies that all traces made it through to MaxLFQ
    assert_eq!(protein.quant.peptide_count, 4);
    assert!(protein.quant.sample_coverage[0]);
    assert!(protein.quant.sample_coverage[1]);
}

#[test]
fn high_q_charge_states_are_excluded_from_maxlfq() {
    // Regression test for bug where high-q precursor traces could bypass FDR filtering
    // when a peptide had at least one passing charge state.
    //
    // Scenario: Peptide A has two charge states:
    //   - Charge +2: q=0.002 (PASSING with threshold 0.01)
    //   - Charge +3: q=0.2   (FAILING with threshold 0.01)
    //
    // Expected: Only the passing trace should reach MaxLFQ (peptide_count = 1)
    // Bug behavior: Both traces reach MaxLFQ (peptide_count = 2)

    let run_count = 2;
    let max_precursor_q = 0.01;
    let db = fake_database(vec![fake_peptide("PEPA", &["P12345"], false)]);

    let traces = vec![
        // Passing charge state
        fake_trace_with_precursor(
            PrecursorId::Charged((PeptideIx(0), 2)),
            0,
            false,
            0.002,
            &[100.0, 50.0],
        ),
        // Failing charge state for the same peptide
        fake_trace_with_precursor(
            PrecursorId::Charged((PeptideIx(0), 3)),
            0,
            false,
            0.2,
            &[400.0, 200.0],
        ),
    ];

    let groups = ProteinQuantTrace::group_by_accession(&db, &traces, run_count, max_precursor_q);

    let results = quantify_protein_groups(
        &traces,
        &groups,
        max_precursor_q,
        MaxLfqConfig {
            min_peptides_per_ratio: 1,
            min_samples_for_protein: 1,
            use_global_normalization: false,
            reference_sample: None,
        },
    )
    .expect("protein quantification should succeed");

    assert_eq!(results.len(), 1);
    let protein = &results[0];

    // Critical assertion: verify high-q trace was excluded
    assert_eq!(
        protein.quant.peptide_count, 1,
        "High-q traces must be excluded from MaxLFQ; only the passing trace should contribute"
    );

    // Additional validation
    assert_eq!(protein.total_peptide_count, 2, "Should count both traces");
    assert_eq!(
        protein.passing_peptide_count, 1,
        "Only one peptide passes q-value threshold"
    );
    assert!(protein.quant.sample_coverage.iter().all(|covered| *covered));
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

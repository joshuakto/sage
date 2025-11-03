//! Ensure that we exhaustively visit all fragment ions matching tolerances

use quickcheck_macros::quickcheck;
use sage_core::database::{Builder, IndexedDatabase, PeptideIx};
use sage_core::fasta::Fasta;
use sage_core::lfq::{build_feature_map, LfqSettings, PeptideQuantTrace, PrecursorId};
use sage_core::mass::Tolerance;
use sage_core::ml::retention_alignment::Alignment;
use sage_core::scoring::Feature;
use sage_core::spectrum::{MS1Spectra, Peak, ProcessedSpectrum};

const FASTA: &'static str = r#"
>sp|Q99536|VAT1_HUMAN Synaptic vesicle membrane protein VAT-1 homolog OS=Homo sapiens OX=9606 GN=VAT1 PE=1 SV=2
MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYD
KVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIA
VGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVL
FDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDY
HTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMA
LARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDS
VWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN
"#;

fn mk_database(bucket_size: usize) -> IndexedDatabase {
    let builder = Builder {
        bucket_size: Some(bucket_size),
        fasta: Some("static".into()),
        ..Default::default()
    };
    let fasta = Fasta::parse(FASTA.into(), "rev_", false);

    builder.make_parameters().build(fasta)
}

#[quickcheck]
fn check_all_ions_visited(target_fragment_mz: f32, bucket_size: usize) {
    let database = mk_database(bucket_size.clamp(1, 8192));

    // Map PeptideIx -> number of fragments between 500 & 700 m/z
    // We want to make sure that IndexedDatabase::query hits *all* of them
    let mut expected = vec![0usize; database.peptides.len()];

    let fragment_tol = Tolerance::Da(-100.0, 100.0);
    let (frag_lo, frag_hi) = fragment_tol.bounds(target_fragment_mz);

    for (chunk_idx, chunk) in database.fragments.chunks(database.bucket_size).enumerate() {
        // Check for total ordering by PeptideIx within a chunk
        let mut last = PeptideIx(0);
        for frag in chunk {
            assert!(frag.peptide_index >= last);
            assert!(frag.fragment_mz >= database.buckets()[chunk_idx]);
            if chunk_idx + 1 < database.buckets().len() {
                assert!(frag.fragment_mz <= database.buckets()[chunk_idx + 1]);
            }

            if frag.fragment_mz >= frag_lo && frag.fragment_mz <= frag_hi {
                expected[frag.peptide_index.0 as usize] += 1;
            }
            last = frag.peptide_index;
        }
    }

    let mut visited = vec![0usize; database.peptides.len()];

    // Hit all peptides in database, track how many of the 500-700 fragment m/z's
    // are returned to us by searching the database.
    let query = database.query(1000.0, Tolerance::Da(-5000.0, 5000.0), fragment_tol);

    for fragment in query.page_search(target_fragment_mz) {
        visited[fragment.peptide_index.0 as usize] += 1;
    }

    // Make sure we visited every possible fragment
    assert_eq!(expected, visited);
}

#[test]
fn quantify_feature_map_emits_target_and_decoy_traces() {
    let fasta = ">test|P00001\nPEPTIDERTESTSEQKPEPTIDERTESTSEQK\n".to_string();

    let fasta = Fasta::parse(fasta, "rev_", false);
    let builder = Builder {
        bucket_size: Some(16),
        fasta: Some("inline".into()),
        peptide_min_mass: Some(100.0),
        peptide_max_mass: Some(5000.0),
        generate_decoys: Some(false),
        ..Default::default()
    };

    let database = builder.make_parameters().build(fasta);
    let (target_ix, target_peptide) = database
        .peptides
        .iter()
        .enumerate()
        .find(|(_, peptide)| !peptide.decoy)
        .map(|(idx, peptide)| (PeptideIx(idx as u32), peptide))
        .expect("expected at least one target peptide");

    let neutral_mass = target_peptide.monoisotopic;
    let peptide_len = target_peptide.sequence.len();

    let mut settings = LfqSettings::default();
    settings.spectral_angle = 0.0;

    let features = vec![Feature {
        peptide_idx: target_ix,
        psm_id: 1,
        peptide_len,
        spec_id: "spec-1".into(),
        file_id: 0,
        rank: 1,
        label: 1,
        expmass: neutral_mass,
        calcmass: neutral_mass,
        charge: 1,
        rt: 5.0,
        aligned_rt: 0.5,
        predicted_rt: 0.5,
        delta_rt_model: 0.0,
        ims: 0.0,
        predicted_ims: 0.0,
        delta_ims_model: 0.0,
        delta_mass: 0.0,
        isotope_error: 0.0,
        average_ppm: 0.0,
        hyperscore: 100.0,
        delta_next: 0.0,
        delta_best: 0.0,
        matched_peaks: 10,
        longest_b: 5,
        longest_y: 5,
        longest_y_pct: 0.5,
        missed_cleavages: 0,
        matched_intensity_pct: 1.0,
        scored_candidates: 1,
        poisson: 0.0,
        discriminant_score: 1.0,
        posterior_error: 0.0,
        spectrum_q: 0.001,
        peptide_q: 0.001,
        protein_q: 0.001,
        ms2_intensity: 1000.0,
        fragments: None,
    }];

    let feature_map = build_feature_map(settings, (1, 1), &features);

    let alignments = vec![
        Alignment {
            file_id: 0,
            max_rt: 10.0,
            slope: 1.0,
            intercept: 0.0,
        },
        Alignment {
            file_id: 1,
            max_rt: 10.0,
            slope: 1.0,
            intercept: 0.0,
        },
    ];

    let spectra = MS1Spectra::NoMobility(vec![
        ProcessedSpectrum {
            level: 1,
            id: "quiet".into(),
            file_id: 0,
            scan_start_time: 5.0,
            ion_injection_time: 0.0,
            precursors: vec![],
            peaks: vec![],
            total_ion_current: 0.0,
            base_peak_intensity: 1.0,
            median_peak_intensity: 1.0,
        },
        ProcessedSpectrum {
            level: 1,
            id: "signal-1".into(),
            file_id: 1,
            scan_start_time: 4.95,
            ion_injection_time: 0.0,
            precursors: vec![],
            peaks: vec![
                Peak {
                    intensity: 500.0,
                    mass: neutral_mass,
                },
                Peak {
                    intensity: 375.0,
                    mass: neutral_mass + 11.06,
                },
            ],
            total_ion_current: 875.0,
            base_peak_intensity: 500.0,
            median_peak_intensity: 500.0,
        },
        ProcessedSpectrum {
            level: 1,
            id: "signal-2".into(),
            file_id: 1,
            scan_start_time: 5.0,
            ion_injection_time: 0.0,
            precursors: vec![],
            peaks: vec![
                Peak {
                    intensity: 1_000.0,
                    mass: neutral_mass,
                },
                Peak {
                    intensity: 750.0,
                    mass: neutral_mass + 11.06,
                },
            ],
            total_ion_current: 1_750.0,
            base_peak_intensity: 1_000.0,
            median_peak_intensity: 1_000.0,
        },
        ProcessedSpectrum {
            level: 1,
            id: "signal-3".into(),
            file_id: 1,
            scan_start_time: 5.05,
            ion_injection_time: 0.0,
            precursors: vec![],
            peaks: vec![
                Peak {
                    intensity: 600.0,
                    mass: neutral_mass,
                },
                Peak {
                    intensity: 450.0,
                    mass: neutral_mass + 11.06,
                },
            ],
            total_ion_current: 1_050.0,
            base_peak_intensity: 600.0,
            median_peak_intensity: 600.0,
        },
    ]);

    let quant = feature_map.quantify(&database, &spectra, &alignments);

    let target_key = (PrecursorId::Combined(target_ix), false);
    let decoy_key = (PrecursorId::Combined(target_ix), true);

    let keys: Vec<_> = quant.keys().cloned().collect();

    assert!(
        quant.contains_key(&target_key),
        "missing target trace: {:?}",
        keys
    );
    assert!(
        quant.contains_key(&decoy_key),
        "missing decoy trace: {:?}",
        keys
    );
    assert_eq!(quant.len(), 2, "unexpected number of traces: {:?}", keys);

    let target_trace = quant.get(&target_key).unwrap();
    let decoy_trace = quant.get(&decoy_key).unwrap();

    let verify_trace = |decoy_flag: bool, trace: &PeptideQuantTrace| {
        assert_eq!(trace.peptide, target_ix);
        assert_eq!(trace.decoy, decoy_flag);
        assert_eq!(trace.intensities.len(), alignments.len());
        assert_eq!(trace.reference_file_id, features[0].file_id);
        assert!(trace.intensities.iter().any(|&i| i > 0.0));
        assert!(trace.dot_product.rows > 0 && trace.dot_product.cols > 0);
        assert!(trace.spectral_angle.rows > 0 && trace.spectral_angle.cols > 0);
        assert!(trace.isotope_traces.rows > 0 && trace.isotope_traces.cols > 0);
        assert!(trace.raw_isotope_traces.rows > 0 && trace.raw_isotope_traces.cols > 0);
        assert_eq!(trace.time_warps.len(), alignments.len());
    };

    verify_trace(false, target_trace);
    verify_trace(true, decoy_trace);
}

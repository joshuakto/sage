use sage_core::database::{Builder, PeptideIx};
use sage_core::fasta::Fasta;
use sage_core::lfq::{build_feature_map, LfqSettings, PrecursorId};
use sage_core::ml::matrix::Matrix;
use sage_core::ml::retention_alignment::Alignment;
use sage_core::scoring::Feature;
use sage_core::spectrum::{MS1Spectra, Peak, ProcessedSpectrum};

fn build_database() -> sage_core::database::IndexedDatabase {
    let builder = Builder {
        bucket_size: Some(2),
        peptide_min_mass: Some(0.0),
        peptide_max_mass: Some(2000.0),
        fasta: Some("inline".into()),
        ..Default::default()
    };

    let fasta = Fasta::parse(
        ">TGT\nPEPTIDE\n".into(),
        "rev_",
        false,
    );

    builder.make_parameters().build(fasta)
}

fn feature_for(peptide_idx: PeptideIx, calcmass: f32, file_id: usize, peptide_len: usize) -> Feature {
    Feature {
        peptide_idx,
        psm_id: 1,
        peptide_len,
        spec_id: "scan".into(),
        file_id,
        rank: 1,
        label: 1,
        expmass: calcmass,
        calcmass,
        charge: 1,
        rt: 0.5,
        aligned_rt: 0.5,
        predicted_rt: 0.5,
        delta_rt_model: 0.0,
        ims: 0.0,
        predicted_ims: 0.0,
        delta_ims_model: 0.0,
        delta_mass: 0.0,
        isotope_error: 0.0,
        average_ppm: 0.0,
        hyperscore: 10.0,
        delta_next: 0.0,
        delta_best: 0.0,
        matched_peaks: 0,
        longest_b: 0,
        longest_y: 0,
        longest_y_pct: 0.0,
        missed_cleavages: 0,
        matched_intensity_pct: 0.0,
        scored_candidates: 1,
        poisson: 0.0,
        discriminant_score: 0.0,
        posterior_error: 0.0,
        spectrum_q: 0.0,
        peptide_q: 0.0,
        protein_q: 0.0,
        ms2_intensity: 0.0,
        fragments: None,
    }
}

fn assert_non_empty(matrix: &Matrix) {
    assert!(
        !matrix.data.is_empty(),
        "matrix should contain quantification data"
    );
    assert!(
        matrix.data.iter().any(|value| *value > 0.0),
        "matrix entries should include non-zero values"
    );
}

#[test]
fn quantify_emits_target_and_decoy_traces() {
    let database = build_database();
    assert_eq!(database.peptides.len(), 2, "expected 1 target and 1 decoy peptide");

    let (target_ix, target) = database
        .peptides
        .iter()
        .enumerate()
        .find(|(_, peptide)| !peptide.decoy)
        .map(|(idx, peptide)| (PeptideIx(idx as u32), peptide))
        .expect("missing target peptide");

    let feature = feature_for(
        target_ix,
        target.monoisotopic,
        0,
        target.sequence.len(),
    );

    let settings = LfqSettings {
        combine_charge_states: true,
        ..Default::default()
    };

    let feature_map = build_feature_map(settings, (1, 1), &[feature.clone()]);

    let decoy_mass = target.monoisotopic + 11.06;

    let spectrum_rt = feature.aligned_rt - 0.005;

    let spectrum = ProcessedSpectrum {
        level: 1,
        id: "scan".into(),
        file_id: 0,
        scan_start_time: spectrum_rt,
        ion_injection_time: 0.0,
        precursors: Vec::new(),
        peaks: vec![
            Peak {
                intensity: 120.0,
                mass: target.monoisotopic,
            },
            Peak {
                intensity: 80.0,
                mass: decoy_mass,
            },
        ],
        total_ion_current: 200.0,
    };

    let spectra = MS1Spectra::NoMobility(vec![spectrum]);
    let alignments = vec![Alignment {
        file_id: 0,
        max_rt: 1.0,
        slope: 1.0,
        intercept: 0.0,
    }];

    let quant = feature_map.quantify(&database, &spectra, &alignments);

    let precursor = PrecursorId::Combined(target_ix);

    let target_trace = quant
        .get(&(precursor, false))
        .expect("missing target quant trace");
    assert_eq!(target_trace.peptide, feature.peptide_idx);
    assert!(!target_trace.decoy);
    assert_eq!(target_trace.intensities.len(), alignments.len());
    assert_non_empty(&target_trace.dot_product);
    assert_non_empty(&target_trace.spectral_angle);
    assert_non_empty(&target_trace.isotope_traces);
    assert_eq!(target_trace.reference_file_id, feature.file_id);

    let decoy_trace = quant
        .get(&(precursor, true))
        .expect("missing decoy quant trace");
    assert_eq!(decoy_trace.peptide, feature.peptide_idx);
    assert!(decoy_trace.decoy);
    assert_eq!(decoy_trace.intensities.len(), alignments.len());
    assert_non_empty(&decoy_trace.dot_product);
    assert_non_empty(&decoy_trace.spectral_angle);
    assert_non_empty(&decoy_trace.isotope_traces);
    assert_eq!(decoy_trace.reference_file_id, feature.file_id);
}

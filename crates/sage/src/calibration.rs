use crate::scoring::Feature;
use serde::{Deserialize, Serialize};

/// Quality assessment of mass calibration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CalibrationQuality {
    /// Median mass error < 2 ppm (typical for well-calibrated instruments)
    Good,
    /// Median mass error 2-4 ppm (acceptable but could be improved)
    Fair,
    /// Median mass error > 4 ppm (poor calibration, correction recommended)
    Poor,
}

/// Mode for handling mass calibration issues
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CalibrationMode {
    /// Auto-detect and correct with warnings (default, transparent)
    Auto,
    /// Never correct, expose calibration issues (strict mode for QC)
    Strict,
    /// Quietly correct like MaxQuant (production pipelines)
    Adaptive,
}

impl Default for CalibrationMode {
    fn default() -> Self {
        CalibrationMode::Auto
    }
}

impl CalibrationMode {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "auto" => Some(CalibrationMode::Auto),
            "strict" => Some(CalibrationMode::Strict),
            "adaptive" => Some(CalibrationMode::Adaptive),
            _ => None,
        }
    }
}

/// Comprehensive report of mass calibration assessment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalibrationReport {
    /// Median mass error in ppm (robust to outliers)
    pub median_error_ppm: f32,
    /// Mean mass error in ppm
    pub mean_error_ppm: f32,
    /// Standard deviation of mass errors in ppm
    pub std_dev_ppm: f32,
    /// Median absolute deviation (MAD) in ppm
    pub mad_ppm: f32,
    /// Confidence in calibration assessment (0-1)
    pub confidence: f32,
    /// Quality classification
    pub quality: CalibrationQuality,
    /// Number of high-confidence PSMs used for assessment
    pub num_psms_used: usize,
    /// Systematic offset to be corrected (same as median_error_ppm)
    pub systematic_offset: f32,
    /// Percentage of PSMs within 1 ppm
    pub within_1ppm_pct: f32,
    /// Percentage of PSMs within 3 ppm
    pub within_3ppm_pct: f32,
    /// Percentage of PSMs within 5 ppm
    pub within_5ppm_pct: f32,
    /// Whether correction should be applied
    pub should_correct: bool,
}

/// Apply systematic mass offset correction to features
///
/// This function corrects the delta_mass (ppm error) by subtracting the systematic offset.
/// This is equivalent to recalibrating the precursor m/z values.
///
/// # Arguments
/// * `features` - Mutable slice of features to update with corrected delta_mass
/// * `offset_ppm` - Systematic offset in ppm to correct (typically median error)
///
/// # Returns
/// Number of features corrected
pub fn apply_mz_correction(features: &mut [Feature], offset_ppm: f32) -> usize {
    if offset_ppm.abs() < 0.01 {
        return 0; // No meaningful correction needed
    }

    let mut corrected = 0;
    for feature in features.iter_mut() {
        // Simply subtract the systematic offset from the delta_mass
        // delta_mass is already in ppm, so we just correct it
        feature.delta_mass -= offset_ppm;
        corrected += 1;
    }

    corrected
}

/// Assess mass calibration quality from high-confidence PSMs
///
/// This function analyzes the precursor mass errors from confident PSMs
/// to detect systematic mass offset that might indicate poor instrument calibration.
///
/// # Arguments
/// * `features` - Array of PSM features from database search
/// * `fdr_threshold` - FDR threshold for selecting high-confidence PSMs (typically 0.001)
///
/// # Returns
/// A `CalibrationReport` containing calibration statistics and recommendations
pub fn assess_calibration(features: &[Feature], fdr_threshold: f32) -> CalibrationReport {
    // Filter to high-confidence target PSMs
    let mut mass_errors: Vec<f32> = features
        .iter()
        .filter(|f| f.label == 1 && f.spectrum_q <= fdr_threshold)
        .map(|f| f.delta_mass)
        .collect();

    let num_psms = mass_errors.len();

    if num_psms == 0 {
        // No PSMs available for calibration assessment
        return CalibrationReport {
            median_error_ppm: 0.0,
            mean_error_ppm: 0.0,
            std_dev_ppm: 0.0,
            mad_ppm: 0.0,
            confidence: 0.0,
            quality: CalibrationQuality::Good,
            num_psms_used: 0,
            systematic_offset: 0.0,
            within_1ppm_pct: 0.0,
            within_3ppm_pct: 0.0,
            within_5ppm_pct: 0.0,
            should_correct: false,
        };
    }

    // Sort for median calculation
    mass_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Calculate median (robust to outliers)
    let median = if num_psms % 2 == 0 {
        (mass_errors[num_psms / 2 - 1] + mass_errors[num_psms / 2]) / 2.0
    } else {
        mass_errors[num_psms / 2]
    };

    // Calculate mean
    let mean: f32 = mass_errors.iter().sum::<f32>() / num_psms as f32;

    // Calculate standard deviation
    let variance: f32 = mass_errors
        .iter()
        .map(|x| {
            let diff = x - mean;
            diff * diff
        })
        .sum::<f32>()
        / num_psms as f32;
    let std_dev = variance.sqrt();

    // Calculate MAD (Median Absolute Deviation)
    let mut abs_deviations: Vec<f32> = mass_errors
        .iter()
        .map(|x| (x - median).abs())
        .collect();
    abs_deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad = if num_psms % 2 == 0 {
        (abs_deviations[num_psms / 2 - 1] + abs_deviations[num_psms / 2]) / 2.0
    } else {
        abs_deviations[num_psms / 2]
    };

    // Calculate percentage within tolerance windows
    let within_1ppm = mass_errors.iter().filter(|&&x| x.abs() < 1.0).count();
    let within_3ppm = mass_errors.iter().filter(|&&x| x.abs() < 3.0).count();
    let within_5ppm = mass_errors.iter().filter(|&&x| x.abs() < 5.0).count();

    let within_1ppm_pct = (within_1ppm as f32 / num_psms as f32) * 100.0;
    let within_3ppm_pct = (within_3ppm as f32 / num_psms as f32) * 100.0;
    let within_5ppm_pct = (within_5ppm as f32 / num_psms as f32) * 100.0;

    // Classify quality based on median absolute error
    let quality = match median.abs() {
        x if x < 2.0 => CalibrationQuality::Good,
        x if x < 4.0 => CalibrationQuality::Fair,
        _ => CalibrationQuality::Poor,
    };

    // Calculate confidence based on:
    // 1. Number of PSMs (more = better)
    // 2. Consistency (low MAD relative to median = better)
    let n_confidence = (num_psms as f32 / (num_psms as f32 + 100.0)).min(1.0);
    let consistency = if median.abs() > 0.1 {
        (1.0 - (mad / median.abs())).max(0.0)
    } else {
        0.5 // Low median, moderate confidence
    };
    let confidence = (n_confidence * 0.5 + consistency * 0.5).clamp(0.0, 1.0);

    // Decide if correction should be applied:
    // - Systematic offset > 2 ppm
    // - High confidence (> 0.8)
    // - Enough PSMs (> 100)
    let should_correct = median.abs() > 2.0 && confidence > 0.8 && num_psms >= 100;

    CalibrationReport {
        median_error_ppm: median,
        mean_error_ppm: mean,
        std_dev_ppm: std_dev,
        mad_ppm: mad,
        confidence,
        quality,
        num_psms_used: num_psms,
        systematic_offset: median,
        within_1ppm_pct,
        within_3ppm_pct,
        within_5ppm_pct,
        should_correct,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::database::PeptideIx;
    use crate::scoring::Feature;

    fn make_test_feature(delta_mass: f32, spectrum_q: f32, label: i32) -> Feature {
        Feature {
            peptide_idx: PeptideIx(0),
            psm_id: 0,
            peptide_len: 10,
            spec_id: "test".to_string(),
            file_id: 0,
            rank: 1,
            label,
            expmass: 1000.0,
            calcmass: 1000.0,
            charge: 2,
            rt: 0.0,
            aligned_rt: 0.0,
            predicted_rt: 0.0,
            delta_rt_model: 0.0,
            ims: 0.0,
            predicted_ims: 0.0,
            delta_ims_model: 0.0,
            delta_mass,
            isotope_error: 0.0,
            average_ppm: 0.0,
            hyperscore: 100.0,
            delta_next: 10.0,
            delta_best: 0.0,
            matched_peaks: 10,
            longest_b: 5,
            longest_y: 5,
            longest_y_pct: 0.5,
            missed_cleavages: 0,
            matched_intensity_pct: 0.5,
            scored_candidates: 100,
            poisson: 0.01,
            discriminant_score: 5.0,
            posterior_error: 0.01,
            spectrum_q,
            peptide_q: 0.01,
            protein_q: 0.01,
            ms2_intensity: 1000.0,
            fragments: None,
        }
    }

    #[test]
    fn test_good_calibration() {
        let features: Vec<Feature> = (0..1000)
            .map(|i| make_test_feature((i % 20) as f32 * 0.1 - 1.0, 0.0001, 1))
            .collect();

        let report = assess_calibration(&features, 0.001);

        assert_eq!(report.quality, CalibrationQuality::Good);
        assert!(report.median_error_ppm.abs() < 2.0);
        assert!(report.confidence > 0.0);
    }

    #[test]
    fn test_poor_calibration() {
        let features: Vec<Feature> = (0..1000)
            .map(|i| make_test_feature(5.5 + (i % 20) as f32 * 0.1, 0.0001, 1))
            .collect();

        let report = assess_calibration(&features, 0.001);

        assert_eq!(report.quality, CalibrationQuality::Poor);
        assert!(report.median_error_ppm > 4.0);
        assert!(report.should_correct);
    }

    #[test]
    fn test_no_psms() {
        let features: Vec<Feature> = vec![];
        let report = assess_calibration(&features, 0.001);

        assert_eq!(report.num_psms_used, 0);
        assert_eq!(report.quality, CalibrationQuality::Good);
        assert!(!report.should_correct);
    }

    #[test]
    fn test_filters_decoys() {
        let mut features: Vec<Feature> = (0..500)
            .map(|i| make_test_feature((i % 20) as f32 * 0.1 - 1.0, 0.0001, 1))
            .collect();

        // Add decoys with different distribution
        features.extend(
            (0..500).map(|i| make_test_feature(10.0 + (i % 20) as f32 * 0.1, 0.0001, -1)),
        );

        let report = assess_calibration(&features, 0.001);

        // Should only use targets, so median should be close to 0
        assert!(report.median_error_ppm.abs() < 2.0);
        assert_eq!(report.num_psms_used, 500);
    }

    #[test]
    fn test_calibration_mode_from_str() {
        assert_eq!(
            CalibrationMode::from_str("auto"),
            Some(CalibrationMode::Auto)
        );
        assert_eq!(
            CalibrationMode::from_str("strict"),
            Some(CalibrationMode::Strict)
        );
        assert_eq!(
            CalibrationMode::from_str("adaptive"),
            Some(CalibrationMode::Adaptive)
        );
        assert_eq!(CalibrationMode::from_str("invalid"), None);
    }

    #[test]
    fn test_apply_mz_correction() {
        let mut features: Vec<Feature> = vec![
            make_test_feature(5.0, 0.0001, 1),  // 5 ppm offset
            make_test_feature(5.5, 0.0001, 1),  // 5.5 ppm offset
            make_test_feature(4.5, 0.0001, 1),  // 4.5 ppm offset
        ];

        // Apply 5 ppm correction
        let corrected = apply_mz_correction(&mut features, 5.0);

        assert_eq!(corrected, 3);
        // After correction, delta_mass should be close to 0
        for feature in &features {
            assert!(feature.delta_mass.abs() < 1.0, "Expected corrected delta_mass < 1 ppm, got {}", feature.delta_mass);
        }
    }

    #[test]
    fn test_apply_mz_correction_no_op() {
        let mut features: Vec<Feature> = vec![make_test_feature(0.5, 0.0001, 1)];
        let original_delta = features[0].delta_mass;

        // Apply negligible correction
        let corrected = apply_mz_correction(&mut features, 0.005);

        assert_eq!(corrected, 0);
        assert_eq!(features[0].delta_mass, original_delta);
    }

    #[test]
    fn test_apply_mz_correction_negative_offset() {
        let mut features: Vec<Feature> = vec![make_test_feature(-5.0, 0.0001, 1)];

        // Apply -5 ppm correction
        let corrected = apply_mz_correction(&mut features, -5.0);

        assert_eq!(corrected, 1);
        assert!(features[0].delta_mass.abs() < 1.0);
    }
}

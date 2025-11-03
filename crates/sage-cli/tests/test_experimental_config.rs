/// Test experimental configuration parsing
use sage_cli::input::{
    DeisotopingConfig, ExperimentalConfig, IntensityNormalization, MatchingConfig,
    PeakSelectionConfig, PeakSelectionMode, ScoringConfig,
};

#[test]
fn test_default_experimental_config() {
    let config = ExperimentalConfig::default();

    assert!(!config.enable_experimental_features);
    assert_eq!(config.peak_selection.mode, PeakSelectionMode::Global);
    assert_eq!(config.peak_selection.global_max_peaks, 150);
    assert_eq!(config.scoring.intensity_normalization, IntensityNormalization::None);
    assert_eq!(config.matching.min_matched_peaks, 4);
    assert!(!config.matching.charge_dependent_threshold);
}

#[test]
fn test_binned_peak_selection_config() {
    let json = r#"{
        "enable_experimental_features": true,
        "peak_selection": {
            "mode": "binned",
            "binned_peaks_per_bin": 10,
            "binned_bin_width_da": 100.0
        }
    }"#;

    let config: ExperimentalConfig = serde_json::from_str(json).unwrap();

    assert!(config.enable_experimental_features);
    assert_eq!(config.peak_selection.mode, PeakSelectionMode::Binned);
    assert_eq!(config.peak_selection.binned_peaks_per_bin, 10);
    assert_eq!(config.peak_selection.binned_bin_width_da, 100.0);
}

#[test]
fn test_intensity_normalization_config() {
    let json = r#"{
        "enable_experimental_features": true,
        "scoring": {
            "intensity_normalization": "base_peak"
        }
    }"#;

    let config: ExperimentalConfig = serde_json::from_str(json).unwrap();

    assert!(config.enable_experimental_features);
    assert_eq!(config.scoring.intensity_normalization, IntensityNormalization::BasePeak);
}

#[test]
fn test_combined_features_config() {
    let json = r#"{
        "enable_experimental_features": true,
        "peak_selection": {
            "mode": "binned",
            "binned_peaks_per_bin": 10,
            "binned_bin_width_da": 100.0
        },
        "scoring": {
            "intensity_normalization": "base_peak"
        },
        "matching": {
            "min_matched_peaks": 3,
            "charge_dependent_threshold": true
        }
    }"#;

    let config: ExperimentalConfig = serde_json::from_str(json).unwrap();

    assert!(config.enable_experimental_features);
    assert_eq!(config.peak_selection.mode, PeakSelectionMode::Binned);
    assert_eq!(config.scoring.intensity_normalization, IntensityNormalization::BasePeak);
    assert_eq!(config.matching.min_matched_peaks, 3);
    assert!(config.matching.charge_dependent_threshold);
}

#[test]
fn test_all_normalization_modes() {
    let modes = vec![
        ("none", IntensityNormalization::None),
        ("base_peak", IntensityNormalization::BasePeak),
        ("tic", IntensityNormalization::Tic),
        ("median", IntensityNormalization::Median),
    ];

    for (mode_str, expected) in modes {
        let json = format!(
            r#"{{
                "scoring": {{
                    "intensity_normalization": "{}"
                }}
            }}"#,
            mode_str
        );

        let config: ExperimentalConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config.scoring.intensity_normalization, expected);
    }
}

#[test]
fn test_all_peak_selection_modes() {
    let modes = vec![
        ("global", PeakSelectionMode::Global),
        ("binned", PeakSelectionMode::Binned),
        ("hybrid", PeakSelectionMode::Hybrid),
    ];

    for (mode_str, expected) in modes {
        let json = format!(
            r#"{{
                "peak_selection": {{
                    "mode": "{}"
                }}
            }}"#,
            mode_str
        );

        let config: ExperimentalConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config.peak_selection.mode, expected);
    }
}

#[test]
fn test_partial_config_uses_defaults() {
    let json = r#"{
        "enable_experimental_features": true
    }"#;

    let config: ExperimentalConfig = serde_json::from_str(json).unwrap();

    // Should use default values for unspecified fields
    assert_eq!(config.peak_selection.mode, PeakSelectionMode::Global);
    assert_eq!(config.peak_selection.global_max_peaks, 150);
    assert_eq!(config.scoring.intensity_normalization, IntensityNormalization::None);
}

#[test]
fn test_empty_config_uses_all_defaults() {
    let json = r#"{}"#;

    let config: ExperimentalConfig = serde_json::from_str(json).unwrap();

    assert!(!config.enable_experimental_features);
    assert_eq!(config.peak_selection.mode, PeakSelectionMode::Global);
    assert_eq!(config.matching.min_matched_peaks, 4);
}

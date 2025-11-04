# Experimental Features Configuration Examples

This directory contains example configurations for testing experimental features in Sage.

## Quick Start

1. **Baseline** (`baseline.json`): Current Sage behavior (no experiments)
2. **Experiment 1** (`exp1_binned_peaks.json`): Binned peak selection
3. **Experiment 2** (`exp2_intensity_normalization.json`): Intensity normalization
4. **Experiment 3** (`exp3_combined.json`): Both features combined (recommended)
5. **Experiment 4** (`exp4_lower_threshold.json`): Lower match threshold

## Running Experiments

```bash
# Update database path and mzML paths in each config first

# Run baseline
sage configs/experimental/baseline.json -o results/baseline/

# Run experiments
sage configs/experimental/exp1_binned_peaks.json -o results/exp1/
sage configs/experimental/exp2_intensity_normalization.json -o results/exp2/
sage configs/experimental/exp3_combined.json -o results/exp3/

# Compare results
python3 scripts/compare_experiments.py \
  --baseline results/baseline/results.sage.tsv \
  --exp1 results/exp1/results.sage.tsv \
  --exp2 results/exp2/results.sage.tsv \
  --exp3 results/exp3/results.sage.tsv \
  --maxquant path/to/maxquant_results.txt \
  --output comparison_report.html
```

## Configuration Options

### Peak Selection

```json
"peak_selection": {
  "mode": "binned",              // "global", "binned", or "hybrid"
  "binned_peaks_per_bin": 10,    // peaks per bin (if binned/hybrid)
  "binned_bin_width_da": 100.0   // bin width in Da
}
```

**Options**:
- `global`: Current behavior (top 150 by intensity)
- `binned`: MaxQuant-style (top N per m/z bin)
- `hybrid`: Take max of global and binned

### Intensity Normalization

```json
"scoring": {
  "intensity_normalization": "base_peak"  // "none", "base_peak", "tic", or "median"
}
```

**Methods**:
- `none`: Current behavior (absolute intensities)
- `base_peak`: Normalize by base peak intensity (recommended)
- `tic`: Normalize by total ion current
- `median`: Normalize by median peak intensity

### Match Threshold

```json
"matching": {
  "min_matched_peaks": 3,           // minimum b+y matches (default: 4)
  "charge_dependent_threshold": false  // adjust by charge state
}
```

## Expected Improvements

Based on investigation findings:

| Experiment | Expected Recovery | Confidence | Notes |
|------------|------------------|------------|-------|
| Binned peaks | +5-10% | Medium | Preserves high m/z diagnostic ions |
| Intensity norm | +30-40% | High | PRIMARY FIX for low-intensity spectra |
| Combined | +40-55% | High | Recommended first test |
| Lower threshold | +5-8% | Medium | May increase false positives |

## Validation Checklist

After each experiment:

- [ ] FDR ≤ 1% (CRITICAL!)
- [ ] Peptide recovery vs MaxQuant measured
- [ ] Score distributions look reasonable
- [ ] Runtime acceptable (<20% increase)
- [ ] No obvious false positives

## Parameter Tuning

### Binned Selection Variations

```json
// More conservative
{"binned_peaks_per_bin": 8, "binned_bin_width_da": 100.0}

// More permissive
{"binned_peaks_per_bin": 15, "binned_bin_width_da": 100.0}

// Finer bins
{"binned_peaks_per_bin": 10, "binned_bin_width_da": 50.0}

// Coarser bins
{"binned_peaks_per_bin": 10, "binned_bin_width_da": 150.0}
```

### Normalization Variations

```json
// Test different methods
{"intensity_normalization": "tic"}
{"intensity_normalization": "median"}
```

## Notes

- All experimental features are **opt-in** via configuration
- Default behavior unchanged (experimental.enable_experimental_features = false)
- Safe to test - can always revert to baseline
- Maintain FDR control (<1%) in all experiments

## Resources

- Implementation plan: `EXPERIMENTAL_FEATURES_IMPLEMENTATION_PLAN.md`
- Investigation findings: `COMPLETE_INVESTIGATION_SYNTHESIS.md`
- Quick start guide: `EXPERIMENTAL_FEATURES_QUICKSTART.md`

# Progress

## Current Task
- [x] Extend LFQ quantification outputs to return enriched `PeptideQuantTrace` objects carrying peptide indices, peak metadata, intensities, and per-isotope trace matrices. This enables downstream protein-level aggregation work from the implementation plan.

## Notes
- Updated CLI, FDR, and parquet serializers to consume the new struct without changing the user-facing TSV schema.
- Stored smoothed isotope trace matrices alongside integrated intensities to support future MaxLFQ ratio graph construction.

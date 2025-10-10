# Progress

## 2025-02-14
- Completed Plan Step 1: Added the `PeptideQuantTrace` struct in `lfq.rs` to retain peptide identifiers, decoy flags, MaxLFQ trace matrices, and per-run intensities returned from `FeatureMap::quantify`.
- Updated FDR assignment, CLI writers, and parquet serialization to consume the richer LFQ traces while keeping existing TSV output semantics intact.
- Ensured time-warp parameters and isotope traces are preserved for downstream reuse, preparing the groundwork for upcoming protein-level LFQ aggregation.
- Ran `cargo test` to validate the changes; current tree builds with existing warnings unrelated to this change.

## Next Steps
- Build protein-level grouping utilities that operate on `PeptideQuantTrace` records.
- Extend CLI outputs and parquet schemas to emit protein-level LFQ summaries.
- Address lingering warnings in auxiliary crates when tackling subsequent milestones.

# Progress

## 2025-02-14
- Completed Plan Step 1: Added the `PeptideQuantTrace` struct in `lfq.rs` to retain peptide identifiers, decoy flags, MaxLFQ trace matrices, and per-run intensities returned from `FeatureMap::quantify`.
- Updated FDR assignment, CLI writers, and parquet serialization to consume the richer LFQ traces while keeping existing TSV output semantics intact.
- Ensured time-warp parameters and isotope traces are preserved for downstream reuse, preparing the groundwork for upcoming protein-level LFQ aggregation.
- Ran `cargo test` to validate the changes; current tree builds with existing warnings unrelated to this change.

## Milestones
1. **Step 1 – Peptide trace retention** — ✅ Completed on 2025-02-14
   - Baseline plumbing for richer peptide-level quantitation is merged and verified by the tests above.

2. **Step 2 – Protein grouping helper** — ⏳ Pending
   - Finalize the `ProteinQuantTrace` data model that aggregates multiple `PeptideQuantTrace` entries per protein group.
   - Implement grouping helpers in `lfq.rs` that respect decoy flags and propagate quality metrics.
   - Validate grouping behavior against representative LFQ runs and document any edge-case heuristics.

3. **Step 3 – CLI outputs** — ⏳ Pending
   - Extend the CLI writers to emit protein-level LFQ summaries alongside existing peptide tables.
   - Update parquet schemas and TSV layouts to include aggregated protein metrics and coverage stats.
   - Ensure backward compatibility toggles are available for users who need peptide-only exports.

4. **Step 4 – TMT roll-up** — ⏳ Pending
   - Implement TMT channel roll-up logic that reuses the grouping helper while honoring reporter ion normalization.
   - Add configuration hooks for labeling strategies and missing-channel imputation.
   - Benchmark roll-up performance on existing TMT integration tests to confirm expected runtime envelopes.

5. **Step 5 – Documentation and tests** — ⏳ Pending
   - Refresh CLI, API, and data model documentation to reflect the new protein-level features.
   - Author targeted integration tests covering protein grouping, CLI exports, and TMT roll-up scenarios.
   - Capture release notes and migration guidance once the full LFQ pipeline has stabilized.

## Next Steps
- Execute Step 2 by landing the protein grouping helper, then iterate through Steps 3–5 to deliver the complete protein LFQ pipeline.
- Address lingering warnings in auxiliary crates when tackling subsequent milestones.

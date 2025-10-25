# Progress

## 2025-02-28
- Completed Milestone Step 2 by enabling the protein grouping acceptance tests (`crates/sage/src/lfq/protein_tests.rs`) on every
  `cargo test` run and verifying they pass without requiring the `protein-quant-prototype` feature gate.
- Confirmed `ProteinQuantTrace::group_by_accession` and its helpers remain available through the main `lfq` module so the new
  always-on tests compile cleanly.
- Captured a full `cargo test` run to validate the broader workspace still passes with the previously gated suite enabled.

## 2025-02-21
- Captured current regression coverage showing `PeptideQuantTrace` is exercised through unit tests in `crates/sage/tests/lfq.rs` and `crates/sage/tests/integration.rs`, CLI end-to-end verification in `crates/sage-cli/tests/lfq.rs`, and parquet serialization round-trips in `crates/sage-cloudpath/src/parquet.rs`, giving us confidence that trace data is consistent across core, CLI, and parquet layers.
- Highlighted that `crates/sage/src/lfq/protein_tests.rs` hosts the ignored acceptance tests for `ProteinQuantTrace::group_by_accession`, which align with Milestone Step 2 and will define our readiness for protein-level aggregation once the helper ships.
- Reaffirmed our roadmap: Step 2 (grouping helper) through Step 5 (documentation/tests) remain to unlock a fully featured, well-tested protein quantification module; clearing the ignored acceptance tests is the next gating task before progressing to downstream milestones.

## 2025-02-14
- Completed Plan Step 1: Added the `PeptideQuantTrace` struct in `lfq.rs` to retain peptide identifiers, decoy flags, MaxLFQ trace matrices, and per-run intensities returned from `FeatureMap::quantify`.
- Updated FDR assignment, CLI writers, and parquet serialization to consume the richer LFQ traces while keeping existing TSV output semantics intact.
- Ensured time-warp parameters and isotope traces are preserved for downstream reuse, preparing the groundwork for upcoming protein-level LFQ aggregation.
- Ran `cargo test` to validate the changes; current tree builds with existing warnings unrelated to this change.

## Milestones
1. **Step 1 – Peptide trace retention** — ✅ Completed on 2025-02-14
   - Baseline plumbing for richer peptide-level quantitation is merged and verified by the tests above.

2. **Step 2 – Protein grouping helper** — ✅ Completed on 2025-02-28
   - `ProteinQuantTrace::group_by_accession` now drives the always-on acceptance tests in `crates/sage/src/lfq/protein_tests.rs`,
     confirming grouping semantics across targets, decoys, q-value filtering, and run coverage without feature gates.

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
- Tackle Step 3 deliverables by extending CLI outputs and storage formats with protein-level summaries.
- Continue iterating through Steps 3–5 while addressing outstanding warnings in auxiliary crates during the rollout.

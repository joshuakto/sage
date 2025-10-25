# `sage` crate

This crate hosts the core search and quantification engine that powers Sage. The
`src/lfq.rs` module contains the label-free quantification (LFQ) logic that
rolls peptide-level chromatographic traces up to protein-level summaries.

## Protein grouping and decoy namespaces

During LFQ aggregation we construct deterministic protein groups by
canonicalising the list of accessions associated with each peptide and using
that canonical list as the key inside a `BTreeMap`. The key property of this
approach is that every peptide whose accessions normalise to the same
collection of identifiers is folded into the same protein group accumulator.
This produces stable output across repeated runs, but it also means that we
must take special care to keep different kinds of decoy evidence separated from
one another.

Sage produces two distinct sources of decoy signal during quantification:

* **FASTA-supplied decoys** – sequences that were already labelled as decoys in
  the searched database (for example when the FASTA ship with a built-in decoy
  prefix or when the picked-protein FDR pipeline injects decoy entries).
* **Synthetic trace decoys** – chromatographic traces that Sage synthesises on
  the fly by applying a retention-time shift and an 11 Da mass offset to each
  confidently scored target peptide. These traces are used to estimate and
  control the false-integration rate of the LFQ pipeline.

Both flavours of decoy traces originate from the same underlying peptide
indices, so without any additional namespacing their canonicalised accession
lists would collide. This collision manifested as a shared `BTreeMap` key in
`ProteinQuantTrace::group_by_accession`, merging FASTA decoy peptides and the
synthetic trace decoys that were spawned from their target counterparts. Once
merged, downstream consumers could no longer determine the provenance of the
aggregated decoy protein entry, which made it impossible to perform the
separate decoy accounting required for integration FDR diagnostics.

To prevent this collapse we mint a dedicated namespace for synthetic trace
peptides by appending the `"#TRACE"` suffix after the database-wide decoy tag.
For example, if the database uses `"REV_"` as its decoy prefix, the accessions
for a target peptide such as `sp|P12345|KINASE` would expand as follows:

| Peptide source            | Aggregated accession |
|---------------------------|----------------------|
| Target                    | `sp|P12345|KINASE`   |
| FASTA-provided decoy      | `REV_sp|P12345|KINASE` |
| Synthetic trace-level decoy | `REV_sp|P12345|KINASE#TRACE` |

The extra suffix guarantees that the synthetic traces do not share a key with
FASTA decoys even when both begin with the same prefix. It also keeps the
relationship between a trace decoy and its source target peptide visible at a
glance, which helps when reasoning about why a particular decoy entry appeared
in the protein table. Should we ever need to inspect or filter trace decoys
separately from FASTA decoys, the namespace now allows us to do so with a
simple string match instead of additional book-keeping structures.

## Related code

* `src/lfq.rs` – includes the `ProteinQuantTrace::group_by_accession` helper
  that applies the namespace logic described above.
* `src/lfq/protein_tests.rs` – contains regression tests that assert FASTA and
  trace decoys stay disjoint while we refactor the LFQ implementation.


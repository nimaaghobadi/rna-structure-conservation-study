# Data provenance

For column names, denominators, and interpretation of every public result table, see the [public data dictionary](data_dictionary.md).

The core analysis uses Rfam 15.1 seed alignments RF00001 (712 sequences, 230 columns) and RF00005 (954 sequences, 118 columns), retrieved 2026-10-30. Official sources: https://rfam.org/downloads, https://rfam.org/family/RF00001, and https://rfam.org/family/RF00005.

SHA-256 RF00001: `18f91da6b0d9a668caf3d50849db437584140923e5477f7b21ed59453100f77c`

SHA-256 RF00005: `d1a66517066e6d4fbad3bd382934406df3018ff9bdaf72f22732f4c9740c9505`

Raw Stockholm alignments are not duplicated here. Retrieve the exact seed release from Rfam and verify both checksums before running the pipeline. The paired-versus-unpaired results are descriptive, not phylogenetically corrected. The restricted 16-taxon mitochondrial tRNA-Phe(GAA) follow-up does not establish orthology across the full RF00005 seed.

## Core result tables

- `tables/qc_summary.csv`: alignment and annotation quality-control counts.
- `tables/per_column_metrics.csv`: per-column conservation, Shannon entropy, valid-base fractions, and counts.
- `tables/base_pair_metrics.csv`: per-consensus-pair counts and valid denominators.
- `tables/conservation_by_structure.csv`: per-family paired/unpaired conservation summaries.
- `tables/entropy_by_structure.csv`: per-family paired/unpaired entropy summaries (bits).
- `tables/base_pair_summary.csv`: valid denominators and observed canonical/wobble pair composition.

These tables support the primary paired-versus-unpaired comparison and the public conservation/entropy and base-pair composition figures.

## Sensitivity and retrieval provenance

- `provenance/rfam_seed_download_manifest.json`: retrieval timestamp, official Rfam family endpoints, alignment dimensions, validation metadata, and SHA-256 checksums.
- `tables/sensitivity_by_validity.csv`: paired-minus-unpaired conservation summaries after minimum valid-base coverage thresholds of 0.50, 0.70, and 0.90.
- `tables/sequence_redundancy_sensitivity.csv`: full-seed summaries compared with deterministic 95% and 90% identity-reduced representative sets. This is not phylogenetic correction.

## Supplementary data

The following machine-readable outputs are under [supplementary/tables/](../supplementary/tables/) and are indexed in [the supplementary README](../supplementary/README.md):

- Restricted mitochondrial tRNA-Phe(GAA) locus: 16-record validation, per-column metrics, pair metrics, summary, and Bovidae composition sensitivity.
- R-scape/CaCoFold/9ZTL/FR3D: primary and CaCoFold `.cov`/`.power`/`.helixcov` outputs, the family-level R-scape summary, alignment-block/tRNA-stem summary, the 82-pair Stockholm-to-RF coordinate bridge, complete RF-model-to-PDB residue map, insertion/local-end records, and pair-level FR3D evidence.
- ViennaRNA: full/95%/90% summary statistics and pair-level probabilities for both families.

Raw primary and CaCoFold R-scape `.cov`/`.power`/`.helixcov` outputs are retained under `supplementary/rscape/`; FR3D geometric annotation records are retained under `supplementary/structural/`. Run settings and source checksums are recorded under `supplementary/provenance/`.

## Third-party data boundary

The Rfam seed alignments and 9ZTL mmCIF coordinate file are not redistributed. Their official source locations and checksums are recorded in the public provenance files. NCBI accession-version identifiers for the restricted tRNA-Phe validation records are provided in the supplementary validation table; the underlying nucleotide records remain hosted by NCBI.

Selected scientific and software references are listed in [REFERENCES.md](../REFERENCES.md).

# Data provenance

The analysis uses Rfam 15.1 seed alignments RF00001 (712 sequences, 230 columns) and RF00005 (954 sequences, 118 columns), retrieved 2026-10-30. Sources: https://rfam.org/downloads, https://rfam.org/family/RF00001, and https://rfam.org/family/RF00005.

SHA-256 RF00001: `18f91da6b0d9a668caf3d50849db437584140923e5477f7b21ed59453100f77c`

SHA-256 RF00005: `d1a66517066e6d4fbad3bd382934406df3018ff9bdaf72f22732f4c9740c9505`

Raw Stockholm alignments are not duplicated here. The paired-versus-unpaired results are descriptive, not phylogenetically corrected. The restricted 16-taxon mitochondrial tRNA-Phe(GAA) follow-up does not establish orthology across the full RF00005 seed.

## Curated result tables

- `tables/qc_summary.csv`: alignment and annotation quality-control counts.
- `tables/conservation_by_structure.csv`: per-family paired/unpaired conservation summaries.
- `tables/entropy_by_structure.csv`: per-family paired/unpaired entropy summaries (bits).
- `tables/base_pair_summary.csv`: valid denominators and observed canonical/wobble pair composition.

These are the exact summary tables used for the main descriptive results. The restricted mitochondrial tRNA-Phe(GAA) subset metrics are in `../supplementary/tables/`.

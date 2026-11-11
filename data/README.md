# Data provenance

The primary analysis uses the exact Rfam 15.1 seed alignments for RF00001 (5S rRNA) and RF00005 (tRNA). Both Stockholm files are included in [`data/raw/`](raw/) so readers can inspect and rerun the core analysis without relying on a live database endpoint.

| Family | Rfam accession | Sequences | Alignment columns | File | SHA-256 |
|---|---|---:|---:|---|---|
| 5S rRNA | RF00001 | 712 | 230 | `rf00001_seed.sto` | `18f91da6b0d9a668caf3d50849db437584140923e5477f7b21ed59453100f77c` |
| tRNA | RF00005 | 954 | 118 | `rf00005_seed.sto` | `d1a66517066e6d4fbad3bd382934406df3018ff9bdaf72f22732f4c9740c9505` |

The files were retrieved on 30 October 2026 from the official Rfam family-alignment endpoints. Rfam 15.1 was the release identified on the official downloads page at that time. The [machine-readable download manifest](provenance/rfam_seed_download_manifest.json) records file sizes, checksums, Stockholm validation, and source URLs. The [release-specific reference note](provenance/rfam_reference_snapshot.md) explains the Rfam statistics that are relevant to interpreting this study.

A Rfam seed alignment is a curated representative sample of a family; it is not the complete set of all known family members. The analysis therefore describes these two selected seed alignments rather than making a population-wide claim about all tRNAs or all 5S rRNAs.

## Core result tables

- `tables/qc_summary.csv`: sequence counts, alignment dimensions, gap and ambiguity summaries, and structural-state counts.
- `tables/per_column_metrics.csv`: conservation, Shannon entropy, valid-base fractions, and missing-data counts for each original Stockholm column.
- `tables/conservation_by_structure.csv` and `tables/entropy_by_structure.csv`: paired and explicitly unpaired summaries for each family.
- `tables/base_pair_metrics.csv`: observed nucleotide-pair counts and denominators for each annotated consensus pair.
- `tables/base_pair_summary.csv`: family-level summaries of observed canonical and wobble pair composition.
- `tables/sensitivity_by_validity.csv`: sensitivity to minimum per-column valid-base coverage.
- `tables/coordinate_overlap_audit.csv` and `tables/coordinate_overlap_sensitivity.csv`: the nine RF00005 sequence IDs flagged by the inclusive coordinate-overlap screen and the effect of excluding them in a separate sensitivity analysis. All nine remain in the primary analysis.
- `tables/sequence_redundancy_sensitivity.csv`: comparison with deterministic 95% and 90% identity-reduced representative sets.

The [data dictionary](data_dictionary.md) defines each field and denominator. The [family contrast table](../supplementary/tables/family_contrast_summary.csv) records the descriptive paired-minus-unpaired contrasts and the difference between the two selected seeds.

## Supplementary data

The machine-readable supplementary outputs are indexed in [supplementary/README.md](../supplementary/README.md). They include:

- the restricted 16-record mitochondrial tRNA-Phe(GAA) validation set, its per-column and pair metrics, and taxon-composition sensitivity;
- R-scape primary and CaCoFold `.cov`, `.power`, and `.helixcov` outputs, the coordinate bridge, and the FR3D evidence tables;
- the 9ZTL covariance-model-to-residue map, chain-insertion records, and FR3D base-pair/stacking annotations;
- ViennaRNA/RNAalifold summary statistics and pair-level probabilities.

## External data and reuse

The Rfam seed alignments are included under the CC0 1.0 terms stated by Rfam. The 9ZTL mmCIF and underlying NCBI nucleotide records are not copied into this repository; their official identifiers, source locations, and relevant checksums or accession versions are provided in the supplementary provenance materials. Selected references are listed in [REFERENCES.md](../REFERENCES.md).

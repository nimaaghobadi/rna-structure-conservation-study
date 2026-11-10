# Public data dictionary

This guide defines the fields and denominators used in the committed result tables. The core analysis uses the Rfam 15.1 RF00001 (5S rRNA) and RF00005 (tRNA) seed alignments. Raw Stockholm alignments are not redistributed; retrieve the official files and verify the checksums in [the data provenance guide](README.md).

## Core tables: `data/tables/`

### `qc_summary.csv`

One row per family. `n_sequences` and `alignment_length` describe the input alignment; `ss_cons_length` checks the consensus annotation length; `rf_annotation_present` records whether the Stockholm reference annotation is available. Gap and ambiguous counts/fractions are calculated over all sequence-by-column characters. `n_paired_columns`, `n_unpaired_columns`, `n_insertion_columns`, and `n_unknown_structure_columns` describe the structure-state classification. `n_coordinate_overlap_sequences` counts IDs flagged by the coordinate-overlap audit; these sequences remain in the primary analysis.

### `per_column_metrics.csv`

One row per original Stockholm alignment column.

- `column_1based`: original 1-based Stockholm column.
- `structure_state`: `paired`, `unpaired`, `insertion`, or `unknown`. WUSS paired positions are paired; explicitly unpaired WUSS positions are unpaired; RF-annotation dots identify insertion columns.
- `conservation`: frequency of the most common valid A/C/G/U nucleotide among valid observations in that column.
- `entropy_bits`: Shannon entropy in bits over valid A/C/G/U observations.
- `n_valid_bases`, `valid_fraction`: valid nucleotide count and its fraction of all aligned sequences.
- `gap_count`, `gap_fraction`: gap/dot count and fraction.
- `ambiguous_count`, `ambiguous_fraction`: ambiguous-symbol count and fraction.

Gaps and ambiguous symbols are excluded from the conservation and entropy denominators, but their counts are retained. If no valid nucleotide is present, conservation and entropy are undefined.

### `conservation_by_structure.csv` and `entropy_by_structure.csv`

One row per family and eligible structure state. `count` is the number of columns contributing to the mean, median, and standard deviation. Only paired and explicitly unpaired columns are included; insertion and unknown columns are excluded. These are descriptive summaries across alignment columns, not estimates from independent evolutionary replicates.

### `base_pair_metrics.csv`

One row per consensus pair, identified by `column_i_1based` and `column_j_1based` in the original Stockholm alignment.

- `n_sequences`: all aligned sequences.
- `n_both_valid`: sequences with valid A/C/G/U nucleotides at both partner columns.
- `valid_pair_fraction`: `n_both_valid / n_sequences`.
- `count_AU`, `count_UA`, `count_GC`, `count_CG`, `count_GU`, `count_UG`: counts among valid pair observations.
- `count_noncanonical_or_other`: valid observed pairs outside those six classes.
- `canonical_or_wobble_fraction`: six canonical/wobble classes divided by `n_both_valid`.
- `canonical_or_wobble_fraction_all_sequences`: six classes divided by all aligned sequences; gaps and ambiguous pairs remain in this denominator.

Pairs touching RF-annotation insertion columns are excluded.

### `base_pair_summary.csv`

One row per family, pooling the pair observations over all annotated consensus pairs. `n_annotated_base_pairs` is the number of included consensus pairs; `n_valid_pair_observations` and `n_possible_pair_observations` are the pooled valid and total sequence-by-pair denominators. The six pair-class fractions and `canonical_or_wobble_fraction` use the pooled valid denominator. `canonical_or_wobble_fraction_all_possible` uses all possible sequence-by-pair observations. These pooled observations are not independent evolutionary events.

### Sensitivity tables

- `sensitivity_by_validity.csv`: paired-minus-unpaired conservation after requiring minimum per-column valid fractions of 0.50, 0.70, or 0.90. The paired/unpaired column counts may change as the threshold changes.
- `sequence_redundancy_sensitivity.csv`: full-seed summaries compared with deterministic greedy representative sets at 95% and 90% identity. Pairwise identity uses shared valid A/C/G/U columns and requires at least 50% overlap of the shorter valid sequence. This is a sensitivity analysis, not phylogenetic correction.

## Supplementary tables

### R-scape and FR3D

- `supplementary/tables/rscape_summary.csv`: one row per family. Reports R-scape v2.6.21 input/filtered alignment dimensions, annotated and detected pairs, additional significant pairs outside the supplied annotation, R-scape sensitivity/positive predictive value/F score relative to the supplied annotation, and power summaries. The final three CaCoFold fields report alternative-structure predicted/significant pair counts. These are method outputs, not confirmation of biological truth or experimental validation.
- `rscape_alignment_block_summary.csv`: R-scape covariation counts grouped by nested alignment blocks and tRNA stems; 5S entries are not assigned named helices.
- `rscape_coordinate_mapping.csv`: 82-row bridge from original Stockholm pair coordinates to Rfam RF-model coordinates and R-scape pair statistics.
- `rfam_to_pdb_9ztl_cmalign.csv`: covariance-model position to 9ZTL residue mapping. PDB label-chain and author-chain identifiers are distinct fields.
- `9ztl_chain_insertions_and_local_ends.csv`: chain insertions and local terminal regions outside the RF-model mapping.
- `rscape_fr3d_structural_evidence.csv` and `rscape_fr3d_evidence_summary.csv`: pair-level and grouped comparison of R-scape candidates with FR3D annotations in one 9ZTL model. FR3D non-detection is not proof of absence.

Primary and CaCoFold `.cov`, `.power`, and `.helixcov` files are archived in `supplementary/rscape/`. Alternative CaCoFold Stockholm files are not redistributed because they contain the full Rfam seed sequences.

### ViennaRNA / RNAalifold

- `viennarna_summary.csv`: one row per family and alignment-redundancy subset (full, 95%, 90%). Includes subset sizes, MFE/ensemble statistics, SCI, pair counts, probabilities, and ViennaRNA version.
- `viennarna_pair_probabilities.csv`: pair-level predicted probabilities and their alignment/model coordinates. Predictions are computational outputs, not posterior probabilities that a pair is biologically real and not experimental validation.

### Restricted mitochondrial tRNA-Phe(GAA) locus

- `trna_phe_locus_validation.csv`: accession-versioned validation evidence for the restricted 16-record locus set.
- `trna_phe_column_metrics.csv`, `trna_phe_pair_metrics.csv`: per-column and pair-level descriptive metrics for that locus.
- `trna_phe_locus_summary.json`: subset size, descriptive results, and declared scope.
- `trna_phe_bovid_sensitivity.csv`: sensitivity to selecting one Bovidae representative. This is not a phylogenetic correction.

All paired/unpaired summaries are descriptive for the selected curated seed alignments. They do not establish compensatory evolution, causation, or universal family-level rules.

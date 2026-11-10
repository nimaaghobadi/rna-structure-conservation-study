# Supplementary resources index

This directory contains the machine-readable and method-level materials needed to inspect the optional locus and structural analyses without exposing internal development notes or the research manuscript.

## Restricted mitochondrial tRNA-Phe(GAA) locus

- [Supplementary methods](supplementary_methods.md): scope, validation and limitations for the 16-taxon locus subset.
- [Locus validation table](tables/trna_phe_locus_validation.csv): accession-versioned public-record checks and gene-neighborhood validation.
- [Per-column metrics](tables/trna_phe_column_metrics.csv) and [pair metrics](tables/trna_phe_pair_metrics.csv).
- [Locus summary](tables/trna_phe_locus_summary.json) and [Bovidae representative sensitivity](tables/trna_phe_bovid_sensitivity.csv).
- [Sensitivity figure](figures/trna_phe_locus_sensitivity.svg).

The locus-specific subset is not generalized to all RF00005 tRNAs. Its paired-minus-unpaired conservation contrast is sensitive to coverage and taxon composition and is not phylogenetically corrected.

## R-scape, 9ZTL and FR3D structural context

- [Methods and reproduction instructions](rscape_fr3d_methods.md).
- [R-scape run provenance](provenance/rscape_run_manifest.md), the [family-level R-scape summary](tables/rscape_summary.csv), primary and CaCoFold `.cov`/`.power`/`.helixcov` outputs under [rscape/](rscape/), and the [alignment-block/tRNA-stem summary](tables/rscape_alignment_block_summary.csv).
- [R-scape coordinate bridge](tables/rscape_coordinate_mapping.csv).
- [Complete Rfam covariance-model to 9ZTL residue map](tables/rfam_to_pdb_9ztl_cmalign.csv), [chain insertions/local ends](tables/9ztl_chain_insertions_and_local_ends.csv), and [mapping manifest](provenance/9ztl_cmalign_mapping_manifest.json).
- [FR3D source annotations](structural/) and [FR3D run manifest](provenance/fr3d_9ztl_manifest.json).
- [Pair-level R-scape/FR3D evidence](tables/rscape_fr3d_structural_evidence.csv) and [summary](tables/rscape_fr3d_evidence_summary.csv).

The deposited 9ZTL mmCIF is not redistributed. Its official source and checksum are recorded in the mapping manifest. FR3D annotations refer to one deposited structure model, not a universal RNA-family structure. CaCoFold outputs are computational alternative-structure predictions, not experimental validation; the alternative `.sto` files are not redistributed because they contain the complete Rfam seed sequences.

## ViennaRNA / RNAalifold

- [Methods and results](viennarna_methods.md).
- [Summary table](tables/viennarna_summary.csv), [pair-probability table](tables/viennarna_pair_probabilities.csv), and [run manifest](provenance/viennarna_manifest.json).

ViennaRNA is an optional supplementary analysis. It is not independent of R-scape's sequence-covariation evidence and is not experimental validation.

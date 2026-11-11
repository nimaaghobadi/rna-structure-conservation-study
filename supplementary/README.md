# Supplementary resources

This directory contains the supporting tables, methods, and provenance records for the optional locus and structural analyses. These materials are separate from the primary paired/unpaired conservation comparison.

## Restricted mitochondrial tRNA-Phe(GAA) locus

- [Supplementary methods](supplementary_methods.md): selection criteria, validation evidence, and limitations for the 16-taxon locus.
- [Orthology acceptance record](provenance/trna_phe_orthology_acceptance.json): machine-readable summary of the evidence and its restricted scope.
- [Locus validation table](tables/trna_phe_locus_validation.csv): accession-versioned checks of sequence identity, genome annotation, tRNA-Phe copy number, and local gene order.
- [Per-column metrics](tables/trna_phe_column_metrics.csv) and [pair metrics](tables/trna_phe_pair_metrics.csv).
- [Locus summary](tables/trna_phe_locus_summary.json) and [Bovidae representative sensitivity](tables/trna_phe_bovid_sensitivity.csv).
- [Sensitivity figure](figures/trna_phe_locus_sensitivity.svg).

The locus set is not generalized to all RF00005 tRNAs. Its paired-minus-unpaired conservation contrast changes under coverage and taxon-composition sensitivity checks; no phylogenetically corrected result is claimed.

## Descriptive comparison between the two Rfam seeds

The [family contrast table](tables/family_contrast_summary.csv) reports the paired-minus-unpaired conservation and entropy contrasts for each selected seed, together with the tRNA-minus-5S difference. These are descriptive arithmetic contrasts, not a formal family-level test.

## R-scape, 9ZTL, and FR3D

- [Methods and reproduction instructions](rscape_fr3d_methods.md).
- [R-scape run provenance](provenance/rscape_run_manifest.md) and the [family-level summary](tables/rscape_summary.csv).
- Primary and CaCoFold `.cov`, `.power`, and `.helixcov` outputs under [rscape/](rscape/).
- [R-scape coordinate bridge](tables/rscape_coordinate_mapping.csv) and [alignment-block/tRNA-stem summary](tables/rscape_alignment_block_summary.csv).
- [Rfam covariance-model to 9ZTL residue map](tables/rfam_to_pdb_9ztl_cmalign.csv), [chain insertions/local ends](tables/9ztl_chain_insertions_and_local_ends.csv), and [mapping manifest](provenance/9ztl_cmalign_mapping_manifest.json).
- [FR3D base-pair and stacking annotations](structural/), [FR3D manifest](provenance/fr3d_9ztl_manifest.json), and the [pair-level R-scape/FR3D comparison](tables/rscape_fr3d_structural_evidence.csv) with its [summary](tables/rscape_fr3d_evidence_summary.csv).
- [PyMOL rendering script](../src/render_9ztl_rna_structure.py) and [render provenance](provenance/9ztl_pymol_render_manifest.json) for the optional molecular-context panels.

The 9ZTL mmCIF is retrieved from PDBe and is not redistributed. FR3D annotations describe one deposited structure model; they are not universal interaction assignments for either RNA family. The optional PyMOL panels are regenerated locally rather than stored as binary assets in this repository. CaCoFold structures are exploratory predictions, not experimental validation. The alternative CaCoFold Stockholm outputs are not stored separately because the corresponding seed alignments are already included; regenerate them with the documented R-scape command when needed.

## ViennaRNA / RNAalifold

- [Methods and results](viennarna_methods.md).
- [Summary table](tables/viennarna_summary.csv), [pair-probability table](tables/viennarna_pair_probabilities.csv), and [run manifest](provenance/viennarna_manifest.json).

ViennaRNA provides comparative-folding predictions. It does not provide an independent experimental test of the R-scape results or of the paired/unpaired conservation contrast.

# Sequence Conservation and Structural Constraints in tRNA and 5S rRNA

A descriptive RNA biology study comparing paired and explicitly unpaired consensus positions in two public Rfam 15.1 seed alignments: 5S ribosomal RNA (RF00001) and transfer RNA (RF00005).

## Research snapshot

The analysis summarizes single-site nucleotide conservation, Shannon entropy, and observed canonical/wobble base-pair composition. In the selected seed alignments, paired columns have lower mean single-site conservation and higher mean entropy than explicitly unpaired columns. The contrasts are descriptive; they are not inferential tests across RNA families or independent evolutionary events.

The work does **not** establish compensatory evolution, causation, a universal family-wide rule, or a new statistical method. A separate 16-species vertebrate mitochondrial tRNA-Phe(GAA) follow-up is restricted to that locus and is sensitive to taxon composition and coverage.

## Read the study

- [Main research report](manuscript/manuscript.md)
- [Data provenance and checksums](data/README.md)
- [Supplementary locus-audit methods](supplementary/supplementary_methods.md)

## Repository contents

- `manuscript/`: reader-facing report and selected main figures.
- `supplementary/`: restricted-locus follow-up figure and supporting methods.
- `data/processed/`: machine-readable summaries for the restricted tRNA-Phe(GAA) locus.

This is a curated reader-facing research snapshot. The raw Rfam seed files are not duplicated here; the source, release, retrieval date, and checksums are documented in the data guide. The analysis scripts and fully automated reproduction environment are still being curated and are not yet included, so this repository is not yet a turnkey reproduction package.

## Author and citation

Nima Ghobadi · [ORCID: 0009-0009-6146-3045](https://orcid.org/0009-0009-6146-3045)

Please cite the manuscript and Rfam 15.1 source data when using the results or alignments. A formal publication DOI will be added only after a release has been reviewed and approved.

## Scope and contact

This is an independent research project in RNA Biology and RNA Bioinformatics. Questions, corrections, and reproducibility issues can be raised through [GitHub Issues](https://github.com/nimaaghobadi/rna-structure-conservation-study/issues).

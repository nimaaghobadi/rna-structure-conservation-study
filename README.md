# RNA Structure Conservation — public research resources

This repository contains the data, analysis code, figures, and supplementary evidence for a descriptive comparison of paired and explicitly unpaired positions in two Rfam 15.1 RNA-family seed alignments: RF00001 (5S rRNA) and RF00005 (tRNA).

The research manuscript is maintained separately. This repository is its reader- and reviewer-facing companion: it provides the materials needed to inspect the analyses and follow their data provenance without exposing internal development records.

## What is here

- `data/raw/`: the exact RF00001 and RF00005 Stockholm seed alignments used in the study.
- `data/tables/`: alignment QC, per-column conservation and entropy, consensus-pair composition, and sensitivity summaries.
- `data/provenance/`: the seed-download manifest and a release-specific reference note, including SHA-256 checksums.
- `figures/`: the main paired/unpaired comparison, observed base-pair composition, family consensus diagrams, R-scape/FR3D evidence, and ViennaRNA comparison.
- `src/`: the core downloader and analysis code, plus optional ViennaRNA and 9ZTL/FR3D workflows.
- `scripts/`: scripts that regenerate the committed figures from the public result tables.
- `supplementary/`: the restricted mitochondrial tRNA-Phe(GAA) locus analysis, R-scape outputs, 9ZTL coordinate mapping, FR3D records, ViennaRNA results, and supporting methods.
- `tests/`: focused tests for the core calculations, Stockholm parsing, downloading, and figure generation.
- `REPRODUCIBILITY.md`: installation and reproduction instructions.
- `REFERENCES.md`: selected literature, database, software, and structure references.
- `CITATION.cff`: citation metadata for this repository.

## Scientific scope

The main comparison is descriptive. It summarizes single-column conservation and Shannon entropy at paired and explicitly unpaired consensus positions, and the observed composition of annotated base pairs. The two seed alignments differ in taxonomic composition, diversity, length, and curation; alignment columns are not independent evolutionary replicates. The results are not phylogenetically corrected and do not establish a general difference between all tRNAs and all 5S rRNAs.

The 16-taxon mitochondrial tRNA-Phe(GAA) analysis applies only to that restricted locus. It is not an orthology claim for the full RF00005 family, and its single-site conservation contrast is sensitive to coverage and taxon composition.

R-scape covariation, ViennaRNA predictions, and FR3D geometric annotations describe different kinds of evidence. They are not interchangeable and do not establish compensatory evolution or a causal mechanism. The TN93 inferential workflow is excluded because its convergence acceptance criteria were not met.

## Start with these files

1. Read the [data provenance note](data/README.md) and [reproduction guide](REPRODUCIBILITY.md).
2. Inspect the [paired/unpaired conservation and entropy figure](figures/conservation_entropy.svg) and its source tables.
3. Use the [supplementary index](supplementary/README.md) to inspect the locus-specific, R-scape/FR3D, and ViennaRNA materials.
4. Consult the [data dictionary](data/data_dictionary.md) before interpreting table columns or denominators.

## Figure index

- [Paired/unpaired conservation and entropy](figures/conservation_entropy.svg)
- [Observed canonical/wobble base-pair composition](figures/base_pair_composition.svg)
- [5S rRNA consensus diagram](figures/5s_rrna_consensus.svg) and [tRNA consensus diagram](figures/trna_consensus.svg)
- [R-scape/FR3D structural evidence](figures/rscape_fr3d_evidence.svg)
- [ViennaRNA/RNAalifold comparison](figures/viennarna_consensus_comparison.svg)
- [Restricted mitochondrial tRNA-Phe(GAA) sensitivity](supplementary/figures/trna_phe_locus_sensitivity.svg)

## Data and software licensing

The two Rfam 15.1 seed alignments in `data/raw/` are made available by Rfam under the CC0 1.0 Universal public-domain dedication. The software in `src/`, `scripts/`, and `tests/` is covered by the MIT license in this repository. That MIT license does not automatically apply to figures, derived research tables, manuscript text, PDB/mmCIF data, or other third-party material. See [LICENSE](LICENSE), the [Rfam reference note](data/provenance/rfam_reference_snapshot.md), and [REFERENCES.md](REFERENCES.md) for the relevant sources and attribution.

## Author and contact

Nima Ghobadi · [ORCID: 0009-0009-6146-3045](https://orcid.org/0009-0009-6146-3045)

Please report reproducibility problems or corrections through [GitHub Issues](https://github.com/nimaaghobadi/rna-structure-conservation-study/issues).

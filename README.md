# RNA Structure Conservation: Public Research Resources

This repository provides reader-facing data, code, figures, provenance, and supplementary methods for a descriptive RNA biology study of Rfam 15.1 seed alignments RF00001 (5S rRNA) and RF00005 (tRNA).

**The research manuscript is intentionally not included.** It is maintained separately and will be prepared for the selected preprint or publication venue after scientific and editorial review.

## Repository contents

- `figures/`: primary paired/unpaired conservation and entropy comparison, observed base-pair composition, consensus RNA diagrams, R-scape/FR3D structural context, and ViennaRNA comparison.
- `data/tables/`: curated quality-control, per-column, conservation, entropy, base-pair, coverage-sensitivity, and sequence-redundancy result tables.
- [`data/data_dictionary.md`](data/data_dictionary.md): definitions, denominators, and interpretation notes for the core and supplementary tables.
- `data/provenance/`: machine-readable Rfam seed retrieval manifest and checksums.
- `supplementary/`: restricted mitochondrial tRNA-Phe(GAA) locus results; R-scape raw outputs; 9ZTL coordinate mapping and FR3D records; ViennaRNA result tables; run manifests and reader-facing methods. Start at [the supplementary index](supplementary/README.md).
- `src/`: core Rfam download/analysis code and optional ViennaRNA, Infernal/9ZTL mapping, and FR3D comparison workflows.
- `scripts/`: figure-generation scripts that read the committed public tables.
- `tests/`: focused tests for core metrics, Stockholm parsing, and the committed figure-generation scripts.
- `requirements.txt`: core Python dependencies.
- `requirements-optional.txt`: optional ViennaRNA 2.7.2 dependency for the comparative-folding supplement.
- `REPRODUCIBILITY.md`: core and supplementary reproduction instructions, dependencies, and explicit workflow boundaries.
- `REFERENCES.md`: selected literature, database, software, and structural-record references.
- `CITATION.cff`: machine-readable citation metadata.

## Scientific scope

The paired-versus-unpaired summaries are descriptive for the two selected Rfam seed alignments. They are not phylogenetically corrected tests across RNA families or independent evolutionary events. The analyses do not establish compensatory evolution, causation, a universal family-wide rule, or a novel statistical method.

The 16-taxon mitochondrial tRNA-Phe(GAA) follow-up applies only to that restricted locus, not to the full RF00005 family. Its paired-minus-unpaired conservation contrast is sensitive to coverage and taxon composition.

The TN93 inferential workflow is excluded from the public results because its convergence acceptance gate was not met.

## Start here

1. Read [the data provenance guide](data/README.md) and [reproduction guide](REPRODUCIBILITY.md).
2. Inspect the primary [conservation/entropy figure](figures/conservation_entropy.svg) and the corresponding summary tables.
3. Use [the supplementary index](supplementary/README.md) for locus-specific, R-scape/CaCoFold/FR3D, and ViennaRNA evidence.

## Figure index

- [Primary paired/unpaired conservation and entropy](figures/conservation_entropy.svg)
- [Observed canonical/wobble base-pair composition](figures/base_pair_composition.svg)
- [5S rRNA consensus diagram](figures/5s_rrna_consensus.svg) and [tRNA consensus diagram](figures/trna_consensus.svg)
- [R-scape/FR3D structural evidence](figures/rscape_fr3d_evidence.svg)
- [ViennaRNA/RNAalifold consensus comparison](figures/viennarna_consensus_comparison.svg)
- [Restricted mitochondrial tRNA-Phe(GAA) locus sensitivity](supplementary/figures/trna_phe_locus_sensitivity.svg)

Raw Rfam seed alignments and the 9ZTL mmCIF are not duplicated here. Their official source locations and checksums are documented so readers can retrieve the exact upstream resources. CaCoFold alternative Stockholm files are also not redistributed because they contain the complete Rfam seed sequences; the numeric `.cov`, `.power`, and `.helixcov` outputs are archived.

## Reuse and attribution

The MIT license applies to software source files under `src/`, `scripts/`, and `tests/` only. It does not license manuscript text, figures, Rfam data, PDB/mmCIF data, or other third-party resources. See the scope statement in [LICENSE](LICENSE) and the attributions in [REFERENCES.md](REFERENCES.md).

## Author and contact

Nima Ghobadi · [ORCID: 0009-0009-6146-3045](https://orcid.org/0009-0009-6146-3045)

Corrections and reproducibility issues can be raised through [GitHub Issues](https://github.com/nimaaghobadi/rna-structure-conservation-study/issues).

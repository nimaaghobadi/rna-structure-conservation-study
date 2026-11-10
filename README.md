# RNA Structure Conservation: Public Research Resources

This repository contains the reader-facing resources needed to inspect and reuse selected outputs from a descriptive RNA biology study of Rfam 15.1 seed alignments RF00001 (5S rRNA) and RF00005 (tRNA).

**The research manuscript is intentionally not included.** It will be prepared separately in the required format for the selected preprint or publication venue after the scientific and editorial review is complete.

## What is here

- `figures/`: selected consensus, base-pair composition, covariation/structure-context, and comparative-folding figures.
- `data/tables/`: curated quality-control, per-column, conservation, entropy, and base-pair result tables.
- `supplementary/`: restricted mitochondrial tRNA-Phe(GAA) locus methods, figures, and machine-readable supporting tables, plus the pair-level R-scape/FR3D structural-context audit.
- `src/`: core Rfam download, Stockholm parsing, and descriptive-analysis code.
- `scripts/`: selected figure-generation scripts.
- `tests/`: focused tests for the core metrics and Stockholm parsing.
- `REPRODUCIBILITY.md`: execution instructions, dependencies, and explicit boundaries of what the public pipeline does and does not regenerate.
- `data/README.md`: Rfam source details, release, retrieval date, checksums, and descriptions of the public tables.
- `data/provenance/`: machine-readable Rfam seed download manifest.
- `REFERENCES.md`: selected literature, database, software, and structural-record references supporting the public resources.
- `CITATION.cff`: machine-readable citation metadata for the public repository.

## Scientific scope

The paired-versus-unpaired summaries are descriptive for the selected seed alignments. They are not phylogenetically corrected tests across RNA families or independent evolutionary events. The analyses do not establish compensatory evolution, causation, a universal family-wide rule, or a novel statistical method. The 16-taxon mitochondrial tRNA-Phe(GAA) follow-up applies to that restricted locus only, not to the full RF00005 family.

The TN93 inferential workflow is not part of the public results because its convergence acceptance gate was not met.

## Data and reuse

Raw Rfam seed alignments are not duplicated here; use the official Rfam sources and verify against the checksums in the data guide. The MIT license applies to software source files under `src/`, `scripts/`, and `tests/` only. It does not license manuscript text, figures, Rfam data, PDB/mmCIF data, or third-party resources; see the scope statement in [LICENSE](LICENSE).

## Author and contact

Nima Ghobadi · [ORCID: 0009-0009-6146-3045](https://orcid.org/0009-0009-6146-3045)

Corrections and reproducibility issues can be raised through [GitHub Issues](https://github.com/nimaaghobadi/rna-structure-conservation-study/issues).

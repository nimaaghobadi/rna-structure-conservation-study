# Rfam 15.1 reference and input snapshot

## Release and input files

The two alignments used in this study are the Rfam 15.1 seed alignments for 5S rRNA (RF00001) and tRNA (RF00005). Rfam 15.1 was released in January 2026. The alignments were retrieved from the official Rfam family endpoints on 30 October 2026 and are included unchanged under `data/raw/`.

| Family | Accession | Sequences | Alignment length | File size (bytes) | SHA-256 |
|---|---|---:|---:|---:|---|
| 5S rRNA | RF00001 | 712 | 230 | 203401 | `18f91da6b0d9a668caf3d50849db437584140923e5477f7b21ed59453100f77c` |
| tRNA | RF00005 | 954 | 118 | 141273 | `d1a66517066e6d4fbad3bd382934406df3018ff9bdaf72f22732f4c9740c9505` |

The checksums are repeated in the [download manifest](rfam_seed_download_manifest.json). They identify the exact input files; a later Rfam release or a curator-updated alignment must not be substituted without treating it as a different dataset.

Rfam describes a seed alignment as a curated representative set used to build a family covariance model. It is not the complete set of known family members. The current Rfam release and archived releases are available from the [official downloads page](https://rfam.org/downloads).

## Definitions relevant to this analysis

Rfam provides several sequence- and structure-related statistics. They should not be treated as interchangeable:

- **Sequence conservation (SC):** the frequency of the most common nucleotide in an alignment column.
- **Base-pair conservation (BPcons):** for each annotated consensus pair, the fraction of all aligned sequences whose two residues form a canonical pair. Rfam's canonical set includes AU, CG, and GU; the reverse orientations are equivalent pair representations.
- **Covariation:** a pairwise statistic that evaluates changes across two alignment columns and uses sequence weighting. It is distinct from counting how often an observed pair is canonical.
- **Relative entropy:** a divergence from background nucleotide frequencies. It is not the Shannon entropy calculated in this project.

The core analysis reports the maximum valid A/C/G/U frequency for each column and Shannon entropy in bits. For each annotated pair, it separately reports the six oriented canonical/wobble observations (AU, UA, GC, CG, GU, UG), the fraction among valid pair observations, and the fraction among all aligned sequences. Gaps and ambiguous residues are excluded from the valid-pair denominator but remain in the all-sequence denominator. The [data dictionary](../data_dictionary.md) gives the exact definitions.

Rfam notes that these statistics are used primarily for family curation and quality control, not as stand-alone evidence of a new biological mechanism. The project therefore uses them as descriptive summaries and does not equate them with evolutionary covariation or compensatory substitutions.

## Data reuse

Rfam's [official downloads page](https://rfam.org/downloads) states that its database downloads are available under the CC0 1.0 Universal public-domain dedication. The two seed files are included here with their original content and SHA-256 values. Cite the Rfam database and the exact release when reusing them.

## References

- Ontiveros-Palacios N, Cooke E, Nawrocki EP, et al. (2025). Rfam 15: RNA families database in 2025. *Nucleic Acids Research*, 53(D1), D258–D267. https://doi.org/10.1093/nar/gkae1023
- Rfam. [Downloads and release archive](https://rfam.org/downloads).
- Rfam. [Structure statistics definitions](https://rfam.org/static/documents/structure_statistics.pdf).

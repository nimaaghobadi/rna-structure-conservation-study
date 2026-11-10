# R-scape run provenance

- Tool: R-scape **2.6.21**
- Release: https://github.com/EddyRivasLab/R-scape/releases/tag/v2.6.21
- Verified source archive SHA-256: `5912d0e2861ca8759824597c8a0e8e396bee0d20e15055efa61c048dfb907ede`
- Random seed: `20261030`
- Statistical mode: two-set test (`-s`), separating the supplied consensus pairs from all other possible pairs
- Secondary exploration: `--cacofold`; this is a computational alternative-structure output, not experimental validation
- Target E-value: R-scape default (0.05 in the archived primary `.cov` outputs)

## Reproduction command

After downloading the two Rfam seed Stockholm alignments from their official family endpoints and verifying the SHA-256 values below:

```bash
R-scape -s --cacofold --seed 20261030 --outdir results/rscape/RF00001 data/raw/rf00001_seed.sto
R-scape -s --cacofold --seed 20261030 --outdir results/rscape/RF00005 data/raw/rf00005_seed.sto
```

## Input alignments

| Family | Rfam accession | Sequences × columns | SHA-256 |
|---|---|---:|---|
| 5S rRNA | RF00001 | 712 × 230 | `18f91da6b0d9a668caf3d50849db437584140923e5477f7b21ed59453100f77c` |
| tRNA | RF00005 | 954 × 118 | `d1a66517066e6d4fbad3bd382934406df3018ff9bdaf72f22732f4c9740c9505` |

The raw seed alignments are not redistributed in this repository; download them from [Rfam](https://rfam.org/downloads) and verify these checksums.

## Archived primary outputs

| File | SHA-256 |
|---|---|
| `RF00001_5S_rRNA.cov` | `95a88e7e6bd030a884704626a88f540ea968209a5c38185e940c11c8cba26e0b` |
| `RF00001_5S_rRNA.power` | `1f647861080d27bec80ef4db4f04d830c701d254ae5f72379612a41b96ecd590` |
| `RF00005_tRNA.cov` | `30dea1605234f8efc5d58472f3d52d209e054afb7b019bf93c15937544b58952` |
| `RF00005_tRNA.power` | `2263defebbf8a984faf1c59e3adb94914aa23e5a19fa09594edc832f5a152a57` |
| `RF00001_5S_rRNA.helixcov` | `6bc9e4363846e5878ae3cc878913a05080096d148fbe621d8c2cb49d3e75f617` |
| `RF00005_tRNA.helixcov` | `6589f052023874d146904f04f50ff8e09d18fb0bdfbebdb3e3bd2b67ff8361e8` |

The six primary R-scape files (`.cov`, `.power`, and `.helixcov` for each family) are stored under [supplementary/rscape/](../rscape/). The alignment-block/tRNA-stem aggregation is available as [rscape_alignment_block_summary.csv](../tables/rscape_alignment_block_summary.csv). The pair-coordinate table derived from these outputs is [rscape_coordinate_mapping.csv](../tables/rscape_coordinate_mapping.csv). It preserves the original Stockholm alignment coordinates and the Rfam RF-model coordinates; the full FR3D-reconciled records are in [rscape_fr3d_structural_evidence.csv](../tables/rscape_fr3d_structural_evidence.csv).

## Interpretation limits

R-scape tests covariation relative to a phylogenetic null expectation. It does not test the paired-versus-unpaired difference in mean single-site conservation. Estimated power is conditional on the method's simulation/model assumptions, not a posterior probability that a pair is real. A non-detected pair is not evidence that the pair is absent.

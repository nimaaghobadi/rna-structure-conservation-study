# Supplementary methods: ViennaRNA / RNAalifold comparative folding

## Purpose

ViennaRNA 2.7.2 was used as a complementary computational analysis of consensus secondary-structure folding for the Rfam 15.1 RF00001 (5S rRNA) and RF00005 (tRNA) seed alignments. The predictions were compared post hoc with the Rfam `SS_cons` annotation and the mapped R-scape/FR3D evidence. This is not experimental validation and is not a phylogenetic test of paired-versus-unpaired single-site conservation.

## Method

The Python ViennaRNA binding (2.7.2) was used for comparative MFE folding and partition-function/base-pair probabilities. Predictions were unconstrained: the Rfam `SS_cons` structure was not passed as a folding constraint and was used only for post hoc comparison. Each family was analysed at three redundancy levels: the full seed and deterministic greedy representative alignments at 95% and 90% pairwise identity. Representative selection uses `src.metrics.greedy_identity_representatives`, lexicographic sequence-ID ordering, and a 50% minimum overlap.

IUPAC ambiguity symbols other than A/C/G/U were converted to N only in the ViennaRNA input. Original Rfam Stockholm alignments were not modified. Individual-sequence MFE folds used for the structure conservation index (SCI) exclude sequences containing ambiguity symbols. The exact software version, subset sizes, model settings, and interpretation limits are recorded in [the ViennaRNA manifest](provenance/viennarna_manifest.json).

## Results

| Family | Subset | Sequences | MFE pairs | MFE pairs matching noncrossing SS_cons | SCI | Pairs with probability ≥0.90 |
|---|---|---:|---:|---:|---:|---:|
| 5S rRNA | Full | 712 | 25 | 21/34 | 0.602 | 24 |
| 5S rRNA | 95% identity | 524 | 25 | 21/34 | 0.592 | 22 |
| 5S rRNA | 90% identity | 397 | 25 | 21/34 | 0.572 | 20 |
| tRNA | Full | 954 | 20 | 20/21 | 1.294 | 20 |
| tRNA | 95% identity | 855 | 20 | 20/21 | 1.299 | 20 |
| tRNA | 90% identity | 719 | 20 | 20/21 | 1.303 | 20 |

No crossing pairs were present in either `SS_cons`; the direct comparison therefore includes all 34 5S pairs and all 21 tRNA pairs. The 5S MFE recovered 21/34 annotated pairs and contained four additional pairs (RF model positions 16–63, 79–99, 82–96, and 83–95). These four are also R-scape covariation candidates and have direct FR3D base-pair annotations in the selected 9ZTL model. The tRNA MFE recovered 20/21 annotated pairs and had no additional pair outside `SS_cons`; its omitted pair is RF 13–21, which has no direct FR3D base-pair annotation in this P-site tRNA model.

The summary and pair-level probability outputs are available as [viennarna_summary.csv](tables/viennarna_summary.csv) and [viennarna_pair_probabilities.csv](tables/viennarna_pair_probabilities.csv). The comparison figure is [viennarna_consensus_comparison.svg](../figures/viennarna_consensus_comparison.svg).

## Reproduction

Install the optional dependency pinned in [requirements-optional.txt](../requirements-optional.txt), download and verify the Rfam seed alignments, and run:

```bash
python -m src.download_rfam
python -m src.run_viennarna_analysis
python scripts/plot_viennarna_pair_overlay.py
```

The analysis writes the summary and pair-probability tables to `supplementary/tables/` and a run manifest to `supplementary/provenance/`. The SVG generation script reads only the archived public tables and writes to `figures/`.

## Interpretation limits

- RNAalifold predicts nested secondary structures and does not model WUSS pseudoknots.
- No GSC sequence weighting is applied; full and identity-reduced subsets are shown separately.
- RNAalifold's covariance scoring uses the same alignment variation that underlies R-scape; these are complementary analyses, not independent replications.
- SCI and base-pair probabilities are model outputs, not posterior probabilities that a pair is biologically real and not experimental validation.
- This analysis does not test the paired-versus-unpaired difference in mean single-site conservation and does not establish compensatory evolution.

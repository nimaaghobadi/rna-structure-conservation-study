# Reproduction guide

## Scope and reproducibility boundary

The public repository supports reproduction of the core descriptive analysis and its main figures from the exact Rfam 15.1 seed alignments. It also archives reader-facing intermediate data and selected workflows for the optional ViennaRNA and 9ZTL/FR3D structural-context analyses.

It is not a single-command rebuild of every historical analysis or every rendered asset. In particular, the repository does not redistribute the Rfam seed alignments or the 9ZTL mmCIF, and it does not claim that all external-tool analyses or 3D rendering are reproduced by the core pipeline.

## Core environment

Python 3.10 or newer is recommended. Install the dependencies in `requirements.txt`:

```bash
python -m venv .venv
# activate the environment for your operating system
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Core analysis

From the repository root:

```bash
python -m src.download_rfam
python -m src.run_analysis
pytest -q
```

The downloader retrieves the Rfam seed alignments into `data/raw/`. The analysis writes generated outputs under `results/`. The public repository intentionally does not commit duplicate raw seed alignments; retrieve them from the official Rfam endpoints and verify the SHA-256 checksums in [data provenance](data/README.md).

The curated tables used for the public results are committed under `data/tables/`. These include per-column metrics, paired/unpaired summaries, pair composition, QC, coverage sensitivity, and sequence-redundancy sensitivity.

## Rebuild public figures

The figure scripts use the committed public tables and write to `figures/` or `supplementary/figures/`:

```bash
python scripts/plot_core_conservation_entropy.py
python scripts/plot_rna_biology_figures.py
python scripts/plot_rscape_fr3d_evidence.py
python scripts/plot_viennarna_pair_overlay.py
python scripts/plot_tRNA_Phe_locus_sensitivity.py
```

The core conservation/entropy and base-pair composition figures require only the core Python environment. The R-scape/FR3D and ViennaRNA figures use the committed supplementary tables; they do not rerun the external analyses.

## Optional ViennaRNA / RNAalifold analysis

The archived comparative-folding analysis used ViennaRNA 2.7.2. Install the optional dependency:

```bash
python -m pip install -r requirements-optional.txt
python -m src.download_rfam
python -m src.run_viennarna_analysis
python scripts/plot_viennarna_pair_overlay.py
```

The script runs the full seed and deterministic 95% and 90% identity-reduced alignments for both families. It writes the summary and pair-probability tables to `supplementary/tables/` and the run manifest to `supplementary/provenance/`. Exact settings and interpretation limits are in [the ViennaRNA supplementary methods](supplementary/viennarna_methods.md).

ViennaRNA covariance scoring uses the same alignment variation that underlies R-scape; these are complementary analyses, not independent replications. The predictions are not experimental validation.

## Optional 9ZTL coordinate mapping and FR3D comparison

The structural-context workflow requires Infernal 1.1.5 (`cmalign`) on `PATH`, the core Python dependencies, and the external FR3D Python repository pinned to commit `494a81ccb024abda18fb4c62115a77c6df608472`.

First retrieve the Rfam alignments, then run the official covariance-model mapping:

```bash
python -m src.download_rfam
python -m src.map_rfam_cm_to_pdb
python -m src.update_rscape_structure_coordinates
```

The mapping script downloads the 9ZTL updated mmCIF from PDBe and the official RF00001/RF00005 covariance models from Rfam. The mmCIF is not redistributed; its source and SHA-256 are recorded in `supplementary/provenance/9ztl_cmalign_mapping_manifest.json`.

To regenerate FR3D interaction records, clone the pinned FR3D repository and run:

```bash
python -m src.run_fr3d_annotation --fr3d-root /path/to/fr3d-python
python -m src.compare_rscape_fr3d
```

The source FR3D base-pair and stacking records are retained in `supplementary/structural/`. The pair-coordinate bridge is retained in `supplementary/tables/rscape_coordinate_mapping.csv`; the primary R-scape `.cov` and `.power` outputs and run provenance are under `supplementary/rscape/` and `supplementary/provenance/`. The resulting 82-row FR3D-reconciled table is `supplementary/tables/rscape_fr3d_structural_evidence.csv`.

The R-scape coordinate bridge is an archived derived input built from the original Stockholm coordinates and primary R-scape outputs. The public package preserves it for inspection and downstream mapping; the core two-family pipeline does not regenerate that bridge automatically. Full method, coordinate conventions, source checksums, and interpretation boundaries are documented in [the structural supplementary methods](supplementary/rscape_fr3d_methods.md).

## Restricted mitochondrial tRNA-Phe(GAA) subset

The restricted 16-taxon locus analysis is represented by its accession-versioned validation table, per-column and pair metrics, summary JSON, Bovidae sensitivity table, and figure under `supplementary/`. The source nucleotide records remain hosted by NCBI and are identified in the validation table. This is a locus-specific descriptive subset, not orthology across all RF00005 and not a phylogenetically corrected analysis. See [the locus supplementary methods](supplementary/supplementary_methods.md).

## What is not included or claimed

- Raw Rfam seed alignments and the 9ZTL mmCIF are not redistributed; official sources and checksums are provided.
- The public repository does not include the research manuscript, private status/decision logs, internal audits, failed experiments, or the retired TN93 inferential workflow.
- The two 9ZTL molecular-context PNGs and the full PyMOL rendering workflow are not included in this public release. The public structural materials focus on residue mapping and FR3D evidence records rather than distributing a third-party coordinate file or making a 3D rendering pipeline claim.
- Not every SVG can be regenerated from only `requirements.txt`; optional analyses require the additional tools listed above.
- The paired-versus-unpaired results are descriptive. Alignment columns are not independent evolutionary replicates, and the public results do not establish compensatory evolution, causation, or a universal family-level rule.

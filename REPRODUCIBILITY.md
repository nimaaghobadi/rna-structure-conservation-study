# Reproduction guide

## What can be reproduced

The repository contains the exact Rfam 15.1 seed alignments, the core analysis code, the result tables used in the study, and scripts for regenerating the committed figures. The optional sections describe the additional ViennaRNA and 9ZTL/FR3D analyses.

These instructions cover the analyses represented in this repository, not every historical exploratory run. The 9ZTL mmCIF is retrieved from PDBe rather than redistributed. The tRNA-Phe subset is provided as accession-versioned validation and result tables; the full NCBI candidate-search and curation process is not rerun here.

## Core environment

Python 3.12 is the reference environment used for the archived core analysis. Install the dependencies from the repository root:

```bash
python -m venv .venv
# Activate the environment for your operating system.
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Verify and run the core analysis

The two exact Stockholm inputs are included in `data/raw/`. First verify their checksums:

```bash
sha256sum data/raw/rf00001_seed.sto data/raw/rf00005_seed.sto
```

Expected values are recorded in [data provenance](data/README.md). Then validate the files and run the descriptive analysis:

```bash
python -m src.download_rfam
python -m src.run_analysis
pytest -q
```

The downloader checks the Stockholm format and the expected Rfam 15.1 SHA-256 values. If the family endpoint begins serving a different release, the script stops rather than silently replacing the study inputs. If the family endpoint changes, retrieve the Rfam 15.1 seed files from the [official Rfam downloads page](https://rfam.org/downloads) and verify them against the recorded checksums.

The analysis writes fresh outputs under the ignored local `results/` directory. The curated result tables used by this public release are committed under `data/tables/`; compare regenerated values with those files rather than treating a new local run as an automatic update to the public tables.

## Regenerate the public figures

The figure scripts read the committed tables and write their outputs to `figures/` or `supplementary/figures/`:

```bash
python scripts/plot_core_conservation_entropy.py
python scripts/plot_rna_biology_figures.py
python scripts/plot_rscape_fr3d_evidence.py
python scripts/plot_viennarna_pair_overlay.py
python scripts/plot_tRNA_Phe_locus_sensitivity.py
```

The core conservation/entropy and base-pair composition figures require only the core Python dependencies. The R-scape/FR3D and ViennaRNA figures use archived supplementary tables; they do not rerun the external analyses.

The 5S rRNA and tRNA consensus diagrams are R2R vector outputs from the pinned R-scape run. Re-running the two R-scape commands in the [run manifest](supplementary/provenance/rscape_run_manifest.md) regenerates these diagrams from the included Stockholm inputs; the committed SVG files and their SHA-256 values are listed in that manifest.

## Optional ViennaRNA / RNAalifold analysis

The comparative-folding analysis used ViennaRNA 2.7.2. Install the optional dependency:

```bash
python -m pip install -r requirements-optional.txt
python -m src.run_viennarna_analysis
python scripts/plot_viennarna_pair_overlay.py
```

The analysis uses the full seed and deterministic 95% and 90% identity-reduced alignments. It writes the summary and pair-probability tables under `supplementary/tables/` and a manifest under `supplementary/provenance/`. See [the ViennaRNA methods](supplementary/viennarna_methods.md) for settings and interpretation limits.

RNAalifold and R-scape use variation in the same alignments. Their results are complementary computational evidence, not independent replications or experimental validation.

## Optional 9ZTL mapping and FR3D comparison

The residue-mapping workflow requires Infernal 1.1.5 (`cmalign`) on `PATH`. It retrieves the 9ZTL updated mmCIF and the official Rfam covariance models, verifies their archived SHA-256 checksums, and aligns the relevant deposited RNA chains to the RF00001 and RF00005 covariance models:

```bash
python -m src.map_rfam_cm_to_pdb
python -m src.update_rscape_structure_coordinates
```

The mmCIF is not redistributed. Its official PDBe source and SHA-256 checksum are recorded in [the mapping manifest](supplementary/provenance/9ztl_cmalign_mapping_manifest.json). The mapping uses separate label-chain and author-chain identifiers; see [the structural methods](supplementary/rscape_fr3d_methods.md).

To regenerate the FR3D geometric annotations, check out the [FR3D Python repository](https://github.com/BGSU-RNA/fr3d-python) at commit `494a81ccb024abda18fb4c62115a77c6df608472`, then run:

```bash
python -m src.run_fr3d_annotation --fr3d-root /path/to/fr3d-python
python -m src.compare_rscape_fr3d
```

The archived FR3D base-pair and stacking records, R-scape outputs, coordinate bridge, and pair-level comparison are provided under `supplementary/`. The R-scape run settings and checksums are in [the run manifest](supplementary/provenance/rscape_run_manifest.md). The alternative CaCoFold Stockholm files are not stored separately, to avoid keeping a second copy of the included seed sequences. Regenerate them with the documented R-scape command when their predicted structures need inspection.

### Optional PyMOL structural-context figures

The two 9ZTL molecular-context panels can be regenerated separately from the residue map and FR3D evidence. This step requires PyMOL 2.5.0 (installed separately for your platform) and Pillow 9.0.1. On Linux, the DejaVu Sans font files used by the reference rendering should also be available. First run the 9ZTL mapping workflow above so that the checksum-verified mmCIF is present under the ignored local `data/structures/` directory. Then run from the repository root:

```bash
python -m pip install Pillow==9.0.1
pymol -cq -r src/render_9ztl_rna_structure.py
```

The script checks the mmCIF and FR3D evidence-table checksums before rendering. It writes the PNGs under `supplementary/figures/` and updates [the rendering manifest](supplementary/provenance/9ztl_pymol_render_manifest.json). The reference PNGs are not committed; their checksums and rendering settings are recorded in that manifest. The renderer marks FR3D interaction endpoints only; it does not draw inferred hydrogen bonds.

The CM mapping, FR3D comparison, and PyMOL panels provide structural context. They do not alter the core paired/unpaired conservation calculations and do not establish compensatory evolution.

## Restricted mitochondrial tRNA-Phe(GAA) subset

The 16-taxon locus follow-up is documented through the accession-versioned validation table, per-column and pair metrics, summary JSON, Bovidae sensitivity table, and [supplementary methods](supplementary/supplementary_methods.md). The public data preserve the exact record identifiers and the validation criteria. This is a locus-specific descriptive subset, not an orthology claim for all RF00005 sequences and not a phylogenetically corrected analysis.

## Interpretation limits

- Rfam seed alignments are curated representative sets, not random samples of all natural family members.
- The paired/unpaired summaries are descriptive across alignment columns; the columns are not independent evolutionary replicates.
- The two families differ in taxonomic composition, diversity, length, and curation.
- The tRNA-Phe follow-up applies only to one restricted mitochondrial locus.
- R-scape covariation, ViennaRNA folding predictions, and FR3D geometric annotations answer different questions.
- The analyses do not establish compensatory evolution, causation, a universal family-wide rule, or a new statistical method.

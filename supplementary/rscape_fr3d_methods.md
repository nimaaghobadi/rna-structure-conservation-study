# Supplementary methods: R-scape, Rfam-to-9ZTL mapping, and FR3D comparison

## Purpose and scope

This supplementary analysis links R-scape covariation results for the Rfam 15.1 seed alignments RF00001 (5S rRNA) and RF00005 (tRNA) to residue coordinates in one deposited ribosome structure, PDB entry 9ZTL. It then compares the mapped residue pairs with FR3D geometric base-pair and stacking annotations. This is a structural-context audit of one model, not an evolutionary test and not evidence that every member of either RNA family has the same interactions.

## R-scape inputs and coordinates

R-scape v2.6.21 was run on the original Rfam seed Stockholm alignments with the two-set test (`-s`), random seed 20261030, and `--cacofold` as a secondary structure-exploration output. The exact command, source archive checksum, input checksums, and checksums of the primary `.cov` and `.power` outputs are in [the R-scape run manifest](provenance/rscape_run_manifest.md). The raw primary outputs are retained under [supplementary/rscape/](rscape/).

The R-scape left/right positions in the primary covariance output are 1-based positions in the original input Stockholm alignment. This was checked against the original `SS_cons` annotations: all significant pairs marked as belonging to the supplied structure match the corresponding pairs at the same Stockholm columns (22/22 for RF00001 and 21/21 for RF00005). The shorter effective/tested alignment length reported in the R-scape log is not a coordinate transform.

Stockholm columns map to Rfam covariance-model positions by counting non-insertion symbols in the `#=GC RF` annotation up to the target column. A dot in `#=GC RF` denotes an insertion column and has no RF model coordinate. The 82-row [R-scape coordinate bridge](tables/rscape_coordinate_mapping.csv) preserves the original alignment positions, RF-model positions, pair class, and R-scape statistics.

## Rfam covariance-model mapping to 9ZTL

The primary residue map uses the official Rfam 15.1 covariance models aligned to the deposited RNA chains with Infernal `cmalign` 1.1.5 in global model-alignment mode (`-g`). The updated 9ZTL mmCIF is retrieved from PDBe; the two Rfam covariance models are retrieved from the official Rfam family endpoints. Their URLs and SHA-256 checksums, the mmCIF checksum, chain identities, and mapping summary are recorded in [the mapping manifest](provenance/9ztl_cmalign_mapping_manifest.json).

Distinct PDB/mmCIF identifiers are retained: label-asym chain, author chain, label sequence number, author residue number, insertion code, and component/base identity. The relevant author chains are `b` for 5S rRNA (label chain `EA`) and `5` for the P-site tRNA (label chain `F`). These identifiers are not interchangeable.

The complete per-model-position map is [rfam_to_pdb_9ztl_cmalign.csv](tables/rfam_to_pdb_9ztl_cmalign.csv). Chain insertions and residues outside the Rfam-mapped region are listed separately in [9ztl_chain_insertions_and_local_ends.csv](tables/9ztl_chain_insertions_and_local_ends.csv). The map contains 116/120 mapped 5S model positions (four model deletions relative to this chain) and 71/71 mapped tRNA model positions. All mapped positions have model-1 atoms. The deposited tRNA chain contains an insertion relative to the 71-position RF model; the original component identifiers for modified nucleotides are retained.

## FR3D geometric annotations

FR3D Python was pinned to commit `494a81ccb024abda18fb4c62115a77c6df608472`. It was used to annotate base-pair and stacking interactions in 9ZTL for author chains `b` and `5`. FR3D unit identifiers use author-chain and author-residue numbering; residue pairs were compared using the mapped author identifiers. Reverse-direction duplicates were treated as one unordered residue pair. Leontis–Westhof interaction codes are retained in the machine-readable table.

The FR3D source outputs are retained in [supplementary/structural/](structural/), and the run metadata are in [the FR3D manifest](provenance/fr3d_9ztl_manifest.json). The 82-row [R-scape/FR3D evidence table](tables/rscape_fr3d_structural_evidence.csv) contains the original Stockholm coordinates, RF-model coordinates, R-scape statistics where applicable, 9ZTL chain/residue identifiers, mapping status, and FR3D base-pair/stacking observations.

## Summary of the structural comparison

| Family and R-scape category | Total pairs | Direct FR3D base-pair annotation | Stacking only | No direct pair/stacking annotation | Not assessable |
|---|---:|---:|---:|---:|---:|
| 5S rRNA: consensus-annotated, R-scape significant | 22 | 19 | 0 | 3 | 0 |
| 5S rRNA: consensus-annotated, R-scape not detected | 12 | 4 | 0 | 8 | 0 |
| 5S rRNA: additional covariation candidates | 9 | 7 | 0 | 1 | 1 |
| tRNA: consensus-annotated, R-scape significant | 21 | 20 | 0 | 1 | 0 |
| tRNA: additional covariation candidates | 18 | 1 | 4 | 13 | 0 |

These categories describe annotations in the selected 9ZTL model. A pair not detected by R-scape can still have a direct FR3D annotation in this structure. Conversely, a covariation candidate is not automatically a geometrically annotated pair.

## Reproduction commands

The Rfam seed alignments are first downloaded using the core downloader. Install Infernal 1.1.5 (`cmalign`) and make it available on `PATH`. Then run:

```bash
python -m src.download_rfam
python -m src.map_rfam_cm_to_pdb
python -m src.update_rscape_structure_coordinates
```

For FR3D, clone the [BGSU-RNA FR3D Python repository](https://github.com/BGSU-RNA/fr3d-python) at the pinned commit and pass its local directory:

```bash
python -m src.run_fr3d_annotation --fr3d-root /path/to/fr3d-python
python -m src.compare_rscape_fr3d
```

The coordinate bridge table and primary R-scape `.cov`/`.power` outputs are archived as data inputs. The mapping and FR3D scripts regenerate the residue-level map and geometric evidence table from the downloaded structure/model resources and these archived pair coordinates.

## Interpretation limits

- FR3D labels are geometric interaction annotations in one deposited structure model; they do not establish a universal interaction across an RNA family.
- Stacking-only is not direct base-pair support.
- Absence of a direct FR3D annotation is not proof that an interaction never occurs; conformation, model quality, and annotation criteria can affect detection.
- The tRNA comparison concerns one P-site tRNA chain in 9ZTL, not all tRNAs represented in RF00005.
- A CM alignment is a model-based homology map, not experimental validation of a base pair.
- This structural comparison does not demonstrate compensatory evolution, causation, or a new biological mechanism.

## Structural source

The deposited structural reference is [PDB 9ZTL](https://doi.org/10.2210/pdb9ZTL/pdb). The mmCIF coordinate file is not redistributed here; its official source and checksum are recorded in the mapping manifest. FR3D methodology is described by Sarver et al. (2008), and Rfam/Infernal/R-scape references are listed in [REFERENCES.md](../REFERENCES.md).

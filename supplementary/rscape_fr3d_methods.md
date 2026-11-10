# Supplementary methods: R-scape coordinate mapping and FR3D structural comparison

## Purpose and scope

This supplementary analysis links R-scape covariation results for the Rfam 15.1 seed alignments RF00001 (5S rRNA) and RF00005 (tRNA) to residue coordinates in one deposited ribosome structure, PDB entry 9ZTL. It then compares the mapped residue pairs with FR3D geometric base-pair and stacking annotations. This is a structural-context audit of one model, not an evolutionary test and not evidence that every member of either RNA family has the same interactions.

## Coordinate mapping

R-scape positions in the primary covariance output are 1-based positions in the original input Stockholm alignment. This was checked against the original `SS_cons` annotations for both families. Stockholm columns map to Rfam covariance-model positions by counting non-insertion symbols in the `#=GC RF` annotation up to the target column; insertion columns have no covariance-model coordinate.

The Rfam covariance models were aligned to the 9ZTL RNA chains with Infernal `cmalign` 1.1.5. The mapping preserves distinct PDB/mmCIF identifiers: label-asym chain, author chain, label sequence number, author residue number, insertion code, and component/base identity. The relevant author chains are `b` for 5S rRNA (label chain `EA`) and `5` for the P-site tRNA (label chain `F`). These chain identifiers are not interchangeable.

## FR3D comparison

FR3D Python was used to annotate base-pair and stacking interactions in 9ZTL. The implementation was pinned to commit `494a81ccb024abda18fb4c62115a77c6df608472`. FR3D unit identifiers use author-chain and author-residue numbering; residue pairs were therefore compared using the mapped author identifiers. Reverse-direction duplicates were treated as one unordered residue pair. Reported Leontis–Westhof interaction codes are retained in the machine-readable table.

The resulting table contains 82 mapped R-scape pair records and includes the original Stockholm coordinates, Rfam model coordinates, R-scape statistics where applicable, 9ZTL chain/residue identifiers, mapping status, and FR3D base-pair/stacking observations.

## Summary of the structural comparison

| Family and R-scape category | Total pairs | Direct FR3D base-pair annotation | Stacking only | No direct pair/stacking annotation | Not assessable |
|---|---:|---:|---:|---:|---:|
| 5S rRNA: consensus-annotated, R-scape significant | 22 | 19 | 0 | 3 | 0 |
| 5S rRNA: consensus-annotated, R-scape not detected | 12 | 4 | 0 | 8 | 0 |
| 5S rRNA: additional covariation candidates | 9 | 7 | 0 | 1 | 1 |
| tRNA: consensus-annotated, R-scape significant | 21 | 20 | 0 | 1 | 0 |
| tRNA: additional covariation candidates | 18 | 1 | 4 | 13 | 0 |

These categories describe annotations in the selected 9ZTL model. A pair not detected by R-scape can still have a direct FR3D annotation in this structure. Conversely, a covariation candidate is not automatically a geometrically annotated pair.

## Interpretation limits

- FR3D labels are geometric interaction annotations in one deposited structure model; they do not establish a universal interaction across an RNA family.
- Stacking-only is not direct base-pair support.
- Absence of a direct FR3D annotation is not proof that an interaction never occurs; conformation, model quality, and annotation criteria can affect detection.
- The tRNA comparison concerns one P-site tRNA chain in 9ZTL, not all tRNAs represented in RF00005.
- The structural comparison does not demonstrate compensatory evolution, causation, or a new biological mechanism.

## Data and references

The pair-level results are provided in [rscape_fr3d_structural_evidence.csv](tables/rscape_fr3d_structural_evidence.csv). The deposited structural reference is [PDB 9ZTL](https://doi.org/10.2210/pdb9ZTL/pdb). FR3D methodology is described by Sarver et al. (2008), cited in [REFERENCES.md](../REFERENCES.md). The FR3D implementation is available from the [BGSU-RNA FR3D Python repository](https://github.com/BGSU-RNA/fr3d-python); its pinned source revision is recorded above. Rfam-to-structure mapping uses Infernal 1.1.5, also listed in [REFERENCES.md](../REFERENCES.md).

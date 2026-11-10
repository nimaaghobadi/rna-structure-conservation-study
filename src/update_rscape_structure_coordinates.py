"""Populate R-scape pair endpoints with validated Infernal RF-to-9ZTL residue IDs."""
from __future__ import annotations

import csv
from pathlib import Path

from src.map_rfam_cm_to_pdb import ROOT

PAIR_TABLE = ROOT / "supplementary" / "tables" / "rscape_coordinate_mapping.csv"
RESIDUE_MAP = ROOT / "supplementary" / "tables" / "rfam_to_pdb_9ztl_cmalign.csv"
MANIFEST = ROOT / "supplementary" / "provenance" / "9ztl_cmalign_mapping_manifest.json"

ENDPOINT_FIELDS = [
    "label_seq_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_comp_id",
    "normalized_base", "has_model1_atoms", "mapping_status",
]


def main() -> None:
    import json

    with RESIDUE_MAP.open(encoding="utf-8", newline="") as handle:
        residue_rows = list(csv.DictReader(handle))
    mapping = {
        (row["family_accession"], int(row["rf_model_position"])): row
        for row in residue_rows
    }
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    summaries = {item["accession"]: item for item in manifest["summaries"]}

    with PAIR_TABLE.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        original_fields = list(reader.fieldnames or [])
        pair_rows = list(reader)

    extra_fields = []
    for side in ("left", "right"):
        for field in ENDPOINT_FIELDS:
            name = f"structure_{side}_{field}"
            if name not in original_fields:
                extra_fields.append(name)
    for field in (
        "pdb_rf_model_positions_mapped", "pdb_rf_model_length",
        "pdb_model1_atom_residues", "pdb_cif_sha256",
    ):
        if field not in original_fields:
            extra_fields.append(field)
    fields = original_fields + extra_fields
    if len(pair_rows) != 82:
        raise ValueError(f"Expected 82 archived R-scape pair rows, found {len(pair_rows)}")

    for pair in pair_rows:
        family = pair["family"]
        accession = family if family.startswith("RF") else (
            "RF00001" if family.lower().startswith("5s") else "RF00005"
        )
        summary = summaries[accession]
        endpoints = {}
        for side in ("left", "right"):
            position = int(pair[f"rf_{side}_model_pos"])
            residue = mapping.get((accession, position))
            endpoints[side] = residue
            for field in ENDPOINT_FIELDS:
                target = f"structure_{side}_{field}"
                if residue is None or residue.get("label_seq_id", "") == "":
                    pair[target] = ""
                else:
                    pair[target] = residue.get(field, "")
            pair[f"structure_{side}_residue"] = (
                residue.get("auth_seq_id", "")
                if residue and residue.get("label_seq_id", "") else ""
            )

        mapped = [
            residue for residue in endpoints.values()
            if residue and residue.get("label_seq_id", "")
        ]
        all_mapped = len(mapped) == 2
        atoms_present = all(
            residue.get("has_model1_atoms", "").lower() == "true"
            for residue in mapped
        ) if all_mapped else False
        if all_mapped and atoms_present:
            pair["coordinate_status"] = (
                "ALIGNMENT_AND_RF_MAPPED;Rfam_COVARIANCE_MODEL_TO_PDB_MAPPED;"
                "BOTH_RESIDUES_HAVE_MODEL1_ATOMS"
            )
        elif all_mapped:
            pair["coordinate_status"] = (
                "ALIGNMENT_AND_RF_MAPPED;Rfam_COVARIANCE_MODEL_TO_PDB_MAPPED;"
                "ONE_OR_BOTH_RESIDUES_UNRESOLVED"
            )
        else:
            pair["coordinate_status"] = (
                "ALIGNMENT_AND_RF_MAPPED;Rfam_COVARIANCE_MODEL_TO_PDB_PARTIAL;"
                "ONE_OR_BOTH_RF_MODEL_POSITIONS_ABSENT_FROM_9ZTL_CHAIN"
            )
        pair["pdb_mapping_status"] = (
            f"Infernal {manifest['infernal_version']} cmalign against Rfam "
            f"{manifest['rfam_release']} CM; {summary['mapped_model_positions']}/"
            f"{summary['cm_model_length']} RF model positions mapped"
        )
        pair["pdb_rf_model_positions_mapped"] = summary["mapped_model_positions"]
        pair["pdb_rf_model_length"] = summary["cm_model_length"]
        pair["pdb_model1_atom_residues"] = summary["model1_atom_residues"]
        pair["pdb_modelled_residues"] = summary["model1_atom_residues"]
        pair["pdb_cif_sha256"] = manifest["cif_sha256"]

    with PAIR_TABLE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(pair_rows)
    print(f"Updated {len(pair_rows)} R-scape pair rows with RF-to-PDB residue identifiers.")


if __name__ == "__main__":
    main()

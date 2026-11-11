"""Align deposited 9ZTL RNA chains to the official Rfam covariance models.

This is the primary RF-model-to-PDB mapping path. Infernal cmalign aligns each
deposited chain sequence to the family CM; the cmalign Stockholm RF annotation
distinguishes model match positions from insertions and local ends.

Requires Infernal (cmalign) on PATH. The 9ZTL mmCIF and chain extraction helpers
are provided by src.structure_mapping_utils.

Run:
    python -m src.map_rfam_cm_to_pdb
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import urllib.request
from pathlib import Path

import pandas as pd
from Bio import AlignIO

from src.structure_mapping_utils import (
    ROOT, FAMILIES, extract_chain, load_cif, obtain_cif,
    best_seed_alignment,
)

RFAM_RELEASE = "15.1"
EXPECTED_CIF_SHA256 = "4401a250c6881fee8bc1f688d7a559b724e7ab7b18bda7d057cf0ce80f14480e"
EXPECTED_CM_SHA256 = {
    "RF00001": "390d717535fedd7d8c6a35e4a1c9f923fa5daa20b9926ced948fdc71bb8f9122",
    "RF00005": "7698e524c9ae5b121186fe1e6608d7ab92bff31ef51e886186d7585a66090865",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def get_cif() -> tuple[Path, str]:
    path = ROOT / "data" / "structures" / "9ztl_updated.cif"
    if path.is_file() and path.stat().st_size > 1_000_000:
        source = "cached PDBe updated mmCIF"
    else:
        path, source = obtain_cif(None)
    observed = sha256(path)
    if observed != EXPECTED_CIF_SHA256:
        raise RuntimeError(
            f"9ZTL mmCIF checksum mismatch: {observed}; "
            f"expected {EXPECTED_CIF_SHA256}. Retrieve the pinned PDBe updated file."
        )
    return path, source


def download_cm(accession: str) -> tuple[Path, str]:
    directory = ROOT / "data" / "structures" / "rfam_cm"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{accession}.cm"
    url = f"https://rfam.org/family/{accession}/cm"
    if not path.exists():
        request = urllib.request.Request(
            url, headers={"User-Agent": "rna-structure-conservation/1.0"}
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
        if b"INFERNAL1" not in payload[:200]:
            raise ValueError(f"Rfam CM endpoint did not return an Infernal model for {accession}")
        path.write_bytes(payload)
    observed = sha256(path)
    expected = EXPECTED_CM_SHA256[accession]
    if observed != expected:
        raise RuntimeError(
            f"{accession} covariance-model checksum mismatch: {observed}; "
            f"expected Rfam {RFAM_RELEASE} checksum {expected}."
        )
    return path, url


def infernal_version() -> str:
    result = subprocess.run(
        ["cmalign", "-h"], capture_output=True, text=True, check=False
    )
    text = result.stdout + "\n" + result.stderr
    match = re.search(r"INFERNAL\s+([0-9.]+)", text)
    if not match:
        raise RuntimeError("Could not determine Infernal version from cmalign -h")
    return match.group(1)


def parse_rf_annotation(path: Path) -> str:
    chunks = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#=GC RF "):
                fields = line.split()
                if len(fields) >= 3:
                    chunks.append(fields[2])
    if not chunks:
        raise ValueError(f"No #=GC RF annotation in {path}")
    return "".join(chunks)


def write_fasta(path: Path, name: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{name}\n{sequence}\n", encoding="utf-8")


def align_chain_to_cm(
    accession: str,
    spec: dict,
    cif: dict,
    cm_path: Path,
    cm_url: str,
    infernal: str,
) -> tuple[list[dict], list[dict], dict]:
    rows = extract_chain(cif, spec)
    chain_sequence = "".join(row["normalized_base"] for row in rows)
    if not chain_sequence or not set(chain_sequence).issubset(set("ACGU")):
        raise ValueError(f"Invalid normalized RNA sequence for {accession}")
    region_rows = []
    for row in rows:
        try:
            auth_number = int(row["auth_seq_id"])
        except (TypeError, ValueError):
            continue
        if spec["rfam_auth_start"] <= auth_number <= spec["rfam_auth_end"]:
            region_rows.append(row)
    expected_region_length = spec["rfam_auth_end"] - spec["rfam_auth_start"] + 1
    if len(region_rows) != expected_region_length:
        raise ValueError(
            f"Rfam mapped region has {len(region_rows)} residues, expected "
            f"{expected_region_length} for {accession}"
        )
    region_sequence = "".join(row["normalized_base"] for row in region_rows)

    fasta_path = ROOT / "data" / "processed" / f"9ztl_{accession}_chain.fasta"
    alignment_path = ROOT / "data" / "processed" / f"9ztl_{accession}_cmalign.sto"
    log_path = ROOT / "data" / "processed" / f"9ztl_{accession}_cmalign.log"
    write_fasta(fasta_path, f"9ZTL_{spec['label_asym_id']}", region_sequence)

    command = [
        "cmalign", "-g", "--informat", "FASTA", "--noprob",
        "-o", str(alignment_path), str(cm_path), str(fasta_path),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    log_path.write_text(
        "COMMAND: " + " ".join(command) + "\n\nSTDOUT:\n" + result.stdout
        + "\nSTDERR:\n" + result.stderr,
        encoding="utf-8",
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"cmalign failed for {accession} (exit {result.returncode}); see {log_path}"
        )

    alignment = AlignIO.read(str(alignment_path), "stockholm")
    if len(alignment) != 1:
        raise ValueError(f"Expected one aligned chain sequence for {accession}")
    aligned_sequence = str(alignment[0].seq).upper().replace("T", "U")
    rf = parse_rf_annotation(alignment_path)
    if len(rf) != len(aligned_sequence):
        raise ValueError("cmalign RF annotation and aligned sequence lengths differ")

    seed_path = ROOT / "data" / "raw" / spec["file"]
    seed_best, seed_top = best_seed_alignment(seed_path, chain_sequence)
    seed_rf = parse_rf_annotation(seed_path)
    if len(seed_rf) != AlignIO.read(str(seed_path), "stockholm").get_alignment_length():
        raise ValueError("Original seed RF annotation length mismatch")
    seed_rf_columns = [i + 1 for i, symbol in enumerate(seed_rf) if symbol != "."]

    sequence_index = 0
    model_position = 0
    mapping_rows = []
    non_model_rows = []
    for column, (rf_symbol, seq_symbol) in enumerate(zip(rf, aligned_sequence), start=1):
        has_residue = seq_symbol not in "-.~"
        residue = None
        if has_residue:
            sequence_index += 1
            if sequence_index > len(region_rows):
                raise ValueError("cmalign output contains more residues than the Rfam-mapped region")
            residue = region_rows[sequence_index - 1]

        if rf_symbol in {".", "~"}:
            if residue:
                non_model_rows.append({
                    "family_accession": accession,
                    "family": spec["family"],
                    "cm_alignment_column": column,
                    "rf_symbol": rf_symbol,
                    "label_asym_id": spec["label_asym_id"],
                    "auth_asym_id": spec["auth_asym_id"],
                    "label_seq_id": residue["label_seq_id"],
                    "auth_seq_id": residue["auth_seq_id"],
                    "pdbx_PDB_ins_code": residue["pdbx_PDB_ins_code"],
                    "label_comp_id": residue["label_comp_id"],
                    "normalized_base": residue["normalized_base"],
                    "has_model1_atoms": residue["has_model1_atoms"],
                    "reason": "CM_INSERTION" if rf_symbol == "." else "LOCAL_END_OUTSIDE_MODEL",
                })
            continue

        model_position += 1
        seed_alignment_column = (
            seed_rf_columns[model_position - 1]
            if model_position <= len(seed_rf_columns) else None
        )
        seed_target_position = (
            seed_best["mapping"].get(seed_alignment_column)
            if seed_alignment_column is not None else None
        )
        seed_residue = rows[seed_target_position - 1] if seed_target_position else None

        if residue:
            try:
                in_rfam_region = (
                    spec["rfam_auth_start"] <= int(residue["auth_seq_id"])
                    <= spec["rfam_auth_end"]
                )
            except (TypeError, ValueError):
                in_rfam_region = False
            status = (
                "MAPPED_MODELLED" if residue["has_model1_atoms"]
                else "MAPPED_POLYMER_UNRESOLVED"
            )
        else:
            in_rfam_region = False
            status = "MODEL_DELETION_IN_9ZTL_CHAIN"

        mapping_rows.append({
            "family_accession": accession,
            "family": spec["family"],
            "rf_model_position": model_position,
            "cm_alignment_column": column,
            "cm_rf_symbol": rf_symbol,
            "chain_sequence_index": residue["label_seq_id"] if residue else "",
            "chain_base": seq_symbol if residue else "",
            "label_asym_id": spec["label_asym_id"],
            "auth_asym_id": spec["auth_asym_id"],
            "label_seq_id": residue["label_seq_id"] if residue else "",
            "auth_seq_id": residue["auth_seq_id"] if residue else "",
            "pdbx_PDB_ins_code": residue["pdbx_PDB_ins_code"] if residue else "",
            "label_comp_id": residue["label_comp_id"] if residue else "",
            "normalized_base": residue["normalized_base"] if residue else "",
            "has_model1_atoms": residue["has_model1_atoms"] if residue else False,
            "is_modified_or_noncanonical_monomer": (
                residue["is_modified_or_noncanonical_monomer"] if residue else False
            ),
            "within_rfam_mapped_region": in_rfam_region,
            "rfam_mapped_auth_region": (
                f"{spec['rfam_auth_start']}-{spec['rfam_auth_end']}"
            ),
            "best_seed_sequence_id": seed_best["sequence_id"],
            "best_seed_identity": seed_best["identity"],
            "best_seed_query_coverage": seed_best["coverage"],
            "best_seed_label_seq_id": seed_target_position or "",
            "best_seed_auth_seq_id": seed_residue["auth_seq_id"] if seed_residue else "",
            "cm_vs_best_seed_agreement": (
                bool(residue and seed_target_position == residue["label_seq_id"])
            ),
            "mapping_status": status,
        })

    if sequence_index != len(region_rows):
        raise ValueError(
            f"cmalign consumed {sequence_index} of {len(region_rows)} Rfam-region residues "
            f"for {accession}"
        )
    mapped_region_ids = {row["label_seq_id"] for row in region_rows}
    for row in rows:
        if row["label_seq_id"] not in mapped_region_ids:
            non_model_rows.append({
                "family_accession": accession,
                "family": spec["family"],
                "cm_alignment_column": "",
                "rf_symbol": "",
                "label_asym_id": spec["label_asym_id"],
                "auth_asym_id": spec["auth_asym_id"],
                "label_seq_id": row["label_seq_id"],
                "auth_seq_id": row["auth_seq_id"],
                "pdbx_PDB_ins_code": row["pdbx_PDB_ins_code"],
                "label_comp_id": row["label_comp_id"],
                "normalized_base": row["normalized_base"],
                "has_model1_atoms": row["has_model1_atoms"],
                "reason": "OUTSIDE_Rfam_STRUCTURE_REGION",
            })
    if model_position != len(seed_rf_columns):
        raise ValueError(
            f"CM model length {model_position} differs from original seed RF length "
            f"{len(seed_rf_columns)} for {accession}"
        )

    mapped_auth = [
        int(row["auth_seq_id"]) for row in mapping_rows
        if row["auth_seq_id"] not in {"", ".", "?"}
    ]
    model_length = model_position
    summary = {
        "accession": accession,
        "family": spec["family"],
        "label_asym_id": spec["label_asym_id"],
        "auth_asym_id": spec["auth_asym_id"],
        "entity_id": spec["entity_id"],
        "polymer_length": len(rows),
        "rfam_region_polymer_length": len(region_rows),
        "model1_atom_residues": sum(row["has_model1_atoms"] for row in rows),
        "cm_model_length": model_length,
        "mapped_model_positions": sum(bool(row["label_seq_id"]) for row in mapping_rows),
        "model_deletions": sum(row["mapping_status"] == "MODEL_DELETION_IN_9ZTL_CHAIN" for row in mapping_rows),
        "non_model_chain_residues": len(non_model_rows),
        "mapped_auth_min": min(mapped_auth) if mapped_auth else None,
        "mapped_auth_max": max(mapped_auth) if mapped_auth else None,
        "rfam_expected_auth_region": f"{spec['rfam_auth_start']}-{spec['rfam_auth_end']}",
        "within_rfam_region": sum(row["within_rfam_mapped_region"] for row in mapping_rows),
        "seed_best_hit": seed_best["sequence_id"],
        "seed_best_identity": seed_best["identity"],
        "seed_best_coverage": seed_best["coverage"],
        "seed_top10_best_identity": seed_top[0]["identity"],
    }
    if summary["mapped_model_positions"] < 0.70 * model_length:
        raise ValueError(f"Less than 70% of {accession} CM model positions mapped to the chain")
    if summary["within_rfam_region"] != summary["mapped_model_positions"]:
        raise ValueError(
            f"Mapped {accession} positions fall outside the Rfam structure-region interval"
        )
    return mapping_rows, non_model_rows, summary


def main() -> None:
    if not shutil_which("cmalign"):
        raise RuntimeError("Infernal cmalign executable not found on PATH")
    infernal = infernal_version()
    cif_path, cif_source = get_cif()
    cif = load_cif(cif_path)
    all_mapping, all_non_model, summaries, cm_records = [], [], [], []
    for accession, spec in FAMILIES.items():
        cm_path, cm_url = download_cm(accession)
        mapping, non_model, summary = align_chain_to_cm(
            accession, spec, cif, cm_path, cm_url, infernal
        )
        all_mapping.extend(mapping)
        all_non_model.extend(non_model)
        summaries.append(summary)
        cm_records.append({
            "accession": accession,
            "source": cm_url,
            "file": str(cm_path.relative_to(ROOT)),
            "bytes": cm_path.stat().st_size,
            "sha256": sha256(cm_path),
        })

    tables = ROOT / "supplementary" / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_mapping).to_csv(
        tables / "rfam_to_pdb_9ztl_cmalign.csv", index=False
    )
    pd.DataFrame(all_non_model).to_csv(
        tables / "9ztl_chain_insertions_and_local_ends.csv", index=False
    )
    manifest = {
        "pdb_id": "9ZTL",
        "rfam_release": RFAM_RELEASE,
        "rfam_structure_map_source": "https://rfam.org/structure/9ztl/av",
        "cif_source": cif_source,
        "cif_file": str(cif_path.relative_to(ROOT)),
        "cif_bytes": cif_path.stat().st_size,
        "cif_sha256": sha256(cif_path),
        "infernal_version": infernal,
        "cm_files": cm_records,
        "summaries": summaries,
        "mapping_method": (
            "Infernal cmalign of the deposited polymer subsequence bounded by the "
            "Rfam structure-region author-numbering interval to the official Rfam 15.1 "
            "covariance model, using global model alignment (-g). Match-state RF columns "
            "define RF model positions; RF-dot columns are insertions relative to the model. "
            "The best seed-sequence pairwise alignment is retained only as a cross-check."
        ),
        "limitations": [
            "The 9ZTL P-site tRNA is one biological tRNA species and does not represent all RF00005 members.",
            "Polymer-scheme residues without model-1 atoms remain flagged as unresolved.",
            "A CM alignment is a model-based homology mapping, not experimental validation of a base pair.",
        ],
    }
    manifest_path = ROOT / "supplementary" / "provenance" / "9ztl_cmalign_mapping_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("CM-based RF-to-PDB mapping completed.")
    for summary in summaries:
        print(json.dumps(summary))


def shutil_which(command: str) -> str | None:
    import shutil
    return shutil.which(command)


if __name__ == "__main__":
    main()

"""Compare R-scape pair candidates with FR3D geometric interactions in 9ZTL."""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAIR_TABLE = ROOT / "supplementary" / "tables" / "rscape_coordinate_mapping.csv"
BASEPAIR_FILE = ROOT / "supplementary" / "structural" / "9ZTL_basepair.txt"
STACKING_FILE = ROOT / "supplementary" / "structural" / "9ZTL_stacking.txt"
OUTPUT = ROOT / "supplementary" / "tables" / "rscape_fr3d_structural_evidence.csv"
SUMMARY = ROOT / "supplementary" / "tables" / "rscape_fr3d_evidence_summary.csv"


def parse_unit(unit_id: str) -> tuple[str, str, int]:
    fields = unit_id.split("|")
    if len(fields) < 5:
        raise ValueError(f"Unexpected FR3D unit ID: {unit_id}")
    # BGSU Unit IDs use auth_asym_id and auth_seq_id (mmCIF author numbering).
    return fields[2], fields[3], int(fields[4])


def read_interactions(path: Path, stacking: bool = False):
    interactions = defaultdict(set)
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 3:
            raise ValueError(f"Malformed FR3D line {line_number} in {path}")
        unit1, code, unit2 = fields[:3]
        chain1, comp1, seq1 = parse_unit(unit1)
        chain2, comp2, seq2 = parse_unit(unit2)
        if chain1 != chain2 or chain1 not in {"b", "5"}:
            continue
        if stacking and (not code.startswith("s") or code.startswith("ns")):
            continue
        key = (chain1, min(seq1, seq2), max(seq1, seq2))
        interactions[key].add(code)
    return interactions


def main() -> None:
    basepairs = read_interactions(BASEPAIR_FILE, stacking=False)
    stackings = read_interactions(STACKING_FILE, stacking=True)

    with PAIR_TABLE.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if len(rows) != 82:
        raise ValueError(f"Expected 82 R-scape pairs, found {len(rows)}")
    for row in rows:
        if row["structure_chain"] not in {"b", "5"}:
            raise ValueError(f"Unexpected structure author chain: {row['structure_chain']}")
        if row.get("structure_left_residue", "") != row.get("structure_left_auth_seq_id", ""):
            raise ValueError("Left author residue ID is inconsistent across mapping columns")
        if row.get("structure_right_residue", "") != row.get("structure_right_auth_seq_id", ""):
            raise ValueError("Right author residue ID is inconsistent across mapping columns")

    additions = [
        "fr3d_basepair_observed", "fr3d_basepair_codes",
        "fr3d_stacking_observed", "fr3d_stacking_codes",
        "fr3d_has_cWW_geometry", "fr3d_structural_evidence",
    ]
    for name in additions:
        if name not in fields:
            fields.append(name)

    summary_counts = Counter()
    summary_rows = []
    for row in rows:
        chain = row["structure_chain"]
        left = row.get("structure_left_residue", "")
        right = row.get("structure_right_residue", "")
        if not left or not right:
            bp_codes, st_codes = set(), set()
            evidence = "NOT_ASSESSABLE_RF_MODEL_POSITION_ABSENT_FROM_9ZTL"
        else:
            left_num, right_num = int(left), int(right)
            key = (chain, min(left_num, right_num), max(left_num, right_num))
            bp_codes = basepairs.get(key, set())
            st_codes = stackings.get(key, set())
            if bp_codes and st_codes:
                evidence = "DIRECT_BASEPAIR_AND_STACKING_ANNOTATED"
            elif bp_codes:
                evidence = "DIRECT_BASEPAIR_ANNOTATED"
            elif st_codes:
                evidence = "STACKING_ONLY_NO_DIRECT_BASEPAIR_ANNOTATION"
            else:
                evidence = "NO_DIRECT_BASEPAIR_OR_STACKING_ANNOTATION"

        row["fr3d_basepair_observed"] = bool(bp_codes)
        row["fr3d_basepair_codes"] = ";".join(sorted(bp_codes))
        row["fr3d_stacking_observed"] = bool(st_codes)
        row["fr3d_stacking_codes"] = ";".join(sorted(st_codes))
        row["fr3d_has_cWW_geometry"] = "cWW" in bp_codes
        row["fr3d_structural_evidence"] = evidence

        key = (row["family"], row["pair_class"], evidence)
        summary_counts[key] += 1

    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    with SUMMARY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["family", "rscape_pair_class", "fr3d_structural_evidence", "n_pairs"],
        )
        writer.writeheader()
        for (family, pair_class, evidence), count in sorted(summary_counts.items()):
            writer.writerow({
                "family": family,
                "rscape_pair_class": pair_class,
                "fr3d_structural_evidence": evidence,
                "n_pairs": count,
            })

    print(f"Compared {len(rows)} R-scape pairs with FR3D annotations.")
    print(f"FR3D direct base-pair keys: {len(basepairs)}; stacking keys: {len(stackings)}")
    for key, count in sorted(summary_counts.items()):
        print("\t".join([*key, str(count)]))


if __name__ == "__main__":
    main()

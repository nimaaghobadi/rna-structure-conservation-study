"""Build the descriptive contrast-of-contrasts table from the committed core summaries."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ("5S rRNA", "tRNA")
STATES = ("paired", "unpaired")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def build_rows(
    qc_rows: list[dict[str, str]],
    conservation_rows: list[dict[str, str]],
    entropy_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    sequence_counts = {row["family"]: int(row["n_sequences"]) for row in qc_rows}
    if set(sequence_counts) != set(FAMILIES):
        raise ValueError(f"Expected exactly these families: {FAMILIES}")

    def means(rows: list[dict[str, str]]) -> dict[tuple[str, str], float]:
        result: dict[tuple[str, str], float] = {}
        for row in rows:
            key = (row["family"], row["structure_state"])
            if key in result:
                raise ValueError(f"Duplicate summary row: {key}")
            result[key] = float(row["mean"])
        expected = {(family, state) for family in FAMILIES for state in STATES}
        if set(result) != expected:
            raise ValueError(f"Summary rows do not match expected family/state pairs: {set(result)}")
        return result

    conservation = means(conservation_rows)
    entropy = means(entropy_rows)
    counts: dict[tuple[str, str], int] = {}
    for row in conservation_rows:
        key = (row["family"], row["structure_state"])
        counts[key] = int(row["count"])

    output = []
    contrasts = {}
    for family in FAMILIES:
        conservation_delta = (
            conservation[(family, "paired")] - conservation[(family, "unpaired")]
        )
        entropy_delta = entropy[(family, "paired")] - entropy[(family, "unpaired")]
        contrasts[family] = (conservation_delta, entropy_delta)
        output.append({
            "family": family,
            "n_sequences": str(sequence_counts[family]),
            "n_paired_columns": str(counts[(family, "paired")]),
            "n_unpaired_columns": str(counts[(family, "unpaired")]),
            "paired_minus_unpaired_mean_conservation": repr(conservation_delta),
            "paired_minus_unpaired_mean_entropy_bits": repr(entropy_delta),
        })

    output.append({
        "family": "tRNA_minus_5S",
        "n_sequences": "NA",
        "n_paired_columns": "NA",
        "n_unpaired_columns": "NA",
        "paired_minus_unpaired_mean_conservation": repr(
            contrasts["tRNA"][0] - contrasts["5S rRNA"][0]
        ),
        "paired_minus_unpaired_mean_entropy_bits": repr(
            contrasts["tRNA"][1] - contrasts["5S rRNA"][1]
        ),
    })
    return output


def main() -> None:
    rows = build_rows(
        read_csv(ROOT / "data" / "tables" / "qc_summary.csv"),
        read_csv(ROOT / "data" / "tables" / "conservation_by_structure.csv"),
        read_csv(ROOT / "data" / "tables" / "entropy_by_structure.csv"),
    )
    output = ROOT / "supplementary" / "tables" / "family_contrast_summary.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "family",
        "n_sequences",
        "n_paired_columns",
        "n_unpaired_columns",
        "paired_minus_unpaired_mean_conservation",
        "paired_minus_unpaired_mean_entropy_bits",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

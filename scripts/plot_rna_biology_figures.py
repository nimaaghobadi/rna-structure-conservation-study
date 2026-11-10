"""Regenerate the project base-pair composition SVG from the archived CSV.

Run from the repository root:
    python analysis/plot_rna_biology_figures.py

Bars show the composition among valid sequence-by-pair observations. Missing or
ambiguous observations are excluded from each bar and reported as coverage.
The script uses only the Python standard library and writes deterministic SVG.
"""
from __future__ import annotations

import csv
import math
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "tables" / "base_pair_summary.csv"
OUTPUT = ROOT / "figures" / "base_pair_composition.svg"
ORDER = ("5S rRNA", "tRNA")
PAIR_TYPES = ("AU", "UA", "GC", "CG", "GU", "UG")
CATEGORIES = (
    ("count_AU", "AU", "#4477AA"),
    ("count_UA", "UA", "#66AADD"),
    ("count_GC", "GC", "#228833"),
    ("count_CG", "CG", "#77AA55"),
    ("count_GU", "GU", "#CCBB44"),
    ("count_UG", "UG", "#EECC66"),
    ("count_noncanonical_or_other", "Other / noncanonical", "#AA4499"),
)


def load_and_validate(path: Path) -> dict[str, dict[str, int | float]]:
    """Read the two expected family rows and verify every plotted denominator."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "family",
            "n_valid_pair_observations",
            "n_possible_pair_observations",
            "valid_pair_observation_fraction",
            "canonical_or_wobble_fraction",
        }
        required.update(column for column, _, _ in CATEGORIES)
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Input table is missing required columns: {sorted(missing)}")
        rows = list(reader)

    if len(rows) != len(ORDER):
        raise ValueError(f"Expected exactly {len(ORDER)} family rows; found {len(rows)}")
    if len({row["family"] for row in rows}) != len(rows):
        raise ValueError("Duplicate family names in input table")
    if {row["family"] for row in rows} != set(ORDER):
        raise ValueError(f"Expected family rows {ORDER}; found {[r['family'] for r in rows]}")

    result: dict[str, dict[str, int | float]] = {}
    for raw in rows:
        family = raw["family"]
        valid = int(raw["n_valid_pair_observations"])
        possible = int(raw["n_possible_pair_observations"])
        counts = {column: int(raw[column]) for column, _, _ in CATEGORIES}
        if valid <= 0 or possible < valid:
            raise ValueError(f"Invalid observation denominators for {family}: {valid}/{possible}")
        if any(value < 0 for value in counts.values()):
            raise ValueError(f"Negative pair count for {family}")
        if sum(counts.values()) != valid:
            raise ValueError(f"Pair counts do not sum to valid observations for {family}")

        fractions = {column: counts[column] / valid for column, _, _ in CATEGORIES}
        if not math.isclose(sum(fractions.values()), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"Plotted categories do not sum to 1 for {family}")
        recorded_coverage = float(raw["valid_pair_observation_fraction"])
        if not math.isclose(recorded_coverage, valid / possible, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"Recorded coverage disagrees with denominators for {family}")
        recorded_canonical = float(raw["canonical_or_wobble_fraction"])
        calculated_canonical = sum(counts[f"count_{pair}"] for pair in PAIR_TYPES) / valid
        if not math.isclose(recorded_canonical, calculated_canonical, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"Canonical/wobble fraction disagrees with pair counts for {family}")
        result[family] = {
            "valid": valid,
            "possible": possible,
            "coverage": valid / possible,
            **counts,
            **{
                f"fraction_{column.removeprefix('count_')}": value
                for column, value in fractions.items()
            },
            "other_fraction": fractions["count_noncanonical_or_other"],
        }
    return result


def render_svg(data: dict[str, dict[str, int | float]]) -> str:
    """Render a fixed, accessible SVG whose geometry is driven by the CSV."""
    x0, bar_width, bar_height = 270.0, 720.0, 66.0
    y_positions = {"5S rRNA": 145.0, "tRNA": 295.0}
    out = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="500" viewBox="0 0 1200 500" role="img" aria-labelledby="title desc">',
        '<title id="title">Observed base-pair composition in 5S rRNA and tRNA</title>',
        '<desc id="desc">Two horizontal stacked bars compare AU, UA, GC, CG, GU, UG, and other or noncanonical base-pair observations among valid pairs in the two Rfam seed alignments.</desc>',
        '<rect width="1200" height="500" fill="#ffffff"/>',
        '<text x="62" y="56" font-family="Arial, sans-serif" font-size="27" font-weight="bold" fill="#172b3a">Observed base-pair composition in two RNA families</text>',
        '<text x="62" y="88" font-family="Arial, sans-serif" font-size="16" fill="#536675">Share of valid paired-base observations · canonical and wobble pairs shown separately</text>',
        '<text x="62" y="178" font-family="Arial, sans-serif" font-size="23" font-weight="bold" fill="#172b3a">5S rRNA</text>',
        '<text x="62" y="328" font-family="Arial, sans-serif" font-size="23" font-weight="bold" fill="#172b3a">tRNA</text>',
        '<g stroke="#dce3e8" stroke-width="1">',
    ]
    for tick in (0, 25, 50, 75, 100):
        x = x0 + bar_width * tick / 100
        out.append(f'<line x1="{x:.2f}" y1="135" x2="{x:.2f}" y2="385"/>')
    out.append("</g>")

    for family in ORDER:
        row = data[family]
        y = y_positions[family]
        out.append("<g>")
        left = x0
        labels = []
        for column, label, color in CATEGORIES:
            count = int(row[column])
            value = 100.0 * count / int(row["valid"])
            width = bar_width * value / 100.0
            out.append(
                f'<rect x="{left:.2f}" y="{y:.0f}" width="{width:.2f}" '
                f'height="{bar_height:.0f}" fill="{color}"/>'
            )
            if value >= 6.0:
                labels.append(
                    f'<text x="{left + width / 2:.1f}" y="{y + 39:.0f}">{value:.1f}</text>'
                )
            left += width
        out.append(
            '<g font-family="Arial, sans-serif" font-size="14" font-weight="bold" '
            'text-anchor="middle" fill="#10212d">'
        )
        out.extend(labels)
        out.append("</g></g>")

    out.append(
        '<g font-family="Arial, sans-serif" font-size="13" fill="#536675" text-anchor="middle">'
    )
    for tick, label in zip((0, 25, 50, 75, 100), ("0", "25", "50", "75", "100%")):
        x = x0 + bar_width * tick / 100
        out.append(f'<text x="{x:.0f}" y="410">{label}</text>')
    out.append("</g>")

    out.append('<g font-family="Arial, sans-serif" font-size="13" fill="#263b49">')
    legend_x = (62, 142, 222, 302, 382, 462, 542)
    for (_, label, color), x in zip(CATEGORIES, legend_x):
        out.append(f'<rect x="{x}" y="440" width="13" height="13" fill="{color}"/>')
        out.append(f'<text x="{x + 20}" y="452">{escape(label)}</text>')
    out.append("</g>")

    five = data["5S rRNA"]
    trna = data["tRNA"]
    footer = (
        f"Valid-pair coverage: 5S rRNA {int(five['valid']):,}/{int(five['possible']):,} "
        f"({100 * float(five['coverage']):.1f}%); tRNA {int(trna['valid']):,}/"
        f"{int(trna['possible']):,} ({100 * float(trna['coverage']):.1f}%). "
        "Percentages within each bar sum to 100%."
    )
    out.append(
        f'<text x="62" y="482" font-family="Arial, sans-serif" font-size="12" '
        f'fill="#536675">{escape(footer)}</text>'
    )
    out.append("</svg>")
    return "\n".join(out) + "\n"


def main() -> None:
    data = load_and_validate(SOURCE)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(render_svg(data), encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT)} from {SOURCE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

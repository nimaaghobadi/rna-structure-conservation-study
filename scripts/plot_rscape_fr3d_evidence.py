"""Generate an original RNA-focused SVG summary from the R-scape/FR3D evidence table."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "supplementary" / "tables" / "rscape_fr3d_evidence_summary.csv"
OUTPUT = ROOT / "figures" / "rscape_fr3d_evidence.svg"

CATEGORIES = [
    ("DIRECT_BASEPAIR_ANNOTATED", "Direct base-pair annotation", "#2A7F62"),
    ("STACKING_ONLY_NO_DIRECT_BASEPAIR_ANNOTATION", "Stacking only", "#D69E36"),
    ("NO_DIRECT_BASEPAIR_OR_STACKING_ANNOTATION", "No direct pair annotation", "#89939D"),
    ("NOT_ASSESSABLE_RF_MODEL_POSITION_ABSENT_FROM_9ZTL", "Not assessable", "#B65353"),
]
GROUPS = [
    ("RF00001", "ANNOTATED_PAIR_SIGNIFICANT", "5S rRNA · annotated, R-scape significant"),
    ("RF00001", "ANNOTATED_PAIR_NOT_DETECTED", "5S rRNA · annotated, R-scape not detected"),
    ("RF00001", "ADDITIONAL_COVARIATION_CANDIDATE_OUTSIDE_ORIGINAL_SS_cons", "5S rRNA · additional covariation candidates"),
    ("RF00005", "ANNOTATED_PAIR_SIGNIFICANT", "tRNA · annotated, R-scape significant"),
    ("RF00005", "ADDITIONAL_COVARIATION_CANDIDATE_OUTSIDE_ORIGINAL_SS_cons", "tRNA · additional covariation candidates"),
]


def main() -> None:
    counts = defaultdict(int)
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            counts[(row["family"], row["rscape_pair_class"], row["fr3d_structural_evidence"])] += int(row["n_pairs"])

    groups = []
    for family, pair_class, label in GROUPS:
        values = [counts[(family, pair_class, category)] for category, _, _ in CATEGORIES]
        groups.append((label, values, sum(values)))
    if len(groups) != 5:
        raise ValueError(f"Expected five RNA-family evidence groups, found {len(groups)}")
    max_total = max(total for _, _, total in groups)
    if max_total != 22:
        raise ValueError(f"Unexpected maximum group size: {max_total}; inspect evidence summary")

    width, height = 1280, 660
    bar_x, bar_max_width = 390, 720
    scale = bar_max_width / max_total
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">R-scape RNA pair candidates compared with FR3D interactions in 9ZTL</title>',
        '<desc id="desc">Five horizontal stacked bars show direct base-pair annotations, stacking-only interactions, no direct pair annotation, and pairs not assessable because a model position is absent from the deposited chain.</desc>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="64" y="52" font-family="Arial,sans-serif" font-size="27" font-weight="bold" fill="#172b3a">R-scape pair candidates mapped to 9ZTL structural interactions</text>',
        '<text x="64" y="82" font-family="Arial,sans-serif" font-size="16" fill="#536675">FR3D geometric annotations · one deposited E. coli ribosome model · RF-model coordinate system</text>',
    ]
    legend_x = 64
    for _, label, color in CATEGORIES:
        out.append(f'<rect x="{legend_x}" y="112" width="14" height="14" rx="2" fill="{color}"/>')
        out.append(f'<text x="{legend_x+21}" y="124" font-family="Arial,sans-serif" font-size="13" fill="#334653">{label}</text>')
        legend_x += [235, 135, 225, 140][CATEGORIES.index((_, label, color))]
    y0, gap, bar_h = 185, 84, 34
    for i, (label, values, total) in enumerate(groups):
        y = y0 + i * gap
        out.append(f'<text x="64" y="{y+22}" font-family="Arial,sans-serif" font-size="16" fill="#203543">{label}</text>')
        x = bar_x
        for value, (_, _, color) in zip(values, CATEGORIES):
            w = value * scale
            if value:
                out.append(f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{bar_h}" fill="{color}"/>')
                if w >= 24:
                    text_color = "#ffffff" if color in {"#2A7F62", "#89939D", "#B65353"} else "#172b3a"
                    out.append(f'<text x="{x+w/2:.1f}" y="{y+22}" text-anchor="middle" font-family="Arial,sans-serif" font-size="14" font-weight="bold" fill="{text_color}">{value}</text>')
                x += w
        out.append(f'<text x="{x+12:.1f}" y="{y+22}" font-family="Arial,sans-serif" font-size="15" font-weight="bold" fill="#203543">n={total}</text>')
        out.append(f'<line x1="{bar_x}" y1="{y+bar_h+13}" x2="{bar_x+bar_max_width}" y2="{y+bar_h+13}" stroke="#e7ecef" stroke-width="1"/>')
    out.extend([
        '<text x="64" y="612" font-family="Arial,sans-serif" font-size="13" fill="#536675">“No direct pair annotation” is not proof of absence. Stacking-only is not base-pair support.</text>',
        '<text x="64" y="634" font-family="Arial,sans-serif" font-size="13" fill="#536675">The tRNA result applies to the single P-site tRNA chain in 9ZTL, not every RF00005 family member.</text>',
        '</svg>',
    ])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

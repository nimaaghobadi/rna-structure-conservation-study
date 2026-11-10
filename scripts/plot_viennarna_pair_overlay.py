"""Create original RF-coordinate arc diagrams comparing RNAalifold MFE with Rfam SS_cons."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "supplementary" / "tables" / "rscape_fr3d_structural_evidence.csv"
VIENNA = ROOT / "supplementary" / "tables" / "viennarna_pair_probabilities.csv"
OUTPUT = ROOT / "figures" / "viennarna_consensus_comparison.svg"

FAMILIES = [
    ("RF00001", "5S rRNA", 120, 430, 165),
    ("RF00005", "tRNA", 71, 820, 555),
]
COLORS = {
    "recovered": "#2A7F62",
    "not_recovered": "#89939D",
    "extra": "#2878A5",
}


def read_rows(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    evidence = read_rows(EVIDENCE)
    vienna = read_rows(VIENNA)
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="1000" viewBox="0 0 1400 1000" role="img" aria-labelledby="title desc">',
        '<title id="title">RNAalifold MFE pair recovery compared with Rfam consensus structures</title>',
        '<desc id="desc">RF-coordinate arc diagrams show Rfam consensus pairs recovered or not recovered by RNAalifold and additional MFE pairs outside the consensus annotation.</desc>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="70" y="54" font-family="Arial,sans-serif" font-size="28" font-weight="bold" fill="#172b3a">RNAalifold consensus-folding agreement with Rfam SS_cons</text>',
        '<text x="70" y="84" font-family="Arial,sans-serif" font-size="16" fill="#536675">Arc endpoints are RF model positions; all structures are unconstrained predictions or reference annotations</text>',
    ]
    legend = [
        ("#2A7F62", "SS_cons pair recovered in MFE"),
        ("#89939D", "SS_cons pair not recovered in MFE"),
        ("#2878A5", "Additional RNAalifold MFE pair"),
    ]
    lx = 70
    for color, label in legend:
        svg.append(f'<line x1="{lx}" y1="120" x2="{lx+28}" y2="120" stroke="{color}" stroke-width="4"/>')
        svg.append(f'<text x="{lx+38}" y="125" font-family="Arial,sans-serif" font-size="13" fill="#334653">{label}</text>')
        lx += [300, 330, 300][legend.index((color, label))]

    for accession, family, length, baseline, title_y in FAMILIES:
        ref_rows = [
            row for row in evidence
            if row["family"] == accession
            and row["pair_class"] in {
                "ANNOTATED_PAIR_SIGNIFICANT",
                "ANNOTATED_PAIR_NOT_DETECTED",
            }
        ]
        ref_pairs = {
            tuple(sorted((int(row["rf_left_model_pos"]), int(row["rf_right_model_pos"]))))
            for row in ref_rows
        }
        mfe_rows = [
            row for row in vienna
            if row["family_accession"] == accession
            and row["subset"] == "full"
            and row["in_mfe_structure"].strip().lower() == "true"
        ]
        mfe_pairs = {
            tuple(sorted((int(row["rf_model_position_left"]), int(row["rf_model_position_right"]))))
            for row in mfe_rows
            if row["rf_model_position_left"] and row["rf_model_position_right"]
        }
        extra = mfe_pairs - ref_pairs
        expected_extra = 4 if accession == "RF00001" else 0
        if len(extra) != expected_extra:
            raise ValueError(
                f"Expected {expected_extra} extra MFE pairs for {accession}, found {len(extra)}"
            )
        recovered = mfe_pairs & ref_pairs
        not_recovered = ref_pairs - mfe_pairs
        x0, plot_width = 90, 850
        scale = plot_width / (length - 1)
        svg.append(f'<text x="70" y="{title_y}" font-family="Arial,sans-serif" font-size="22" font-weight="bold" fill="#172b3a">{family} · {len(recovered)}/{len(ref_pairs)} SS_cons pairs recovered</text>')
        svg.append(f'<line x1="{x0}" y1="{baseline}" x2="{x0+plot_width}" y2="{baseline}" stroke="#263b49" stroke-width="2"/>')

        def draw_pair(pair, color, dash=None, width=2.0):
            left, right = pair
            x1 = x0 + (left - 1) * scale
            x2 = x0 + (right - 1) * scale
            span = right - left
            height = 24 + (span / max(1, length - 1)) * 130
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            svg.append(
                f'<path d="M {x1:.1f} {baseline} Q {(x1+x2)/2:.1f} {baseline-height:.1f} {x2:.1f} {baseline}" '
                f'fill="none" stroke="{color}" stroke-width="{width}"{dash_attr}/>'
            )

        for pair in sorted(not_recovered, key=lambda p: p[1]-p[0], reverse=True):
            draw_pair(pair, COLORS["not_recovered"], dash="5 4", width=1.5)
        for pair in sorted(recovered, key=lambda p: p[1]-p[0], reverse=True):
            draw_pair(pair, COLORS["recovered"], width=2.2)
        for pair in sorted(extra, key=lambda p: p[1]-p[0], reverse=True):
            draw_pair(pair, COLORS["extra"], width=3.2)

        tick_step = 20 if length > 100 else 10
        ticks = list(range(1, length + 1, tick_step))
        if ticks[-1] != length:
            ticks.append(length)
        for tick in ticks:
            x = x0 + (tick - 1) * scale
            svg.append(f'<line x1="{x:.1f}" y1="{baseline}" x2="{x:.1f}" y2="{baseline+7}" stroke="#263b49" stroke-width="1"/>')
            svg.append(f'<text x="{x:.1f}" y="{baseline+25}" text-anchor="middle" font-family="Arial,sans-serif" font-size="11" fill="#536675">{tick}</text>')
        svg.append(f'<text x="{x0+plot_width/2:.1f}" y="{baseline+46}" text-anchor="middle" font-family="Arial,sans-serif" font-size="12" fill="#536675">RF model position</text>')

        note_x = 990
        if accession == "RF00001":
            svg.append(f'<text x="{note_x}" y="{baseline-175}" font-family="Arial,sans-serif" font-size="16" font-weight="bold" fill="#172b3a">Additional MFE pairs</text>')
            for j, row in enumerate(sorted(
                [r for r in mfe_rows if tuple(sorted((int(r["rf_model_position_left"]), int(r["rf_model_position_right"])))) in extra],
                key=lambda r: int(r["rf_model_position_left"])
            )):
                y = baseline - 145 + j * 31
                p = float(row["base_pair_probability"])
                svg.append(f'<text x="{note_x}" y="{y}" font-family="Arial,sans-serif" font-size="14" fill="#2878A5">RF {row["rf_model_position_left"]}–{row["rf_model_position_right"]} · p={p:.2f}</text>')
            svg.append(f'<text x="{note_x}" y="{baseline-12}" font-family="Arial,sans-serif" font-size="12" fill="#536675">All four are R-scape candidates</text>')
            svg.append(f'<text x="{note_x}" y="{baseline+7}" font-family="Arial,sans-serif" font-size="12" fill="#536675">with direct FR3D base-pair</text>')
            svg.append(f'<text x="{note_x}" y="{baseline+26}" font-family="Arial,sans-serif" font-size="12" fill="#536675">annotations in 9ZTL.</text>')
        else:
            svg.append(f'<text x="{note_x}" y="{baseline-150}" font-family="Arial,sans-serif" font-size="16" font-weight="bold" fill="#172b3a">Prediction summary</text>')
            svg.append(f'<text x="{note_x}" y="{baseline-118}" font-family="Arial,sans-serif" font-size="14" fill="#334653">20 of 21 SS_cons pairs</text>')
            svg.append(f'<text x="{note_x}" y="{baseline-96}" font-family="Arial,sans-serif" font-size="14" fill="#334653">recovered in the MFE.</text>')
            svg.append(f'<text x="{note_x}" y="{baseline-55}" font-family="Arial,sans-serif" font-size="13" fill="#89939D">Omitted: RF 13–21</text>')
            svg.append(f'<text x="{note_x}" y="{baseline-34}" font-family="Arial,sans-serif" font-size="12" fill="#536675">No direct FR3D pair</text>')
            svg.append(f'<text x="{note_x}" y="{baseline-15}" font-family="Arial,sans-serif" font-size="12" fill="#536675">annotation in this P-site tRNA.</text>')

    svg.extend([
        '<text x="70" y="960" font-family="Arial,sans-serif" font-size="13" fill="#536675">The arc diagram is an original alignment-derived comparison, not an atomic structure. RNAalifold is unconstrained and predicts pseudoknot-free secondary structures.</text>',
        '</svg>',
    ])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(svg) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

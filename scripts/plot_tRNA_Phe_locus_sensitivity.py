#!/usr/bin/env python3
"""Render an original SVG summarizing tRNA-Phe locus context and sampling sensitivity."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "supplementary/tables/trna_phe_locus_summary.json"
CONSERVATION = ROOT / "data/tables/conservation_by_structure.csv"
ENTROPY = ROOT / "data/tables/entropy_by_structure.csv"
BOVID = ROOT / "supplementary/tables/trna_phe_bovid_sensitivity.csv"
OUTPUT = ROOT / "supplementary/figures/trna_phe_locus_sensitivity.svg"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def broad_delta(rows: list[dict[str, str]], family: str) -> float:
    values = {row["structure_state"]: float(row["mean"])
              for row in rows if row["family"] == family}
    if set(values) != {"paired", "unpaired"}:
        raise ValueError(f"Expected paired and unpaired values for {family}.")
    return values["paired"] - values["unpaired"]


def build_points(summary: dict[str, Any], conservation_rows: list[dict[str, str]],
                 entropy_rows: list[dict[str, str]], bovid_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    thresholds = summary["valid_fraction_threshold_sensitivity"]
    return [
        {"label": "5S rRNA (n=712)", "conservation": broad_delta(conservation_rows, "5S rRNA"),
         "entropy": broad_delta(entropy_rows, "5S rRNA"), "kind": "broad"},
        {"label": "RF00005 tRNA (n=954)", "conservation": broad_delta(conservation_rows, "tRNA"),
         "entropy": broad_delta(entropy_rows, "tRNA"), "kind": "broad"},
        {"label": "tRNA-Phe, all 16 taxa",
         "conservation": summary["paired_minus_unpaired_conservation"],
         "entropy": summary["paired_minus_unpaired_entropy_bits"], "kind": "locus"},
        {"label": "tRNA-Phe, coverage ≥50%",
         "conservation": thresholds["0.5"]["paired_minus_unpaired_conservation"],
         "entropy": thresholds["0.5"]["paired_minus_unpaired_entropy_bits"], "kind": "sensitivity"},
        {"label": "tRNA-Phe, coverage ≥70%",
         "conservation": thresholds["0.7"]["paired_minus_unpaired_conservation"],
         "entropy": thresholds["0.7"]["paired_minus_unpaired_entropy_bits"], "kind": "sensitivity"},
        {"label": "tRNA-Phe, coverage ≥90%",
         "conservation": thresholds["0.9"]["paired_minus_unpaired_conservation"],
         "entropy": thresholds["0.9"]["paired_minus_unpaired_entropy_bits"], "kind": "sensitivity"},
        {"label": "One-Bovidae sensitivity (range)",
         "conservation_range": (
             min(float(row["paired_minus_unpaired_conservation"]) for row in bovid_rows),
             max(float(row["paired_minus_unpaired_conservation"]) for row in bovid_rows),
         ),
         "entropy_range": (
             min(float(row["paired_minus_unpaired_entropy_bits"]) for row in bovid_rows),
             max(float(row["paired_minus_unpaired_entropy_bits"]) for row in bovid_rows),
         ), "kind": "range"},
    ]


def xmap(value: float, xmin: float, xmax: float, x0: float, width: float) -> float:
    return x0 + (value - xmin) / (xmax - xmin) * width


def panel_svg(points: list[dict[str, Any]], metric: str, title: str,
              x0: float, y0: float, width: float, row_height: float,
              xmin: float, xmax: float, color: str, panel_letter: str) -> str:
    out = [
        f'<line x1="{x0}" y1="{y0+18}" x2="{x0+width}" y2="{y0+18}" class="axis"/>',
    ]
    chart_y = y0 + 52
    zero_x = xmap(0.0, xmin, xmax, x0, width)
    out.append(f'<line x1="{zero_x:.1f}" y1="{chart_y-18}" x2="{zero_x:.1f}" y2="{chart_y+row_height*len(points)-8}" class="zero"/>')
    for i, point in enumerate(points):
        y = chart_y + i * row_height
        out.append(f'<text x="{x0-16}" y="{y+5:.1f}" text-anchor="end" class="row-label">{point["label"]}</text>')
        if "conservation_range" in point and metric == "conservation":
            lo, hi = point["conservation_range"]
            left, right = xmap(lo, xmin, xmax, x0, width), xmap(hi, xmin, xmax, x0, width)
            out.extend([
                f'<line x1="{left:.1f}" y1="{y}" x2="{right:.1f}" y2="{y}" stroke="{color}" stroke-width="5" stroke-linecap="round"/>',
                f'<line x1="{left:.1f}" y1="{y-7}" x2="{left:.1f}" y2="{y+7}" stroke="{color}" stroke-width="2"/>',
                f'<line x1="{right:.1f}" y1="{y-7}" x2="{right:.1f}" y2="{y+7}" stroke="{color}" stroke-width="2"/>',
                f'<text x="{right+9:.1f}" y="{y+5:.1f}" class="value-label">{lo:+.3f} to {hi:+.3f}</text>',
            ])
        elif "entropy_range" in point and metric == "entropy":
            lo, hi = point["entropy_range"]
            left, right = xmap(lo, xmin, xmax, x0, width), xmap(hi, xmin, xmax, x0, width)
            out.extend([
                f'<line x1="{left:.1f}" y1="{y}" x2="{right:.1f}" y2="{y}" stroke="{color}" stroke-width="5" stroke-linecap="round"/>',
                f'<line x1="{left:.1f}" y1="{y-7}" x2="{left:.1f}" y2="{y+7}" stroke="{color}" stroke-width="2"/>',
                f'<line x1="{right:.1f}" y1="{y-7}" x2="{right:.1f}" y2="{y+7}" stroke="{color}" stroke-width="2"/>',
                f'<text x="{right+9:.1f}" y="{y+5:.1f}" class="value-label">{lo:+.3f} to {hi:+.3f}</text>',
            ])
        else:
            value = float(point[metric])
            x = xmap(value, xmin, xmax, x0, width)
            out.append(f'<circle cx="{x:.1f}" cy="{y}" r="6.5" fill="{color}" stroke="#ffffff" stroke-width="1.5"/>')
            # Keep labels near the zero reference line from crossing it.
            if value < 0 and abs(value) < 0.012:
                dx, anchor = -9, "end"
            else:
                dx, anchor = 9, "start"
            out.append(f'<text x="{x+dx:.1f}" y="{y+5:.1f}" text-anchor="{anchor}" class="value-label">{value:+.3f}</text>')
    ticks = [xmin, (xmin + xmax) / 2, 0.0, xmax]
    for tick in sorted(set(round(value, 6) for value in ticks)):
        x = xmap(tick, xmin, xmax, x0, width)
        out.append(f'<line x1="{x:.1f}" y1="{chart_y+row_height*len(points)-7}" x2="{x:.1f}" y2="{chart_y+row_height*len(points)-2}" class="axis"/>')
        out.append(f'<text x="{x:.1f}" y="{chart_y+row_height*len(points)+17}" text-anchor="middle" class="tick">{tick:+.2f}</text>')
    out.append(f'<text x="{x0+width/2:.1f}" y="{chart_y+row_height*len(points)+42}" text-anchor="middle" class="axis-label">Δ = paired − unpaired</text>')
    return "\n".join(out)


def render_svg(summary: dict[str, Any], conservation_rows: list[dict[str, str]],
               entropy_rows: list[dict[str, str]], bovid_rows: list[dict[str, str]]) -> str:
    points = build_points(summary, conservation_rows, entropy_rows, bovid_rows)
    taxa = int(summary["n_taxa"])
    valid = int(summary["valid_pair_observations"])
    possible = int(summary["possible_pair_observations"])
    canonical = float(summary["canonical_or_wobble_pair_fraction"]) * 100
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1320" height="1120" viewBox="0 0 1320 1120" role="img" aria-labelledby="title desc">
<title id="title">Mitochondrial tRNA-Phe(GAA) locus context and paired/unpaired sensitivity</title>
<desc id="desc">Schematic conserved mitochondrial gene neighborhood in 16 vertebrate mitogenomes, followed by paired-minus-unpaired conservation and entropy differences for 5S rRNA, broad tRNA, the restricted tRNA-Phe locus, coverage thresholds, and one-Bovidae representative sensitivity.</desc>
<style>
  text {{ font-family: Arial, Helvetica, sans-serif; fill: #202733; }}
  .main-title {{ font-size: 25px; font-weight: 700; }}
  .subtitle {{ font-size: 14px; fill: #4b5563; }}
  .panel-title {{ font-size: 18px; font-weight: 700; }}
  .gene-title {{ font-size: 16px; font-weight: 700; }}
  .gene-sub {{ font-size: 12px; }}
  .evidence {{ font-size: 13px; font-weight: 600; }}
  .row-label {{ font-size: 13px; }}
  .value-label {{ font-size: 12px; font-weight: 700; }}
  .tick {{ font-size: 11px; fill: #4b5563; }}
  .axis-label {{ font-size: 12px; fill: #4b5563; }}
  .axis {{ stroke: #9ca3af; stroke-width: 1; }}
  .zero {{ stroke: #111827; stroke-width: 1.6; stroke-dasharray: 5 4; }}
  .note {{ font-size: 12px; fill: #4b5563; }}
  .gene {{ stroke: #334155; stroke-width: 1.2; rx: 8; }}
</style>
<rect width="1320" height="1050" fill="#ffffff"/>
<text x="48" y="48" class="main-title">One validated mitochondrial tRNA-Phe(GAA) locus</text>
<text x="48" y="74" class="subtitle">Genomic context and the sensitivity of paired-versus-unpaired conservation estimates</text>

<text x="48" y="116" class="panel-title">A. Conserved local gene neighborhood</text>
<text x="1270" y="116" text-anchor="end" class="note">Observed in {taxa}/{taxa} complete mitogenomes · schematic, not to scale</text>
<line x1="130" y1="188" x2="1190" y2="188" stroke="#64748b" stroke-width="2"/>
<rect x="105" y="145" width="220" height="86" class="gene" fill="#d9f0e3"/>
<text x="215" y="178" text-anchor="middle" class="gene-title">tRNA-Phe</text>
<text x="215" y="199" text-anchor="middle" class="gene-sub">anticodon GAA</text>
<path d="M330 188 L365 188 M355 180 L365 188 L355 196" stroke="#64748b" stroke-width="2" fill="none"/>
<rect x="370" y="145" width="220" height="86" class="gene" fill="#dbeafe"/>
<text x="480" y="178" text-anchor="middle" class="gene-title">12S rRNA</text>
<text x="480" y="199" text-anchor="middle" class="gene-sub">mitochondrial small subunit rRNA</text>
<path d="M595 188 L630 188 M620 180 L630 188 L620 196" stroke="#64748b" stroke-width="2" fill="none"/>
<rect x="635" y="145" width="220" height="86" class="gene" fill="#d9f0e3"/>
<text x="745" y="178" text-anchor="middle" class="gene-title">tRNA-Val</text>
<text x="745" y="199" text-anchor="middle" class="gene-sub">adjacent tRNA gene</text>
<path d="M860 188 L895 188 M885 180 L895 188 L885 196" stroke="#64748b" stroke-width="2" fill="none"/>
<rect x="900" y="145" width="220" height="86" class="gene" fill="#dbeafe"/>
<text x="1010" y="178" text-anchor="middle" class="gene-title">16S rRNA</text>
<text x="1010" y="199" text-anchor="middle" class="gene-sub">mitochondrial large subunit rRNA</text>

<rect x="48" y="250" width="290" height="58" rx="7" fill="#f8fafc" stroke="#cbd5e1"/>
<text x="193" y="274" text-anchor="middle" class="evidence">{taxa}/{taxa} genomes</text>
<text x="193" y="293" text-anchor="middle" class="note">exactly one tRNA-Phe feature</text>
<rect x="355" y="250" width="290" height="58" rx="7" fill="#f8fafc" stroke="#cbd5e1"/>
<text x="500" y="274" text-anchor="middle" class="evidence">22 tRNA features</text>
<text x="500" y="293" text-anchor="middle" class="note">in every complete mitogenome</text>
<rect x="662" y="250" width="290" height="58" rx="7" fill="#f8fafc" stroke="#cbd5e1"/>
<text x="807" y="274" text-anchor="middle" class="evidence">GAA mapped in {taxa}/{taxa}</text>
<text x="807" y="293" text-anchor="middle" class="note">Rfam model 33–35 / Stockholm 44–46</text>
<rect x="969" y="250" width="303" height="58" rx="7" fill="#f8fafc" stroke="#cbd5e1"/>
<text x="1120" y="274" text-anchor="middle" class="evidence">{valid}/{possible} valid pair observations</text>
<text x="1120" y="293" text-anchor="middle" class="note">{canonical:.1f}% canonical or wobble</text>

<text x="48" y="356" class="panel-title">B. Paired-minus-unpaired mean conservation</text>
<text x="1270" y="356" text-anchor="end" class="note">Positive values: paired columns have higher mean conservation</text>
{panel_svg(points, "conservation", "Paired minus unpaired mean conservation", 330, 378, 720, 30, -0.22, 0.06, "#2563a6", "B")}
<text x="48" y="720" class="panel-title">C. Paired-minus-unpaired mean entropy</text>
<text x="1270" y="720" text-anchor="end" class="note">Positive values: paired columns have higher mean entropy</text>
{panel_svg(points, "entropy", "Paired minus unpaired mean entropy (bits)", 330, 720, 720, 30, -0.12, 0.55, "#c46b25", "C")}

<text x="48" y="1080" class="note">Δ = paired − unpaired. Descriptive means only; no p-values or phylogenetic weighting. Nine of the 16 locus taxa are Bovidae.</text>
<text x="48" y="1100" class="note">One-Bovidae range: nine possible representatives, each combined with all seven non-Bovidae taxa; this is not a phylogenetic correction.</text>
</svg>
'''


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    parser.add_argument("--conservation", type=Path, default=CONSERVATION)
    parser.add_argument("--entropy", type=Path, default=ENTROPY)
    parser.add_argument("--bovid", type=Path, default=BOVID)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    svg = render_svg(summary, read_csv(args.conservation), read_csv(args.entropy), read_csv(args.bovid))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(svg, encoding="utf-8")
    print(json.dumps({"output": str(args.output.relative_to(ROOT)), "bytes": args.output.stat().st_size}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

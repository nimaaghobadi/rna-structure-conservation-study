#!/usr/bin/env python3
"""Create the main manuscript SVG for paired/unpaired conservation and entropy."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
CONSERVATION = ROOT / "results/tables/conservation_by_structure.csv"
ENTROPY = ROOT / "results/tables/entropy_by_structure.csv"
OUTPUT = ROOT / "results/figures/figure_1_conservation_entropy.svg"

FAMILIES = ("5S rRNA", "tRNA")
STATES = ("paired", "unpaired")
COLORS = {"paired": "#2563a6", "unpaired": "#c46b25"}
COLUMN_COUNTS = {
    ("5S rRNA", "paired"): 68, ("5S rRNA", "unpaired"): 52,
    ("tRNA", "paired"): 42, ("tRNA", "unpaired"): 29,
}
SEQUENCE_COUNTS = {"5S rRNA": 712, "tRNA": 954}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def values_by_group(rows: list[dict[str, str]]) -> dict[tuple[str, str], float]:
    result: dict[tuple[str, str], float] = {}
    for row in rows:
        key = (row["family"], row["structure_state"])
        if key in result:
            raise ValueError(f"Duplicate metric row for {key}.")
        result[key] = float(row["mean"])
    expected = {(family, state) for family in FAMILIES for state in STATES}
    if set(result) != expected:
        raise ValueError(f"Expected exactly four family/state metric rows; observed {sorted(result)}.")
    return result


def panel_svg(values: dict[tuple[str, str], float], title: str, panel_x: int,
              max_value: float, tick_step: float, panel_letter: str,
              unit: str) -> str:
    plot_x = panel_x + 68
    plot_y = 160
    plot_width = 430
    plot_height = 330
    baseline = plot_y + plot_height
    out = [
        f'<text x="{panel_x}" y="118" class="panel-title">{panel_letter}. {title}</text>',
        f'<line x1="{plot_x}" y1="{plot_y}" x2="{plot_x}" y2="{baseline}" class="axis"/>',
        f'<line x1="{plot_x}" y1="{baseline}" x2="{plot_x+plot_width}" y2="{baseline}" class="axis"/>',
    ]
    tick = 0.0
    while tick <= max_value + 1e-9:
        y = baseline - tick / max_value * plot_height
        out.append(f'<line x1="{plot_x}" y1="{y:.1f}" x2="{plot_x+plot_width}" y2="{y:.1f}" class="grid"/>')
        out.append(f'<text x="{plot_x-12}" y="{y+4:.1f}" text-anchor="end" class="tick">{tick:g}</text>')
        tick += tick_step

    group_centres = {"5S rRNA": plot_x + 112, "tRNA": plot_x + 322}
    bar_width = 62
    bar_gap = 14
    for family in FAMILIES:
        centre = group_centres[family]
        for j, state in enumerate(STATES):
            value = values[(family, state)]
            x = centre + (j - 0.5) * (bar_width + bar_gap) - bar_width / 2
            height = value / max_value * plot_height
            y = baseline - height
            out.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width}" height="{height:.1f}" rx="5" fill="{COLORS[state]}"/>')
            out.append(f'<text x="{x+bar_width/2:.1f}" y="{y-10:.1f}" text-anchor="middle" class="value">{value:.3f}</text>')
            out.append(f'<text x="{x+bar_width/2:.1f}" y="{baseline+22}" text-anchor="middle" class="bar-label">{state}</text>')
            out.append(f'<text x="{x+bar_width/2:.1f}" y="{baseline+39}" text-anchor="middle" class="column-label">n={COLUMN_COUNTS[(family,state)]} columns</text>')
        out.append(f'<text x="{centre:.1f}" y="{baseline+68}" text-anchor="middle" class="family-label">{family}</text>')
        out.append(f'<text x="{centre:.1f}" y="{baseline+87}" text-anchor="middle" class="sequence-label">{SEQUENCE_COUNTS[family]} sequences</text>')
    out.append(f'<text x="{plot_x+plot_width/2}" y="{baseline+118}" text-anchor="middle" class="axis-label">{unit}</text>')
    return "\n".join(out)


def render_svg(conservation_rows: list[dict[str, str]],
               entropy_rows: list[dict[str, str]]) -> str:
    conservation = values_by_group(conservation_rows)
    entropy = values_by_group(entropy_rows)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="680" viewBox="0 0 1200 680" role="img" aria-labelledby="title desc">
<title id="title">Paired and unpaired sequence conservation and entropy in 5S rRNA and tRNA</title>
<desc id="desc">Two panels compare mean per-column conservation and Shannon entropy for paired and unpaired columns in Rfam 5S rRNA and tRNA seed alignments. Values are descriptive means across columns.</desc>
<style>
  text {{ font-family: Arial, Helvetica, sans-serif; fill: #202733; }}
  .main-title {{ font-size: 25px; font-weight: 700; }}
  .subtitle {{ font-size: 14px; fill: #4b5563; }}
  .panel-title {{ font-size: 19px; font-weight: 700; }}
  .axis {{ stroke: #64748b; stroke-width: 1.4; }}
  .grid {{ stroke: #e2e8f0; stroke-width: 1; }}
  .tick {{ font-size: 12px; fill: #64748b; }}
  .value {{ font-size: 13px; font-weight: 700; }}
  .bar-label {{ font-size: 12px; }}
  .column-label {{ font-size: 11px; fill: #64748b; }}
  .family-label {{ font-size: 15px; font-weight: 700; }}
  .sequence-label {{ font-size: 11px; fill: #64748b; }}
  .axis-label {{ font-size: 12px; fill: #4b5563; }}
  .legend {{ font-size: 13px; }}
  .note {{ font-size: 12px; fill: #4b5563; }}
</style>
<rect width="1200" height="680" fill="#ffffff"/>
<text x="48" y="48" class="main-title">Paired and unpaired positions show different sequence patterns</text>
<text x="48" y="75" class="subtitle">Rfam 15.1 seed alignments · descriptive column means · paired status follows the consensus secondary-structure annotation</text>
<rect x="820" y="92" width="16" height="16" rx="3" fill="{COLORS['paired']}"/>
<text x="844" y="105" class="legend">Paired</text>
<rect x="925" y="92" width="16" height="16" rx="3" fill="{COLORS['unpaired']}"/>
<text x="949" y="105" class="legend">Unpaired</text>
{panel_svg(conservation, "Mean single-site conservation", 55, 1.0, 0.25, "A", "Mean conservation (0–1)")}
{panel_svg(entropy, "Mean Shannon entropy", 650, 2.0, 0.5, "B", "Entropy (bits)")}
<text x="48" y="660" class="note">Each bar is the mean across eligible alignment columns (column counts shown below). Columns are not independent replicates; no p-values or inferential error bars are shown.</text>
</svg>
'''


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conservation", type=Path, default=CONSERVATION)
    parser.add_argument("--entropy", type=Path, default=ENTROPY)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    svg = render_svg(read_rows(args.conservation), read_rows(args.entropy))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(svg, encoding="utf-8")
    print(f"Wrote {args.output.relative_to(ROOT)} ({args.output.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

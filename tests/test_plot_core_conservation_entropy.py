from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import plot_core_conservation_entropy as FIG


def rows(metric_values):
    return [
        {"family": family, "structure_state": state, "mean": str(value)}
        for (family, state), value in metric_values.items()
    ]


def test_values_by_group_requires_all_four_family_state_cells():
    data = rows({
        ("5S rRNA", "paired"): 0.649,
        ("5S rRNA", "unpaired"): 0.720,
        ("tRNA", "paired"): 0.474,
        ("tRNA", "unpaired"): 0.668,
    })
    result = FIG.values_by_group(data)
    assert result[("5S rRNA", "paired")] == 0.649
    assert result[("tRNA", "unpaired")] == 0.668


def test_core_figure_contains_means_column_counts_and_noninferential_note():
    conservation = rows({
        ("5S rRNA", "paired"): 0.649,
        ("5S rRNA", "unpaired"): 0.720,
        ("tRNA", "paired"): 0.474,
        ("tRNA", "unpaired"): 0.668,
    })
    entropy = rows({
        ("5S rRNA", "paired"): 1.201,
        ("5S rRNA", "unpaired"): 1.033,
        ("tRNA", "paired"): 1.656,
        ("tRNA", "unpaired"): 1.167,
    })
    svg = FIG.render_svg(conservation, entropy)
    assert "Mean single-site conservation" in svg
    assert "Mean Shannon entropy" in svg
    assert "n=68 columns" in svg and "n=29 columns" in svg
    assert "712 sequences" in svg and "954 sequences" in svg
    assert "no p-values or inferential error bars" in svg

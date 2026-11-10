from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import plot_tRNA_Phe_locus_sensitivity as FIG


def fixture_data():
    summary = {
        "n_taxa": 16, "valid_pair_observations": 324, "possible_pair_observations": 336,
        "canonical_or_wobble_pair_fraction": 305 / 324,
        "paired_minus_unpaired_conservation": 0.0363574,
        "paired_minus_unpaired_entropy_bits": -0.0950199,
        "valid_fraction_threshold_sensitivity": {
            "0.5": {"paired_minus_unpaired_conservation": 0.0346797,
                    "paired_minus_unpaired_entropy_bits": -0.0957614},
            "0.7": {"paired_minus_unpaired_conservation": 0.0212908,
                    "paired_minus_unpaired_entropy_bits": -0.0686114},
            "0.9": {"paired_minus_unpaired_conservation": -0.0057292,
                    "paired_minus_unpaired_entropy_bits": 0.0159537},
        },
    }
    conservation = [
        {"family": "5S rRNA", "structure_state": "paired", "mean": "0.649385"},
        {"family": "5S rRNA", "structure_state": "unpaired", "mean": "0.720457"},
        {"family": "tRNA", "structure_state": "paired", "mean": "0.474460"},
        {"family": "tRNA", "structure_state": "unpaired", "mean": "0.667855"},
    ]
    entropy = [
        {"family": "5S rRNA", "structure_state": "paired", "mean": "1.201186"},
        {"family": "5S rRNA", "structure_state": "unpaired", "mean": "1.033100"},
        {"family": "tRNA", "structure_state": "paired", "mean": "1.655722"},
        {"family": "tRNA", "structure_state": "unpaired", "mean": "1.167336"},
    ]
    bovid = [
        {"paired_minus_unpaired_conservation": "-0.01564",
         "paired_minus_unpaired_entropy_bits": "0.05686"},
        {"paired_minus_unpaired_conservation": "0.00324",
         "paired_minus_unpaired_entropy_bits": "-0.03937"},
    ]
    return summary, conservation, entropy, bovid


def test_build_points_preserves_broad_family_and_locus_specific_differences():
    summary, conservation, entropy, bovid = fixture_data()
    points = FIG.build_points(summary, conservation, entropy, bovid)
    assert len(points) == 7
    assert points[0]["conservation"] < 0
    assert points[1]["conservation"] < 0
    assert points[2]["conservation"] > 0
    assert points[5]["conservation"] < 0
    assert points[6]["conservation_range"] == (-0.01564, 0.00324)


def test_svg_contains_genomic_context_and_noninferential_sensitivity_plot():
    summary, conservation, entropy, bovid = fixture_data()
    svg = FIG.render_svg(summary, conservation, entropy, bovid)
    assert "tRNA-Phe" in svg and "12S rRNA" in svg and "16S rRNA" in svg
    assert "94.1% canonical or wobble" in svg
    assert "One-Bovidae sensitivity (range)" in svg
    assert "no p-values or phylogenetic weighting" in svg
    assert "0.000" not in svg

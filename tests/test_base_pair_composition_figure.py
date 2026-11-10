"""Integrity tests for the reproducible base-pair composition figure."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "plot_rna_biology_figures.py"
SPEC = importlib.util.spec_from_file_location("base_pair_figure", SCRIPT)
FIGURE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FIGURE)

SOURCE = ROOT / "data" / "tables" / "base_pair_summary.csv"
TEST_CSV = """family,n_valid_pair_observations,n_possible_pair_observations,valid_pair_observation_fraction,canonical_or_wobble_fraction,count_AU,count_UA,count_GC,count_CG,count_GU,count_UG,count_noncanonical_or_other
5S rRNA,23823,24208,0.9840961665565102,0.917978424211896,2551,3019,6056,6944,1287,2012,1954
tRNA,19973,20034,0.9969551762004593,0.9641516046662995,3639,3678,6664,3972,834,470,716
"""


def test_archived_metrics_validate():
    data = FIGURE.load_and_validate(SOURCE)
    assert data["5S rRNA"]["valid"] == 23823
    assert data["tRNA"]["possible"] == 20034
    assert data["5S rRNA"]["fraction_GU"] == pytest.approx(1287 / 23823)


def test_svg_is_deterministic_and_contains_denominators():
    data = FIGURE.load_and_validate(SOURCE)
    first = FIGURE.render_svg(data)
    assert first == FIGURE.render_svg(data)
    committed_svg = ROOT / "figures" / "base_pair_composition.svg"
    assert committed_svg.read_text(encoding="utf-8") == first
    assert "23,823/24,208 (98.4%)" in first
    assert "19,973/20,034 (99.7%)" in first
    assert "Other / noncanonical" in first


def test_rejects_inconsistent_pair_counts(tmp_path):
    path = tmp_path / "bad_counts.csv"
    path.write_text(
        TEST_CSV.replace(
            "0.917978424211896,2551",
            "0.917978424211896,2552",
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="do not sum"):
        FIGURE.load_and_validate(path)


def test_rejects_inconsistent_coverage(tmp_path):
    path = tmp_path / "bad_coverage.csv"
    path.write_text(TEST_CSV.replace("0.9840961665565102", "0.9"), encoding="utf-8")
    with pytest.raises(ValueError, match="coverage disagrees"):
        FIGURE.load_and_validate(path)

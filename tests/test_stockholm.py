from pathlib import Path

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from src.stockholm import (
    classify_structure,
    find_overlapping_sequence_ids,
    read_ss_cons,
    validate_alignment,
    validate_wuss,
)


def test_wuss_pairing_and_dot_is_unpaired():
    assert classify_structure("((..))") == [
        "paired", "paired", "unpaired", "unpaired", "paired", "paired"
    ]


def test_explicit_single_strand_wuss_symbols():
    assert classify_structure("-,:_~") == [
        "unpaired", "unpaired", "unpaired", "unpaired", "unpaired"
    ]


def test_rf_annotation_marks_insertions_separately():
    assert classify_structure("....", "AC..") == [
        "unpaired", "unpaired", "insertion", "insertion"
    ]


def test_rf_annotation_length_must_match_structure():
    with pytest.raises(ValueError, match="RF reference annotation"):
        classify_structure("..", ".")


def test_wuss_letter_pseudoknot_pairing():
    assert classify_structure("Aa") == ["paired", "paired"]


@pytest.mark.parametrize("structure", ["(", ")", "A", "a", "(x)", "?."])
def test_invalid_or_unbalanced_wuss_is_rejected(structure):
    with pytest.raises(ValueError):
        validate_wuss(structure)


def test_missing_ss_cons_raises(tmp_path: Path):
    p = tmp_path / "no_ss.sto"
    p.write_text("# STOCKHOLM 1.0\nseq1 AC-G\n//\n", encoding="utf-8")
    with pytest.raises(ValueError):
        read_ss_cons(p)


def test_alignment_and_ss_length_validation(tmp_path: Path):
    p = tmp_path / "ok.sto"
    p.write_text(
        "# STOCKHOLM 1.0\nseq1 AC-G\nseq2 ACUG\n"
        "#=GC SS_cons (..)\n//\n",
        encoding="utf-8",
    )
    assert validate_alignment(p) == {
        "n_sequences": 2,
        "alignment_length": 4,
        "ss_cons_length": 4,
        "rf_annotation_present": False,
    }


def test_rf_length_validation(tmp_path: Path):
    p = tmp_path / "bad_rf.sto"
    p.write_text(
        "# STOCKHOLM 1.0\nseq1 AC-G\nseq2 ACUG\n"
        "#=GC SS_cons (..)\n#=GC RF xxx\n//\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="RF.*length"):
        validate_alignment(p)


def test_mismatched_ss_length_raises(tmp_path: Path):
    p = tmp_path / "bad.sto"
    p.write_text(
        "# STOCKHOLM 1.0\nseq1 AC-G\nseq2 ACUG\n"
        "#=GC SS_cons (...)\n//\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        validate_alignment(p)


def test_current_rfam_trna_annotation_fixture():
    # Short annotation fixture copied from the official RF00005 Stockholm seed:
    # https://rfam.org/family/RF00005/alignment (accessed 2026-10-30).
    ss_cons = (
        "(((((((,.,.<<<<___...___..._>>>>,.<.<<<<_.______>>>>.>,,"
        "........................,...,<<<.<<_.___.__.._>>>>.>.))))))).:"
    )
    reference = (
        "ggggaugU.A.GCuuAaU...GGU...AaaGCa.u.ugGaC.UuauAAuCca.aAg"
        "........................g...cgcg.GGU.UCg.Aa..UCCcg.c.caucccc.A"
    )
    qc = validate_wuss(ss_cons)
    states = classify_structure(ss_cons, reference)
    assert qc["length"] == len(reference)
    assert len(states) == len(ss_cons)
    assert states.count("insertion") == reference.count(".")
    assert states.count("paired") > 0
    assert states.count("unpaired") > 0


def test_coordinate_overlap_audit_flags_only_overlapping_sequence():
    alignment = MultipleSeqAlignment([
        SeqRecord(Seq("ACG"), id="ACC/10-1"),
        SeqRecord(Seq("ACG"), id="ACC/1-5"),
        SeqRecord(Seq("ACG"), id="ACC/20-15"),
        SeqRecord(Seq("ACG"), id="OTHER/5-1"),
    ])
    original_ids = [record.id for record in alignment]
    assert find_overlapping_sequence_ids(alignment) == {"ACC/1-5"}
    assert [record.id for record in alignment] == original_ids

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest
from src.metrics import column_metrics, base_pair_metrics, base_pair_family_summary, greedy_identity_representatives

def make_alignment(sequences):
    return MultipleSeqAlignment([SeqRecord(Seq(s), id=f"s{i}") for i, s in enumerate(sequences)])

def test_conservation_and_entropy_ignore_gaps_and_ambiguity():
    aln = make_alignment(["A-", "AC", "AN", "AG"])
    df = column_metrics(aln, ["unpaired", "unpaired"], "demo")
    assert df.loc[0, "conservation"] == pytest.approx(1.0)
    assert df.loc[0, "valid_fraction"] == pytest.approx(1.0)
    assert df.loc[1, "n_valid_bases"] == 2
    assert df.loc[1, "gap_count"] == 1
    assert df.loc[1, "ambiguous_count"] == 1

def test_pair_metrics_has_explicit_valid_denominator():
    aln = make_alignment(["AU", "GC", "GN", "A-"])
    df = base_pair_metrics(aln, "()", "demo")
    assert len(df) == 1
    assert df.loc[0, "n_both_valid"] == 2
    assert df.loc[0, "valid_pair_fraction"] == pytest.approx(0.5)
    assert df.loc[0, "canonical_or_wobble_fraction"] == pytest.approx(1.0)
    # Rfam BPcons: canonical-pair indicator summed over all aligned sequences.
    assert df.loc[0, "canonical_or_wobble_fraction_all_sequences"] == pytest.approx(2 / 4)


def test_pair_metrics_excludes_rf_insertion_columns():
    aln = make_alignment(["AU", "GC"])
    df = base_pair_metrics(aln, "()", "demo", reference="..")
    assert df.empty


def test_pair_metrics_checks_reference_length():
    aln = make_alignment(["AU", "GC"])
    with pytest.raises(ValueError, match="RF reference annotation"):
        base_pair_metrics(aln, "()", "demo", reference=".")


def test_pair_metrics_reports_noncanonical_valid_pairs():
    aln = make_alignment(["AA", "GC", "GN", "A-"])
    df = base_pair_metrics(aln, "()", "demo")
    assert df.loc[0, "n_both_valid"] == 2
    assert df.loc[0, "count_noncanonical_or_other"] == 1
    assert df.loc[0, "canonical_or_wobble_fraction"] == pytest.approx(0.5)


def test_family_pair_summary_uses_pooled_valid_denominator():
    aln = make_alignment(["AU", "GC", "AA", "GN"])
    pair_df = base_pair_metrics(aln, "()", "demo")
    summary = base_pair_family_summary(pair_df)
    row = summary.iloc[0]
    assert row["n_annotated_base_pairs"] == 1
    assert row["n_valid_pair_observations"] == 3
    assert row["n_possible_pair_observations"] == 4
    assert row["valid_pair_observation_fraction"] == pytest.approx(0.75)
    assert row["canonical_or_wobble_fraction"] == pytest.approx(2 / 3)
    assert row["canonical_or_wobble_fraction_all_possible"] == pytest.approx(2 / 4)
    # Family pooled all-sequence fraction equals the mean per-pair BPcons value.
    assert row["canonical_or_wobble_fraction_all_possible"] == pytest.approx(
        pair_df["canonical_or_wobble_fraction_all_sequences"].mean()
    )
    assert row["noncanonical_or_other_fraction"] == pytest.approx(1 / 3)
    assert row["fraction_AU"] == pytest.approx(1 / 3)


def test_greedy_identity_filter_is_deterministic():
    aln = make_alignment(["ACGU", "ACGU", "ACGA"])
    assert greedy_identity_representatives(aln, threshold=0.90) == {"s0", "s2"}

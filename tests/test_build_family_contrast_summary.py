from scripts.build_family_contrast_summary import build_rows


def test_family_contrast_summary_uses_full_precision_means():
    qc = [
        {"family": "5S rRNA", "n_sequences": "712"},
        {"family": "tRNA", "n_sequences": "954"},
    ]
    conservation = [
        {"family": "5S rRNA", "structure_state": "paired", "mean": "0.6493851318991198", "count": "68"},
        {"family": "5S rRNA", "structure_state": "unpaired", "mean": "0.7204566551928051", "count": "52"},
        {"family": "tRNA", "structure_state": "paired", "mean": "0.4744601951577073", "count": "42"},
        {"family": "tRNA", "structure_state": "unpaired", "mean": "0.6678548319279626", "count": "29"},
    ]
    entropy = [
        {"family": "5S rRNA", "structure_state": "paired", "mean": "1.201185780024612"},
        {"family": "5S rRNA", "structure_state": "unpaired", "mean": "1.0330997147836134"},
        {"family": "tRNA", "structure_state": "paired", "mean": "1.6557222540353518"},
        {"family": "tRNA", "structure_state": "unpaired", "mean": "1.1673364703850397"},
    ]

    rows = build_rows(qc, conservation, entropy)

    assert rows[0]["paired_minus_unpaired_mean_conservation"] == "-0.0710715232936853"
    assert rows[1]["paired_minus_unpaired_mean_conservation"] == "-0.1933946367702553"
    assert rows[2]["paired_minus_unpaired_mean_conservation"] == "-0.12232311347657"
    assert rows[2]["paired_minus_unpaired_mean_entropy_bits"] == "0.320299718409314"

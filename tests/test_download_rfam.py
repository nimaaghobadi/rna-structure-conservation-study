import pytest

import json
from src.download_rfam import download, sha256, validate_stockholm_bytes


def test_valid_stockholm_alignment_passes_qc():
    data = (
        "# STOCKHOLM 1.0\n"
        "#=GF ID demo\n"
        "seq1 AC-G\n"
        "seq2 ACUG\n"
        "#=GC SS_cons (..)\n"
        "#=GC RF ACGU\n"
        "//\n"
    ).encode()
    assert validate_stockholm_bytes(data, "RFTEST") == {
        "n_sequences": 2,
        "alignment_length": 4,
        "ss_cons_length": 4,
        "rf_annotation_present": True,
        "wuss_validation": {
            "length": 4,
            "paired_positions": 2,
            "unpaired_positions": 2,
        },
    }


@pytest.mark.parametrize(
    "data",
    [
        b"<html>rate limited</html>",
        b"# STOCKHOLM 1.0\nseq1 AC-G\n",
        (
            "# STOCKHOLM 1.0\n"
            "seq1 AC-G\n"
            "seq2 ACUG\n"
            "#=GC SS_cons (...)\n"
            "//\n"
        ).encode(),
        (
            "# STOCKHOLM 1.0\n"
            "seq1 AC-G\n"
            "seq2 ACUG\n"
            "#=GC SS_cons (..x\n"
            "//\n"
        ).encode(),
        (
            "# STOCKHOLM 1.0\n"
            "seq1 AC-G\n"
            "seq2 ACUG\n"
            "#=GC SS_cons (..)\n"
            "#=GC RF xxx\n"
            "//\n"
        ).encode(),
    ],
)
def test_invalid_download_is_rejected(data):
    with pytest.raises((RuntimeError, ValueError)):
        validate_stockholm_bytes(data, "RFTEST")



def test_reused_seed_files_do_not_rewrite_unchanged_manifest(tmp_path, monkeypatch):
    import src.download_rfam as downloader

    raw = tmp_path / "raw"
    processed = tmp_path / "processed"
    raw.mkdir()
    processed.mkdir()
    monkeypatch.setattr(downloader, "RAW", raw)
    monkeypatch.setattr(downloader, "PROCESSED", processed)
    monkeypatch.setattr(
        downloader,
        "FAMILIES",
        {"RFTEST": {"name": "test family", "file": "test_seed.sto"}},
    )

    data = (
        "# STOCKHOLM 1.0\n"
        "#=GF ID demo\n"
        "seq1 AC-G\n"
        "seq2 ACUG\n"
        "#=GC SS_cons (..)\n"
        "#=GC RF ACGU\n"
        "//\n"
    ).encode()
    (raw / "test_seed.sto").write_bytes(data)

    previous = {
        "retrieved_at_utc": "2026-01-01T00:00:00+00:00",
        "source": "Official Rfam family alignment endpoint",
        "families": {
            "RFTEST": {
                "name": "test family",
                "file": "data/raw/test_seed.sto",
                "url": "https://rfam.org/family/RFTEST/alignment",
                "sha256": sha256(data),
                "bytes": len(data),
                "action": "downloaded",
                "stockholm_header": data[:200].decode(),
                "response_content_type": "text/plain",
                "final_url": "https://rfam.org/family/RFTEST/alignment",
                "validation": validate_stockholm_bytes(data, "RFTEST"),
            }
        },
    }
    manifest_path = processed / "download_manifest.json"
    original = json.dumps(previous, indent=2, ensure_ascii=False) + "\n"
    manifest_path.write_text(original, encoding="utf-8")

    result = download()

    assert result == previous
    assert manifest_path.read_text(encoding="utf-8") == original

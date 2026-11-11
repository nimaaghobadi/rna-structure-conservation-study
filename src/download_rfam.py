"""Download official Rfam seed alignments and record immutable provenance."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

from Bio import AlignIO

from src.stockholm import validate_wuss

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
RFAM_RELEASE = "15.1"
RFAM_ARCHIVE_URL = "https://rfam.org/downloads"
FAMILIES = {
    "RF00001": {"name": "5S rRNA", "file": "rf00001_seed.sto", "sha256": "18f91da6b0d9a668caf3d50849db437584140923e5477f7b21ed59453100f77c"},
    "RF00005": {"name": "tRNA", "file": "rf00005_seed.sto", "sha256": "d1a66517066e6d4fbad3bd382934406df3018ff9bdaf72f22732f4c9740c9505"},
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def validate_stockholm_bytes(data: bytes, accession: str) -> dict:
    """Reject HTML, truncated, malformed, or annotation-mismatched downloads."""
    if not data.strip():
        raise RuntimeError(f"{accession}: empty response")
    text = data.decode("utf-8")
    if not text.startswith("# STOCKHOLM 1.0"):
        raise RuntimeError(f"{accession}: response lacks the Stockholm 1.0 header")
    if not any(line.strip() == "//" for line in text.splitlines()):
        raise RuntimeError(f"{accession}: Stockholm terminator // is missing")

    alignment = AlignIO.read(io.StringIO(text), "stockholm")
    if len(alignment) == 0:
        raise RuntimeError(f"{accession}: no sequences found")
    lengths = {len(record.seq) for record in alignment}
    if len(lengths) != 1:
        raise RuntimeError(f"{accession}: aligned sequences have unequal lengths")

    ss_cons_chunks = []
    for line in text.splitlines():
        if line.startswith("#=GC SS_cons"):
            fields = line.split()
            if len(fields) >= 3:
                ss_cons_chunks.append(fields[2])
    if not ss_cons_chunks:
        raise RuntimeError(f"{accession}: #=GC SS_cons annotation is missing")
    ss_cons = "".join(ss_cons_chunks)
    alignment_length = next(iter(lengths))
    if len(ss_cons) != alignment_length:
        raise RuntimeError(
            f"{accession}: SS_cons length {len(ss_cons)} does not match "
            f"alignment length {alignment_length}"
        )
    try:
        wuss_qc = validate_wuss(ss_cons)
    except ValueError as exc:
        raise RuntimeError(f"{accession}: invalid SS_cons WUSS annotation: {exc}") from exc

    rf_chunks = []
    for line in text.splitlines():
        if line.startswith("#=GC RF"):
            fields = line.split()
            if len(fields) >= 3 and fields[1] == "RF":
                rf_chunks.append(fields[2])
    rf = "".join(rf_chunks) if rf_chunks else None
    if rf is not None and len(rf) != alignment_length:
        raise RuntimeError(
            f"{accession}: RF annotation length {len(rf)} does not match "
            f"alignment length {alignment_length}"
        )
    return {
        "n_sequences": len(alignment),
        "alignment_length": alignment_length,
        "ss_cons_length": len(ss_cons),
        "rf_annotation_present": rf is not None,
        "wuss_validation": wuss_qc,
    }


def download(force: bool = False) -> dict:
    RAW.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    manifest_path = PROCESSED / "download_manifest.json"

    previous_manifest = None
    if manifest_path.exists():
        try:
            previous_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous_manifest = None

    previous_families = (
        previous_manifest.get("families", {})
        if isinstance(previous_manifest, dict) else {}
    )
    retrieved_at = (
        previous_manifest.get("retrieved_at_utc")
        if isinstance(previous_manifest, dict) else None
    )
    manifest = {
        "retrieved_at_utc": retrieved_at or datetime.now(timezone.utc).isoformat(),
        "rfam_release": RFAM_RELEASE,
        "source": "Official Rfam 15.1 family alignment endpoints",
        "families": {},
    }
    downloaded_any = False
    all_reused_unchanged = (
        previous_manifest is not None
        and not force
        and previous_manifest.get("rfam_release", RFAM_RELEASE) == RFAM_RELEASE
        and set(previous_families) == set(FAMILIES)
    )

    for accession, info in FAMILIES.items():
        path = RAW / info["file"]
        content_type, final_url = None, None
        if path.exists() and not force:
            data, action = path.read_bytes(), "existing_file_reused"
        else:
            url = f"https://rfam.org/family/{accession}/alignment"
            request = Request(
                url,
                headers={
                    "User-Agent": (
                        "rna-structure-conservation/0.1 "
                        "reproducible academic project"
                    )
                },
            )
            with urlopen(request, timeout=60) as response:
                data = response.read()
                content_type = response.headers.get("Content-Type", "")
                final_url = response.geturl()
            action = "downloaded"
            downloaded_any = True

        # Validate both new downloads and reused local files before accepting them.
        qc = validate_stockholm_bytes(data, accession)
        checksum = sha256(data)
        expected_checksum = info.get("sha256")
        if expected_checksum and checksum != expected_checksum:
            raise RuntimeError(
                f"{accession}: SHA-256 {checksum} does not match the archived "
                f"Rfam {RFAM_RELEASE} seed ({expected_checksum}). The family endpoint "
                f"may have moved to a different release; retrieve the pinned release "
                f"from the Rfam 15.1 release archive linked there and verify its SHA-256 value."
            )
        previous_family = previous_families.get(accession, {})
        if (
            action == "existing_file_reused"
            and previous_family.get("sha256") == checksum
        ):
            # Preserve the exact provenance record. A validation-only run must
            # not create a new manifest timestamp or an unnecessary Git commit.
            manifest["families"][accession] = previous_family
        else:
            all_reused_unchanged = False
            manifest["families"][accession] = {
                "name": info["name"],
                "file": str(path.relative_to(ROOT)).replace("\\", "/"),
                "url": f"https://rfam.org/family/{accession}/alignment",
                "sha256": checksum,
                "bytes": len(data),
                "action": action,
                "stockholm_header": data[:200].decode("utf-8", errors="replace"),
                "response_content_type": content_type,
                "final_url": final_url,
                "validation": qc,
            }
        if action == "downloaded":
            path.write_bytes(data)

    if all_reused_unchanged:
        return previous_manifest

    if downloaded_any:
        manifest["retrieved_at_utc"] = datetime.now(timezone.utc).isoformat()
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace local raw alignments with a fresh official response",
    )
    args = parser.parse_args()
    print(json.dumps(download(force=args.force), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

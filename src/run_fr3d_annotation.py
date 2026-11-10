"""Run pinned FR3D base-pair and stacking annotations on 9ZTL RNA chains."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FR3D_COMMIT = "494a81ccb024abda18fb4c62115a77c6df608472"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fr3d-root", type=Path, required=True)
    args = parser.parse_args()

    cif_source = ROOT / "data" / "structures" / "9ztl_updated.cif"
    cif_input = ROOT / "data" / "structures" / "9ZTL.cif"
    if not cif_source.is_file() or cif_source.stat().st_size < 1_000_000:
        raise FileNotFoundError("Official 9ZTL updated mmCIF is missing")
    shutil.copy2(cif_source, cif_input)

    commit_check = subprocess.run(
        ["git", "-C", str(args.fr3d_root), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=False,
    )
    actual_commit = commit_check.stdout.strip()
    if commit_check.returncode != 0 or actual_commit != FR3D_COMMIT:
        raise RuntimeError(
            f"FR3D checkout must be pinned to {FR3D_COMMIT}; found {actual_commit or 'unknown'}"
        )
    annotator = args.fr3d_root / "fr3d" / "classifiers" / "NA_pairwise_interactions.py"
    if not annotator.is_file():
        raise FileNotFoundError(annotator)

    output = ROOT / "supplementary" / "structural"
    output.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable, str(annotator),
        "--input", str(cif_input.parent),
        "--output", str(output),
        "-c", "basepair,stacking",
        "--chain", "b,5",
        "9ZTL",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    log = ROOT / "supplementary" / "provenance" / "fr3d_9ztl_run.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    command_template = (
        "python fr3d/classifiers/NA_pairwise_interactions.py "
        "--input data/structures --output supplementary/structural "
        "-c basepair,stacking --chain b,5 9ZTL"
    )
    log.write_text(
        "COMMAND TEMPLATE: " + command_template + "\n\nSTDOUT:\n" + result.stdout
        + "\nSTDERR:\n" + result.stderr,
        encoding="utf-8",
    )
    if result.returncode != 0:
        raise RuntimeError(f"FR3D annotation failed; inspect {log}")

    outputs = sorted(
        str(path.relative_to(ROOT))
        for path in output.iterdir() if path.is_file()
    )
    if not outputs:
        raise RuntimeError("FR3D completed but produced no annotation files")
    manifest = {
        "pdb_id": "9ZTL",
        "pdb_cif": str(cif_source.relative_to(ROOT)),
        "pdb_sha256": sha256(cif_source),
        "fr3d_repository": "https://github.com/BGSU-RNA/fr3d-python",
        "fr3d_commit": FR3D_COMMIT,
        "author_chains": {"5S_rRNA": "b", "P_site_tRNA": "5"},
        "annotation_categories": ["basepair", "stacking"],
        "command_template": [
            "python fr3d/classifiers/NA_pairwise_interactions.py",
            "--input data/structures",
            "--output supplementary/structural",
            "-c basepair,stacking",
            "--chain b,5",
            "9ZTL",
        ],
        "outputs": outputs,
        "note": (
            "FR3D output is a geometric interaction annotation of one deposited "
            "structure model, not experimental re-refinement or evolutionary evidence."
        ),
    }
    manifest_path = ROOT / "supplementary" / "provenance" / "fr3d_9ztl_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

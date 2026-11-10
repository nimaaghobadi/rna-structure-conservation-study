"""Run pinned ViennaRNA consensus-folding analyses on Rfam seed alignments."""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import statistics
from pathlib import Path

from Bio import AlignIO

from src.metrics import greedy_identity_representatives
from src.stockholm import read_rf_annotation, read_ss_cons, _make_pairs

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = {
    "RF00001": {"name": "5S rRNA", "file": "rf00001_seed.sto"},
    "RF00005": {"name": "tRNA", "file": "rf00005_seed.sto"},
}
SUBSETS = [("full", None), ("identity95", 0.95), ("identity90", 0.90)]
PROBABILITY_CUTOFF = 0.05


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_subset_stockholm(path: Path, records, ss_cons: str, rf: str, family: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# STOCKHOLM 1.0\n")
        handle.write(f"#=GF ID {family}\n")
        for record in records:
            sequence = str(record.seq).upper().replace("T", "U")
            handle.write(f"{record.id} {sequence}\n")
        handle.write(f"#=GC SS_cons {ss_cons}\n")
        handle.write(f"#=GC RF {rf}\n//\n")


def dotbracket_pairs(structure: str) -> set[tuple[int, int]]:
    stack = []
    pairs = set()
    for position, symbol in enumerate(structure, start=1):
        if symbol == "(":
            stack.append(position)
        elif symbol == ")":
            if not stack:
                raise ValueError("Unbalanced ViennaRNA dot-bracket structure")
            left = stack.pop()
            pairs.add((left, position))
    if stack:
        raise ValueError("Unbalanced ViennaRNA dot-bracket structure")
    return pairs


def rf_model_index(rf: str) -> list[int | None]:
    result = []
    model_position = 0
    for symbol in rf:
        if symbol == ".":
            result.append(None)
        else:
            model_position += 1
            result.append(model_position)
    return result


def read_rscape_pair_lookup():
    path = ROOT / "supplementary" / "tables" / "rscape_fr3d_structural_evidence.csv"
    lookup = {}
    if not path.exists():
        return lookup
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            left = int(row["rf_left_model_pos"])
            right = int(row["rf_right_model_pos"])
            lookup[(row["family"], min(left, right), max(left, right))] = row
    return lookup


def run_one_alignment(RNA, accession, family_name, subset, records, ss_cons, rf, rscape_lookup):
    aligned = []
    ambiguous_normalized = 0
    for record in records:
        source = str(record.seq).upper().replace("T", "U").replace(".", "-")
        normalized = []
        for symbol in source:
            if symbol in "ACGU-":
                normalized.append(symbol)
            else:
                normalized.append("N")
                ambiguous_normalized += 1
        aligned.append("".join(normalized))
    if not aligned or len({len(sequence) for sequence in aligned}) != 1:
        raise ValueError(f"Invalid aligned sequence set for {accession}/{subset}")

    fc = RNA.fold_compound(aligned)
    mfe_structure, mfe_energy = fc.mfe()
    fc.exp_params_rescale(mfe_energy)
    pf_structure, ensemble_energy = fc.pf()

    probabilities = {}
    plist = fc.plist_from_probs(PROBABILITY_CUTOFF)
    for item in plist:
        if int(item.i) == 0 and int(item.j) == 0:
            break
        left, right = sorted((int(item.i), int(item.j)))
        if left > 0 and right > left:
            probabilities[(left, right)] = float(item.p)

    mfe_pairs = dotbracket_pairs(mfe_structure)
    predicted_pairs = set(probabilities) | mfe_pairs
    column_to_rf = rf_model_index(rf)

    # RNAalifold is pseudoknot-free, so compare against all WUSS pairs that do not
    # cross any other SS_cons pair, irrespective of the bracket alphabet used.
    wuss_map = _make_pairs(ss_cons)
    wuss_pairs = {
        (min(left, right), max(left, right))
        for left, right in wuss_map.items() if left < right
    }
    def crosses(first, second):
        a, b = first
        c, d = second
        return (a < c < b < d) or (c < a < d < b)
    crossing_pairs = {
        pair for pair in wuss_pairs
        if any(crosses(pair, other) for other in wuss_pairs if other != pair)
    }
    noncrossing_pairs = wuss_pairs - crossing_pairs
    ss_noncrossing = set()
    ss_crossing = set()
    for left0, right0 in noncrossing_pairs:
        if column_to_rf[left0] and column_to_rf[right0]:
            ss_noncrossing.add(tuple(sorted((column_to_rf[left0], column_to_rf[right0]))))
    for left0, right0 in crossing_pairs:
        if column_to_rf[left0] and column_to_rf[right0]:
            ss_crossing.add(tuple(sorted((column_to_rf[left0], column_to_rf[right0]))))

    rscape_pairs = {
        (left, right) for (family, left, right), row in rscape_lookup.items()
        if family == accession and row["pair_class"] != "ANNOTATED_PAIR_NOT_DETECTED"
    }
    rscape_class_lookup = {
        (left, right): row["pair_class"]
        for (family, left, right), row in rscape_lookup.items()
        if family == accession
    }
    fr3d_lookup = {
        (left, right): row["fr3d_structural_evidence"]
        for (family, left, right), row in rscape_lookup.items()
        if family == accession
    }

    pair_rows = []
    for left, right in sorted(predicted_pairs):
        left_rf = column_to_rf[left - 1]
        right_rf = column_to_rf[right - 1]
        rf_pair = (
            tuple(sorted((left_rf, right_rf)))
            if left_rf is not None and right_rf is not None else None
        )
        pair_meta = rscape_lookup.get((accession, *rf_pair)) if rf_pair else None
        probability = probabilities.get((left, right), 0.0)
        pair_rows.append({
            "family_accession": accession,
            "family": family_name,
            "subset": subset,
            "n_sequences": len(records),
            "alignment_column_left": left,
            "alignment_column_right": right,
            "rf_model_position_left": left_rf if left_rf is not None else "",
            "rf_model_position_right": right_rf if right_rf is not None else "",
            "base_pair_probability": probability,
            "in_mfe_structure": (left, right) in mfe_pairs,
            "in_rfam_ss_cons_noncrossing_pairs": bool(rf_pair and rf_pair in ss_noncrossing),
            "in_rfam_ss_cons_crossing_pairs": bool(rf_pair and rf_pair in ss_crossing),
            "rscape_pair_class": pair_meta["pair_class"] if pair_meta else "",
            "rscape_detected_or_candidate": bool(rf_pair and rf_pair in rscape_pairs),
            "fr3d_9ztl_evidence": pair_meta["fr3d_structural_evidence"] if pair_meta else "",
            "fr3d_direct_basepair_annotated": (
                pair_meta["fr3d_basepair_observed"].lower() == "true"
                if pair_meta else False
            ),
        })

    single_energies = []
    for record in records:
        sequence = str(record.seq).upper().replace("T", "U")
        sequence = sequence.replace("-", "").replace(".", "")
        if len(sequence) < 20 or set(sequence) - set("ACGU"):
            continue
        _, energy = RNA.fold(sequence)
        single_energies.append(float(energy))
    mean_single = statistics.mean(single_energies) if single_energies else None
    sci = float(mfe_energy) / mean_single if mean_single not in (None, 0) else None

    try:
        centroid_structure, centroid_distance = fc.centroid()
    except Exception:
        centroid_structure, centroid_distance = "", None
    try:
        mea_structure, mea_score = fc.MEA()
    except Exception:
        mea_structure, mea_score = "", None
    try:
        mfe_probability = float(fc.pr_structure(mfe_structure))
    except Exception:
        mfe_probability = None
    try:
        ensemble_diversity = float(fc.mean_bp_distance())
    except Exception:
        ensemble_diversity = None

    summary = {
        "family_accession": accession,
        "family": family_name,
        "subset": subset,
        "n_sequences": len(records),
        "alignment_length": len(aligned[0]),
        "n_iupac_ambiguity_symbols_normalized_to_N": ambiguous_normalized,
        "rf_model_length": sum(symbol != "." for symbol in rf),
        "mfe_energy_kcal_mol": float(mfe_energy),
        "ensemble_energy_kcal_mol": float(ensemble_energy),
        "mean_single_sequence_mfe_kcal_mol": mean_single,
        "n_single_sequence_folds_for_sci": len(single_energies),
        "structure_conservation_index": sci,
        "mfe_pair_count": len(mfe_pairs),
        "mfe_pairs_matching_rfam_noncrossing_ss_cons": sum(
            1 for left, right in mfe_pairs
            if column_to_rf[left - 1] and column_to_rf[right - 1]
            and tuple(sorted((column_to_rf[left - 1], column_to_rf[right - 1]))) in ss_noncrossing
        ),
        "probability_pairs_at_or_above_0_05": len(probabilities),
        "probability_pairs_at_or_above_0_50": sum(p >= 0.50 for p in probabilities.values()),
        "probability_pairs_at_or_above_0_90": sum(p >= 0.90 for p in probabilities.values()),
        "mfe_structure_probability": mfe_probability,
        "ensemble_diversity_mean_bp_distance": ensemble_diversity,
        "centroid_distance": float(centroid_distance) if centroid_distance is not None else "",
        "mea_score": float(mea_score) if mea_score is not None else "",
        "mfe_structure": mfe_structure,
        "centroid_structure": centroid_structure,
        "mea_structure": mea_structure,
        "n_rfam_noncrossing_pairs": len(ss_noncrossing),
        "n_rfam_crossing_pairs_excluded_from_direct_comparison": len(ss_crossing),
    }
    return summary, pair_rows


def main() -> None:
    import RNA

    version = importlib.metadata.version("ViennaRNA")
    if version != "2.7.2":
        raise RuntimeError(f"Expected ViennaRNA 2.7.2, found {version}")

    data_dir = ROOT / "data" / "processed" / "viennarna"
    data_dir.mkdir(parents=True, exist_ok=True)
    rscape_lookup = read_rscape_pair_lookup()
    summaries, pair_rows, inputs = [], [], []

    for accession, spec in FAMILIES.items():
        seed_path = ROOT / "data" / "raw" / spec["file"]
        alignment = AlignIO.read(str(seed_path), "stockholm")
        ss_cons = read_ss_cons(seed_path)
        rf = read_rf_annotation(seed_path)
        if rf is None or len(rf) != len(ss_cons):
            raise ValueError(f"Missing or invalid RF annotation for {accession}")

        all_records = list(alignment)
        representative_sets = {
            "full": {record.id for record in all_records},
            "identity95": greedy_identity_representatives(alignment, 0.95),
            "identity90": greedy_identity_representatives(alignment, 0.90),
        }

        for subset, threshold in SUBSETS:
            selected_ids = representative_sets[subset]
            records = sorted(
                [record for record in all_records if record.id in selected_ids],
                key=lambda record: record.id,
            )
            if not records:
                raise ValueError(f"Empty representative set for {accession}/{subset}")
            stockholm_path = data_dir / f"{accession}_{subset}.sto"
            ids_path = data_dir / f"{accession}_{subset}_sequence_ids.txt"
            write_subset_stockholm(stockholm_path, records, ss_cons, rf, accession)
            ids_path.write_text("\n".join(record.id for record in records) + "\n", encoding="utf-8")
            inputs.extend([
                {"file": str(stockholm_path.relative_to(ROOT)), "sha256": sha256(stockholm_path)},
                {"file": str(ids_path.relative_to(ROOT)), "sha256": sha256(ids_path)},
            ])
            summary, pairs = run_one_alignment(
                RNA, accession, spec["name"], subset, records, ss_cons, rf, rscape_lookup
            )
            summary["viennarna_version"] = version
            summary["identity_threshold"] = threshold if threshold is not None else "full"
            summaries.append(summary)
            pair_rows.extend(pairs)
            print(
                f"{accession} {subset}: n={len(records)}, "
                f"MFE={summary['mfe_energy_kcal_mol']:.2f}, "
                f"SCI={summary['structure_conservation_index']}"
            )

    if len(summaries) != 6 or not pair_rows:
        raise RuntimeError(
            f"Expected six family/subset runs and non-empty pair output; " ,
            f"got {len(summaries)} summaries and {len(pair_rows)} pairs" ,
        )
    results_dir = ROOT / "supplementary" / "tables"
    results_dir.mkdir(parents=True, exist_ok=True)
    with (results_dir / "viennarna_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    with (results_dir / "viennarna_pair_probabilities.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)

    manifest = {
        "viennarna_version": version,
        "python_binding": "ViennaRNA PyPI 2.7.2",
        "model": "default ViennaRNA Turner parameter set at 37 C",
        "alignment_method": "RNA.fold_compound(alignment), comparative MFE and partition function",
        "constraint_policy": "unconstrained; Rfam SS_cons is used only for post hoc comparison",
        "subsets": "full seed, deterministic greedy 95% and 90% identity representatives",
        "representative_selection": "src.metrics.greedy_identity_representatives; lexicographic deterministic order; 50% minimum overlap",
        "probability_cutoff": PROBABILITY_CUTOFF,
        "ambiguity_policy": "IUPAC symbols other than A/C/G/U are converted to N only in ViennaRNA input; source Stockholm files are preserved unchanged.",
        "inputs": inputs,
        "interpretation_limits": [
            "RNAalifold predicts nested secondary structures and does not model WUSS pseudoknots.",
            "No GSC weighting is applied; full and identity-reduced subsets are reported separately.",
            "RNAalifold covariance scoring is not independent of sequence covariation and is not a substitute for R-scape.",
            "Predicted structures and SCI are computational support, not experimental validation.",
        ],
    }
    manifest_path = ROOT / "supplementary" / "provenance" / "viennarna_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"ViennaRNA {version} analysis complete.")


if __name__ == "__main__":
    main()

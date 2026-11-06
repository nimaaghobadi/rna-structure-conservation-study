"""Run QC and descriptive analyses on downloaded Rfam seed alignments."""
from __future__ import annotations
import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path
from Bio.Align import MultipleSeqAlignment
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import Bio
import numpy as np
import pandas as pd
from src.download_rfam import FAMILIES
from src.stockholm import read_stockholm, read_ss_cons, read_rf_annotation, validate_alignment, classify_structure, find_overlapping_sequence_ids
from src.metrics import column_metrics, base_pair_metrics, base_pair_family_summary, greedy_identity_representatives

ROOT = Path(__file__).resolve().parents[1]
RAW, TABLES, FIGURES = ROOT/"data"/"raw", ROOT/"results"/"tables", ROOT/"results"/"figures"
LABELS = {"RF00001": "5S rRNA", "RF00005": "tRNA"}

def run() -> dict:
    TABLES.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    qc, columns, pairs, overlap_audit, overlap_sensitivity, redundancy_sensitivity = [], [], [], [], [], []
    for accession, info in FAMILIES.items():
        path = RAW/info["file"]
        if not path.exists():
            raise FileNotFoundError(f"Missing {path}; run python -m src.download_rfam first.")
        summary = validate_alignment(path)
        alignment, structure = read_stockholm(path), read_ss_cons(path)
        reference, family = read_rf_annotation(path), LABELS[accession]
        states = classify_structure(structure, reference)
        overlap_ids = find_overlapping_sequence_ids(alignment)
        filtered_alignment = MultipleSeqAlignment(
            [record for record in alignment if record.id not in overlap_ids]
        )
        for sequence_id in sorted(overlap_ids):
            overlap_audit.append({
                "accession": accession, "family": family,
                "sequence_id": sequence_id,
                "excluded_in_overlap_sensitivity": True,
            })
        seqs = [str(rec.seq).upper() for rec in alignment]
        total = sum(map(len, seqs))
        gaps = sum(seq.count("-")+seq.count(".") for seq in seqs)
        ambiguous = sum(sum(c not in "ACGUT-." for c in seq) for seq in seqs)
        summary.update({"accession": accession, "family": family, "raw_file": str(path.relative_to(ROOT)),
            "gap_count_total": gaps, "gap_fraction_total": gaps/total if total else None,
            "ambiguous_count_total": ambiguous, "ambiguous_fraction_total": ambiguous/total if total else None,
            "n_paired_columns": states.count("paired"), "n_unpaired_columns": states.count("unpaired"),
            "n_insertion_columns": states.count("insertion"), "n_unknown_structure_columns": states.count("unknown"),
            "n_coordinate_overlap_sequences": len(overlap_ids),
            "n_sequences_after_overlap_filter": len(filtered_alignment)})
        full_column_metrics = column_metrics(alignment, states, family)
        full_pair_metrics = base_pair_metrics(alignment, structure, family, reference)
        filtered_column_metrics = column_metrics(filtered_alignment, states, family)
        filtered_pair_metrics = base_pair_metrics(filtered_alignment, structure, family, reference)
        qc.append(summary)
        columns.append(full_column_metrics)
        pairs.append(full_pair_metrics)

        def state_mean(frame, state):
            values = frame.loc[frame.structure_state == state, "conservation"].dropna()
            return float(values.mean()) if not values.empty else None

        full_paired = state_mean(full_column_metrics, "paired")
        full_unpaired = state_mean(full_column_metrics, "unpaired")
        filtered_paired = state_mean(filtered_column_metrics, "paired")
        filtered_unpaired = state_mean(filtered_column_metrics, "unpaired")
        full_pair_summary = base_pair_family_summary(full_pair_metrics)
        filtered_pair_summary = base_pair_family_summary(filtered_pair_metrics)
        full_pair_fraction = (
            float(full_pair_summary.iloc[0]["canonical_or_wobble_fraction"])
            if not full_pair_summary.empty else None
        )
        filtered_pair_fraction = (
            float(filtered_pair_summary.iloc[0]["canonical_or_wobble_fraction"])
            if not filtered_pair_summary.empty else None
        )
        overlap_sensitivity.append({
            "family": family,
            "n_sequences_full": len(alignment),
            "n_coordinate_overlap_sequences_flagged": len(overlap_ids),
            "n_sequences_after_overlap_filter": len(filtered_alignment),
            "full_paired_mean_conservation": full_paired,
            "full_unpaired_mean_conservation": full_unpaired,
            "full_paired_minus_unpaired_mean": (
                full_paired - full_unpaired
                if full_paired is not None and full_unpaired is not None else None
            ),
            "filtered_paired_mean_conservation": filtered_paired,
            "filtered_unpaired_mean_conservation": filtered_unpaired,
            "filtered_paired_minus_unpaired_mean": (
                filtered_paired - filtered_unpaired
                if filtered_paired is not None and filtered_unpaired is not None else None
            ),
            "full_canonical_or_wobble_fraction": full_pair_fraction,
            "filtered_canonical_or_wobble_fraction": filtered_pair_fraction,
        })
    for identity_threshold in (0.95, 0.90):
        for accession, info in FAMILIES.items():
            path = RAW/info["file"]
            alignment = read_stockholm(path)
            structure = read_ss_cons(path)
            reference = read_rf_annotation(path)
            family = LABELS[accession]
            states = classify_structure(structure, reference)
            representative_ids = greedy_identity_representatives(
                alignment, threshold=identity_threshold, min_overlap_fraction=0.50
            )
            filtered_alignment = MultipleSeqAlignment([
                record for record in alignment if record.id in representative_ids
            ])
            full_columns = column_metrics(alignment, states, family)
            filtered_columns = column_metrics(filtered_alignment, states, family)
            full_pairs = base_pair_metrics(alignment, structure, family, reference)
            filtered_pairs = base_pair_metrics(filtered_alignment, structure, family, reference)
            full_pair_summary = base_pair_family_summary(full_pairs)
            filtered_pair_summary = base_pair_family_summary(filtered_pairs)

            def mean_metric(frame, state, metric):
                values = frame.loc[frame.structure_state == state, metric].dropna()
                return float(values.mean()) if not values.empty else None

            full_paired = mean_metric(full_columns, "paired", "conservation")
            full_unpaired = mean_metric(full_columns, "unpaired", "conservation")
            filtered_paired = mean_metric(filtered_columns, "paired", "conservation")
            filtered_unpaired = mean_metric(filtered_columns, "unpaired", "conservation")
            full_pair_fraction = (
                float(full_pair_summary.iloc[0]["canonical_or_wobble_fraction"])
                if not full_pair_summary.empty else None
            )
            filtered_pair_fraction = (
                float(filtered_pair_summary.iloc[0]["canonical_or_wobble_fraction"])
                if not filtered_pair_summary.empty else None
            )
            redundancy_sensitivity.append({
                "family": family,
                "identity_threshold": identity_threshold,
                "minimum_overlap_fraction": 0.50,
                "n_sequences_full": len(alignment),
                "n_sequences_retained": len(filtered_alignment),
                "n_sequences_removed": len(alignment) - len(filtered_alignment),
                "full_paired_mean_conservation": full_paired,
                "full_unpaired_mean_conservation": full_unpaired,
                "full_paired_minus_unpaired_mean": (
                    full_paired - full_unpaired
                    if full_paired is not None and full_unpaired is not None else None
                ),
                "filtered_paired_mean_conservation": filtered_paired,
                "filtered_unpaired_mean_conservation": filtered_unpaired,
                "filtered_paired_minus_unpaired_mean": (
                    filtered_paired - filtered_unpaired
                    if filtered_paired is not None and filtered_unpaired is not None else None
                ),
                "full_paired_mean_entropy_bits": mean_metric(full_columns, "paired", "entropy_bits"),
                "full_unpaired_mean_entropy_bits": mean_metric(full_columns, "unpaired", "entropy_bits"),
                "filtered_paired_mean_entropy_bits": mean_metric(filtered_columns, "paired", "entropy_bits"),
                "filtered_unpaired_mean_entropy_bits": mean_metric(filtered_columns, "unpaired", "entropy_bits"),
                "full_canonical_or_wobble_fraction": full_pair_fraction,
                "filtered_canonical_or_wobble_fraction": filtered_pair_fraction,
            })

    qc_df, col_df, pair_df = pd.DataFrame(qc), pd.concat(columns, ignore_index=True), pd.concat(pairs, ignore_index=True)
    pair_summary_df = base_pair_family_summary(pair_df)
    qc_df.to_csv(TABLES/"qc_summary.csv", index=False)
    col_df.to_csv(TABLES/"per_column_metrics.csv", index=False)
    pair_df.to_csv(TABLES/"base_pair_metrics.csv", index=False)
    pair_summary_df.to_csv(TABLES/"base_pair_summary.csv", index=False)
    pd.DataFrame(overlap_audit, columns=[
        "accession", "family", "sequence_id", "excluded_in_overlap_sensitivity"
    ]).to_csv(TABLES/"coordinate_overlap_audit.csv", index=False)
    pd.DataFrame(overlap_sensitivity).to_csv(
        TABLES/"coordinate_overlap_sensitivity.csv", index=False
    )
    pd.DataFrame(redundancy_sensitivity).to_csv(
        TABLES/"sequence_redundancy_sensitivity.csv", index=False
    )
    grouped = col_df[col_df.structure_state.isin(["paired","unpaired"])].groupby(["family","structure_state"])["conservation"].agg(["count","mean","median","std"]).reset_index()
    grouped.to_csv(TABLES/"conservation_by_structure.csv", index=False)
    entropy_grouped = col_df[
        col_df.structure_state.isin(["paired", "unpaired"])
    ].groupby(["family", "structure_state"])["entropy_bits"].agg(
        ["count", "mean", "median", "std"]
    ).reset_index()
    entropy_grouped.to_csv(TABLES/"entropy_by_structure.csv", index=False)
    sensitivity = []
    for threshold in (0.50, 0.70, 0.90):
        eligible = col_df[
            col_df.structure_state.isin(["paired", "unpaired"])
            & col_df.conservation.notna()
            & (col_df.valid_fraction >= threshold)
        ]
        for family in LABELS.values():
            sub = eligible[eligible.family == family]
            paired = sub.loc[sub.structure_state == "paired", "conservation"]
            unpaired = sub.loc[sub.structure_state == "unpaired", "conservation"]
            sensitivity.append({
                "family": family, "minimum_valid_fraction": threshold,
                "n_paired_columns": int(paired.count()), "n_unpaired_columns": int(unpaired.count()),
                "paired_mean_conservation": float(paired.mean()) if paired.count() else None,
                "unpaired_mean_conservation": float(unpaired.mean()) if unpaired.count() else None,
                "paired_minus_unpaired_mean": float(paired.mean()-unpaired.mean()) if paired.count() and unpaired.count() else None
            })
    pd.DataFrame(sensitivity).to_csv(TABLES/"sensitivity_by_validity.csv", index=False)
    plot_df = col_df[col_df.structure_state.isin(["paired","unpaired"]) & col_df.conservation.notna()]
    groups, labels = [], []
    for family in LABELS.values():
        for state in ("paired","unpaired"):
            vals = plot_df[(plot_df.family==family)&(plot_df.structure_state==state)].conservation.tolist()
            groups.append(vals); labels.append(f"{family}\n{state} (n={len(vals)})")
    fig, ax = plt.subplots(figsize=(9,5)); ax.boxplot(groups, showfliers=False); ax.set_xticks(range(1, len(labels) + 1)); ax.set_xticklabels(labels)
    ax.tick_params(axis="x", labelsize=9)
    ax.set_ylabel("Per-column conservation (maximum valid-base frequency)")
    ax.set_title("Seed-alignment conservation by consensus structure state")
    fig.tight_layout(); fig.savefig(FIGURES/"conservation_by_structure.png", dpi=180); plt.close(fig)

    entropy_df = col_df[
        col_df.structure_state.isin(["paired", "unpaired"]) & col_df.entropy_bits.notna()
    ]
    entropy_groups, entropy_labels = [], []
    for family in LABELS.values():
        for state in ("paired", "unpaired"):
            vals = entropy_df[
                (entropy_df.family == family) & (entropy_df.structure_state == state)
            ].entropy_bits.tolist()
            entropy_groups.append(vals)
            entropy_labels.append(f"{family}\n{state} (n={len(vals)})")
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.boxplot(entropy_groups, showfliers=False); ax.set_xticks(range(1, len(entropy_labels) + 1)); ax.set_xticklabels(entropy_labels)
    ax.tick_params(axis="x", labelsize=9)
    ax.set_ylabel("Per-column Shannon entropy (bits; valid A/C/G/U only)")
    ax.set_title("Seed-alignment entropy by consensus structure state")
    fig.tight_layout()
    fig.savefig(FIGURES/"entropy_by_structure.png", dpi=180)
    plt.close(fig)

    if not pair_summary_df.empty:
        pair_types = ("AU", "UA", "GC", "CG", "GU", "UG", "Other")
        positions = list(range(len(pair_types)))
        width = 0.36
        fig, ax = plt.subplots(figsize=(9, 4.8))
        for family, offset in (("5S rRNA", -width / 2), ("tRNA", width / 2)):
            row = pair_summary_df.loc[pair_summary_df.family == family].iloc[0]
            values = [float(row[f"fraction_{pair}"]) for pair in pair_types[:-1]]
            values.append(float(row["noncanonical_or_other_fraction"]))
            ax.bar([position + offset for position in positions], values, width, label=family)
        ax.set_xticks(positions)
        ax.set_xticklabels(pair_types)
        ax.tick_params(axis="x", labelsize=9)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Fraction of valid paired-base observations")
        ax.set_title("Observed base-pair composition by RNA family")
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIGURES/"base_pair_composition_by_family.png", dpi=180)
        plt.close(fig)

    manifest = {
        "analysis_status": "completed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_revision": os.environ.get("GITHUB_SHA"),
        "software": {
            "python": platform.python_version(),
            "biopython": Bio.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "note": "Descriptive summaries of exact local seed alignments; see download_manifest.json for provenance.",
        "outputs": [
            "results/tables/qc_summary.csv",
            "results/tables/per_column_metrics.csv",
            "results/tables/conservation_by_structure.csv",
            "results/tables/entropy_by_structure.csv",
            "results/tables/sensitivity_by_validity.csv",
            "results/tables/base_pair_metrics.csv",
            "results/tables/base_pair_summary.csv",
            "results/tables/coordinate_overlap_audit.csv",
            "results/tables/coordinate_overlap_sensitivity.csv",
            "results/tables/sequence_redundancy_sensitivity.csv",
            "results/figures/conservation_by_structure.png",
            "results/figures/entropy_by_structure.png",
            "results/figures/base_pair_composition_by_family.png",
        ],
        "families": qc,
    }
    (ROOT/"data"/"processed"/"analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return {
        "qc": qc, "n_columns": len(col_df), "n_pair_columns": len(pair_df),
        "n_coordinate_overlap_sequences": len(overlap_audit),
        "n_redundancy_sensitivity_rows": len(redundancy_sensitivity),
    }

if __name__ == "__main__":
    print(json.dumps(run(), indent=2, ensure_ascii=False))

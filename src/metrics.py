"""Conservation, entropy, and base-pair summaries for aligned RNA families."""
from __future__ import annotations
from collections import Counter
import math
import numpy as np
import pandas as pd
from src.stockholm import _make_pairs, validate_wuss

VALID = set("ACGU")
CANONICAL = {"AU", "UA", "GC", "CG", "GU", "UG"}

def column_metrics(alignment, structure_classes: list[str], family: str) -> pd.DataFrame:
    if alignment.get_alignment_length() != len(structure_classes):
        raise ValueError("Structure classes must match alignment length")
    rows = []
    for j, state in enumerate(structure_classes):
        chars = [str(rec.seq[j]).upper().replace("T", "U") for rec in alignment]
        valid = [c for c in chars if c in VALID]
        counts = Counter(valid)
        n_valid = len(valid)
        conservation = max(counts.values()) / n_valid if n_valid else float("nan")
        entropy = -sum((n/n_valid)*math.log2(n/n_valid) for n in counts.values()) if n_valid else float("nan")
        gaps = sum(c in {"-", "."} for c in chars)
        ambiguous = sum(c not in VALID and c not in {"-", "."} for c in chars)
        rows.append({"family": family, "column_1based": j+1, "structure_state": state,
            "conservation": conservation, "entropy_bits": entropy, "n_sequences": len(chars),
            "n_valid_bases": n_valid, "valid_fraction": n_valid/len(chars) if chars else float("nan"),
            "gap_count": gaps, "gap_fraction": gaps/len(chars) if chars else float("nan"),
            "ambiguous_count": ambiguous, "ambiguous_fraction": ambiguous/len(chars) if chars else float("nan")})
    return pd.DataFrame(rows)

def base_pair_metrics(
    alignment, structure: str, family: str, reference: str | None = None
) -> pd.DataFrame:
    if len(structure) != alignment.get_alignment_length():
        raise ValueError("SS_cons must match alignment length")
    if reference is not None and len(reference) != len(structure):
        raise ValueError("RF reference annotation must match SS_cons length")
    validate_wuss(structure)
    pairs = _make_pairs(structure)
    seen, rows = set(), []
    for i, j in sorted(pairs.items()):
        if i in seen:
            continue
        seen.update((i, j))
        if reference is not None and (reference[i] == "." or reference[j] == "."):
            continue
        counts, n_both_valid = Counter(), 0
        for record in alignment:
            a = str(record.seq[i]).upper().replace("T", "U")
            b = str(record.seq[j]).upper().replace("T", "U")
            if a in VALID and b in VALID:
                n_both_valid += 1
                counts[a+b] += 1
        rows.append({"family": family, "column_i_1based": i+1, "column_j_1based": j+1,
            "n_sequences": len(alignment), "n_both_valid": n_both_valid,
            "valid_pair_fraction": n_both_valid/len(alignment) if len(alignment) else float("nan"),
            **{f"count_{p}": counts[p] for p in sorted(CANONICAL)},
            "count_noncanonical_or_other": n_both_valid - sum(counts[p] for p in CANONICAL),
            "canonical_or_wobble_fraction": sum(counts[p] for p in CANONICAL)/n_both_valid if n_both_valid else float("nan"),
            "canonical_or_wobble_fraction_all_sequences": sum(counts[p] for p in CANONICAL)/len(alignment) if len(alignment) else float("nan")})
    return pd.DataFrame(rows)


def base_pair_family_summary(pair_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate paired-base observations by family using explicit denominators."""
    pair_types = ("AU", "UA", "GC", "CG", "GU", "UG")
    columns = [
        "family", "n_annotated_base_pairs", "n_sequences",
        "n_valid_pair_observations", "n_possible_pair_observations",
        "valid_pair_observation_fraction",
        *[f"count_{pair}" for pair in pair_types],
        "count_noncanonical_or_other",
        *[f"fraction_{pair}" for pair in pair_types],
        "canonical_or_wobble_fraction", "canonical_or_wobble_fraction_all_possible",
        "noncanonical_or_other_fraction",
    ]
    if pair_df.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    for family, group in pair_df.groupby("family", sort=False):
        n_pairs = len(group)
        n_sequences = int(group["n_sequences"].iloc[0])
        n_valid = int(group["n_both_valid"].sum())
        n_possible = int(group["n_sequences"].sum())
        counts = {pair: int(group[f"count_{pair}"].sum()) for pair in pair_types}
        n_other = int(group["count_noncanonical_or_other"].sum())
        n_canonical = sum(counts.values())
        rows.append({
            "family": family,
            "n_annotated_base_pairs": n_pairs,
            "n_sequences": n_sequences,
            "n_valid_pair_observations": n_valid,
            "n_possible_pair_observations": n_possible,
            "valid_pair_observation_fraction": n_valid / n_possible if n_possible else float("nan"),
            **{f"count_{pair}": counts[pair] for pair in pair_types},
            "count_noncanonical_or_other": n_other,
            **{f"fraction_{pair}": counts[pair] / n_valid if n_valid else float("nan")
               for pair in pair_types},
            "canonical_or_wobble_fraction": n_canonical / n_valid if n_valid else float("nan"),
            "canonical_or_wobble_fraction_all_possible": n_canonical / n_possible if n_possible else float("nan"),
            "noncanonical_or_other_fraction": n_other / n_valid if n_valid else float("nan"),
        })
    return pd.DataFrame(rows, columns=columns)


def greedy_identity_representatives(
    alignment, threshold: float = 0.90, min_overlap_fraction: float = 0.50
) -> set[str]:
    """Choose deterministic greedy representatives below a pairwise identity threshold.

    Sequence identity is matches divided by columns where both sequences have
    valid A/C/G/U bases. A comparison is eligible only when those shared valid
    columns cover at least min_overlap_fraction of the shorter valid sequence.
    IDs are processed lexicographically to make the representative choice
    reproducible. This is a sensitivity filter, not a phylogenetic correction.
    """
    if not 0 < threshold <= 1:
        raise ValueError("Identity threshold must be in (0, 1]")
    if not 0 < min_overlap_fraction <= 1:
        raise ValueError("Minimum overlap fraction must be in (0, 1]")
    records = sorted(alignment, key=lambda record: record.id)
    if not records:
        return set()

    matrix = np.asarray([
        list(str(record.seq).upper().replace("T", "U"))
        for record in records
    ], dtype="U1")
    valid = np.isin(matrix, np.asarray(sorted(VALID)))
    valid_counts = valid.sum(axis=1)
    kept_indices = []

    for index in range(len(records)):
        if not kept_indices:
            kept_indices.append(index)
            continue
        kept = np.asarray(kept_indices, dtype=int)
        comparable = valid[kept] & valid[index]
        denominators = comparable.sum(axis=1)
        shorter_lengths = np.minimum(valid_counts[kept], valid_counts[index])
        overlap_fraction = np.divide(
            denominators, shorter_lengths,
            out=np.zeros(denominators.shape, dtype=float),
            where=shorter_lengths > 0,
        )
        matches = ((matrix[kept] == matrix[index]) & comparable).sum(axis=1)
        identities = np.divide(
            matches, denominators,
            out=np.zeros(denominators.shape, dtype=float),
            where=denominators > 0,
        )
        is_redundant = np.any(
            (denominators > 0)
            & (overlap_fraction >= min_overlap_fraction)
            & (identities >= threshold)
        )
        if not is_redundant:
            kept_indices.append(index)

    return {records[index].id for index in kept_indices}

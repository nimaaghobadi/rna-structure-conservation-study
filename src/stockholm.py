"""Helpers for validating Stockholm alignments and WUSS SS_cons."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re

from Bio import AlignIO

WUSS_BRACKETS = {"(": ")", "<": ">", "[": "]", "{": "}"}
WUSS_CLOSERS = {v: k for k, v in WUSS_BRACKETS.items()}
WUSS_UNPAIRED = set("._-,:~")
WUSS_LETTERS = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
WUSS_ALLOWED = WUSS_UNPAIRED | set(WUSS_BRACKETS) | set(WUSS_CLOSERS) | WUSS_LETTERS


def read_stockholm(path: str | Path):
    alignment = AlignIO.read(str(path), "stockholm")
    if len(alignment) == 0:
        raise ValueError(f"No sequences found in {path}")
    return alignment


def _read_gc_annotation(path: str | Path, annotation: str) -> str | None:
    chunks = []
    prefix = f"#=GC {annotation}"
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(prefix):
                fields = line.rstrip("\r\n").split()
                if len(fields) >= 3 and fields[1] == annotation:
                    chunks.append(fields[2])
    return "".join(chunks) if chunks else None


def find_overlapping_sequence_ids(alignment) -> set[str]:
    """Return sequence IDs flagged by a greedy inclusive coordinate-overlap audit.

    This mirrors the read-only overlap rule used by Rfam's seed-QC tool:
    sequences are grouped by accession, ordered by the lower coordinate, and
    a sequence is flagged if its coordinate interval overlaps one already kept.
    Raw alignments are never modified by this function.
    """
    accession_groups = defaultdict(list)
    for record in alignment:
        sequence_id = record.id
        if "/" not in sequence_id:
            continue
        accession, coordinates = sequence_id.split("/", 1)
        match = re.match(r"(\d+)-(\d+)", coordinates)
        if not match:
            continue
        start, end = int(match.group(1)), int(match.group(2))
        accession_groups[accession].append((sequence_id, start, end))

    flagged = set()
    for entries in accession_groups.values():
        kept = []
        for sequence_id, start, end in sorted(entries, key=lambda item: min(item[1], item[2])):
            lower, upper = min(start, end), max(start, end)
            overlaps = any(
                lower <= max(kept_start, kept_end)
                and min(kept_start, kept_end) <= upper
                for _, kept_start, kept_end in kept
            )
            if overlaps:
                flagged.add(sequence_id)
            else:
                kept.append((sequence_id, start, end))
    return flagged


def read_ss_cons(path: str | Path) -> str:
    ss = _read_gc_annotation(path, "SS_cons")
    if ss is None:
        raise ValueError(f"No #=GC SS_cons annotation found in {path}")
    return ss


def read_rf_annotation(path: str | Path) -> str | None:
    """Read the Stockholm #=GC RF reference annotation, if present."""
    return _read_gc_annotation(path, "RF")


def validate_wuss(structure: str) -> dict:
    """Raise on invalid WUSS characters or unmatched bracket/pseudoknot symbols."""
    invalid = sorted(set(structure) - WUSS_ALLOWED)
    if invalid:
        raise ValueError(f"Invalid WUSS character(s): {invalid}")

    pairs = _make_pairs(structure)
    paired_positions = {
        i for i, char in enumerate(structure)
        if char in WUSS_BRACKETS or char in WUSS_CLOSERS
        or char in WUSS_LETTERS
    }
    if len(pairs) != len(paired_positions):
        unmatched = [
            f"{i + 1}:{char}" for i, char in enumerate(structure)
            if i not in pairs and i in paired_positions
        ]
        raise ValueError("Unmatched WUSS pairing symbol(s) at " + ", ".join(unmatched))
    return {
        "length": len(structure),
        "paired_positions": len(pairs),
        "unpaired_positions": sum(char in WUSS_UNPAIRED for char in structure),
    }


def validate_alignment(path: str | Path) -> dict:
    alignment = read_stockholm(path)
    lengths = {len(record.seq) for record in alignment}
    if len(lengths) != 1:
        raise ValueError(f"Unequal sequence lengths in {path}: {sorted(lengths)}")
    ss = read_ss_cons(path)
    length = next(iter(lengths))
    if len(ss) != length:
        raise ValueError(f"SS_cons length {len(ss)} does not match alignment length {length}")
    validate_wuss(ss)
    rf = read_rf_annotation(path)
    if rf is not None and len(rf) != length:
        raise ValueError(f"RF annotation length {len(rf)} does not match alignment length {length}")
    return {
        "n_sequences": len(alignment),
        "alignment_length": length,
        "ss_cons_length": len(ss),
        "rf_annotation_present": rf is not None,
    }


def _make_pairs(structure: str) -> dict[int, int]:
    """Return paired indices for bracket and letter pseudoknot symbols."""
    stacks = defaultdict(list)
    pairs = {}
    for i, char in enumerate(structure):
        if char in WUSS_BRACKETS:
            stacks[char].append(i)
        elif char in WUSS_CLOSERS:
            opener = WUSS_CLOSERS[char]
            if stacks[opener]:
                j = stacks[opener].pop()
                pairs[i], pairs[j] = j, i
        elif char.isupper() and char.isalpha():
            stacks[char].append(i)
        elif char.islower() and char.isalpha():
            opener = char.upper()
            if stacks[opener]:
                j = stacks[opener].pop()
                pairs[i], pairs[j] = j, i
    return pairs


def classify_structure(structure: str, reference: str | None = None) -> list[str]:
    """Classify SS_cons columns; use #=GC RF to identify insertion columns.

    In SS_cons, dot is treated as an unpaired WUSS position. A dot in the
    separate RF reference annotation marks an insertion column relative to
    the family consensus. If RF is absent, insertion status cannot be inferred
    from SS_cons alone and such dots remain unpaired.
    """
    if reference is not None and len(reference) != len(structure):
        raise ValueError("RF reference annotation must match SS_cons length")
    validate_wuss(structure)
    pairs = _make_pairs(structure)
    result = []
    for i, char in enumerate(structure):
        if reference is not None and reference[i] == ".":
            result.append("insertion")
        elif i in pairs:
            result.append("paired")
        elif char in WUSS_UNPAIRED:
            result.append("unpaired")
        else:
            result.append("unknown")
    return result

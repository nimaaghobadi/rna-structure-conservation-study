"""Render project-derived RNA interaction evidence in the deposited 9ZTL ribosome.

Run in the dedicated optional workflow with:
    pymol -cq -r src/render_9ztl_rna_structure.py

The script loads the archived updated mmCIF, detects whether PyMOL exposes
author or label chain identifiers, and uses the matching residue-number column
from the validated mapping/evidence table. Highlighted residues are endpoints
of FR3D-annotated pairs; no hydrogen-bond geometry is inferred or drawn here.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path

from pymol import cmd
from PIL import Image, ImageDraw, ImageFont
import PIL

ROOT = Path.cwd().resolve()
CIF = ROOT / "data" / "structures" / "9ztl_updated.cif"
EVIDENCE = ROOT / "supplementary" / "tables" / "rscape_fr3d_structural_evidence.csv"
OUT = ROOT / "supplementary" / "figures"
MANIFEST = ROOT / "supplementary" / "provenance" / "9ztl_pymol_render_manifest.json"
CIF_SOURCE = "https://www.ebi.ac.uk/pdbe/static/entry/9ztl_updated.cif"
EXPECTED_CIF_SHA256 = "4401a250c6881fee8bc1f688d7a559b724e7ab7b18bda7d057cf0ce80f14480e"
EXPECTED_EVIDENCE_SHA256 = "5a5373fc5a4c39e185b6db71e7a639fa17b0bc926b8cacf939532f3fabb4b944"

FAMILIES = {
    "RF00001": {
        "name": "5S rRNA",
        "auth_chain": "b",
        "label_chain": "EA",
        "png": "9ztl_5s_rRNA_structural_context.png",
    },
    "RF00005": {
        "name": "P-site tRNA",
        "auth_chain": "5",
        "label_chain": "F",
        "png": "9ztl_tRNA_structural_context.png",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_evidence():
    with EVIDENCE.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def load_font(size: int, bold: bool = False):
    filename = (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
        if bold
        else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
    )
    return ImageFont.truetype(filename, size)


def add_legend(path: Path, family: str, accession: str, chain: str,
               has_stacking: bool) -> None:
    """Add a reproducible, in-panel key without altering molecular coordinates."""
    image = Image.open(path).convert("RGBA")
    overlay = Image.new("RGBA", image.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(overlay)

    x, y = 58, 58
    panel_width = 900
    row_height = 42
    entries = [
        ("forest", "Target RNA"),
        ("orange", "Direct consensus-pair endpoints"),
        ("magenta", "Additional candidates outside SS_cons"),
    ]
    if has_stacking:
        entries.append(("cyan", "Stacking-only endpoints"))

    panel_height = 112 + row_height * len(entries) + 42
    draw.rounded_rectangle(
        (x, y, x + panel_width, y + panel_height),
        radius=22,
        fill=(255, 255, 255, 232),
        outline=(95, 105, 115, 210),
        width=2,
    )
    title_font = load_font(31, bold=True)
    item_font = load_font(23)
    note_font = load_font(17)
    draw.text(
        (x + 28, y + 18),
        f"{family} — {accession} | 9ZTL chain {chain}",
        font=title_font,
        fill=(25, 35, 45, 255),
    )
    draw.text(
        (x + 28, y + 58),
        "FR3D-annotated residue endpoints",
        font=load_font(19),
        fill=(70, 80, 90, 255),
    )

    palette = {
        "forest": (34, 139, 34, 255),
        "orange": (235, 135, 25, 255),
        "magenta": (205, 45, 155, 255),
        "cyan": (0, 160, 190, 255),
    }
    for index, (color, label) in enumerate(entries):
        cy = y + 112 + index * row_height + 12
        draw.ellipse((x + 30, cy, x + 48, cy + 18), fill=palette[color])
        draw.text((x + 64, cy - 5), label, font=item_font, fill=(35, 45, 55, 255))

    note_y = y + 112 + len(entries) * row_height + 10
    draw.text(
        (x + 28, note_y),
        "Endpoint annotations only; no hydrogen bonds are inferred.",
        font=note_font,
        fill=(85, 85, 85, 255),
    )
    image = Image.alpha_composite(image, overlay).convert("RGB")
    image.save(path, format="PNG", optimize=True)


def main():
    missing = [str(path) for path in (CIF, EVIDENCE) if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            f"Missing render inputs: {missing}. Run the public 9ZTL mapping workflow first."
        )
    cif_sha = sha256(CIF)
    evidence_sha = sha256(EVIDENCE)
    if cif_sha != EXPECTED_CIF_SHA256:
        raise ValueError(f"9ZTL mmCIF checksum mismatch: {cif_sha}")
    if evidence_sha != EXPECTED_EVIDENCE_SHA256:
        raise ValueError(f"FR3D evidence table checksum mismatch: {evidence_sha}")
    OUT.mkdir(parents=True, exist_ok=True)
    cmd.reinitialize()
    cmd.set("cif_use_auth", 1)
    cmd.load(str(CIF), "9ztl")
    chains = set(cmd.get_chains("9ztl"))
    rows = read_evidence()
    manifest = {
        "pdb_id": "9ZTL",
        "cif_source": CIF_SOURCE,
        "cif_sha256": sha256(CIF),
        "evidence_table": "supplementary/tables/rscape_fr3d_structural_evidence.csv",
        "evidence_table_sha256": sha256(EVIDENCE),
        "pymol_version": cmd.get_version()[0],
        "pillow_version": PIL.__version__,
        "renderer_sha256": os.environ.get("RNA_RENDERER_SHA256") or sha256(Path(__file__).resolve()),
        "rendering": "Target RNA cartoon with translucent nearby ribosomal context, restrained C1-prime endpoint markers, and an in-panel legend",
        "rendering_settings": {
            "protein_context": "complete protein chains with atoms within 14 angstroms of the target RNA",
            "rna_context": "single strongest neighboring RNA chain by local atom count within 12 angstroms of target; local-residue fallback if none",
            "context_appearance": "muted gray cartoon with transparency",
            "camera": "oriented on target RNA; 3.5 angstrom zoom buffer",
            "endpoint_markers": "C1-prime atom spheres; original disjoint FR3D residue sets retained",
            "output_size_px": [2400, 1800],
        },
        "interpretation_limit": "Highlighted residues are endpoints of FR3D geometric annotations. The rendering does not depict inferred hydrogen bonds or establish interactions absent from FR3D output.",
        "reference_outputs": [
            {
                "family": "5S rRNA",
                "path": "supplementary/figures/9ztl_5s_rRNA_structural_context.png",
                "sha256": "9e1fca6548d16a99a94b820e967c5cc16a39badf5fc4f07fbdbd8c770a051a1e",
            },
            {
                "family": "P-site tRNA",
                "path": "supplementary/figures/9ztl_tRNA_structural_context.png",
                "sha256": "9ef90294ad9a341ef77ce1fac18393c7acc2063948342245d3b74c686ea42357",
            },
        ],
        "distribution_note": "Reference PNGs are supplied with the manuscript rather than duplicated here; regenerate them with this script if needed. Their original checksums are listed above.",
        "figures": [],
    }
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("cartoon_nucleic_acid_mode", 1)
    cmd.set("stick_radius", 0.18)
    cmd.set("label_size", 18)
    cmd.set("orthoscopic", 1)

    for accession, spec in FAMILIES.items():
        if spec["auth_chain"] in chains:
            chain = spec["auth_chain"]
            residue_field = "structure_left_auth_seq_id"
            right_field = "structure_right_auth_seq_id"
            numbering = "auth_seq_id"
        elif spec["label_chain"] in chains:
            chain = spec["label_chain"]
            residue_field = "structure_left_label_seq_id"
            right_field = "structure_right_label_seq_id"
            numbering = "label_seq_id"
        else:
            raise RuntimeError(
                f"Neither author chain {spec['auth_chain']} nor label chain "
                f"{spec['label_chain']} is present in PyMOL chains: {sorted(chains)}"
            )

        target = f"9ztl and chain {chain}"
        protein_context = f"protein_context_{accession}"
        rna_context = f"rna_context_{accession}"
        nearby_rna_atoms = f"nearby_rna_atoms_{accession}"
        cmd.hide("everything", "all")
        # Keep complete neighboring protein chains so the local ribosomal
        # environment reads as molecular context rather than isolated fragments.
        cmd.select(
            protein_context,
            f"bychain ((polymer.protein) within 14 of ({target})) and not ({target})",
        )
        # Identify local RNA residues first, then select only the strongest
        # neighboring RNA chain for a continuous context view. A transparent
        # local-residue fallback is retained if no RNA chain is within 12 Å.
        cmd.select(
            nearby_rna_atoms,
            f"byres ((polymer.nucleic) within 12 of ({target})) and not ({target})",
        )
        nearby_model = cmd.get_model(nearby_rna_atoms)
        candidate_chains = sorted({atom.chain for atom in nearby_model.atom if atom.chain != chain})
        rna_chain_atom_counts = {}
        for candidate_chain in candidate_chains:
            candidate_local = (
                f"(9ztl and chain {candidate_chain} and polymer.nucleic) "
                f"within 12 of ({target})"
            )
            local_count = cmd.count_atoms(candidate_local)
            if local_count:
                rna_chain_atom_counts[candidate_chain] = local_count
        selected_rna_chain = (
            max(rna_chain_atom_counts, key=rna_chain_atom_counts.get)
            if rna_chain_atom_counts
            else None
        )
        print(
            f"{accession}: nearby RNA-chain atom counts={rna_chain_atom_counts}; "
            f"selected context chain={selected_rna_chain!r}"
        )
        if selected_rna_chain is not None:
            cmd.select(
                rna_context,
                f"bychain (9ztl and chain {selected_rna_chain} and polymer.nucleic)",
            )
        else:
            cmd.select(rna_context, nearby_rna_atoms)
        cmd.show("cartoon", protein_context)
        cmd.color("gray85", protein_context)
        cmd.set("cartoon_transparency", 0.42, protein_context)
        cmd.show("cartoon", rna_context)
        cmd.color("gray75", rna_context)
        cmd.set("cartoon_transparency", 0.28, rna_context)
        cmd.show("cartoon", target)
        cmd.color("forest", target)

        direct_consensus = set()
        direct_additional = set()
        stacking_only = set()
        for row in rows:
            if row["family"] != accession:
                continue
            evidence = row["fr3d_structural_evidence"]
            try:
                left = int(row[residue_field])
                right = int(row[right_field])
            except (KeyError, TypeError, ValueError):
                continue
            endpoints = {left, right}
            if evidence in {
                "DIRECT_BASEPAIR_ANNOTATED",
                "DIRECT_BASEPAIR_AND_STACKING_ANNOTATED",
            }:
                if row["pair_class"] == "ADDITIONAL_COVARIATION_CANDIDATE_OUTSIDE_ORIGINAL_SS_cons":
                    direct_additional.update(endpoints)
                else:
                    direct_consensus.update(endpoints)
            elif evidence == "STACKING_ONLY_NO_DIRECT_BASEPAIR_ANNOTATION":
                stacking_only.update(endpoints)

        # Keep categories disjoint for an unambiguous visual legend.
        direct_additional -= direct_consensus
        stacking_only -= (direct_consensus | direct_additional)
        def selection(residues):
            if not residues:
                return None
            residue_list = "+".join(str(value) for value in sorted(residues))
            return f"9ztl and chain {chain} and resi {residue_list}"

        sel_consensus = selection(direct_consensus)
        sel_additional = selection(direct_additional)
        sel_stacking = selection(stacking_only)
        marker_counts = {}
        for category, residues, color in (
            ("consensus", direct_consensus, "orange"),
            ("additional", direct_additional, "magenta"),
            ("stacking", stacking_only, "cyan"),
        ):
            residue_selection = selection(residues)
            if residue_selection:
                marker_selection = f"({residue_selection}) and name C1'"
                marker_name = f"{category}_{accession}"
                cmd.select(marker_name, marker_selection)
                count = cmd.count_atoms(marker_name)
                if count == 0:
                    raise RuntimeError(
                        f"No C1-prime marker atoms found for {category} endpoints "
                        f"in {accession}; refusing to silently omit evidence residues."
                    )
                marker_counts[category] = count
                cmd.show("spheres", marker_name)
                cmd.set("sphere_scale", 0.42, marker_name)
                cmd.color(color, marker_name)
            else:
                marker_counts[category] = 0

        cmd.orient(target)
        cmd.zoom(target, 3.5)
        cmd.set("ray_trace_mode", 1)
        cmd.set("ray_trace_fog", 0)
        output = OUT / spec["png"]
        cmd.png(str(output), width=2400, height=1800, dpi=300, ray=1)
        if not output.is_file() or output.stat().st_size < 10000:
            raise RuntimeError(f"PyMOL did not create a valid rendering: {output}")
        add_legend(
            output,
            spec["name"],
            accession,
            chain,
            has_stacking=bool(stacking_only),
        )
        manifest["figures"].append({
            "accession": accession,
            "family": spec["name"],
            "chain_used": chain,
            "rna_context_chain_used": selected_rna_chain,
            "rna_context_candidate_local_atom_counts": rna_chain_atom_counts,
            "residue_numbering_used": numbering,
            "direct_consensus_pair_endpoint_residues": sorted(direct_consensus),
            "direct_additional_candidate_endpoint_residues": sorted(direct_additional),
            "stacking_only_endpoint_residues": sorted(stacking_only),
            "endpoint_marker_atom_counts": marker_counts,
            "legend_categories": [
                "target RNA",
                "direct consensus-pair endpoints",
                "additional candidates outside SS_cons",
            ] + (["stacking-only endpoints"] if stacking_only else []),
            "output": str(output.relative_to(ROOT)),
            "output_sha256": sha256(output),
        })
        for name in (
            protein_context,
            rna_context,
            nearby_rna_atoms,
            f"consensus_{accession}",
            f"additional_{accession}",
            f"stacking_{accession}",
        ):
            cmd.delete(name)

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    cmd.quit()


if __name__ == "__main__":
    main()

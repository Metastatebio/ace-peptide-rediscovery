#!/usr/bin/env python3
"""Run a small ACE docking probe for positives and hard decoys."""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
import time
from datetime import date
from pathlib import Path
from statistics import mean
from typing import Mapping

import requests
from Bio.PDB import PDBIO
import PeptideBuilder

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from backend.app.analysis.pdbqt_conversion import convert_pdb_to_pdbqt_python  # noqa: E402
from scripts.ace_rediscovery_common import PAPER_DIR, auroc_from_scores, read_csv_rows, save_json, write_csv_rows  # noqa: E402

DEFAULT_RUN_DIR = PAPER_DIR / "experiments" / "2026-08-29-expanded-literature-panel-001"
PRIMARY_ANCHORS = {"VPP", "IPP", "LKPNM", "LKP"}
RESULT_FIELDS = (
    "sequence",
    "label",
    "source_rank",
    "source_score",
    "length",
    "panel_role",
    "vina_affinity_kcal_mol",
    "vina_score_positive",
    "affinity_per_residue",
    "status",
    "seconds",
    "notes",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default=str(DEFAULT_RUN_DIR))
    parser.add_argument("--out-dir", default="")
    parser.add_argument("--pdb-id", default="1O86")
    parser.add_argument("--max-positive", type=int, default=8)
    parser.add_argument("--decoys-per-positive", type=int, default=2)
    parser.add_argument("--max-length", type=int, default=5)
    parser.add_argument("--exhaustiveness", type=int, default=1)
    parser.add_argument("--num-modes", type=int, default=1)
    parser.add_argument("--box-size", type=float, default=24.0)
    parser.add_argument("--seed", type=int, default=20260829)
    return parser.parse_args()


def clean_label(label: str) -> str:
    return label.replace("/", "_").replace(" ", "_")


def download_pdb(pdb_id: str, out_path: Path) -> str:
    if out_path.exists():
        return out_path.read_text(encoding="utf-8")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(response.text, encoding="utf-8")
    return response.text


def prepare_receptor(pdb_text: str, receptor_pdb: Path, receptor_pdbqt: Path) -> tuple[float, float, float]:
    if receptor_pdbqt.exists() and receptor_pdb.exists():
        zn = find_zinc_center(receptor_pdb.read_text(encoding="utf-8"))
        if zn is not None:
            return zn

    filtered: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            filtered.append(line)
        elif line.startswith("HETATM") and line[76:78].strip().upper() in {"ZN", "CL"}:
            filtered.append(line)
    filtered.append("END")
    receptor_pdb.parent.mkdir(parents=True, exist_ok=True)
    receptor_pdb.write_text("\n".join(filtered) + "\n", encoding="utf-8")
    zn = find_zinc_center(receptor_pdb.read_text(encoding="utf-8"))
    if zn is None:
        raise RuntimeError("could not find zinc center in ACE receptor")
    convert_pdb_to_pdbqt_python(receptor_pdb, receptor_pdbqt, mode="receptor")
    return zn


def find_zinc_center(pdb_text: str) -> tuple[float, float, float] | None:
    for line in pdb_text.splitlines():
        if line.startswith("HETATM") and line[76:78].strip().upper() == "ZN":
            return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return None


def write_peptide_pdb(sequence: str, path: Path) -> None:
    structure = PeptideBuilder.make_extended_structure(sequence)
    io = PDBIO()
    io.set_structure(structure)
    path.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(path))
    ensure_c_terminal_oxt(path)


def _xyz_from_pdb_line(line: str) -> tuple[float, float, float]:
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    norm = math.sqrt(sum(value * value for value in vector))
    if norm == 0.0:
        raise ValueError("cannot normalize zero-length vector")
    return tuple(value / norm for value in vector)


def _format_atom_line(
    *,
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    xyz: tuple[float, float, float],
    occupancy: float = 1.0,
    temp_factor: float = 0.0,
    element: str = "O",
) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:>4s} {residue_name:>3s} {chain_id:1s}"
        f"{residue_number:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        f"{occupancy:6.2f}{temp_factor:6.2f}          {element:>2s}  "
    )


def ensure_c_terminal_oxt(path: Path) -> bool:
    """Add a C-terminal OXT atom to PeptideBuilder output if it is missing."""
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if not atom_lines:
        raise ValueError(f"{path} does not contain PDB atoms")

    terminal_residue = max(int(line[22:26]) for line in atom_lines)
    terminal_atoms = {
        line[12:16].strip(): line for line in atom_lines if int(line[22:26]) == terminal_residue
    }
    if "OXT" in terminal_atoms:
        return False
    required = ("CA", "C", "O")
    if any(atom not in terminal_atoms for atom in required):
        missing = ", ".join(atom for atom in required if atom not in terminal_atoms)
        raise ValueError(f"terminal residue {terminal_residue} is missing required atoms: {missing}")

    ca = _xyz_from_pdb_line(terminal_atoms["CA"])
    carbon = _xyz_from_pdb_line(terminal_atoms["C"])
    oxygen = _xyz_from_pdb_line(terminal_atoms["O"])
    ca_direction = _unit(tuple(ca[index] - carbon[index] for index in range(3)))
    o_direction = _unit(tuple(oxygen[index] - carbon[index] for index in range(3)))
    oxt_direction = _unit(tuple(-ca_direction[index] - o_direction[index] for index in range(3)))
    oxt = tuple(carbon[index] + 1.25 * oxt_direction[index] for index in range(3))

    template = terminal_atoms["C"]
    serials = [int(line[6:11]) for line in atom_lines]
    oxt_line = _format_atom_line(
        serial=max(serials) + 1,
        atom_name="OXT",
        residue_name=template[17:20].strip(),
        chain_id=(template[21] or "A"),
        residue_number=terminal_residue,
        xyz=oxt,
    )

    insert_index = len(lines)
    for index, line in enumerate(lines):
        if line.startswith(("TER", "END")):
            insert_index = index
            break
    lines.insert(insert_index, oxt_line)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return True


def select_panel(results_path: Path, max_positive: int, decoys_per_positive: int, max_length: int) -> list[dict[str, str]]:
    rows = read_csv_rows(results_path)
    positives = [
        row
        for row in rows
        if row["label"] in {"known_positive", "digestome_known_positive"} and int(row["length"]) <= max_length
    ]
    primary = [row for row in positives if row["sequence"] in PRIMARY_ANCHORS]
    primary.sort(key=lambda row: (0 if row["sequence"] in PRIMARY_ANCHORS else 1, int(row["rank"])))
    fill = [row for row in positives if row["sequence"] not in {item["sequence"] for item in primary}]
    fill.sort(key=lambda row: int(row["rank"]))
    selected_positives = (primary + fill)[:max_positive]
    selected_positive_sequences = {row["sequence"] for row in selected_positives}
    decoys_needed = len(selected_positives) * decoys_per_positive
    decoys = [
        row
        for row in rows
        if row["label"] == "matched_decoy"
        and int(row["length"]) <= max_length
        and row["sequence"] not in selected_positive_sequences
    ]
    decoys.sort(key=lambda row: int(row["rank"]))
    selected: list[dict[str, str]] = []
    for row in selected_positives:
        selected.append({**row, "panel_role": "positive"})
    for row in decoys[:decoys_needed]:
        selected.append({**row, "panel_role": "hard_matched_decoy"})
    return selected


def parse_vina_affinity(pdbqt_path: Path) -> float | None:
    if not pdbqt_path.exists():
        return None
    for line in pdbqt_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = re.match(r"REMARK VINA RESULT:\s+(-?[0-9.]+)", line)
        if match:
            return float(match.group(1))
    return None


def dock_sequence(
    *,
    sequence: str,
    receptor_pdbqt: Path,
    center: tuple[float, float, float],
    out_dir: Path,
    exhaustiveness: int,
    num_modes: int,
    box_size: float,
    seed: int,
) -> tuple[float | None, str, float, str]:
    start = time.monotonic()
    sequence_dir = out_dir / clean_label(sequence)
    sequence_dir.mkdir(parents=True, exist_ok=True)
    ligand_pdb = sequence_dir / f"{sequence}.pdb"
    ligand_pdbqt = sequence_dir / f"{sequence}.pdbqt"
    poses = sequence_dir / f"{sequence}_vina.pdbqt"
    try:
        if not ligand_pdb.exists():
            write_peptide_pdb(sequence, ligand_pdb)
        if not ligand_pdbqt.exists():
            convert_pdb_to_pdbqt_python(ligand_pdb, ligand_pdbqt, mode="ligand")
        command = [
            "vina",
            "--receptor",
            str(receptor_pdbqt),
            "--ligand",
            str(ligand_pdbqt),
            "--out",
            str(poses),
            "--center_x",
            str(center[0]),
            "--center_y",
            str(center[1]),
            "--center_z",
            str(center[2]),
            "--size_x",
            str(box_size),
            "--size_y",
            str(box_size),
            "--size_z",
            str(box_size),
            "--exhaustiveness",
            str(exhaustiveness),
            "--num_modes",
            str(num_modes),
            "--cpu",
            "1",
            "--seed",
            str(seed),
            "--verbosity",
            "0",
        ]
        subprocess.run(command, check=True, capture_output=True, text=True)
        affinity = parse_vina_affinity(poses)
        return affinity, "ok" if affinity is not None else "missing_affinity", time.monotonic() - start, ""
    except Exception as exc:  # pragma: no cover - integration path
        return None, "failed", time.monotonic() - start, str(exc)


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def mean_or_blank(values: list[float]) -> str:
    return f"{mean(values):.3f}" if values else ""


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir)
    out_dir = Path(args.out_dir) if args.out_dir else run_dir / "docking-probe-001"
    out_dir.mkdir(parents=True, exist_ok=True)
    receptor_dir = out_dir / "receptor"
    pdb_path = receptor_dir / f"{args.pdb_id}.pdb"
    receptor_pdb = receptor_dir / f"{args.pdb_id}_ace_receptor.pdb"
    receptor_pdbqt = receptor_dir / f"{args.pdb_id}_ace_receptor.pdbqt"
    pdb_text = download_pdb(args.pdb_id, pdb_path)
    center = prepare_receptor(pdb_text, receptor_pdb, receptor_pdbqt)
    panel = select_panel(
        run_dir / "ace-rediscovery-results.csv",
        max_positive=args.max_positive,
        decoys_per_positive=args.decoys_per_positive,
        max_length=args.max_length,
    )

    result_rows: list[dict[str, object]] = []
    for index, row in enumerate(panel):
        affinity, status, seconds, notes = dock_sequence(
            sequence=row["sequence"],
            receptor_pdbqt=receptor_pdbqt,
            center=center,
            out_dir=out_dir / "ligands",
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            box_size=args.box_size,
            seed=args.seed + index,
        )
        vina_positive = -affinity if affinity is not None else ""
        per_residue = affinity / int(row["length"]) if affinity is not None else ""
        result_rows.append(
            {
                "sequence": row["sequence"],
                "label": row["label"],
                "source_rank": row["rank"],
                "source_score": row["ace_score"],
                "length": row["length"],
                "panel_role": row["panel_role"],
                "vina_affinity_kcal_mol": affinity if affinity is not None else "",
                "vina_score_positive": vina_positive,
                "affinity_per_residue": per_residue,
                "status": status,
                "seconds": round(seconds, 3),
                "notes": notes,
            }
        )

    write_csv_rows(out_dir / "ace-docking-probe-results.csv", result_rows, RESULT_FIELDS)
    ok_rows = [row for row in result_rows if row["status"] == "ok"]
    positives = [row for row in ok_rows if row["panel_role"] == "positive"]
    decoys = [row for row in ok_rows if row["panel_role"] == "hard_matched_decoy"]
    positive_scores = [float(row["vina_score_positive"]) for row in positives]
    decoy_scores = [float(row["vina_score_positive"]) for row in decoys]
    auc = auroc_from_scores(positive_scores, decoy_scores)
    table_rows = [
        [
            row["sequence"],
            row["panel_role"],
            row["source_rank"],
            row["source_score"],
            row["length"],
            row["vina_affinity_kcal_mol"],
            row["status"],
        ]
        for row in sorted(ok_rows, key=lambda item: float(item["vina_affinity_kcal_mol"]))
    ]
    manifest = {
        "date": date.today().isoformat(),
        "pdb_id": args.pdb_id,
        "pdb_url": f"https://files.rcsb.org/download/{args.pdb_id}.pdb",
        "receptor_source": "Human ACE lisinopril complex; PDB 1O86; zinc-centered box.",
        "ligand_terminal_chemistry": "PeptideBuilder extended conformer with explicit C-terminal OXT atom before PDBQT conversion.",
        "center": center,
        "box_size_angstrom": args.box_size,
        "exhaustiveness": args.exhaustiveness,
        "num_modes": args.num_modes,
        "max_positive": args.max_positive,
        "decoys_per_positive": args.decoys_per_positive,
        "max_length": args.max_length,
        "successful_docks": len(ok_rows),
        "positive_count": len(positives),
        "decoy_count": len(decoys),
        "mean_positive_affinity": mean_or_blank([float(row["vina_affinity_kcal_mol"]) for row in positives]),
        "mean_decoy_affinity": mean_or_blank([float(row["vina_affinity_kcal_mol"]) for row in decoys]),
        "auroc_positive_vs_hard_decoy": auc,
        "claim_boundary": "Post-hoc Vina docking probe. This is target-informed and label-free, but not a validated binding free-energy calculation.",
    }
    save_json(out_dir / "ace-docking-probe-manifest.json", manifest)
    report = f"""# ACE Docking Probe

Date: {date.today().isoformat()}

## Scope

This is a post-hoc molecular-simulation probe, not a trained ML model. It uses human ACE structure `{args.pdb_id}` with a zinc-centered Vina docking box and generated extended peptide conformers.

It is target-informed and label-free: no ACE peptide activity labels are used for scoring. It is not a validated binding free-energy workflow.

## Summary

- Successful docks: `{len(ok_rows)}`
- Positives docked: `{len(positives)}`
- Hard matched decoys docked: `{len(decoys)}`
- Mean positive affinity: `{manifest["mean_positive_affinity"]}` kcal/mol
- Mean hard-decoy affinity: `{manifest["mean_decoy_affinity"]}` kcal/mol
- AUROC positives versus hard decoys: `{auc:.4f}`

## Docking Results

{markdown_table(["Sequence", "Role", "Source rank", "Source score", "Length", "Vina affinity", "Status"], table_rows)}

## Claim Boundary

This probe can support a methods statement that PepLab can add a target-structure-informed physics tier. It cannot support a standalone claim that Vina proves ACE inhibition.
"""
    (out_dir / "ace-docking-probe-report.md").write_text(report, encoding="utf-8")
    print(out_dir / "ace-docking-probe-report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

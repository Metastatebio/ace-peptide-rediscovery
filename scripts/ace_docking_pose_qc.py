#!/usr/bin/env python3
"""Rescore ACE docking poses with fixed, label-free pose-quality metrics."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from statistics import mean
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_docking_probe import clean_label  # noqa: E402
from scripts.ace_rediscovery_common import PAPER_DIR, auroc_from_scores, read_csv_rows, save_json, write_csv_rows  # noqa: E402

DEFAULT_DOCKING_DIR = PAPER_DIR / "experiments" / "2026-08-29-docking-probe-001"
FLOAT_RE = re.compile(r"[-+]?\d+(?:\.\d+)?")
POSE_QC_FIELDS = (
    "sequence",
    "panel_role",
    "label",
    "source_rank",
    "length",
    "vina_affinity_kcal_mol",
    "raw_vina_score",
    "size_normalized_vina_score",
    "min_heavy_atom_to_zinc_a",
    "min_polar_heavy_atom_to_zinc_a",
    "min_oxygen_to_zinc_a",
    "min_terminal_oxygen_to_zinc_a",
    "polar_heavy_atoms_within_3_2a",
    "oxygen_atoms_within_3_2a",
    "terminal_oxygen_atoms_detected",
    "zinc_contact_score",
    "terminal_oxygen_zinc_score",
    "pose_qc_score",
    "status",
    "notes",
)


@dataclass(frozen=True)
class PdbqtAtom:
    serial: int
    name: str
    residue_name: str
    chain_id: str
    residue_number: int | None
    x: float
    y: float
    z: float
    atom_type: str
    element: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docking-dir", default=str(DEFAULT_DOCKING_DIR))
    parser.add_argument("--out-dir", default="")
    return parser.parse_args()


def safe_float(value: object) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def safe_int(value: object) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def infer_element(name: str, atom_type: str) -> str:
    atom_type = atom_type.strip().upper()
    if atom_type in {"OA", "O"}:
        return "O"
    if atom_type in {"NA", "N"}:
        return "N"
    if atom_type in {"SA", "S"}:
        return "S"
    if atom_type in {"HD", "H"}:
        return "H"
    if atom_type in {"A", "C"}:
        return "C"
    stripped_name = name.strip().upper()
    if not stripped_name:
        return ""
    if stripped_name.startswith("CL"):
        return "CL"
    return stripped_name[0]


def parse_xyz(line: str) -> tuple[float, float, float]:
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except ValueError:
        values = FLOAT_RE.findall(line)
        if len(values) < 3:
            raise
        return (float(values[-6]), float(values[-5]), float(values[-4]))


def parse_pdbqt_first_model_atoms(text: str) -> list[PdbqtAtom]:
    atoms: list[PdbqtAtom] = []
    in_first_model = False
    saw_model = False
    for line in text.splitlines():
        if line.startswith("MODEL"):
            if saw_model:
                break
            saw_model = True
            in_first_model = True
            continue
        if line.startswith("ENDMDL") and in_first_model:
            break
        if saw_model and not in_first_model:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        xyz = parse_xyz(line)
        atom_type = line[77:].strip() or line.split()[-1]
        serial = safe_int(line[6:11]) or len(atoms) + 1
        residue_number = safe_int(line[22:26])
        name = line[12:16].strip()
        atoms.append(
            PdbqtAtom(
                serial=serial,
                name=name,
                residue_name=line[17:20].strip(),
                chain_id=line[21:22].strip(),
                residue_number=residue_number,
                x=xyz[0],
                y=xyz[1],
                z=xyz[2],
                atom_type=atom_type,
                element=infer_element(name, atom_type),
            )
        )
    return atoms


def distance(atom: PdbqtAtom, center: tuple[float, float, float]) -> float:
    return math.dist((atom.x, atom.y, atom.z), center)


def minimum_distance(
    atoms: Sequence[PdbqtAtom],
    center: tuple[float, float, float],
    *,
    allowed_elements: set[str] | None = None,
    allowed_names: set[str] | None = None,
) -> float | None:
    filtered = []
    for atom in atoms:
        if allowed_elements is not None and atom.element not in allowed_elements:
            continue
        if allowed_names is not None and atom.name not in allowed_names:
            continue
        filtered.append(distance(atom, center))
    return min(filtered) if filtered else None


def triangular_score(value: float | None, target: float, width: float) -> float:
    if value is None or width <= 0:
        return 0.0
    return max(0.0, min(1.0, 1.0 - abs(value - target) / width))


def zinc_distance_score(value: float | None) -> float:
    return triangular_score(value, target=2.3, width=2.2)


def terminal_oxygen_atoms(atoms: Sequence[PdbqtAtom]) -> list[PdbqtAtom]:
    residue_numbers = [atom.residue_number for atom in atoms if atom.residue_number is not None]
    if not residue_numbers:
        return []
    c_terminal_residue = max(residue_numbers)
    return [
        atom
        for atom in atoms
        if atom.residue_number == c_terminal_residue and atom.element == "O" and atom.name in {"O", "OXT"}
    ]


def minmax_score(value: float | None, values: Sequence[float]) -> float:
    if value is None or not values:
        return 0.0
    low = min(values)
    high = max(values)
    if high == low:
        return 0.5
    return (value - low) / (high - low)


def score_pose_row(
    row: Mapping[str, str],
    atoms: Sequence[PdbqtAtom],
    center: tuple[float, float, float],
    *,
    raw_vina_values: Sequence[float],
    size_normalized_values: Sequence[float],
) -> dict[str, object]:
    affinity = safe_float(row.get("vina_affinity_kcal_mol"))
    length = safe_int(row.get("length")) or len(str(row.get("sequence") or ""))
    raw_vina = -affinity if affinity is not None else None
    size_normalized = raw_vina / length if raw_vina is not None and length else None
    heavy_atoms = [atom for atom in atoms if atom.element != "H"]
    polar_heavy_atoms = [atom for atom in heavy_atoms if atom.element in {"O", "N", "S"}]
    oxygen_atoms = [atom for atom in heavy_atoms if atom.element == "O"]
    terminal_oxygens = terminal_oxygen_atoms(heavy_atoms)
    min_heavy = minimum_distance(heavy_atoms, center)
    min_polar = minimum_distance(polar_heavy_atoms, center)
    min_oxygen = minimum_distance(oxygen_atoms, center)
    min_terminal_oxygen = minimum_distance(terminal_oxygens, center)
    zinc_contact = zinc_distance_score(min_polar)
    terminal_oxygen_score = zinc_distance_score(min_terminal_oxygen)
    raw_vina_component = minmax_score(raw_vina, raw_vina_values)
    size_vina_component = minmax_score(size_normalized, size_normalized_values)
    pose_qc = (
        0.35 * terminal_oxygen_score
        + 0.25 * zinc_contact
        + 0.25 * raw_vina_component
        + 0.15 * size_vina_component
    )
    return {
        "sequence": row.get("sequence", ""),
        "panel_role": row.get("panel_role", ""),
        "label": row.get("label", ""),
        "source_rank": row.get("source_rank", ""),
        "length": length,
        "vina_affinity_kcal_mol": affinity if affinity is not None else "",
        "raw_vina_score": raw_vina if raw_vina is not None else "",
        "size_normalized_vina_score": size_normalized if size_normalized is not None else "",
        "min_heavy_atom_to_zinc_a": min_heavy if min_heavy is not None else "",
        "min_polar_heavy_atom_to_zinc_a": min_polar if min_polar is not None else "",
        "min_oxygen_to_zinc_a": min_oxygen if min_oxygen is not None else "",
        "min_terminal_oxygen_to_zinc_a": min_terminal_oxygen if min_terminal_oxygen is not None else "",
        "polar_heavy_atoms_within_3_2a": sum(1 for atom in polar_heavy_atoms if distance(atom, center) <= 3.2),
        "oxygen_atoms_within_3_2a": sum(1 for atom in oxygen_atoms if distance(atom, center) <= 3.2),
        "terminal_oxygen_atoms_detected": len(terminal_oxygens),
        "zinc_contact_score": zinc_contact,
        "terminal_oxygen_zinc_score": terminal_oxygen_score,
        "pose_qc_score": pose_qc,
        "status": "ok" if atoms else "missing_pose_atoms",
        "notes": "",
    }


def metric_auroc(rows: Sequence[Mapping[str, object]], metric: str) -> float:
    positives = [
        float(row[metric])
        for row in rows
        if row.get("panel_role") == "positive" and row.get(metric) not in ("", None)
    ]
    decoys = [
        float(row[metric])
        for row in rows
        if row.get("panel_role") == "hard_matched_decoy" and row.get(metric) not in ("", None)
    ]
    return auroc_from_scores(positives, decoys)


def metric_mean(rows: Sequence[Mapping[str, object]], metric: str, role: str) -> str:
    values = [
        float(row[metric])
        for row in rows
        if row.get("panel_role") == role and row.get(metric) not in ("", None)
    ]
    return f"{mean(values):.4f}" if values else ""


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def load_zinc_center(docking_dir: Path) -> tuple[float, float, float]:
    manifest_path = docking_dir / "ace-docking-probe-manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    center = manifest.get("center")
    if not isinstance(center, list) or len(center) != 3:
        raise ValueError(f"missing docking center in {manifest_path}")
    return (float(center[0]), float(center[1]), float(center[2]))


def main() -> int:
    args = parse_args()
    docking_dir = Path(args.docking_dir)
    out_dir = Path(args.out_dir) if args.out_dir else docking_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = [row for row in read_csv_rows(docking_dir / "ace-docking-probe-results.csv") if row.get("status") == "ok"]
    center = load_zinc_center(docking_dir)
    raw_vina_values: list[float] = []
    size_normalized_values: list[float] = []
    for row in rows:
        affinity = safe_float(row.get("vina_affinity_kcal_mol"))
        length = safe_int(row.get("length")) or len(row.get("sequence", ""))
        if affinity is None or not length:
            continue
        raw = -affinity
        raw_vina_values.append(raw)
        size_normalized_values.append(raw / length)

    scored_rows: list[dict[str, object]] = []
    for row in rows:
        sequence = row["sequence"]
        pose_path = docking_dir / "ligands" / clean_label(sequence) / f"{sequence}_vina.pdbqt"
        atoms = parse_pdbqt_first_model_atoms(pose_path.read_text(encoding="utf-8", errors="ignore")) if pose_path.exists() else []
        scored_rows.append(
            score_pose_row(
                row,
                atoms,
                center,
                raw_vina_values=raw_vina_values,
                size_normalized_values=size_normalized_values,
            )
        )

    write_csv_rows(out_dir / "ace-docking-pose-qc-results.csv", scored_rows, POSE_QC_FIELDS)
    metrics = [
        "raw_vina_score",
        "size_normalized_vina_score",
        "zinc_contact_score",
        "terminal_oxygen_zinc_score",
        "pose_qc_score",
    ]
    summary_rows = []
    for metric in metrics:
        summary_rows.append(
            {
                "metric": metric,
                "positive_mean": metric_mean(scored_rows, metric, "positive"),
                "hard_decoy_mean": metric_mean(scored_rows, metric, "hard_matched_decoy"),
                "auroc_positive_vs_hard_decoy": metric_auroc(scored_rows, metric),
            }
        )
    save_json(
        out_dir / "ace-docking-pose-qc-summary.json",
        {
            "date": date.today().isoformat(),
            "docking_dir": str(docking_dir),
            "zinc_center": center,
            "successful_pose_qc_rows": len([row for row in scored_rows if row["status"] == "ok"]),
            "metrics": summary_rows,
            "claim_boundary": "Fixed label-free pose QC over Vina outputs. Diagnostic only; not a validated binding free-energy workflow.",
        },
    )

    metric_table = [
        [
            row["metric"],
            row["positive_mean"],
            row["hard_decoy_mean"],
            f'{float(row["auroc_positive_vs_hard_decoy"]):.4f}',
        ]
        for row in summary_rows
    ]
    ranked_table = [
        [
            row["sequence"],
            row["panel_role"],
            f'{float(row["pose_qc_score"]):.4f}',
            f'{float(row["raw_vina_score"]):.3f}',
            f'{float(row["min_polar_heavy_atom_to_zinc_a"]):.2f}' if row["min_polar_heavy_atom_to_zinc_a"] != "" else "",
            f'{float(row["min_terminal_oxygen_to_zinc_a"]):.2f}' if row["min_terminal_oxygen_to_zinc_a"] != "" else "",
        ]
        for row in sorted(scored_rows, key=lambda item: -float(item["pose_qc_score"]))
    ]
    report = f"""# ACE Docking Pose QC

Date: {date.today().isoformat()}

## Scope

This is a fixed, label-free pose-quality pass over the existing Vina outputs. It tests whether zinc-proximal pose chemistry adds discriminative signal beyond raw docking affinity. It is diagnostic and should not be read as validated binding free energy.

## Metric Summary

{markdown_table(["Metric", "Positive mean", "Hard-decoy mean", "AUROC"], metric_table)}

## Pose-QC Ranking

{markdown_table(["Sequence", "Role", "Pose-QC", "Raw Vina score", "Min polar-Zn A", "Min terminal O-Zn A"], ranked_table)}

## Claim Boundary

Use this as a molecular-simulation audit. If the hard-decoy AUROC remains weak, the current docking tier is not publication-grade evidence and should be upgraded before it appears as a main result.
"""
    (out_dir / "ace-docking-pose-qc-report.md").write_text(report, encoding="utf-8")
    print(out_dir / "ace-docking-pose-qc-report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

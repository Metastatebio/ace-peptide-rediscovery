#!/usr/bin/env python3
"""Build a cross-checkpoint consensus table for ACE PLM masked variants."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from datetime import date
from pathlib import Path
from statistics import mean
from typing import Any, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO_ROOT / "docs" / "publication-engine" / "papers" / "01-ace-rediscovery"
DEFAULT_OUT_DIR = PAPER_DIR / "experiments" / "2026-08-30-ace-plm-variant-consensus-001"
DEFAULT_VARIANTS = (
    (
        "esm2_8m",
        PAPER_DIR
        / "experiments"
        / "2026-08-30-ace-protein-lm-fusion-001"
        / "ace_plm_masked_variant_candidates.csv",
    ),
    (
        "esm2_35m",
        PAPER_DIR
        / "experiments"
        / "2026-08-30-ace-protein-lm-fusion-002-esm35m"
        / "ace_plm_masked_variant_candidates.csv",
    ),
    (
        "esm2_150m",
        PAPER_DIR
        / "experiments"
        / "2026-08-30-ace-protein-lm-fusion-003-esm150m"
        / "ace_plm_masked_variant_candidates.csv",
    ),
)

CONSENSUS_FIELDS = (
    "consensus_rank",
    "sequence",
    "models_present",
    "model_ids",
    "mean_model_rank",
    "best_model_rank",
    "mean_plm_mean_log_likelihood",
    "mean_delta_parent_plm_mean",
    "parent_sequences",
    "edit_counts",
    "mutations_by_model",
    "nearest_known_positive",
    "nearest_known_positive_distance",
    "developability_flags",
    "mol_wt",
    "mol_logp",
    "tpsa",
    "rotatable_bonds",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument(
        "--variant",
        action="append",
        default=[],
        help="Variant table in the form label=/path/to/ace_plm_masked_variant_candidates.csv.",
    )
    parser.add_argument("--top-n", type=int, default=80)
    return parser.parse_args()


def clean_id(text: str) -> str:
    return "".join(char if char.isalnum() or char in {"_", "-", "."} else "_" for char in text)


def variant_inputs(raw_variants: Sequence[str]) -> list[tuple[str, Path]]:
    if not raw_variants:
        return [(label, path) for label, path in DEFAULT_VARIANTS]
    parsed: list[tuple[str, Path]] = []
    for item in raw_variants:
        if "=" not in item:
            raise ValueError(f"variant input must use label=path format: {item}")
        label, path = item.split("=", 1)
        parsed.append((clean_id(label.strip()), Path(path)))
    return parsed


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(payload), indent=2, sort_keys=True), encoding="utf-8")


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return default if math.isnan(number) else number


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def first_nonempty(values: Sequence[Any]) -> Any:
    for value in values:
        if str(value) != "":
            return value
    return ""


def build_consensus(inputs: Sequence[tuple[str, Path]], top_n: int) -> list[dict[str, Any]]:
    rows_by_sequence: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for model_id, path in inputs:
        for row in read_csv_rows(path):
            sequence = row.get("sequence", "")
            if not sequence:
                continue
            rows_by_sequence[sequence].append({"model_id": model_id, **row})

    consensus_rows = []
    for sequence, model_rows in rows_by_sequence.items():
        model_ids = sorted({str(row["model_id"]) for row in model_rows})
        ranks = [safe_int(row.get("rank")) for row in model_rows]
        plm_scores = [safe_float(row.get("plm_mean_log_likelihood")) for row in model_rows]
        deltas = [safe_float(row.get("delta_parent_plm_mean")) for row in model_rows]
        nearest_distances = [safe_int(row.get("nearest_known_positive_distance"), 999999) for row in model_rows]
        row = {
            "sequence": sequence,
            "models_present": len(model_ids),
            "model_ids": ";".join(model_ids),
            "mean_model_rank": round(mean(ranks), 3),
            "best_model_rank": min(ranks),
            "mean_plm_mean_log_likelihood": round(mean(plm_scores), 8),
            "mean_delta_parent_plm_mean": round(mean(deltas), 8),
            "parent_sequences": ";".join(sorted({str(item.get("parent_sequence", "")) for item in model_rows})),
            "edit_counts": ";".join(sorted({str(item.get("edit_count", "")) for item in model_rows})),
            "mutations_by_model": ";".join(
                f"{item['model_id']}:{item.get('mutations', '')}" for item in sorted(model_rows, key=lambda r: r["model_id"])
            ),
            "nearest_known_positive": first_nonempty(
                [
                    item.get("nearest_known_positive", "")
                    for item in sorted(
                        model_rows,
                        key=lambda r: safe_int(r.get("nearest_known_positive_distance"), 999999),
                    )
                ]
            ),
            "nearest_known_positive_distance": min(nearest_distances),
            "developability_flags": ";".join(sorted({str(item.get("developability_flag", "")) for item in model_rows})),
            "mol_wt": first_nonempty([item.get("mol_wt", "") for item in model_rows]),
            "mol_logp": first_nonempty([item.get("mol_logp", "") for item in model_rows]),
            "tpsa": first_nonempty([item.get("tpsa", "") for item in model_rows]),
            "rotatable_bonds": first_nonempty([item.get("rotatable_bonds", "") for item in model_rows]),
        }
        consensus_rows.append(row)

    consensus_rows.sort(
        key=lambda row: (
            -safe_int(row["models_present"]),
            safe_float(row["mean_model_rank"]),
            -safe_float(row["mean_delta_parent_plm_mean"]),
            str(row["sequence"]),
        )
    )
    for rank, row in enumerate(consensus_rows, start=1):
        row["consensus_rank"] = rank
    return consensus_rows[:top_n]


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def report_markdown(
    *,
    inputs: Sequence[tuple[str, Path]],
    rows: Sequence[Mapping[str, Any]],
    out_dir: Path,
) -> str:
    top_rows = [
        [
            row["consensus_rank"],
            row["sequence"],
            row["models_present"],
            row["model_ids"],
            row["mean_model_rank"],
            row["parent_sequences"],
            row["nearest_known_positive"],
            row["nearest_known_positive_distance"],
            row["developability_flags"],
        ]
        for row in rows[:20]
    ]
    return f"""# ACE PLM Variant Consensus

Date: {date.today().isoformat()}

## Inputs

{markdown_table(["Model", "Variant table"], [[label, path] for label, path in inputs])}

## Consensus Candidates

Candidates are ranked first by checkpoint agreement, then by mean rank within the model-specific variant tables. These are not validated inhibitors; they are the strongest current computational candidates for a wet-lab or receptor-conditioned follow-up batch.

{markdown_table(["Rank", "Sequence", "Models", "Model IDs", "Mean model rank", "Parents", "Nearest known positive", "Distance", "Developability"], top_rows)}

## Output

- Consensus CSV: `{out_dir / "ace_plm_consensus_variant_candidates.csv"}`
- Manifest: `{out_dir / "ace_plm_variant_consensus_manifest.json"}`
"""


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = variant_inputs(args.variant)
    rows = build_consensus(inputs, args.top_n)
    write_csv(out_dir / "ace_plm_consensus_variant_candidates.csv", rows, CONSENSUS_FIELDS)
    write_json(
        out_dir / "ace_plm_variant_consensus_manifest.json",
        {
            "date": date.today().isoformat(),
            "inputs": [{"model_id": label, "path": str(path)} for label, path in inputs],
            "rows": len(rows),
            "ranking": "checkpoint agreement desc, mean model rank asc, mean parent delta desc",
            "claim_boundary": "Computational candidate prioritization only; no wet-lab potency or receptor-contact claim.",
        },
    )
    (out_dir / "ace_plm_variant_consensus_report.md").write_text(
        report_markdown(inputs=inputs, rows=rows, out_dir=out_dir),
        encoding="utf-8",
    )
    print(out_dir / "ace_plm_variant_consensus_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

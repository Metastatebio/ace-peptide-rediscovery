#!/usr/bin/env python3
"""Apply a frozen ACE two-stage rerank recipe to a scored candidate universe."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rank_fusion_search import (  # noqa: E402
    BASE_FIELDS,
    DEFAULT_BASELINE_SWEEP_DIR,
    DEFAULT_MODEL_DIR,
    POSITIVE_LABELS,
    combined_score,
    length_stratified_components,
    markdown_table,
    normalize_components,
    ranked_metrics,
    read_baseline_scores,
    read_wide_scores,
    short_float,
)
from scripts.ace_two_stage_rerank_search import staged_scores  # noqa: E402
from scripts.ace_rediscovery_common import PAPER_DIR, read_csv_rows, save_json, sha256_text, write_csv_rows  # noqa: E402

DEFAULT_OUT_DIR = PAPER_DIR / "experiments" / "2026-08-29-v4-independent-validation-001" / "frozen-rerank"

RANKED_FIELDS = (
    "candidate_id",
    "rank",
    "sequence",
    "ace_score",
    "length",
    "label",
    "origin",
    "provenance",
    "gate_component",
    "rank_component",
    "gate_percent",
    "gated",
)
SUMMARY_FIELDS = (
    "recipe_id",
    "gate_component",
    "rank_component",
    "gate_percent",
    "gate_count",
    "top_1pct_count",
    "top_1pct_fraction",
    "top_1pct_enrichment_over_random",
    "top_1pct_hypergeom_p",
    "top_10pct_count",
    "top_10pct_fraction",
    "top_100_count",
    "mean_positive_rank",
    "median_positive_rank",
    "best_positive_rank",
    "auroc_vs_all_negatives",
    "auroc_vs_matched_decoys",
    "average_precision_vs_all_negatives",
)
POSITIVE_FIELDS = ("recipe_id", "rank", "sequence", "fusion_score", "label", "origin", "provenance")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", required=True, help="Frozen recipe JSON path.")
    parser.add_argument(
        "--model-dir",
        action="append",
        default=None,
        help="Model score directory containing eval-model-scores.csv. Can be supplied more than once.",
    )
    parser.add_argument("--baseline-sweep-dir", default=str(DEFAULT_BASELINE_SWEEP_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    return parser.parse_args()


def load_recipe(path: Path) -> dict[str, object]:
    recipe = json.loads(path.read_text(encoding="utf-8"))
    required = {"recipe_id", "gate_component", "rank_component", "gate_percent"}
    missing = sorted(required - set(recipe))
    if missing:
        raise ValueError(f"recipe is missing required keys: {missing}")
    return recipe


def add_recipe_fusions(
    *,
    normalized: Mapping[str, np.ndarray],
    fusion_components: Sequence[Mapping[str, object]],
) -> dict[str, np.ndarray]:
    components = dict(normalized)
    for fusion in fusion_components:
        component_id = str(fusion["component_id"])
        weights_raw = fusion["weights"]
        if not isinstance(weights_raw, dict):
            raise ValueError(f"fusion component {component_id} has non-dict weights")
        weights = {str(component): float(weight) for component, weight in weights_raw.items()}
        missing = sorted(set(weights) - set(normalized))
        if missing:
            raise ValueError(f"fusion component {component_id} references missing components: {missing[:8]}")
        components[component_id] = combined_score(weights, normalized)
    return normalize_components(components)


def ranked_rows(
    *,
    rows: Sequence[Mapping[str, str]],
    scores: np.ndarray,
    gate_indices: np.ndarray,
    gate_component: str,
    rank_component: str,
    gate_percent: float,
) -> list[dict[str, object]]:
    gated = set(int(index) for index in gate_indices)
    order = sorted(range(len(rows)), key=lambda idx: (-float(scores[idx]), len(rows[idx]["sequence"]), rows[idx]["sequence"]))
    output_rows: list[dict[str, object]] = []
    for rank, index in enumerate(order, start=1):
        row = rows[index]
        output_rows.append(
            {
                "candidate_id": row["candidate_id"],
                "rank": rank,
                "sequence": row["sequence"],
                "ace_score": round(float(scores[index]), 10),
                "length": len(row["sequence"]),
                "label": row.get("label", ""),
                "origin": row.get("origin", ""),
                "provenance": row.get("provenance", ""),
                "gate_component": gate_component,
                "rank_component": rank_component,
                "gate_percent": gate_percent,
                "gated": int(index in gated),
            }
        )
    return output_rows


def recipe_report(
    *,
    recipe: Mapping[str, object],
    recipe_path: Path,
    model_dirs: Sequence[Path],
    baseline_sweep_dir: Path,
    rows: Sequence[Mapping[str, str]],
    summary_row: Mapping[str, object] | None,
    component_names: Sequence[str],
    out_dir: Path,
) -> str:
    label_counts = Counter(row.get("label", "") for row in rows)
    if summary_row:
        metric_rows = [
            [
                summary_row["top_1pct_count"],
                summary_row["top_10pct_count"],
                summary_row["top_100_count"],
                short_float(summary_row["auroc_vs_matched_decoys"]),
                short_float(summary_row["auroc_vs_all_negatives"]),
                short_float(summary_row["mean_positive_rank"], 1),
            ]
        ]
        metrics_section = markdown_table(
            ["Top 1%", "Top 10%", "Top 100", "AUROC matched", "AUROC all", "Mean positive rank"],
            metric_rows,
        )
    else:
        metrics_section = "No positive labels were available in the scored universe, so only ranks were written."

    return f"""# ACE Frozen Rerank Recipe Application

Date: {date.today().isoformat()}

## Frozen Recipe

- Recipe: `{recipe["recipe_id"]}`
- Recipe file: `{recipe_path}`
- Gate component: `{recipe["gate_component"]}`
- Rank component: `{recipe["rank_component"]}`
- Gate percent: `{recipe["gate_percent"]}`
- Include length-stratified components: `{bool(recipe.get("include_length_stratified_components", False))}`
- Candidate rows: `{len(rows)}`
- Label counts: `{dict(sorted(label_counts.items()))}`

## Inputs

- Model score directories: `{", ".join(str(path) for path in model_dirs)}`
- Baseline sweep directory: `{baseline_sweep_dir}`
- Components available after frozen expansion: `{len(component_names)}`

## Validation Metrics

{metrics_section}

## Interpretation

This run applies a locked recipe and does not search gate components, rank components, gate percentages, or fusion weights on the current panel. Any validation labels are used only after ranking to compute readout metrics.

## Outputs

- Ranked candidates: `{out_dir / "frozen-rerank-ranked-candidates.csv"}`
- Summary: `{out_dir / "frozen-rerank-summary.csv"}`
- Positive ranks: `{out_dir / "frozen-rerank-positive-ranks.csv"}`
"""


def main() -> int:
    args = parse_args()
    recipe_path = Path(args.recipe)
    recipe = load_recipe(recipe_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    model_dirs = [Path(path) for path in args.model_dir] if args.model_dir else [DEFAULT_MODEL_DIR]
    rows, components = read_wide_scores(model_dirs)
    components.update(read_baseline_scores(rows, Path(args.baseline_sweep_dir)))
    if bool(recipe.get("include_length_stratified_components", False)):
        components.update(length_stratified_components(rows, components))

    normalized = normalize_components(components)
    normalized = add_recipe_fusions(
        normalized=normalized,
        fusion_components=recipe.get("fusion_components", []),  # type: ignore[arg-type]
    )

    gate_component = str(recipe["gate_component"])
    rank_component = str(recipe["rank_component"])
    missing_components = sorted({gate_component, rank_component} - set(normalized))
    if missing_components:
        raise ValueError(f"recipe references missing components: {missing_components}")

    gate_percent = float(recipe["gate_percent"])
    if gate_percent <= 0:
        raise ValueError("gate_percent must be positive")
    top_1pct_cutoff = max(1, math.ceil(len(rows) * 0.01))
    gate_count = max(top_1pct_cutoff, math.ceil(len(rows) * gate_percent / 100.0))
    gate_scores = normalized[gate_component]
    if gate_count >= len(gate_scores):
        gate_indices = np.arange(len(gate_scores))
    else:
        gate_indices = np.argpartition(gate_scores, -gate_count)[-gate_count:]

    scores = staged_scores(gate_indices, normalized[rank_component])
    labels = np.array([1 if row.get("label", "") in POSITIVE_LABELS else 0 for row in rows], dtype=int)
    output_rows = ranked_rows(
        rows=rows,
        scores=scores,
        gate_indices=gate_indices,
        gate_component=gate_component,
        rank_component=rank_component,
        gate_percent=gate_percent,
    )
    write_csv_rows(out_dir / "frozen-rerank-ranked-candidates.csv", output_rows, RANKED_FIELDS)

    summary_row: dict[str, object] | None = None
    positive_rows: list[dict[str, object]] = []
    if int(labels.sum()) > 0:
        summary, positives = ranked_metrics(
            fusion_id=str(recipe["recipe_id"]),
            rows=rows,
            labels=labels,
            scores=scores,
        )
        summary_row = {
            "recipe_id": recipe["recipe_id"],
            "gate_component": gate_component,
            "rank_component": rank_component,
            "gate_percent": gate_percent,
            "gate_count": gate_count,
            **{key: value for key, value in summary.items() if key != "fusion_id"},
        }
        positive_rows = [
            {"recipe_id": recipe["recipe_id"], **{key: value for key, value in row.items() if key != "fusion_id"}}
            for row in positives
        ]
        write_csv_rows(out_dir / "frozen-rerank-summary.csv", [summary_row], SUMMARY_FIELDS)
        write_csv_rows(out_dir / "frozen-rerank-positive-ranks.csv", positive_rows, POSITIVE_FIELDS)
    else:
        write_csv_rows(out_dir / "frozen-rerank-summary.csv", [], SUMMARY_FIELDS)
        write_csv_rows(out_dir / "frozen-rerank-positive-ranks.csv", [], POSITIVE_FIELDS)

    component_names = sorted(normalized)
    (out_dir / "frozen-rerank-report.md").write_text(
        recipe_report(
            recipe=recipe,
            recipe_path=recipe_path,
            model_dirs=model_dirs,
            baseline_sweep_dir=Path(args.baseline_sweep_dir),
            rows=rows,
            summary_row=summary_row,
            component_names=component_names,
            out_dir=out_dir,
        ),
        encoding="utf-8",
    )
    save_json(
        out_dir / "frozen-rerank-manifest.json",
        {
            "date": date.today().isoformat(),
            "recipe_path": str(recipe_path),
            "recipe_sha256": sha256_text(json.dumps(recipe, sort_keys=True)),
            "recipe": recipe,
            "model_dirs": [str(path) for path in model_dirs],
            "baseline_sweep_dir": args.baseline_sweep_dir,
            "components": component_names,
            "gate_count": gate_count,
            "label_counts": dict(Counter(row.get("label", "") for row in rows)),
            "ranked_sha256": sha256_text(json.dumps(output_rows, sort_keys=True)),
            "summary_sha256": sha256_text(json.dumps(summary_row or {}, sort_keys=True)),
            "claim_boundary": (
                "Frozen application only. Current labels were not used to select gate/rank components, "
                "gate percent, or fusion weights."
            ),
        },
    )
    print(out_dir / "frozen-rerank-report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Search two-stage ACE rerankers over model, heuristic, and fusion scores."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rank_fusion_search import (  # noqa: E402
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
from scripts.ace_rediscovery_common import PAPER_DIR, read_csv_rows, save_json, sha256_text, write_csv_rows  # noqa: E402

DEFAULT_OUT_DIR = PAPER_DIR / "experiments" / "2026-08-29-model-iteration-001" / "two-stage-rerank-v3-dev-001"

SUMMARY_FIELDS = (
    "rerank_id",
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
POSITIVE_FIELDS = ("rerank_id", "rank", "sequence", "fusion_score", "label", "origin", "provenance")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", action="append", default=None)
    parser.add_argument("--baseline-sweep-dir", default=str(DEFAULT_BASELINE_SWEEP_DIR))
    parser.add_argument("--fusion-dir", action="append", default=[])
    parser.add_argument("--top-fusions-per-dir", type=int, default=12)
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--gate-percent", type=float, action="append", default=None)
    parser.add_argument(
        "--gate-percent-range",
        action="append",
        default=[],
        help="Inclusive floating range as start:end:step, for example 5:10:0.25.",
    )
    parser.add_argument("--gate-component", action="append", default=None)
    parser.add_argument("--rank-component", action="append", default=None)
    parser.add_argument("--keep", type=int, default=80)
    parser.add_argument("--include-length-stratified-components", action="store_true")
    parser.add_argument("--development-label", default="development panel")
    parser.add_argument("--next-holdout-label", default="new independent primary-literature holdout")
    return parser.parse_args()


def expand_gate_percents(args: argparse.Namespace) -> list[float]:
    values = list(args.gate_percent or [])
    for range_text in args.gate_percent_range:
        parts = [float(part) for part in range_text.split(":")]
        if len(parts) != 3:
            raise ValueError(f"invalid gate percent range: {range_text}")
        start, stop, step = parts
        if step <= 0:
            raise ValueError(f"invalid non-positive gate range step: {range_text}")
        current = start
        while current <= stop + step / 10.0:
            values.append(round(current, 6))
            current += step
    if not values:
        values = [2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0]
    return sorted(set(values))


def load_fusion_weights(fusion_dir: Path, keep: int) -> list[tuple[str, dict[str, float]]]:
    summary_rows = read_csv_rows(fusion_dir / "fusion-summary.csv")
    keep_ids = [row["fusion_id"] for row in summary_rows[:keep]]
    keep_set = set(keep_ids)
    by_fusion: dict[str, dict[str, float]] = defaultdict(dict)
    for row in read_csv_rows(fusion_dir / "fusion-weights.csv"):
        if row["fusion_id"] in keep_set:
            by_fusion[row["fusion_id"]][row["component"]] = float(row["weight"])
    return [(fusion_id, by_fusion[fusion_id]) for fusion_id in keep_ids if fusion_id in by_fusion]


def add_fusion_components(
    *,
    normalized: Mapping[str, np.ndarray],
    fusion_dirs: Sequence[str],
    top_fusions_per_dir: int,
) -> dict[str, np.ndarray]:
    components = dict(normalized)
    for fusion_dir_text in fusion_dirs:
        fusion_dir = Path(fusion_dir_text)
        for fusion_id, weights in load_fusion_weights(fusion_dir, top_fusions_per_dir):
            missing = sorted(set(weights) - set(normalized))
            if missing:
                raise ValueError(f"{fusion_dir} fusion {fusion_id} references missing components: {missing[:5]}")
            component_id = f"fusion:{fusion_dir.name}:{fusion_id}"
            components[component_id] = combined_score(weights, normalized)
    return normalize_components(components)


def top_count_from_gate_rank(
    labels: np.ndarray,
    gate_indices: np.ndarray,
    rank_scores: np.ndarray,
    cutoff: int,
) -> int:
    if cutoff <= 0:
        return 0
    candidate_indices = gate_indices
    if len(candidate_indices) > cutoff:
        local = np.argpartition(rank_scores[candidate_indices], -cutoff)[-cutoff:]
        candidate_indices = candidate_indices[local]
    return int(labels[candidate_indices].sum())


def staged_scores(gate_indices: np.ndarray, rank_scores: np.ndarray) -> np.ndarray:
    scores = rank_scores - 2.0
    scores = scores.copy()
    scores[gate_indices] = rank_scores[gate_indices]
    return scores


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    model_dirs = [Path(path) for path in args.model_dir] if args.model_dir else [DEFAULT_MODEL_DIR]
    rows, components = read_wide_scores(model_dirs)
    components.update(read_baseline_scores(rows, Path(args.baseline_sweep_dir)))
    if args.include_length_stratified_components:
        components.update(length_stratified_components(rows, components))
    normalized = normalize_components(components)
    normalized = add_fusion_components(
        normalized=normalized,
        fusion_dirs=args.fusion_dir,
        top_fusions_per_dir=args.top_fusions_per_dir,
    )

    labels = np.array([1 if row["label"] in POSITIVE_LABELS else 0 for row in rows], dtype=int)
    gate_percents = expand_gate_percents(args)
    top_1pct_cutoff = max(1, math.ceil(len(rows) * 0.01))
    top_10pct_cutoff = max(1, math.ceil(len(rows) * 0.10))
    component_names = sorted(normalized)
    gate_component_names = sorted(args.gate_component or component_names)
    rank_component_names = sorted(args.rank_component or component_names)
    missing_gate = sorted(set(gate_component_names) - set(component_names))
    missing_rank = sorted(set(rank_component_names) - set(component_names))
    if missing_gate or missing_rank:
        raise ValueError(f"missing components: gate={missing_gate[:5]} rank={missing_rank[:5]}")

    gates: list[tuple[str, float, np.ndarray]] = []
    for gate_percent in gate_percents:
        gate_count = max(top_1pct_cutoff, math.ceil(len(rows) * gate_percent / 100.0))
        for gate_component in gate_component_names:
            gate_scores = normalized[gate_component]
            if gate_count >= len(gate_scores):
                gate_indices = np.arange(len(gate_scores))
            else:
                gate_indices = np.argpartition(gate_scores, -gate_count)[-gate_count:]
            gates.append((gate_component, gate_percent, gate_indices))

    candidates: list[tuple[int, int, str, str, str, float, int]] = []
    for gate_component, gate_percent, gate_indices in gates:
        gate_count = len(gate_indices)
        for rank_component in rank_component_names:
            rank_scores = normalized[rank_component]
            top_1pct = top_count_from_gate_rank(labels, gate_indices, rank_scores, top_1pct_cutoff)
            top_10pct = top_count_from_gate_rank(labels, gate_indices, rank_scores, min(top_10pct_cutoff, gate_count))
            rerank_id = (
                f"gate{gate_percent:g}__{gate_component.replace(':', '__')}__"
                f"rank__{rank_component.replace(':', '__')}"
            )
            candidates.append((top_1pct, top_10pct, rerank_id, gate_component, rank_component, gate_percent, gate_count))

    candidates.sort(key=lambda item: (-item[0], -item[1], item[2]))
    shortlisted = candidates[: args.keep]
    summary_rows: list[dict[str, object]] = []
    positive_rows: list[dict[str, object]] = []
    for _top_1pct, _top_10pct, rerank_id, gate_component, rank_component, gate_percent, gate_count in shortlisted:
        gate_scores = normalized[gate_component]
        if gate_count >= len(gate_scores):
            gate_indices = np.arange(len(gate_scores))
        else:
            gate_indices = np.argpartition(gate_scores, -gate_count)[-gate_count:]
        scores = staged_scores(gate_indices, normalized[rank_component])
        summary, positives = ranked_metrics(fusion_id=rerank_id, rows=rows, labels=labels, scores=scores)
        summary_rows.append(
            {
                "rerank_id": rerank_id,
                "gate_component": gate_component,
                "rank_component": rank_component,
                "gate_percent": gate_percent,
                "gate_count": gate_count,
                **{key: value for key, value in summary.items() if key != "fusion_id"},
            }
        )
        for row in positives:
            positive_rows.append({"rerank_id": rerank_id, **{key: value for key, value in row.items() if key != "fusion_id"}})

    summary_rows.sort(
        key=lambda row: (
            -int(row["top_1pct_count"]),
            -int(row["top_10pct_count"]),
            -float(row["auroc_vs_matched_decoys"]),
            float(row["mean_positive_rank"]),
            str(row["rerank_id"]),
        )
    )
    write_csv_rows(out_dir / "two-stage-summary.csv", summary_rows, SUMMARY_FIELDS)
    write_csv_rows(out_dir / "two-stage-positive-ranks.csv", positive_rows, POSITIVE_FIELDS)

    table_rows = [
        [
            row["rerank_id"][:72],
            row["top_1pct_count"],
            row["top_10pct_count"],
            short_float(row["auroc_vs_matched_decoys"]),
            short_float(row["auroc_vs_all_negatives"]),
            short_float(row["mean_positive_rank"], 1),
        ]
        for row in summary_rows[:12]
    ]
    best = summary_rows[0] if summary_rows else {}
    report = f"""# ACE Development Two-Stage Rerank Search

Date: {date.today().isoformat()}

## Setup

- Development panel: `{args.development_label}`
- Model score directories: `{", ".join(str(path) for path in model_dirs)}`
- Baseline sweep directory: `{args.baseline_sweep_dir}`
- Fusion directories: `{", ".join(args.fusion_dir)}`
- Components after fusion expansion: `{len(component_names)}`
- Gate components searched: `{len(gate_component_names)}`
- Rank components searched: `{len(rank_component_names)}`
- Gate percents: `{", ".join(str(value) for value in gate_percents)}`
- Candidate rows: `{len(rows)}`
- Positive rows: `{int(labels.sum())}`

## Top Two-Stage Rerankers

{markdown_table(["Reranker", "Top 1%", "Top 10%", "AUROC matched", "AUROC all", "Mean positive rank"], table_rows)}

## Current Best

Best reranker: `{best.get("rerank_id", "")}`

- Gate component: `{best.get("gate_component", "")}`
- Rank component: `{best.get("rank_component", "")}`
- Gate percent: `{best.get("gate_percent", "")}`
- Top 1 percent recovery: `{best.get("top_1pct_count", "")}/{int(labels.sum())}`
- Top 10 percent recovery: `{best.get("top_10pct_count", "")}/{int(labels.sum())}`
- Matched-decoy AUROC: `{short_float(best.get("auroc_vs_matched_decoys", ""))}`

## Interpretation

This is a development search against `{args.development_label}`. A two-stage recipe selected here must be frozen before `{args.next_holdout_label}` is curated.
"""
    (out_dir / "two-stage-rerank-report.md").write_text(report, encoding="utf-8")
    save_json(
        out_dir / "two-stage-rerank-manifest.json",
        {
            "date": date.today().isoformat(),
            "model_dirs": [str(path) for path in model_dirs],
            "baseline_sweep_dir": args.baseline_sweep_dir,
            "fusion_dirs": args.fusion_dir,
            "top_fusions_per_dir": args.top_fusions_per_dir,
            "gate_percents": gate_percents,
            "gate_components": gate_component_names,
            "rank_components": rank_component_names,
            "components": component_names,
            "label_counts": dict(Counter(row["label"] for row in rows)),
            "development_label": args.development_label,
            "next_holdout_label": args.next_holdout_label,
            "claim_boundary": f"Development search only. {args.development_label} has been used for reranker selection and cannot be clean validation.",
            "summary_sha256": sha256_text(json.dumps(summary_rows, sort_keys=True)),
        },
    )
    print(out_dir / "two-stage-rerank-report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Search rank fusions over ACE heuristic and supervised model scores.

This is an exploratory development tool. If ``--labels`` come from a panel used
for this search, that panel is no longer a clean validation set.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import rankdata

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import (  # noqa: E402
    PAPER_DIR,
    average_precision,
    auroc_from_scores,
    read_csv_rows,
    save_json,
    sha256_text,
    write_csv_rows,
)
from scripts.ace_unblind_rediscovery_report import hypergeom_sf  # noqa: E402

DEFAULT_MODEL_DIR = PAPER_DIR / "experiments" / "2026-08-29-model-iteration-001" / "expanded-v2-plus-h12"
DEFAULT_BASELINE_SWEEP_DIR = PAPER_DIR / "experiments" / "2026-08-29-post-freeze-v3-primary-panel-001" / "score-mode-sweep"
DEFAULT_OUT_DIR = PAPER_DIR / "experiments" / "2026-08-29-model-iteration-001" / "fusion-search-v3-dev-001"

BASE_FIELDS = ("candidate_id", "sequence", "label", "origin", "provenance")
POSITIVE_LABELS = {"known_positive", "digestome_known_positive"}

SUMMARY_FIELDS = (
    "fusion_id",
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

WEIGHT_FIELDS = ("fusion_id", "component", "weight")
POSITIVE_FIELDS = ("fusion_id", "rank", "sequence", "fusion_score", "label", "origin", "provenance")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model-dir",
        action="append",
        default=None,
        help="Model score directory containing eval-model-scores.csv. Can be supplied more than once.",
    )
    parser.add_argument("--baseline-sweep-dir", default=str(DEFAULT_BASELINE_SWEEP_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--seed", type=int, default=2026082907)
    parser.add_argument("--trials", type=int, default=2500)
    parser.add_argument("--keep", type=int, default=40)
    parser.add_argument("--development-label", default="development panel")
    parser.add_argument("--next-holdout-label", default="new independent primary-literature holdout")
    parser.add_argument(
        "--include-length-stratified-components",
        action="store_true",
        help="Add within-length percentile ranks for every score component.",
    )
    return parser.parse_args()


def read_wide_scores(model_dirs: Sequence[Path]) -> tuple[list[dict[str, str]], dict[str, np.ndarray]]:
    base_rows: list[dict[str, str]] | None = None
    base_keys: list[tuple[str, str]] | None = None
    components: dict[str, np.ndarray] = {}

    for model_dir in model_dirs:
        rows = read_csv_rows(model_dir / "eval-model-scores.csv")
        if not rows:
            raise ValueError(f"no rows in {model_dir / 'eval-model-scores.csv'}")
        keys = [(row["candidate_id"], row["sequence"]) for row in rows]
        if base_rows is None:
            base_rows = rows
            base_keys = keys
        elif keys != base_keys:
            raise ValueError(f"{model_dir / 'eval-model-scores.csv'} does not align with the first model directory")

        component_names = [name for name in rows[0] if name not in BASE_FIELDS]
        prefix = model_dir.name
        for component_name in component_names:
            component_id = f"ml:{prefix}:{component_name}"
            components[component_id] = np.array([float(row[component_name]) for row in rows], dtype=float)

    if base_rows is None:
        raise ValueError("at least one model directory is required")
    return base_rows, components


def read_baseline_scores(rows: Sequence[Mapping[str, str]], baseline_sweep_dir: Path) -> dict[str, np.ndarray]:
    candidate_ids = [row["candidate_id"] for row in rows]
    components: dict[str, np.ndarray] = {}
    for mode_dir in sorted(path for path in baseline_sweep_dir.iterdir() if path.is_dir()):
        result_path = mode_dir / "ace-rediscovery-results.csv"
        if not result_path.exists():
            continue
        score_by_id = {row["candidate_id"]: float(row["ace_score"]) for row in read_csv_rows(result_path)}
        if set(candidate_ids) - set(score_by_id):
            raise ValueError(f"baseline {mode_dir.name} is missing candidate scores")
        components[f"mode:{mode_dir.name}"] = np.array([score_by_id[candidate_id] for candidate_id in candidate_ids], dtype=float)
    return components


def normalize_components(components: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    normalized: dict[str, np.ndarray] = {}
    for name, scores in components.items():
        if np.nanmax(scores) == np.nanmin(scores):
            normalized[name] = np.zeros_like(scores, dtype=float)
            continue
        normalized[name] = (rankdata(scores, method="average") - 1.0) / max(1.0, len(scores) - 1.0)
    return normalized


def length_stratified_components(
    rows: Sequence[Mapping[str, str]],
    components: Mapping[str, np.ndarray],
) -> dict[str, np.ndarray]:
    by_length: dict[int, list[int]] = {}
    for index, row in enumerate(rows):
        by_length.setdefault(len(row["sequence"]), []).append(index)

    adjusted: dict[str, np.ndarray] = {}
    for name, scores in components.items():
        length_scores = np.zeros_like(scores, dtype=float)
        for indices in by_length.values():
            bucket_scores = scores[indices]
            if len(indices) == 1 or np.nanmax(bucket_scores) == np.nanmin(bucket_scores):
                length_scores[indices] = 0.5
                continue
            bucket_ranks = (rankdata(bucket_scores, method="average") - 1.0) / max(1.0, len(indices) - 1.0)
            length_scores[indices] = bucket_ranks
        adjusted[f"lenrank:{name}"] = length_scores
    return adjusted


def top_count_fast(labels: np.ndarray, scores: np.ndarray, cutoff: int) -> int:
    if cutoff >= len(scores):
        return int(labels.sum())
    top_indices = np.argpartition(scores, -cutoff)[-cutoff:]
    return int(labels[top_indices].sum())


def ranked_metrics(
    *,
    fusion_id: str,
    rows: Sequence[Mapping[str, str]],
    labels: np.ndarray,
    scores: np.ndarray,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    total = len(rows)
    top_1pct_cutoff = max(1, math.ceil(total * 0.01))
    top_10pct_cutoff = max(1, math.ceil(total * 0.10))
    order = sorted(range(total), key=lambda idx: (-float(scores[idx]), len(rows[idx]["sequence"]), rows[idx]["sequence"]))
    ranks = np.empty(total, dtype=int)
    for rank, index in enumerate(order, start=1):
        ranks[index] = rank
    positive_indices = np.flatnonzero(labels == 1)
    negative_indices = np.flatnonzero(labels == 0)
    matched_negative_indices = np.array(
        [idx for idx, row in enumerate(rows) if row["label"] == "matched_decoy"],
        dtype=int,
    )
    positive_ranks = ranks[positive_indices]
    positive_count = len(positive_indices)
    top_1pct_count = int((positive_ranks <= top_1pct_cutoff).sum())
    top_10pct_count = int((positive_ranks <= top_10pct_cutoff).sum())
    top_100_count = int((positive_ranks <= 100).sum())
    expected_top_1pct = positive_count * top_1pct_cutoff / total if total else 0.0
    summary = {
        "fusion_id": fusion_id,
        "top_1pct_count": top_1pct_count,
        "top_1pct_fraction": top_1pct_count / positive_count if positive_count else float("nan"),
        "top_1pct_enrichment_over_random": top_1pct_count / max(1e-12, expected_top_1pct),
        "top_1pct_hypergeom_p": hypergeom_sf(top_1pct_count, total, positive_count, top_1pct_cutoff),
        "top_10pct_count": top_10pct_count,
        "top_10pct_fraction": top_10pct_count / positive_count if positive_count else float("nan"),
        "top_100_count": top_100_count,
        "mean_positive_rank": float(positive_ranks.mean()) if positive_count else float("nan"),
        "median_positive_rank": float(np.median(positive_ranks)) if positive_count else float("nan"),
        "best_positive_rank": int(positive_ranks.min()) if positive_count else "",
        "auroc_vs_all_negatives": auroc_from_scores(
            scores[positive_indices].tolist(),
            scores[negative_indices].tolist(),
        ),
        "auroc_vs_matched_decoys": auroc_from_scores(
            scores[positive_indices].tolist(),
            scores[matched_negative_indices].tolist(),
        ),
        "average_precision_vs_all_negatives": average_precision(
            [(int(label), float(score)) for label, score in zip(labels, scores)]
        ),
    }
    positive_rows = [
        {
            "fusion_id": fusion_id,
            "rank": int(ranks[index]),
            "sequence": rows[index]["sequence"],
            "fusion_score": round(float(scores[index]), 10),
            "label": rows[index]["label"],
            "origin": rows[index].get("origin", ""),
            "provenance": rows[index].get("provenance", ""),
        }
        for index in positive_indices
    ]
    positive_rows.sort(key=lambda row: int(row["rank"]))
    return summary, positive_rows


def candidate_weight_sets(component_names: Sequence[str], seed: int, trials: int) -> list[tuple[str, dict[str, float]]]:
    rng = np.random.default_rng(seed)
    preferred_ml_models = {
        "logreg_l2_c0_05",
        "logreg_l2_c0_2",
        "sgd_hinge_l2_alpha1e_4",
    }
    preferred_modes = {
        "mode:ace_no_exact_proline_motif",
        "mode:ace_pharmacophore_zero_label",
        "mode:oral_stability_zero_shot",
    }
    groups = {
        "all": list(component_names),
        "ml": [name for name in component_names if name.startswith("ml:")],
        "mode": [name for name in component_names if name.startswith("mode:")],
        "best_manual": sorted(
            name
            for name in component_names
            if (name.startswith("ml:") and name.rsplit(":", 1)[-1] in preferred_ml_models) or name in preferred_modes
        ),
    }
    weight_sets: list[tuple[str, dict[str, float]]] = []
    for name in component_names:
        weight_sets.append((f"single__{name.replace(':', '__')}", {name: 1.0}))
    for group_id, names in groups.items():
        if not names:
            continue
        equal = {name: 1.0 / len(names) for name in names}
        weight_sets.append((f"mean__{group_id}", equal))
    for group_id, names in groups.items():
        if len(names) < 2:
            continue
        matrix_size = len(names)
        group_trials = max(1, trials // len(groups))
        for index in range(group_trials):
            alpha = np.ones(matrix_size) * 0.55
            weights = rng.dirichlet(alpha)
            weight_sets.append(
                (
                    f"dirichlet__{group_id}__{index:04d}",
                    {name: float(weight) for name, weight in zip(names, weights) if weight > 1e-5},
                )
            )
    return weight_sets


def combined_score(weights: Mapping[str, float], normalized: Mapping[str, np.ndarray]) -> np.ndarray:
    scores: np.ndarray | None = None
    total_weight = 0.0
    for component, weight in weights.items():
        if weight <= 0:
            continue
        total_weight += weight
        if scores is None:
            scores = normalized[component] * weight
        else:
            scores = scores + normalized[component] * weight
    if scores is None or total_weight <= 0:
        raise ValueError("empty fusion weights")
    return scores / total_weight


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def short_float(value: object, digits: int = 4) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


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
    component_names = sorted(normalized)
    labels = np.array([1 if row["label"] in POSITIVE_LABELS else 0 for row in rows], dtype=int)

    top_1pct_cutoff = max(1, math.ceil(len(rows) * 0.01))
    top_10pct_cutoff = max(1, math.ceil(len(rows) * 0.10))
    candidates: list[tuple[int, int, str, dict[str, float]]] = []
    for fusion_id, weights in candidate_weight_sets(component_names, args.seed, args.trials):
        scores = combined_score(weights, normalized)
        top_1pct = top_count_fast(labels, scores, top_1pct_cutoff)
        top_10pct = top_count_fast(labels, scores, top_10pct_cutoff)
        candidates.append((top_1pct, top_10pct, fusion_id, weights))

    candidates.sort(key=lambda item: (-item[0], -item[1], item[2]))
    shortlisted: list[tuple[str, dict[str, float]]] = []
    seen: set[str] = set()
    for _top_1pct, _top_10pct, fusion_id, weights in candidates[: max(args.keep * 3, args.keep)]:
        signature = json.dumps(weights, sort_keys=True)
        if signature in seen:
            continue
        seen.add(signature)
        shortlisted.append((fusion_id, weights))
        if len(shortlisted) >= args.keep:
            break

    summary_rows: list[dict[str, object]] = []
    weight_rows: list[dict[str, object]] = []
    positive_rows: list[dict[str, object]] = []
    for fusion_id, weights in shortlisted:
        scores = combined_score(weights, normalized)
        summary, positives = ranked_metrics(fusion_id=fusion_id, rows=rows, labels=labels, scores=scores)
        summary_rows.append(summary)
        positive_rows.extend(positives)
        for component, weight in sorted(weights.items(), key=lambda item: -item[1]):
            weight_rows.append({"fusion_id": fusion_id, "component": component, "weight": round(float(weight), 8)})

    summary_rows.sort(
        key=lambda row: (
            -int(row["top_1pct_count"]),
            -int(row["top_10pct_count"]),
            -float(row["auroc_vs_matched_decoys"]),
            float(row["mean_positive_rank"]),
            str(row["fusion_id"]),
        )
    )
    write_csv_rows(out_dir / "fusion-summary.csv", summary_rows, SUMMARY_FIELDS)
    write_csv_rows(out_dir / "fusion-weights.csv", weight_rows, WEIGHT_FIELDS)
    write_csv_rows(out_dir / "fusion-positive-ranks.csv", positive_rows, POSITIVE_FIELDS)

    table_rows = [
        [
            row["fusion_id"],
            row["top_1pct_count"],
            short_float(row["top_1pct_fraction"], 3),
            short_float(row["top_1pct_enrichment_over_random"], 2) + "x",
            row["top_10pct_count"],
            short_float(row["auroc_vs_matched_decoys"]),
            short_float(row["auroc_vs_all_negatives"]),
            short_float(row["mean_positive_rank"], 1),
        ]
        for row in summary_rows[:12]
    ]
    best = summary_rows[0] if summary_rows else {}
    best_weights = [row for row in weight_rows if row["fusion_id"] == best.get("fusion_id")]
    best_weight_rows = [[row["component"], row["weight"]] for row in best_weights]
    report = f"""# ACE Development Rank-Fusion Search

Date: {date.today().isoformat()}

## Setup

- Development panel: `{args.development_label}`
- Model score directories: `{", ".join(str(path) for path in model_dirs)}`
- Baseline sweep directory: `{args.baseline_sweep_dir}`
- Components: `{len(component_names)}`
- Length-stratified components: `{args.include_length_stratified_components}`
- Search trials: `{args.trials}`
- Candidate rows: `{len(rows)}`
- Positive rows: `{int(labels.sum())}`

## Top Fusions

{markdown_table(["Fusion", "Top 1%", "Top 1% frac", "Enrichment", "Top 10%", "AUROC matched", "AUROC all", "Mean positive rank"], table_rows)}

## Best Fusion Weights

{markdown_table(["Component", "Weight"], best_weight_rows)}

## Interpretation

This is a development search against `{args.development_label}`, so it is not a validation result. If a fusion is materially better than the frozen modes, the next step is to freeze its recipe and curate `{args.next_holdout_label}`.
"""
    (out_dir / "fusion-search-report.md").write_text(report, encoding="utf-8")
    save_json(
        out_dir / "fusion-search-manifest.json",
        {
            "date": date.today().isoformat(),
            "model_dirs": [str(path) for path in model_dirs],
            "baseline_sweep_dir": args.baseline_sweep_dir,
            "seed": args.seed,
            "trials": args.trials,
            "include_length_stratified_components": args.include_length_stratified_components,
            "components": component_names,
            "label_counts": dict(Counter(row["label"] for row in rows)),
            "development_label": args.development_label,
            "next_holdout_label": args.next_holdout_label,
            "claim_boundary": f"Development search only. {args.development_label} has been used for fusion selection and cannot be a clean validation holdout.",
            "summary_sha256": sha256_text(json.dumps(summary_rows, sort_keys=True)),
        },
    )
    print(out_dir / "fusion-search-report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

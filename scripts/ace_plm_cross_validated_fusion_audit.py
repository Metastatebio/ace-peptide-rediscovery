#!/usr/bin/env python3
"""Cross-validate ACE ProteinLM fusion scenarios on parent-grouped matched decoys."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import date
from pathlib import Path
from statistics import mean, median, stdev
from typing import Any, Mapping, Sequence

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO_ROOT / "docs" / "publication-engine" / "papers" / "01-ace-rediscovery"
DEFAULT_OUT_DIR = PAPER_DIR / "experiments" / "2026-08-30-ace-plm-cv-fusion-audit-002-with-150m"
DEFAULT_SCENARIOS = (
    (
        "esm2_8m",
        PAPER_DIR
        / "experiments"
        / "2026-08-30-ace-protein-lm-fusion-001"
        / "ace_plm_scored_panel.csv",
    ),
    (
        "esm2_35m",
        PAPER_DIR
        / "experiments"
        / "2026-08-30-ace-protein-lm-fusion-002-esm35m"
        / "ace_plm_scored_panel.csv",
    ),
    (
        "esm2_150m",
        PAPER_DIR
        / "experiments"
        / "2026-08-30-ace-protein-lm-fusion-003-esm150m"
        / "ace_plm_scored_panel.csv",
    ),
)
POSITIVE_LABELS = {"known_positive", "digestome_known_positive"}
NEGATIVE_LABEL = "matched_decoy"

SUMMARY_FIELDS = (
    "scenario_id",
    "plm_feature",
    "repeats",
    "folds",
    "fold_results",
    "mean_selected_ace_weight",
    "median_selected_ace_weight",
    "mean_test_auroc",
    "sd_test_auroc",
    "mean_ace_only_auroc",
    "mean_plm_only_auroc",
    "mean_delta_vs_ace",
    "median_delta_vs_ace",
    "win_count_vs_ace",
    "win_rate_vs_ace",
    "sign_test_p_value_vs_ace",
    "mean_test_average_precision",
    "mean_ace_only_average_precision",
    "mean_plm_only_average_precision",
    "mean_top_10pct_positives",
    "mean_ace_only_top_10pct_positives",
    "label_count_positive",
    "label_count_matched_decoy",
    "group_count",
)
FOLD_FIELDS = (
    "scenario_id",
    "plm_feature",
    "repeat",
    "fold",
    "train_rows",
    "test_rows",
    "train_positive_count",
    "test_positive_count",
    "selected_ace_weight",
    "train_auroc",
    "test_auroc",
    "test_average_precision",
    "test_ace_only_auroc",
    "test_ace_only_average_precision",
    "test_plm_only_auroc",
    "test_plm_only_average_precision",
    "test_top_10pct_positives",
    "test_ace_only_top_10pct_positives",
)
FULL_PANEL_FIELDS = (
    "scenario_id",
    "plm_feature",
    "ace_weight",
    "plm_weight",
    "eval_rows",
    "positive_count",
    "matched_decoy_count",
    "top_1pct_count",
    "top_5pct_count",
    "top_10pct_count",
    "best_positive_rank",
    "median_positive_rank",
    "mean_positive_rank",
    "auroc_vs_matched_decoys",
    "average_precision_vs_matched_decoys",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260830)
    parser.add_argument(
        "--scenario",
        action="append",
        default=[],
        help="Scenario in the form label=/path/to/ace_plm_scored_panel.csv.",
    )
    return parser.parse_args()


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


def clean_id(text: str) -> str:
    return "".join(char if char.isalnum() or char in {"_", "-", "."} else "_" for char in text)


def scenario_inputs(raw_scenarios: Sequence[str]) -> list[tuple[str, Path]]:
    if not raw_scenarios:
        return [(label, path) for label, path in DEFAULT_SCENARIOS]
    parsed: list[tuple[str, Path]] = []
    for item in raw_scenarios:
        if "=" not in item:
            raise ValueError(f"scenario must use label=path format: {item}")
        label, path = item.split("=", 1)
        parsed.append((clean_id(label.strip()), Path(path)))
    return parsed


def group_key(row: Mapping[str, str]) -> str:
    if row.get("label") in POSITIVE_LABELS:
        return str(row.get("sequence", ""))
    origin = str(row.get("origin", ""))
    if origin.startswith("composition_shuffle:") or origin.startswith("length_matched_random:"):
        return origin.split(":", 1)[1]
    return str(row.get("candidate_id", row.get("sequence", "")))


def evaluation_rows(path: Path) -> list[dict[str, str]]:
    rows = [
        row
        for row in read_csv_rows(path)
        if row.get("label") in POSITIVE_LABELS or row.get("label") == NEGATIVE_LABEL
    ]
    rows.sort(key=lambda row: (group_key(row), row.get("label", ""), row.get("sequence", "")))
    return rows


def minmax_fit(values: Sequence[float]) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    return float(np.min(arr)), float(np.max(arr))


def minmax_apply(values: Sequence[float], bounds: tuple[float, float]) -> np.ndarray:
    minimum, maximum = bounds
    arr = np.asarray(values, dtype=float)
    if math.isclose(minimum, maximum):
        return np.full(len(arr), 0.5, dtype=float)
    return (arr - minimum) / (maximum - minimum)


def length_z_fit(rows: Sequence[Mapping[str, str]], values: Sequence[float]) -> dict[int, tuple[float, float]]:
    arr = np.asarray(values, dtype=float)
    global_mean = float(np.mean(arr))
    global_sd = float(np.std(arr)) or 1.0
    stats: dict[int, tuple[float, float]] = {-1: (global_mean, global_sd)}
    lengths = sorted({int(row["length"]) for row in rows})
    for length in lengths:
        indices = [index for index, row in enumerate(rows) if int(row["length"]) == length]
        subset = arr[indices]
        stats[length] = (float(np.mean(subset)), float(np.std(subset)) or global_sd)
    return stats


def length_z_apply(
    rows: Sequence[Mapping[str, str]],
    values: Sequence[float],
    stats: Mapping[int, tuple[float, float]],
) -> np.ndarray:
    fallback = stats[-1]
    output = np.zeros(len(rows), dtype=float)
    for index, row in enumerate(rows):
        value = float(values[index])
        center, scale = stats.get(int(row["length"]), fallback)
        output[index] = (value - center) / scale
    return output


def labels_for(rows: Sequence[Mapping[str, str]]) -> np.ndarray:
    return np.asarray([1 if row.get("label") in POSITIVE_LABELS else 0 for row in rows], dtype=int)


def top_10pct_positives(labels: np.ndarray, scores: np.ndarray) -> int:
    cutoff = max(1, math.ceil(len(labels) * 0.10))
    order = np.argsort(-scores, kind="mergesort")
    return int(np.sum(labels[order[:cutoff]]))


def train_test_features(
    train_rows: Sequence[Mapping[str, str]],
    test_rows: Sequence[Mapping[str, str]],
    *,
    plm_feature: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    train_ace_raw = [safe_float(row["ace_score"]) for row in train_rows]
    test_ace_raw = [safe_float(row["ace_score"]) for row in test_rows]
    ace_bounds = minmax_fit(train_ace_raw)
    train_ace = minmax_apply(train_ace_raw, ace_bounds)
    test_ace = minmax_apply(test_ace_raw, ace_bounds)

    train_plm_raw = [safe_float(row["plm_mean_log_likelihood"]) for row in train_rows]
    test_plm_raw = [safe_float(row["plm_mean_log_likelihood"]) for row in test_rows]
    if plm_feature == "length_z":
        stats = length_z_fit(train_rows, train_plm_raw)
        train_plm_values = length_z_apply(train_rows, train_plm_raw, stats)
        test_plm_values = length_z_apply(test_rows, test_plm_raw, stats)
    elif plm_feature == "raw":
        train_plm_values = np.asarray(train_plm_raw, dtype=float)
        test_plm_values = np.asarray(test_plm_raw, dtype=float)
    else:
        raise ValueError(f"unknown plm_feature: {plm_feature}")
    plm_bounds = minmax_fit(train_plm_values)
    train_plm = minmax_apply(train_plm_values, plm_bounds)
    test_plm = minmax_apply(test_plm_values, plm_bounds)
    return train_ace, test_ace, train_plm, test_plm


def metric_pair(labels: np.ndarray, scores: np.ndarray) -> tuple[float, float]:
    return float(roc_auc_score(labels, scores)), float(average_precision_score(labels, scores))


def cross_validate(
    *,
    scenario_id: str,
    rows: Sequence[Mapping[str, str]],
    plm_feature: str,
    repeats: int,
    folds: int,
    seed: int,
) -> list[dict[str, Any]]:
    labels = labels_for(rows)
    groups = np.asarray([group_key(row) for row in rows])
    weights = [round(index / 20.0, 2) for index in range(21)]
    results: list[dict[str, Any]] = []
    for repeat in range(repeats):
        splitter = StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed + repeat)
        for fold_index, (train_indices, test_indices) in enumerate(
            splitter.split(np.zeros(len(rows)), labels, groups),
            start=1,
        ):
            train_rows = [rows[index] for index in train_indices]
            test_rows = [rows[index] for index in test_indices]
            train_labels = labels[train_indices]
            test_labels = labels[test_indices]
            train_ace, test_ace, train_plm, test_plm = train_test_features(
                train_rows,
                test_rows,
                plm_feature=plm_feature,
            )

            best_weight = 1.0
            best_train_auc = -1.0
            for weight in weights:
                train_scores = (weight * train_ace) + ((1.0 - weight) * train_plm)
                train_auc, _ = metric_pair(train_labels, train_scores)
                if train_auc > best_train_auc:
                    best_train_auc = train_auc
                    best_weight = weight

            test_scores = (best_weight * test_ace) + ((1.0 - best_weight) * test_plm)
            test_auc, test_ap = metric_pair(test_labels, test_scores)
            ace_auc, ace_ap = metric_pair(test_labels, test_ace)
            plm_auc, plm_ap = metric_pair(test_labels, test_plm)
            results.append(
                {
                    "scenario_id": scenario_id,
                    "plm_feature": plm_feature,
                    "repeat": repeat + 1,
                    "fold": fold_index,
                    "train_rows": len(train_rows),
                    "test_rows": len(test_rows),
                    "train_positive_count": int(np.sum(train_labels)),
                    "test_positive_count": int(np.sum(test_labels)),
                    "selected_ace_weight": best_weight,
                    "train_auroc": best_train_auc,
                    "test_auroc": test_auc,
                    "test_average_precision": test_ap,
                    "test_ace_only_auroc": ace_auc,
                    "test_ace_only_average_precision": ace_ap,
                    "test_plm_only_auroc": plm_auc,
                    "test_plm_only_average_precision": plm_ap,
                    "test_top_10pct_positives": top_10pct_positives(test_labels, test_scores),
                    "test_ace_only_top_10pct_positives": top_10pct_positives(test_labels, test_ace),
                }
            )
    return results


def sign_test_p_value(wins: int, trials: int) -> float:
    if trials <= 0:
        return 1.0
    observed = max(wins, trials - wins)
    tail = sum(math.comb(trials, k) for k in range(observed, trials + 1)) / float(2**trials)
    return min(1.0, 2.0 * tail)


def summarize(
    *,
    scenario_id: str,
    plm_feature: str,
    rows: Sequence[Mapping[str, str]],
    fold_rows: Sequence[Mapping[str, Any]],
    repeats: int,
    folds: int,
) -> dict[str, Any]:
    deltas = [
        safe_float(row["test_auroc"]) - safe_float(row["test_ace_only_auroc"])
        for row in fold_rows
    ]
    wins = sum(1 for delta in deltas if delta > 0)
    labels = labels_for(rows)
    weights = [safe_float(row["selected_ace_weight"]) for row in fold_rows]
    test_aurocs = [safe_float(row["test_auroc"]) for row in fold_rows]
    ace_aurocs = [safe_float(row["test_ace_only_auroc"]) for row in fold_rows]
    plm_aurocs = [safe_float(row["test_plm_only_auroc"]) for row in fold_rows]
    return {
        "scenario_id": scenario_id,
        "plm_feature": plm_feature,
        "repeats": repeats,
        "folds": folds,
        "fold_results": len(fold_rows),
        "mean_selected_ace_weight": mean(weights),
        "median_selected_ace_weight": median(weights),
        "mean_test_auroc": mean(test_aurocs),
        "sd_test_auroc": stdev(test_aurocs) if len(test_aurocs) > 1 else 0.0,
        "mean_ace_only_auroc": mean(ace_aurocs),
        "mean_plm_only_auroc": mean(plm_aurocs),
        "mean_delta_vs_ace": mean(deltas),
        "median_delta_vs_ace": median(deltas),
        "win_count_vs_ace": wins,
        "win_rate_vs_ace": wins / float(len(deltas)),
        "sign_test_p_value_vs_ace": sign_test_p_value(wins, len(deltas)),
        "mean_test_average_precision": mean(safe_float(row["test_average_precision"]) for row in fold_rows),
        "mean_ace_only_average_precision": mean(
            safe_float(row["test_ace_only_average_precision"]) for row in fold_rows
        ),
        "mean_plm_only_average_precision": mean(
            safe_float(row["test_plm_only_average_precision"]) for row in fold_rows
        ),
        "mean_top_10pct_positives": mean(safe_float(row["test_top_10pct_positives"]) for row in fold_rows),
        "mean_ace_only_top_10pct_positives": mean(
            safe_float(row["test_ace_only_top_10pct_positives"]) for row in fold_rows
        ),
        "label_count_positive": int(np.sum(labels)),
        "label_count_matched_decoy": int(len(labels) - np.sum(labels)),
        "group_count": len(set(group_key(row) for row in rows)),
    }


def full_panel_metrics(
    *,
    scenario_id: str,
    rows: Sequence[Mapping[str, str]],
    plm_feature: str,
) -> list[dict[str, Any]]:
    labels = labels_for(rows)
    ace = minmax_apply([safe_float(row["ace_score"]) for row in rows], minmax_fit([safe_float(row["ace_score"]) for row in rows]))
    raw_plm = [safe_float(row["plm_mean_log_likelihood"]) for row in rows]
    if plm_feature == "length_z":
        stats = length_z_fit(rows, raw_plm)
        plm_values = length_z_apply(rows, raw_plm, stats)
    else:
        plm_values = np.asarray(raw_plm, dtype=float)
    plm = minmax_apply(plm_values, minmax_fit(plm_values))
    output = []
    for weight in [round(index / 20.0, 2) for index in range(21)]:
        scores = (weight * ace) + ((1.0 - weight) * plm)
        order = sorted(range(len(rows)), key=lambda index: (-float(scores[index]), rows[index]["sequence"]))
        positive_ranks = [rank for rank, index in enumerate(order, start=1) if int(labels[index]) == 1]
        output.append(
            {
                "scenario_id": scenario_id,
                "plm_feature": plm_feature,
                "ace_weight": weight,
                "plm_weight": round(1.0 - weight, 2),
                "eval_rows": len(rows),
                "positive_count": int(np.sum(labels)),
                "matched_decoy_count": int(len(labels) - np.sum(labels)),
                "top_1pct_count": sum(1 for rank in positive_ranks if rank <= max(1, math.ceil(len(rows) * 0.01))),
                "top_5pct_count": sum(1 for rank in positive_ranks if rank <= max(1, math.ceil(len(rows) * 0.05))),
                "top_10pct_count": sum(1 for rank in positive_ranks if rank <= max(1, math.ceil(len(rows) * 0.10))),
                "best_positive_rank": min(positive_ranks),
                "median_positive_rank": median(positive_ranks),
                "mean_positive_rank": mean(positive_ranks),
                "auroc_vs_matched_decoys": float(roc_auc_score(labels, scores)),
                "average_precision_vs_matched_decoys": float(average_precision_score(labels, scores)),
            }
        )
    return output


def fmt(value: Any, digits: int = 4) -> str:
    number = safe_float(value, math.nan)
    if math.isnan(number):
        return ""
    return f"{number:.{digits}f}"


def fmt_p_value(value: Any) -> str:
    number = safe_float(value, math.nan)
    if math.isnan(number):
        return ""
    if number < 0.0001:
        return f"{number:.2e}"
    return f"{number:.4f}"


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
    summaries: Sequence[Mapping[str, Any]],
    full_panel_rows: Sequence[Mapping[str, Any]],
    scenario_paths: Sequence[tuple[str, Path]],
    out_dir: Path,
) -> str:
    best = max(summaries, key=lambda row: safe_float(row["mean_delta_vs_ace"], -999.0))
    summary_table = markdown_table(
        [
            "Scenario",
            "PLM feature",
            "CV AUROC",
            "ACE-only AUROC",
            "Delta",
            "Win rate",
            "Median ACE weight",
            "Sign-test p",
        ],
        [
            [
                row["scenario_id"],
                row["plm_feature"],
                fmt(row["mean_test_auroc"], 4),
                fmt(row["mean_ace_only_auroc"], 4),
                fmt(row["mean_delta_vs_ace"], 4),
                fmt(row["win_rate_vs_ace"], 3),
                fmt(row["median_selected_ace_weight"], 2),
                fmt_p_value(row["sign_test_p_value_vs_ace"]),
            ]
            for row in sorted(summaries, key=lambda item: -safe_float(item["mean_delta_vs_ace"]))
        ],
    )
    full_best_by_scenario = []
    for scenario_id in sorted({row["scenario_id"] for row in full_panel_rows}):
        scenario_rows = [row for row in full_panel_rows if row["scenario_id"] == scenario_id]
        best_full = max(scenario_rows, key=lambda row: safe_float(row["auroc_vs_matched_decoys"], -999.0))
        full_best_by_scenario.append(
            [
                best_full["scenario_id"],
                best_full["plm_feature"],
                fmt(best_full["ace_weight"], 2),
                fmt(best_full["plm_weight"], 2),
                best_full["top_1pct_count"],
                fmt(best_full["auroc_vs_matched_decoys"], 4),
                fmt(best_full["average_precision_vs_matched_decoys"], 4),
                fmt(best_full["median_positive_rank"], 1),
            ]
        )
    full_table = markdown_table(
        [
            "Scenario",
            "PLM feature",
            "ACE wt",
            "PLM wt",
            "Top 1%",
            "AUROC",
            "AP",
            "Median positive rank",
        ],
        full_best_by_scenario,
    )
    verdict = (
        "strengthens"
        if safe_float(best["mean_delta_vs_ace"]) >= 0.015 and safe_float(best["win_rate_vs_ace"]) >= 0.65
        else "supporting_only"
    )
    return f"""# ACE ProteinLM Cross-Validated Fusion Audit

Date: {date.today().isoformat()}

## Inputs

{markdown_table(["Scenario", "Scored panel"], [[label, path] for label, path in scenario_paths])}

## Parent-Grouped Cross-Validation

Parent grouping keeps each literature positive and its composition-shuffled or length-matched decoys in the same fold. Fusion weights are selected on each training fold, then evaluated on held-out parent groups.

{summary_table}

Verdict: `{verdict}`.

Best cross-validated scenario: `{best["scenario_id"]}` with `{best["plm_feature"]}` PLM feature. Mean held-out AUROC `{fmt(best["mean_test_auroc"], 4)}` versus ACE-only `{fmt(best["mean_ace_only_auroc"], 4)}`, delta `{fmt(best["mean_delta_vs_ace"], 4)}`.

## Full-Panel Exploratory Readout

The full-panel table is exploratory because it selects a weight after seeing all labels. It is useful for prioritizing the next preregistered run, not as a stand-alone validation claim.

{full_table}

## Files

- Summary: `{out_dir / "ace_plm_cv_fusion_summary.csv"}`
- Fold-level rows: `{out_dir / "ace_plm_cv_fusion_folds.csv"}`
- Full-panel exploratory metrics: `{out_dir / "ace_plm_full_panel_exploratory_metrics.csv"}`
- Manifest: `{out_dir / "ace_plm_cv_fusion_manifest.json"}`
"""


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    scenarios = scenario_inputs(args.scenario)

    all_fold_rows: list[dict[str, Any]] = []
    all_summary_rows: list[dict[str, Any]] = []
    all_full_panel_rows: list[dict[str, Any]] = []
    manifest_scenarios: list[dict[str, Any]] = []
    for scenario_id, path in scenarios:
        rows = evaluation_rows(path)
        labels = labels_for(rows)
        manifest_scenarios.append(
            {
                "scenario_id": scenario_id,
                "path": str(path),
                "rows": len(rows),
                "label_counts": dict(Counter(row.get("label", "") for row in rows)),
                "group_count": len(set(group_key(row) for row in rows)),
            }
        )
        for plm_feature in ("raw", "length_z"):
            fold_rows = cross_validate(
                scenario_id=scenario_id,
                rows=rows,
                plm_feature=plm_feature,
                repeats=args.repeats,
                folds=args.folds,
                seed=args.seed,
            )
            all_fold_rows.extend(fold_rows)
            all_summary_rows.append(
                summarize(
                    scenario_id=scenario_id,
                    plm_feature=plm_feature,
                    rows=rows,
                    fold_rows=fold_rows,
                    repeats=args.repeats,
                    folds=args.folds,
                )
            )
            all_full_panel_rows.extend(
                full_panel_metrics(scenario_id=scenario_id, rows=rows, plm_feature=plm_feature)
            )

    write_csv(out_dir / "ace_plm_cv_fusion_folds.csv", all_fold_rows, FOLD_FIELDS)
    write_csv(out_dir / "ace_plm_cv_fusion_summary.csv", all_summary_rows, SUMMARY_FIELDS)
    write_csv(
        out_dir / "ace_plm_full_panel_exploratory_metrics.csv",
        all_full_panel_rows,
        FULL_PANEL_FIELDS,
    )
    write_json(
        out_dir / "ace_plm_cv_fusion_manifest.json",
        {
            "date": date.today().isoformat(),
            "repeats": args.repeats,
            "folds": args.folds,
            "seed": args.seed,
            "positive_labels": sorted(POSITIVE_LABELS),
            "negative_label": NEGATIVE_LABEL,
            "scenarios": manifest_scenarios,
            "claim_boundary": (
                "Parent-grouped cross-validation only. Full-panel best weights are exploratory and should be "
                "frozen before the next independent panel."
            ),
        },
    )
    (out_dir / "ace_plm_cv_fusion_report.md").write_text(
        report_markdown(
            summaries=all_summary_rows,
            full_panel_rows=all_full_panel_rows,
            scenario_paths=scenarios,
            out_dir=out_dir,
        ),
        encoding="utf-8",
    )
    print(out_dir / "ace_plm_cv_fusion_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

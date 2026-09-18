#!/usr/bin/env python3
"""Train ACE rank-fusion model iterations and evaluate on a held-out universe.

The intended use is development, not claim-bearing validation. A panel used as
``--eval-*`` here must be treated as contaminated for future publication claims.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import re
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MaxAbsScaler
from sklearn.svm import LinearSVC

try:
    from xgboost import XGBClassifier
except Exception:  # pragma: no cover - optional dependency
    XGBClassifier = None

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import (  # noqa: E402
    AA,
    PAPER_DIR,
    SCORING_MODE_DESCRIPTIONS,
    ace_feature_row,
    average_precision,
    auroc_from_scores,
    read_csv_rows,
    score_features_for_mode,
    save_json,
    sha256_text,
    write_csv_rows,
)
from scripts.ace_unblind_rediscovery_report import hypergeom_sf  # noqa: E402

DEFAULT_OUT_DIR = PAPER_DIR / "experiments" / "2026-08-29-model-iteration-001"
DEFAULT_EVAL_DIR = PAPER_DIR / "experiments" / "2026-08-29-post-freeze-v3-primary-panel-001"

POSITIVE_LABELS = {"known_positive", "digestome_known_positive"}
NEGATIVE_LABELS = {"matched_decoy", "random_background", "digestome_background"}

SUMMARY_FIELDS = (
    "model_id",
    "train_id",
    "positive_train_count",
    "negative_train_count",
    "eval_candidate_count",
    "eval_positive_count",
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

POSITIVE_RANK_FIELDS = (
    "model_id",
    "rank",
    "sequence",
    "model_score",
    "label",
    "origin",
    "provenance",
)

TOP_CANDIDATE_FIELDS = (
    "model_id",
    "rank",
    "sequence",
    "model_score",
    "label",
    "origin",
)

WIDE_SCORE_BASE_FIELDS = (
    "candidate_id",
    "sequence",
    "label",
    "origin",
    "provenance",
)


@dataclass(frozen=True)
class LabeledSequence:
    candidate_id: str
    sequence: str
    label: int
    source_label: str
    origin: str
    provenance: str
    sample_weight: float = 1.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--train-dir",
        action="append",
        required=True,
        help="Experiment directory containing blinded-candidate-universe.csv and blinding-key.private.csv.",
    )
    parser.add_argument("--train-id", default="train-panels")
    parser.add_argument("--eval-dir", default=str(DEFAULT_EVAL_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--seed", type=int, default=2026082904)
    parser.add_argument("--background-per-positive", type=int, default=120)
    parser.add_argument("--top-candidates", type=int, default=200)
    parser.add_argument(
        "--positive-weight-panel",
        action="append",
        default=[],
        help="CSV panel with notes containing parsed IC50 values; positives in this panel receive potency weights.",
    )
    parser.add_argument("--weight-reference-ic50-um", type=float, default=1000.0)
    parser.add_argument("--max-positive-weight", type=float, default=6.0)
    parser.add_argument(
        "--model-id",
        action="append",
        default=None,
        help="Restrict training to selected model IDs. Can be supplied more than once.",
    )
    return parser.parse_args()


def load_key_rows(experiment_dir: Path) -> list[dict[str, str]]:
    key_path = experiment_dir / "blinding-key.private.csv"
    if not key_path.exists():
        raise FileNotFoundError(key_path)
    return read_csv_rows(key_path)


def ic50_weight(ic50_um: float, *, reference_ic50_um: float, max_weight: float) -> float:
    if ic50_um <= 0:
        return max_weight
    raw = 1.0 + max(0.0, math.log10(reference_ic50_um / max(ic50_um, 1e-6)))
    return max(1.0, min(max_weight, raw))


def load_positive_weights(
    panel_paths: Sequence[str],
    *,
    reference_ic50_um: float,
    max_weight: float,
) -> dict[str, float]:
    weights: dict[str, float] = {}
    ic50_pattern = re.compile(r"retained best parsed IC50 ([0-9.eE+-]+) uM")
    for path_text in panel_paths:
        path = Path(path_text)
        if not path.exists():
            raise FileNotFoundError(path)
        for row in read_csv_rows(path):
            sequence = row["sequence"].strip().upper()
            match = ic50_pattern.search(row.get("notes", ""))
            if not match:
                continue
            weight = ic50_weight(
                float(match.group(1)),
                reference_ic50_um=reference_ic50_um,
                max_weight=max_weight,
            )
            weights[sequence] = max(weights.get(sequence, 1.0), weight)
    return weights


def merge_training_sequences(
    train_dirs: Sequence[Path],
    *,
    rng: random.Random,
    background_per_positive: int,
    positive_weights: Mapping[str, float],
) -> list[LabeledSequence]:
    by_sequence: dict[str, LabeledSequence] = {}
    negatives_by_type: dict[str, list[LabeledSequence]] = {label: [] for label in NEGATIVE_LABELS}

    for train_dir in train_dirs:
        for row in load_key_rows(train_dir):
            sequence = row["sequence"]
            source_label = row["label"]
            labeled = LabeledSequence(
                candidate_id=row.get("candidate_id", ""),
                sequence=sequence,
                label=1 if source_label in POSITIVE_LABELS else 0,
                source_label=source_label,
                origin=row.get("origin", ""),
                provenance=row.get("provenance", ""),
                sample_weight=positive_weights.get(sequence, 1.0) if source_label in POSITIVE_LABELS else 1.0,
            )
            existing = by_sequence.get(sequence)
            if (
                existing is None
                or (existing.label == 0 and labeled.label == 1)
                or (existing.label == labeled.label == 1 and labeled.sample_weight > existing.sample_weight)
            ):
                by_sequence[sequence] = labeled
            if source_label in NEGATIVE_LABELS:
                negatives_by_type[source_label].append(labeled)

    positives = [row for row in by_sequence.values() if row.label == 1]
    blocked_positive_sequences = {row.sequence for row in positives}
    matched = [
        row
        for row in {row.sequence: row for row in negatives_by_type["matched_decoy"]}.values()
        if row.sequence not in blocked_positive_sequences
    ]
    digestome = [
        row
        for row in {row.sequence: row for row in negatives_by_type["digestome_background"]}.values()
        if row.sequence not in blocked_positive_sequences
    ]
    random_bg = [
        row
        for row in {row.sequence: row for row in negatives_by_type["random_background"]}.values()
        if row.sequence not in blocked_positive_sequences
    ]

    max_bg = max(0, background_per_positive * max(1, len(positives)))
    rng.shuffle(digestome)
    rng.shuffle(random_bg)
    negatives = matched + digestome[: max_bg // 2] + random_bg[: max_bg]
    return positives + negatives


def load_eval_sequences(eval_dir: Path) -> list[LabeledSequence]:
    rows: list[LabeledSequence] = []
    for row in load_key_rows(eval_dir):
        rows.append(
            LabeledSequence(
                candidate_id=row.get("candidate_id", ""),
                sequence=row["sequence"],
                label=1 if row["label"] in POSITIVE_LABELS else 0,
                source_label=row["label"],
                origin=row.get("origin", ""),
                provenance=row.get("provenance", ""),
            )
        )
    return rows


def terminal_slice(sequence: str, *, start: int, length: int) -> str:
    if len(sequence) < start + length:
        return "_"
    return sequence[start : start + length]


def cterminal_slice(sequence: str, *, end_offset: int, length: int) -> str:
    if len(sequence) < end_offset + length:
        return "_"
    end = len(sequence) - end_offset
    start = end - length
    return sequence[start:end]


def edge_ngram_counts(sequence: str, max_n: int = 3) -> dict[str, float]:
    features: dict[str, float] = {}
    for ngram_len in range(1, max_n + 1):
        if len(sequence) < ngram_len:
            continue
        features[f"n_edge_{ngram_len}={sequence[:ngram_len]}"] = 1.0
        features[f"c_edge_{ngram_len}={sequence[-ngram_len:]}"] = 1.0
    return features


def internal_kmer_features(sequence: str) -> dict[str, float]:
    features: dict[str, float] = {}
    for k in (2, 3):
        if len(sequence) < k:
            continue
        denominator = max(1, len(sequence) - k + 1)
        counts = Counter(sequence[index : index + k] for index in range(denominator))
        for kmer, count in counts.items():
            features[f"kmer_{k}={kmer}"] = count / denominator
    return features


def feature_dict(sequence: str) -> dict[str, float | str]:
    base = ace_feature_row(sequence)
    features: dict[str, float | str] = {}
    for key, value in base.items():
        if key == "sequence":
            continue
        try:
            features[key] = float(value)
        except (TypeError, ValueError):
            pass

    for mode in SCORING_MODE_DESCRIPTIONS:
        features[f"score_mode:{mode}"] = score_features_for_mode(base, mode)

    length = len(sequence)
    features["length_squared"] = float(length * length)
    features["has_terminal_proline"] = 1.0 if sequence.endswith("P") else 0.0
    features["has_terminal_aromatic"] = 1.0 if sequence[-1:] in {"F", "W", "Y"} else 0.0
    features["has_terminal_basic"] = 1.0 if sequence[-1:] in {"K", "R", "H"} else 0.0
    features["has_terminal_acidic"] = 1.0 if sequence[-1:] in {"D", "E"} else 0.0
    features["contains_pp"] = 1.0 if "PP" in sequence else 0.0
    features["contains_gp"] = 1.0 if "GP" in sequence else 0.0
    features["contains_pr"] = 1.0 if "PR" in sequence else 0.0
    features["contains_lrp"] = 1.0 if "LRP" in sequence else 0.0
    features["contains_yy"] = 1.0 if "YY" in sequence else 0.0

    counts = Counter(sequence)
    for aa in AA:
        features[f"aa_count:{aa}"] = float(counts[aa])
        features[f"aa_frac:{aa}"] = counts[aa] / length if length else 0.0
    for index in range(5):
        features[f"n_pos_{index + 1}"] = sequence[index] if index < length else "_"
        features[f"c_pos_{index + 1}"] = sequence[-(index + 1)] if index < length else "_"
    for start in range(3):
        features[f"n2_at_{start + 1}"] = terminal_slice(sequence, start=start, length=2)
        features[f"c2_at_{start + 1}"] = cterminal_slice(sequence, end_offset=start, length=2)
    features.update(edge_ngram_counts(sequence))
    features.update(internal_kmer_features(sequence))
    return features


def rows_to_feature_matrix(rows: Sequence[LabeledSequence]) -> tuple[list[dict[str, float | str]], np.ndarray]:
    return [feature_dict(row.sequence) for row in rows], np.array([row.label for row in rows], dtype=int)


def model_specs(seed: int) -> list[tuple[str, object]]:
    specs: list[tuple[str, object]] = [
        (
            "logreg_l2_c0_05",
            LogisticRegression(
                C=0.05,
                class_weight="balanced",
                max_iter=3000,
                random_state=seed,
                solver="liblinear",
            ),
        ),
        (
            "logreg_l2_c0_2",
            LogisticRegression(
                C=0.2,
                class_weight="balanced",
                max_iter=3000,
                random_state=seed,
                solver="liblinear",
            ),
        ),
        (
            "logreg_l2_c1",
            LogisticRegression(
                C=1.0,
                class_weight="balanced",
                max_iter=3000,
                random_state=seed,
                solver="liblinear",
            ),
        ),
        (
            "linear_svc_c0_05",
            LinearSVC(C=0.05, class_weight="balanced", dual="auto", max_iter=5000, random_state=seed),
        ),
        (
            "linear_svc_c0_2",
            LinearSVC(C=0.2, class_weight="balanced", dual="auto", max_iter=5000, random_state=seed),
        ),
        (
            "sgd_log_l2_alpha1e_4",
            SGDClassifier(
                alpha=1e-4,
                class_weight="balanced",
                loss="log_loss",
                max_iter=3000,
                random_state=seed,
                tol=1e-4,
            ),
        ),
        (
            "sgd_hinge_l2_alpha1e_4",
            SGDClassifier(
                alpha=1e-4,
                class_weight="balanced",
                loss="hinge",
                max_iter=3000,
                random_state=seed,
                tol=1e-4,
            ),
        ),
    ]
    if XGBClassifier is not None:
        specs.extend(
            [
                (
                    "xgb_depth2_spw30_n300",
                    XGBClassifier(
                        objective="binary:logistic",
                        eval_metric="logloss",
                        n_estimators=300,
                        max_depth=2,
                        learning_rate=0.05,
                        min_child_weight=5,
                        subsample=0.85,
                        colsample_bytree=0.75,
                        reg_lambda=6.0,
                        scale_pos_weight=30.0,
                        tree_method="hist",
                        n_jobs=8,
                        random_state=seed,
                    ),
                ),
                (
                    "xgb_depth3_spw30_n350",
                    XGBClassifier(
                        objective="binary:logistic",
                        eval_metric="logloss",
                        n_estimators=350,
                        max_depth=3,
                        learning_rate=0.04,
                        min_child_weight=4,
                        subsample=0.85,
                        colsample_bytree=0.75,
                        reg_lambda=8.0,
                        scale_pos_weight=30.0,
                        tree_method="hist",
                        n_jobs=8,
                        random_state=seed,
                    ),
                ),
                (
                    "xgb_depth3_spw80_n350",
                    XGBClassifier(
                        objective="binary:logistic",
                        eval_metric="logloss",
                        n_estimators=350,
                        max_depth=3,
                        learning_rate=0.04,
                        min_child_weight=4,
                        subsample=0.85,
                        colsample_bytree=0.75,
                        reg_lambda=8.0,
                        scale_pos_weight=80.0,
                        tree_method="hist",
                        n_jobs=8,
                        random_state=seed,
                    ),
                ),
            ]
        )
    return specs


def selected_model_specs(seed: int, selected_model_ids: Sequence[str] | None) -> list[tuple[str, object]]:
    specs = model_specs(seed)
    if not selected_model_ids:
        return specs
    selected = set(selected_model_ids)
    known = {model_id for model_id, _estimator in specs}
    missing = selected - known
    if missing:
        raise ValueError(f"unknown model IDs: {', '.join(sorted(missing))}")
    return [(model_id, estimator) for model_id, estimator in specs if model_id in selected]


def make_pipeline(model: object) -> Pipeline:
    return Pipeline(
        [
            ("vectorizer", DictVectorizer(sparse=True)),
            ("scaler", MaxAbsScaler()),
            ("model", model),
        ]
    )


def model_scores(pipeline: Pipeline, features: Sequence[dict[str, float | str]]) -> np.ndarray:
    estimator = pipeline.named_steps["model"]
    if hasattr(estimator, "predict_proba"):
        return pipeline.predict_proba(features)[:, 1]
    if hasattr(estimator, "decision_function"):
        return pipeline.decision_function(features)
    return pipeline.predict(features)


def pct_rank_threshold(total: int, percent: float) -> int:
    return max(1, math.ceil(total * percent / 100.0))


def safe_auc(labels: Sequence[int], scores: Sequence[float]) -> float:
    if len(set(labels)) < 2:
        return float("nan")
    return float(roc_auc_score(labels, scores))


def safe_average_precision(labels: Sequence[int], scores: Sequence[float]) -> float:
    if sum(labels) == 0:
        return float("nan")
    return float(average_precision_score(labels, scores))


def evaluate_model(
    *,
    model_id: str,
    train_id: str,
    train_rows: Sequence[LabeledSequence],
    eval_rows: Sequence[LabeledSequence],
    scores: Sequence[float],
    top_candidates: int,
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    ranked = sorted(
        zip(eval_rows, scores),
        key=lambda item: (-float(item[1]), len(item[0].sequence), item[0].sequence),
    )
    total = len(ranked)
    top_1pct_cutoff = pct_rank_threshold(total, 1.0)
    top_10pct_cutoff = pct_rank_threshold(total, 10.0)
    positives: list[tuple[int, LabeledSequence, float]] = []
    matched_scores: list[tuple[int, float]] = []
    labels_and_scores: list[tuple[int, float]] = []
    positive_rank_rows: list[dict[str, object]] = []
    top_candidate_rows: list[dict[str, object]] = []

    for rank, (row, score) in enumerate(ranked, start=1):
        score_float = float(score)
        labels_and_scores.append((row.label, score_float))
        if row.source_label == "matched_decoy":
            matched_scores.append((row.label, score_float))
        if row.label == 1:
            positives.append((rank, row, score_float))
            positive_rank_rows.append(
                {
                    "model_id": model_id,
                    "rank": rank,
                    "sequence": row.sequence,
                    "model_score": round(score_float, 8),
                    "label": row.source_label,
                    "origin": row.origin,
                    "provenance": row.provenance,
                }
            )
        if rank <= top_candidates:
            top_candidate_rows.append(
                {
                    "model_id": model_id,
                    "rank": rank,
                    "sequence": row.sequence,
                    "model_score": round(score_float, 8),
                    "label": row.source_label,
                    "origin": row.origin,
                }
            )

    positive_ranks = [rank for rank, _row, _score in positives]
    positive_scores = [score for _rank, _row, score in positives]
    negative_scores = [score for row, score in ranked if row.label == 0]
    matched_negative_scores = [score for label, score in matched_scores if label == 0]
    top_1pct_count = sum(rank <= top_1pct_cutoff for rank in positive_ranks)
    top_10pct_count = sum(rank <= top_10pct_cutoff for rank in positive_ranks)
    top_100_count = sum(rank <= 100 for rank in positive_ranks)
    positive_count = len(positive_ranks)
    expected_top_1pct = positive_count * top_1pct_cutoff / total if total else 0.0
    summary = {
        "model_id": model_id,
        "train_id": train_id,
        "positive_train_count": sum(row.label for row in train_rows),
        "negative_train_count": sum(1 - row.label for row in train_rows),
        "eval_candidate_count": total,
        "eval_positive_count": positive_count,
        "top_1pct_count": top_1pct_count,
        "top_1pct_fraction": top_1pct_count / positive_count if positive_count else float("nan"),
        "top_1pct_enrichment_over_random": top_1pct_count / max(1e-12, expected_top_1pct),
        "top_1pct_hypergeom_p": hypergeom_sf(top_1pct_count, total, positive_count, top_1pct_cutoff),
        "top_10pct_count": top_10pct_count,
        "top_10pct_fraction": top_10pct_count / positive_count if positive_count else float("nan"),
        "top_100_count": top_100_count,
        "mean_positive_rank": sum(positive_ranks) / positive_count if positive_count else float("nan"),
        "median_positive_rank": float(np.median(positive_ranks)) if positive_count else float("nan"),
        "best_positive_rank": min(positive_ranks) if positive_ranks else "",
        "auroc_vs_all_negatives": auroc_from_scores(positive_scores, negative_scores),
        "auroc_vs_matched_decoys": auroc_from_scores(positive_scores, matched_negative_scores),
        "average_precision_vs_all_negatives": average_precision(labels_and_scores),
    }
    return summary, positive_rank_rows, top_candidate_rows


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


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    args = parse_args()
    rng = random.Random(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    train_dirs = [Path(path) for path in args.train_dir]
    eval_dir = Path(args.eval_dir)
    positive_weights = load_positive_weights(
        args.positive_weight_panel,
        reference_ic50_um=args.weight_reference_ic50_um,
        max_weight=args.max_positive_weight,
    )
    specs = selected_model_specs(args.seed, args.model_id)

    train_rows = merge_training_sequences(
        train_dirs,
        rng=rng,
        background_per_positive=args.background_per_positive,
        positive_weights=positive_weights,
    )
    eval_rows = load_eval_sequences(eval_dir)
    train_features, train_labels = rows_to_feature_matrix(train_rows)
    eval_features, _eval_labels = rows_to_feature_matrix(eval_rows)
    summary_rows: list[dict[str, object]] = []
    positive_rank_rows: list[dict[str, object]] = []
    top_candidate_rows: list[dict[str, object]] = []
    wide_score_rows: list[dict[str, object]] = [
        {
            "candidate_id": row.candidate_id,
            "sequence": row.sequence,
            "label": row.source_label,
            "origin": row.origin,
            "provenance": row.provenance,
        }
        for row in eval_rows
    ]

    sample_weights = np.array([row.sample_weight for row in train_rows], dtype=float)
    for model_id, estimator in specs:
        pipeline = make_pipeline(estimator)
        pipeline.fit(train_features, train_labels, model__sample_weight=sample_weights)
        scores = model_scores(pipeline, eval_features)
        for score_row, score in zip(wide_score_rows, scores):
            score_row[model_id] = round(float(score), 10)
        summary, positives, top_rows = evaluate_model(
            model_id=model_id,
            train_id=args.train_id,
            train_rows=train_rows,
            eval_rows=eval_rows,
            scores=scores,
            top_candidates=args.top_candidates,
        )
        summary_rows.append(summary)
        positive_rank_rows.extend(positives)
        top_candidate_rows.extend(top_rows)

    summary_rows.sort(
        key=lambda row: (
            -int(row["top_1pct_count"]),
            -float(row["auroc_vs_matched_decoys"]),
            float(row["mean_positive_rank"]),
            str(row["model_id"]),
        )
    )
    write_csv_rows(out_dir / "model-iteration-summary.csv", summary_rows, SUMMARY_FIELDS)
    write_csv_rows(out_dir / "positive-ranks.csv", positive_rank_rows, POSITIVE_RANK_FIELDS)
    write_csv_rows(out_dir / "top-candidates.csv", top_candidate_rows, TOP_CANDIDATE_FIELDS)
    score_fields = list(WIDE_SCORE_BASE_FIELDS) + [model_id for model_id, _estimator in specs]
    write_csv_rows(out_dir / "eval-model-scores.csv", wide_score_rows, score_fields)

    positive_train_weights = [row.sample_weight for row in train_rows if row.label == 1]
    table_rows = [
        [
            row["model_id"],
            row["top_1pct_count"],
            short_float(row["top_1pct_fraction"], 3),
            short_float(row["top_1pct_enrichment_over_random"], 2) + "x",
            row["top_10pct_count"],
            short_float(row["auroc_vs_matched_decoys"]),
            short_float(row["auroc_vs_all_negatives"]),
            short_float(row["mean_positive_rank"], 1),
        ]
        for row in summary_rows
    ]
    best = summary_rows[0] if summary_rows else {}
    report = f"""# ACE Rank-Fusion Model Iteration

Date: {date.today().isoformat()}

## Setup

- Development status: exploratory model iteration. The evaluation panel used here must not be reused as a clean publication holdout.
- Train ID: `{args.train_id}`
- Train directories: `{", ".join(str(path) for path in train_dirs)}`
- Eval directory: `{eval_dir}`
- Training positives: `{sum(row.label for row in train_rows)}`
- Training negatives: `{sum(1 - row.label for row in train_rows)}`
- Potency-weighted positives: `{sum(1 for weight in positive_train_weights if weight > 1.0)}`
- Positive weight range: `{short_float(min(positive_train_weights), 2) if positive_train_weights else ""}` to `{short_float(max(positive_train_weights), 2) if positive_train_weights else ""}`
- Evaluation candidates: `{len(eval_rows)}`
- Evaluation positives: `{sum(row.label for row in eval_rows)}`

## Results

{markdown_table(["Model", "Top 1%", "Top 1% frac", "Enrichment", "Top 10%", "AUROC matched", "AUROC all", "Mean positive rank"], table_rows)}

## Current Best

Best model by the pre-specified sort is `{best.get("model_id", "")}`.

- Top 1 percent recovery: `{best.get("top_1pct_count", "")}/{best.get("eval_positive_count", "")}`
- Top 1 percent enrichment: `{short_float(best.get("top_1pct_enrichment_over_random", ""), 2)}x`
- Matched-decoy AUROC: `{short_float(best.get("auroc_vs_matched_decoys", ""))}`
- Mean positive rank: `{short_float(best.get("mean_positive_rank", ""), 1)}`

## Interpretation

This is a development run. A useful model must improve top-rank recovery on v3 without destroying matched-decoy AUROC, then be frozen before a new v4 literature panel is curated.
"""
    (out_dir / "model-iteration-report.md").write_text(report, encoding="utf-8")
    manifest = {
        "date": date.today().isoformat(),
        "script": str(Path(__file__)),
        "script_sha256": file_sha256(Path(__file__)),
        "train_id": args.train_id,
        "train_dirs": [str(path) for path in train_dirs],
        "eval_dir": str(eval_dir),
        "seed": args.seed,
        "background_per_positive": args.background_per_positive,
        "positive_weight_panels": args.positive_weight_panel,
        "weight_reference_ic50_um": args.weight_reference_ic50_um,
        "max_positive_weight": args.max_positive_weight,
        "selected_model_ids": [model_id for model_id, _estimator in specs],
        "positive_train_weight_summary": {
            "weighted_positive_count": sum(1 for weight in positive_train_weights if weight > 1.0),
            "min": min(positive_train_weights) if positive_train_weights else None,
            "mean": sum(positive_train_weights) / len(positive_train_weights) if positive_train_weights else None,
            "max": max(positive_train_weights) if positive_train_weights else None,
        },
        "training_label_counts": dict(Counter(row.source_label for row in train_rows)),
        "eval_label_counts": dict(Counter(row.source_label for row in eval_rows)),
        "claim_boundary": "Exploratory rank-fusion iteration. Do not use as final validation because the evaluation panel is now a development set.",
        "feature_families": [
            "ACE heuristic feature rows",
            "all frozen score-mode values",
            "amino-acid composition",
            "N- and C-terminal position categories",
            "edge n-grams",
            "internal 2-mer and 3-mer frequencies",
            "simple motif flags",
        ],
        "summary_sha256": sha256_text(json.dumps(summary_rows, sort_keys=True)),
    }
    save_json(out_dir / "model-iteration-manifest.json", manifest)
    print(out_dir / "model-iteration-report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

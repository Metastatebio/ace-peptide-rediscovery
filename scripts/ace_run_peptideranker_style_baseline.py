#!/usr/bin/env python3
"""Run a transparent PeptideRanker-style baseline on an ACE blinded universe.

This is not the official PeptideRanker neural-network server. It is a local,
deterministic baseline inspired by Mooney et al. 2012: a general bioactivity
score using peptide length and amino-acid composition, without ACE-specific
motifs, ACE labels, or target information.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import (  # noqa: E402
    PAPER_DIR,
    average_precision,
    auroc_from_scores,
    clean_sequence,
    read_csv_rows,
    save_json,
    write_csv_rows,
)
from scripts.ace_unblind_rediscovery_report import hypergeom_sf, pct_rank_threshold  # noqa: E402

DEFAULT_EVAL_DIR = PAPER_DIR / "experiments" / "2026-08-29-v6-independent-validation-001"
POSITIVE_LABELS = {"known_positive", "digestome_known_positive"}

RANKED_FIELDS = (
    "candidate_id",
    "rank",
    "sequence",
    "ace_score",
    "length",
)
SUMMARY_FIELDS = (
    "baseline_id",
    "baseline_family",
    "total_candidates",
    "known_positive_count",
    "top_1pct_count",
    "top_1pct_fraction",
    "top_1pct_enrichment_over_random",
    "top_1pct_hypergeom_p",
    "top_5pct_count",
    "top_5pct_fraction",
    "top_5pct_hypergeom_p",
    "top_10pct_count",
    "top_10pct_fraction",
    "top_10pct_hypergeom_p",
    "top_100_count",
    "mean_positive_rank",
    "median_positive_rank",
    "best_positive_rank",
    "auroc_vs_all_negatives",
    "auroc_vs_matched_decoys",
    "average_precision_vs_all_negatives",
)
POSITIVE_FIELDS = ("baseline_id", "rank", "sequence", "score", "label", "origin", "provenance")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eval-dir", default=str(DEFAULT_EVAL_DIR))
    parser.add_argument("--out-dir", default="")
    return parser.parse_args()


def triangular(value: float, target: float, width: float) -> float:
    if width <= 0:
        return 0.0
    return max(0.0, min(1.0, 1.0 - abs(value - target) / width))


def fraction(sequence: str, residues: set[str]) -> float:
    if not sequence:
        return 0.0
    return sum(1 for residue in sequence if residue in residues) / len(sequence)


def length_window_score(length: int) -> float:
    if 4 <= length <= 20:
        return 1.0
    if 2 <= length <= 3:
        return 0.55 + 0.15 * (length - 2)
    if 21 <= length <= 35:
        return max(0.0, 1.0 - (length - 20) / 15.0)
    return 0.0


def peptideranker_style_score(sequence: str) -> float:
    """Return a general bioactivity score in [0, 1].

    The formula intentionally avoids ACE-specific residues at fixed positions.
    It follows broad short-peptide composition trends reported for
    PeptideRanker rather than reproducing the original N-to-1 neural network.
    """
    sequence = clean_sequence(sequence)
    length = len(sequence)
    if not sequence:
        return 0.0
    f_fraction = fraction(sequence, {"F"})
    g_fraction = fraction(sequence, {"G"})
    c_fraction = fraction(sequence, {"C"})
    e_fraction = fraction(sequence, {"E"})
    t_fraction = fraction(sequence, {"T"})
    hydrophobic_fraction = fraction(sequence, set("AILMFWYV"))
    charged_fraction = fraction(sequence, set("DEKRH"))
    composition_score = (
        0.24 * triangular(f_fraction, 0.073, 0.10)
        + 0.18 * triangular(g_fraction, 0.10, 0.14)
        + 0.10 * triangular(c_fraction, 0.035, 0.07)
        + 0.16 * triangular(hydrophobic_fraction, 0.42, 0.42)
        + 0.14 * triangular(charged_fraction, 0.18, 0.24)
        + 0.10 * (1.0 - min(1.0, e_fraction / 0.12))
        + 0.08 * (1.0 - min(1.0, t_fraction / 0.12))
    )
    score = 0.35 * length_window_score(length) + 0.65 * composition_score
    return round(max(0.0, min(1.0, score)), 10)


def stable_tiebreaker(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def ranked_rows(public_rows: Sequence[Mapping[str, str]]) -> list[dict[str, object]]:
    scored = [
        {
            "candidate_id": row["candidate_id"],
            "sequence": row["sequence"],
            "length": int(row.get("length") or len(row["sequence"])),
            "ace_score": peptideranker_style_score(row["sequence"]),
        }
        for row in public_rows
    ]
    scored.sort(
        key=lambda row: (
            -float(row["ace_score"]),
            int(row["length"]),
            str(row["sequence"]),
            stable_tiebreaker(str(row["sequence"])),
        )
    )
    for rank, row in enumerate(scored, start=1):
        row["rank"] = rank
    return scored


def summarize(
    *,
    baseline_id: str,
    ranked: Sequence[Mapping[str, object]],
    key_rows: Sequence[Mapping[str, str]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    key_by_id = {row["candidate_id"]: row for row in key_rows}
    total = len(ranked)
    top_1pct_cutoff = pct_rank_threshold(total, 1.0)
    top_5pct_cutoff = pct_rank_threshold(total, 5.0)
    top_10pct_cutoff = pct_rank_threshold(total, 10.0)
    positives: list[tuple[int, Mapping[str, object], Mapping[str, str]]] = []
    positive_scores: list[float] = []
    negative_scores: list[float] = []
    matched_scores: list[float] = []
    labels_and_scores: list[tuple[int, float]] = []
    positive_rows: list[dict[str, object]] = []

    for row in ranked:
        key = key_by_id.get(str(row["candidate_id"]), {})
        label = key.get("label", "")
        is_positive = label in POSITIVE_LABELS
        score = float(row["ace_score"])
        labels_and_scores.append((1 if is_positive else 0, score))
        if is_positive:
            rank = int(row["rank"])
            positives.append((rank, row, key))
            positive_scores.append(score)
            positive_rows.append(
                {
                    "baseline_id": baseline_id,
                    "rank": rank,
                    "sequence": row["sequence"],
                    "score": score,
                    "label": label,
                    "origin": key.get("origin", ""),
                    "provenance": key.get("provenance", ""),
                }
            )
        else:
            negative_scores.append(score)
        if label == "matched_decoy":
            matched_scores.append(score)

    positive_ranks = [rank for rank, _row, _key in positives]
    positive_count = len(positive_ranks)
    top_1pct_count = sum(rank <= top_1pct_cutoff for rank in positive_ranks)
    top_5pct_count = sum(rank <= top_5pct_cutoff for rank in positive_ranks)
    top_10pct_count = sum(rank <= top_10pct_cutoff for rank in positive_ranks)
    expected_top_1pct = positive_count * top_1pct_cutoff / total if total else 0.0
    summary = {
        "baseline_id": baseline_id,
        "baseline_family": "external_literature_inspired_general_bioactivity",
        "total_candidates": total,
        "known_positive_count": positive_count,
        "top_1pct_count": top_1pct_count,
        "top_1pct_fraction": top_1pct_count / positive_count if positive_count else float("nan"),
        "top_1pct_enrichment_over_random": top_1pct_count / max(1e-12, expected_top_1pct),
        "top_1pct_hypergeom_p": hypergeom_sf(top_1pct_count, total, positive_count, top_1pct_cutoff),
        "top_5pct_count": top_5pct_count,
        "top_5pct_fraction": top_5pct_count / positive_count if positive_count else float("nan"),
        "top_5pct_hypergeom_p": hypergeom_sf(top_5pct_count, total, positive_count, top_5pct_cutoff),
        "top_10pct_count": top_10pct_count,
        "top_10pct_fraction": top_10pct_count / positive_count if positive_count else float("nan"),
        "top_10pct_hypergeom_p": hypergeom_sf(top_10pct_count, total, positive_count, top_10pct_cutoff),
        "top_100_count": sum(rank <= 100 for rank in positive_ranks),
        "mean_positive_rank": sum(positive_ranks) / positive_count if positive_count else float("nan"),
        "median_positive_rank": float(sorted(positive_ranks)[positive_count // 2]) if positive_count else float("nan"),
        "best_positive_rank": min(positive_ranks) if positive_ranks else "",
        "auroc_vs_all_negatives": auroc_from_scores(positive_scores, negative_scores),
        "auroc_vs_matched_decoys": auroc_from_scores(positive_scores, matched_scores),
        "average_precision_vs_all_negatives": average_precision(labels_and_scores),
    }
    return summary, positive_rows


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def report_text(eval_dir: Path, out_dir: Path, summary: Mapping[str, object]) -> str:
    rows = [
        [
            summary["top_1pct_count"],
            summary["top_5pct_count"],
            summary["top_10pct_count"],
            f"{float(summary['top_1pct_hypergeom_p']):.3g}",
            f"{float(summary['auroc_vs_matched_decoys']):.4f}",
            f"{float(summary['auroc_vs_all_negatives']):.4f}",
        ]
    ]
    return f"""# ACE PeptideRanker-Style Baseline

Date: {date.today().isoformat()}

## Baseline

- Baseline ID: `{summary["baseline_id"]}`
- Family: general bioactive-peptide comparator inspired by Mooney et al. 2012 PeptideRanker.
- Literature reference: https://doi.org/10.1371/journal.pone.0045012
- Evaluation directory: `{eval_dir}`
- Candidate rows: `{summary["total_candidates"]}`
- Known positives: `{summary["known_positive_count"]}`

This is not the official PeptideRanker neural-network server. It is a transparent local comparator using peptide length and amino-acid composition only. It does not use ACE-specific motifs, ACE labels, target structure, or holdout labels during ranking.

## Results

{markdown_table(["Top 1%", "Top 5%", "Top 10%", "Top-1 p", "AUROC matched", "AUROC all"], rows)}

## Outputs

- Ranked candidates: `{out_dir / "peptideranker-style-ranked-candidates.csv"}`
- Positive ranks: `{out_dir / "peptideranker-style-positive-ranks.csv"}`
- Summary: `{out_dir / "peptideranker-style-summary.csv"}`
"""


def main() -> int:
    args = parse_args()
    eval_dir = Path(args.eval_dir)
    out_dir = Path(args.out_dir) if args.out_dir else eval_dir / "external-baselines" / "peptideranker-style"
    out_dir.mkdir(parents=True, exist_ok=True)
    public_rows = read_csv_rows(eval_dir / "blinded-candidate-universe.csv")
    key_rows = read_csv_rows(eval_dir / "blinding-key.private.csv")
    baseline_id = "peptideranker_style_mooney2012_composition_v1"
    ranked = ranked_rows(public_rows)
    summary, positive_rows = summarize(baseline_id=baseline_id, ranked=ranked, key_rows=key_rows)
    write_csv_rows(out_dir / "peptideranker-style-ranked-candidates.csv", ranked, RANKED_FIELDS)
    write_csv_rows(out_dir / "peptideranker-style-positive-ranks.csv", positive_rows, POSITIVE_FIELDS)
    write_csv_rows(out_dir / "peptideranker-style-summary.csv", [summary], SUMMARY_FIELDS)
    save_json(
        out_dir / "peptideranker-style-manifest.json",
        {
            "date": date.today().isoformat(),
            "script": str(Path(__file__)),
            "eval_dir": str(eval_dir),
            "baseline_id": baseline_id,
            "literature_reference": "https://doi.org/10.1371/journal.pone.0045012",
            "claim_boundary": "Local PeptideRanker-style comparator, not official PeptideRanker server output.",
            "feature_policy": "Length and amino-acid composition only; no ACE target, ACE labels, or holdout labels.",
        },
    )
    (out_dir / "peptideranker-style-report.md").write_text(report_text(eval_dir, out_dir, summary), encoding="utf-8")
    print(out_dir / "peptideranker-style-report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

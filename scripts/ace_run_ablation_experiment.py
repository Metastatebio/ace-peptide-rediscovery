#!/usr/bin/env python3
"""Run ACE scoring ablations against the blinded rediscovery universe."""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Callable, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import (  # noqa: E402
    BLINDING_KEY_PATH,
    PAPER_DIR,
    RANKED_CANDIDATES_PATH,
    auroc_from_scores,
    average_precision,
    clamp01,
    read_csv_rows,
    write_csv_rows,
)

ABLATION_RESULTS_PATH = PAPER_DIR / "ace-ablation-results.csv"
ABLATION_REPORT_PATH = PAPER_DIR / "ace-ablation-report.md"


def f(row: Mapping[str, str], column: str) -> float:
    return float(row.get(column, 0.0) or 0.0)


def full_score(row: Mapping[str, str]) -> float:
    return f(row, "ace_score")


def no_prodrug_triad(row: Mapping[str, str]) -> float:
    return clamp01(full_score(row) - 0.12 * f(row, "prodrug_triad"))


def no_cterm_proline_bonus(row: Mapping[str, str]) -> float:
    return clamp01(full_score(row) - 0.12 * f(row, "cterm_proline"))


def no_position_priors(row: Mapping[str, str]) -> float:
    removable = (
        0.16 * f(row, "cterm_anchor")
        + 0.12 * f(row, "cterm_proline")
        + 0.07 * f(row, "cterm_aromatic")
        + 0.12 * f(row, "penultimate_proline")
        + 0.08 * f(row, "antepenult_anchor")
        + 0.09 * f(row, "nterm_hydrophobic")
        + 0.08 * f(row, "second_basic_or_proline")
        + 0.07 * f(row, "third_proline")
        + 0.12 * f(row, "prodrug_triad")
    )
    return clamp01(full_score(row) - removable)


def composition_only(row: Mapping[str, str]) -> float:
    return clamp01(
        0.18 * f(row, "short_score")
        + 0.10 * f(row, "hydrophobic_balance")
        + 0.09 * f(row, "proline_balance")
        + 0.05 * f(row, "charge_balance")
        - 0.12 * f(row, "acidic_penalty")
        - 0.18 * f(row, "cysteine_penalty")
        - 0.10 * f(row, "extreme_charge_penalty")
    )


def random_hash_baseline(row: Mapping[str, str]) -> float:
    digest = hashlib.sha256(str(row["candidate_id"]).encode("utf-8")).hexdigest()
    return int(digest[:12], 16) / float(16**12 - 1)


ABLATIONS: dict[str, Callable[[Mapping[str, str]], float]] = {
    "full": full_score,
    "no_prodrug_triad": no_prodrug_triad,
    "no_cterm_proline_bonus": no_cterm_proline_bonus,
    "no_position_priors": no_position_priors,
    "composition_only": composition_only,
    "random_hash_baseline": random_hash_baseline,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranked", default=str(RANKED_CANDIDATES_PATH))
    parser.add_argument("--key", default=str(BLINDING_KEY_PATH))
    parser.add_argument("--out", default=str(ABLATION_RESULTS_PATH))
    parser.add_argument("--report-out", default=str(ABLATION_REPORT_PATH))
    return parser.parse_args()


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    ranked_rows = read_csv_rows(Path(args.ranked))
    key_rows = read_csv_rows(Path(args.key))
    key_by_id = {row["candidate_id"]: row for row in key_rows}
    total = len(ranked_rows)
    top_1pct_cutoff = max(1, math.ceil(total * 0.01))
    output_rows: list[dict[str, object]] = []
    report_rows: list[list[object]] = []

    for ablation_id, scorer in ABLATIONS.items():
        scored = []
        for row in ranked_rows:
            label = key_by_id[row["candidate_id"]]["label"]
            score = scorer(row)
            scored.append({**row, "ablation_id": ablation_id, "ablation_score": score, "label": label})
        scored.sort(key=lambda row: (-float(row["ablation_score"]), int(row["length"]), str(row["sequence"])))
        for index, row in enumerate(scored, start=1):
            row["ablation_rank"] = index

        positives = [row for row in scored if row["label"] in {"known_positive", "digestome_known_positive"}]
        negatives = [row for row in scored if row["label"] not in {"known_positive", "digestome_known_positive"}]
        matched_decoys = [row for row in scored if row["label"] == "matched_decoy"]
        positive_scores = [float(row["ablation_score"]) for row in positives]
        negative_scores = [float(row["ablation_score"]) for row in negatives]
        matched_scores = [float(row["ablation_score"]) for row in matched_decoys]
        top_counts = Counter()
        for row in positives:
            rank = int(row["ablation_rank"])
            if rank <= top_1pct_cutoff:
                top_counts["top_1pct"] += 1
            if rank <= 10:
                top_counts["top_10"] += 1
            if rank <= 25:
                top_counts["top_25"] += 1
            if rank <= 100:
                top_counts["top_100"] += 1
            output_rows.append(
                {
                    "ablation_id": ablation_id,
                    "sequence": row["sequence"],
                    "label": row["label"],
                    "ablation_rank": rank,
                    "ablation_score": round(float(row["ablation_score"]), 6),
                }
            )
        auc_all = auroc_from_scores(positive_scores, negative_scores)
        auc_matched = auroc_from_scores(positive_scores, matched_scores)
        ap_all = average_precision(
            [(1, float(row["ablation_score"])) for row in positives]
            + [(0, float(row["ablation_score"])) for row in negatives]
        )
        report_rows.append(
            [
                ablation_id,
                top_counts["top_1pct"],
                top_counts["top_10"],
                top_counts["top_25"],
                top_counts["top_100"],
                f"{mean([int(row['ablation_rank']) for row in positives]):.2f}",
                f"{auc_all:.4f}",
                f"{auc_matched:.4f}",
                f"{ap_all:.4f}",
            ]
        )

    write_csv_rows(
        Path(args.out),
        output_rows,
        ("ablation_id", "sequence", "label", "ablation_rank", "ablation_score"),
    )
    report = f"""# ACE Rediscovery Ablation Report

Date: {date.today().isoformat()}

This is a post-hoc sensitivity analysis over the same blinded ranking table. It asks whether known-positive recovery depends on specific ACE-position priors.

Total candidates: `{total}`  
Top 1 percent cutoff rank: `{top_1pct_cutoff}`

{markdown_table(["Ablation", "Top 1%", "Top 10", "Top 25", "Top 100", "Mean positive rank", "AUROC all", "AUROC matched", "AUPRC all"], report_rows)}

Interpretation:

- `full` is the registered pilot score.
- `no_prodrug_triad` removes the N-terminal hydrophobic/basic/proline pattern that helps `LKPNM`.
- `no_cterm_proline_bonus` removes only the terminal-proline bonus.
- `no_position_priors` removes most sequence-position priors and leaves mainly global composition/developability features.
- `composition_only` is a weaker baseline using length, hydrophobicity, proline fraction, charge, and penalties.
- `random_hash_baseline` is a deterministic random-ranking control.
"""
    Path(args.report_out).write_text(report, encoding="utf-8")
    print(f"wrote ablation report for {len(ABLATIONS)} scoring variants")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Unblind ACE rediscovery ranks and write a pilot report."""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import (
    BLINDING_KEY_PATH,
    FIGURE_DIR,
    RANKED_CANDIDATES_PATH,
    REPORT_PATH,
    RESULTS_PATH,
    average_precision,
    auroc_from_scores,
    mean_rank,
    read_csv_rows,
    save_json,
    write_csv_rows,
)


RESULT_COLUMNS = (
    "candidate_id",
    "rank",
    "sequence",
    "ace_score",
    "length",
    "label",
    "origin",
    "provenance",
    "top_1pct",
    "top_10",
    "top_25",
    "top_50",
    "top_100",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranked", default=str(RANKED_CANDIDATES_PATH))
    parser.add_argument("--key", default=str(BLINDING_KEY_PATH))
    parser.add_argument("--out", default=str(RESULTS_PATH))
    parser.add_argument("--report-out", default=str(REPORT_PATH))
    parser.add_argument("--figure-dir", default=str(FIGURE_DIR))
    parser.add_argument("--primary-min-top1pct-count", type=int, default=2)
    parser.add_argument("--primary-min-top1pct-fraction", type=float, default=0.0)
    parser.add_argument("--primary-min-auroc-matched", type=float, default=0.0)
    return parser.parse_args()


def pct_rank_threshold(total: int, percent: float) -> int:
    return max(1, math.ceil(total * percent / 100.0))


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def sparkline(counts: Counter[int], width: int = 80) -> str:
    if not counts:
        return ""
    max_count = max(counts.values())
    parts: list[str] = []
    for rank_bucket in sorted(counts):
        bar_width = max(1, round(counts[rank_bucket] / max_count * width))
        parts.append(f"{rank_bucket:>5}: {'#' * bar_width} {counts[rank_bucket]}")
    return "\n".join(parts)


def hypergeom_sf(k: int, population: int, positives: int, draws: int) -> float:
    """Return P(X >= k) for X following a hypergeometric distribution."""
    if k <= 0:
        return 1.0
    if population <= 0 or positives <= 0 or draws <= 0:
        return float("nan")
    upper = min(positives, draws)
    if k > upper:
        return 0.0
    denominator = (
        math.lgamma(population + 1)
        - math.lgamma(draws + 1)
        - math.lgamma(population - draws + 1)
    )
    terms: list[float] = []
    for observed in range(k, upper + 1):
        if positives - observed > population - draws:
            continue
        if draws - observed > population - positives:
            continue
        log_prob = (
            math.lgamma(positives + 1)
            - math.lgamma(observed + 1)
            - math.lgamma(positives - observed + 1)
            + math.lgamma(population - positives + 1)
            - math.lgamma(draws - observed + 1)
            - math.lgamma(population - positives - draws + observed + 1)
            - denominator
        )
        terms.append(log_prob)
    if not terms:
        return 0.0
    max_log = max(terms)
    return min(1.0, sum(math.exp(term - max_log) for term in terms) * math.exp(max_log))


def main() -> int:
    args = parse_args()
    ranked_rows = read_csv_rows(Path(args.ranked))
    key_rows = read_csv_rows(Path(args.key))
    key_by_id = {row["candidate_id"]: row for row in key_rows}
    total = len(ranked_rows)
    top_1pct_cutoff = pct_rank_threshold(total, 1.0)
    joined: list[dict[str, object]] = []
    for row in ranked_rows:
        key = key_by_id.get(row["candidate_id"], {})
        rank = int(row["rank"])
        joined.append(
            {
                "candidate_id": row["candidate_id"],
                "rank": rank,
                "sequence": row["sequence"],
                "ace_score": row["ace_score"],
                "length": row["length"],
                "label": key.get("label", "missing_key"),
                "origin": key.get("origin", ""),
                "provenance": key.get("provenance", ""),
                "top_1pct": int(rank <= top_1pct_cutoff),
                "top_10": int(rank <= 10),
                "top_25": int(rank <= 25),
                "top_50": int(rank <= 50),
                "top_100": int(rank <= 100),
            }
        )
    write_csv_rows(Path(args.out), joined, RESULT_COLUMNS)

    positives = [row for row in joined if row["label"] in {"known_positive", "digestome_known_positive"}]
    negatives = [row for row in joined if row["label"] not in {"known_positive", "digestome_known_positive"}]
    matched_decoys = [row for row in joined if row["label"] == "matched_decoy"]
    positive_scores = [float(row["ace_score"]) for row in positives]
    negative_scores = [float(row["ace_score"]) for row in negatives]
    matched_scores = [float(row["ace_score"]) for row in matched_decoys]
    top_counts = {
        "top_1pct": sum(int(row["top_1pct"]) for row in positives),
        "top_10": sum(int(row["top_10"]) for row in positives),
        "top_25": sum(int(row["top_25"]) for row in positives),
        "top_50": sum(int(row["top_50"]) for row in positives),
        "top_100": sum(int(row["top_100"]) for row in positives),
    }
    auc_all = auroc_from_scores(positive_scores, negative_scores)
    auc_matched = auroc_from_scores(positive_scores, matched_scores)
    ap_all = average_precision([(1, float(row["ace_score"])) for row in positives] + [(0, float(row["ace_score"])) for row in negatives])
    top_1pct_fraction = top_counts["top_1pct"] / len(positives) if positives else float("nan")
    random_expected_top_1pct = len(positives) * top_1pct_cutoff / total if total else float("nan")
    top_1pct_enrichment = top_counts["top_1pct"] / max(1e-12, random_expected_top_1pct)
    top_1pct_p = hypergeom_sf(top_counts["top_1pct"], total, len(positives), top_1pct_cutoff)
    success = (
        top_counts["top_1pct"] >= args.primary_min_top1pct_count
        and top_1pct_fraction >= args.primary_min_top1pct_fraction
        and (math.isnan(auc_matched) or auc_matched >= args.primary_min_auroc_matched)
    )
    label_counts = Counter(str(row["label"]) for row in joined)

    positive_rows = [
        [
            row["sequence"],
            row["rank"],
            row["ace_score"],
            row["length"],
            row["label"],
            row["origin"],
        ]
        for row in sorted(positives, key=lambda item: int(item["rank"]))
    ]
    top_rows = [
        [
            row["rank"],
            row["sequence"],
            row["ace_score"],
            row["length"],
            row["label"],
        ]
        for row in joined[:25]
    ]

    rank_buckets: Counter[int] = Counter()
    for row in positives:
        bucket = (int(row["rank"]) - 1) // 100 + 1
        rank_buckets[bucket * 100] += 1

    figure_dir = Path(args.figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    (figure_dir / "positive-rank-buckets.txt").write_text(sparkline(rank_buckets) + "\n", encoding="utf-8")
    save_json(
        figure_dir / "ace-rediscovery-summary.json",
        {
            "date": date.today().isoformat(),
            "total_candidates": total,
            "top_1pct_cutoff_rank": top_1pct_cutoff,
            "label_counts": dict(sorted(label_counts.items())),
            "known_positive_count": len(positives),
            "top_counts": top_counts,
            "top_1pct_fraction": top_1pct_fraction,
            "top_1pct_enrichment_over_random": top_1pct_enrichment,
            "top_1pct_hypergeom_p": top_1pct_p,
            "mean_positive_rank": mean_rank([int(row["rank"]) for row in positives]),
            "mean_matched_decoy_rank": mean_rank([int(row["rank"]) for row in matched_decoys]),
            "auroc_vs_all_negatives": auc_all,
            "auroc_vs_matched_decoys": auc_matched,
            "average_precision_vs_all_negatives": ap_all,
            "primary_success": success,
        },
    )

    report = f"""# ACE Rediscovery Pilot Report

Date: {date.today().isoformat()}

## Claim Boundary

This is a computational pilot. It tests whether a sequence-only ACE-specific PepLab score can recover known literature-positive ACE inhibitory peptides after blinding. It does not demonstrate new wet-lab activity or clinical efficacy.

## Primary Endpoint

Success criterion: all configured gates below must pass.

Configured gates:

- Minimum known positives in top 1 percent: `{args.primary_min_top1pct_count}`
- Minimum top 1 percent recovery fraction: `{args.primary_min_top1pct_fraction}`
- Minimum AUROC versus matched decoys: `{args.primary_min_auroc_matched}`

- Total candidates: `{total}`
- Top 1 percent cutoff rank: `{top_1pct_cutoff}`
- Known positives: `{len(positives)}`
- Known positives in top 1 percent: `{top_counts["top_1pct"]}`
- Top 1 percent recovery fraction: `{top_1pct_fraction:.4f}`
- Top 1 percent enrichment over random: `{top_1pct_enrichment:.2f}x`
- Top 1 percent hypergeometric p-value: `{top_1pct_p:.3g}`
- Primary endpoint passed: `{"yes" if success else "no"}`

## Recovery Metrics

- Top 10 recall: `{top_counts["top_10"]}/{len(positives)}`
- Top 25 recall: `{top_counts["top_25"]}/{len(positives)}`
- Top 50 recall: `{top_counts["top_50"]}/{len(positives)}`
- Top 100 recall: `{top_counts["top_100"]}/{len(positives)}`
- Mean positive rank: `{mean_rank([int(row["rank"]) for row in positives]):.2f}`
- Mean matched-decoy rank: `{mean_rank([int(row["rank"]) for row in matched_decoys]):.2f}`
- AUROC versus all negatives: `{auc_all:.4f}`
- AUROC versus matched decoys: `{auc_matched:.4f}`
- Average precision versus all negatives: `{ap_all:.4f}`

## Label Counts

{markdown_table(["Label", "Count"], [[label, count] for label, count in sorted(label_counts.items())])}

## Known Positive Ranks

{markdown_table(["Sequence", "Rank", "ACE score", "Length", "Label", "Origin"], positive_rows)}

## Top 25 After Unblinding

{markdown_table(["Rank", "Sequence", "ACE score", "Length", "Label"], top_rows)}

## Interpretation

If the primary endpoint passed, the defensible claim is:

> PepLab rediscovered literature-validated ACE-inhibitory peptides from a blinded sequence search space in a computational pilot.

Do not claim novel bioactivity. The next gate is a larger universe with stricter source-protein split controls and a no-motif-prior ablation.
"""
    Path(args.report_out).write_text(report, encoding="utf-8")
    print(
        f"unblinded {total} candidates; positives in top 1%: "
        f"{top_counts['top_1pct']}/{len(positives)}; primary success={success}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

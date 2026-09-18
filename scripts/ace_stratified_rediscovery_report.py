#!/usr/bin/env python3
"""Write subset metrics for an ACE rediscovery run."""

from __future__ import annotations

import argparse
import re
import sys
from datetime import date
from pathlib import Path
from statistics import median
from typing import Callable, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import PAPER_DIR, read_csv_rows, write_csv_rows  # noqa: E402
from scripts.ace_unblind_rediscovery_report import hypergeom_sf  # noqa: E402

DEFAULT_PANEL = PAPER_DIR / "known-positive-expanded-panel.csv"
DEFAULT_RESULTS = PAPER_DIR / "ace-rediscovery-results.csv"
DEFAULT_OUT = PAPER_DIR / "ace-stratified-results.csv"
DEFAULT_REPORT = PAPER_DIR / "ace-stratified-report.md"
PRIMARY_ANCHORS = {"VPP", "IPP", "LKPNM", "LKP"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", default=str(DEFAULT_RESULTS))
    parser.add_argument("--panel", default=str(DEFAULT_PANEL))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--report-out", default=str(DEFAULT_REPORT))
    return parser.parse_args()


def reported_ic50_um(panel_row: Mapping[str, str]) -> float | None:
    notes = panel_row.get("notes", "")
    match = re.search(r"IC50\s+([0-9.]+)\s+uM", notes)
    if not match:
        return None
    return float(match.group(1))


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def summarize_subset(
    subset_id: str,
    rows: list[Mapping[str, object]],
    total_candidates: int,
    top_1pct_cutoff: int,
) -> dict[str, object]:
    count = len(rows)
    ranks = [int(row["rank"]) for row in rows]
    top_1pct = sum(int(row["top_1pct"]) for row in rows)
    top_100 = sum(int(row["top_100"]) for row in rows)
    expected = count * top_1pct_cutoff / total_candidates if total_candidates else 0.0
    enrichment = top_1pct / expected if expected else 0.0
    return {
        "subset_id": subset_id,
        "positive_count": count,
        "top_1pct_count": top_1pct,
        "top_1pct_fraction": round(top_1pct / count, 6) if count else "",
        "top_1pct_enrichment_over_random": round(enrichment, 6) if count else "",
        "top_1pct_hypergeom_p": hypergeom_sf(top_1pct, total_candidates, count, top_1pct_cutoff) if count else "",
        "top_100_count": top_100,
        "mean_rank": round(sum(ranks) / count, 3) if count else "",
        "median_rank": median(ranks) if count else "",
        "best_rank": min(ranks) if count else "",
    }


def main() -> int:
    args = parse_args()
    result_rows = read_csv_rows(Path(args.results))
    panel_rows = read_csv_rows(Path(args.panel))
    panel_by_sequence = {row["sequence"]: row for row in panel_rows}
    positives = [
        {**row, "panel": panel_by_sequence.get(row["sequence"], {})}
        for row in result_rows
        if row["label"] in {"known_positive", "digestome_known_positive"}
    ]
    total_candidates = len(result_rows)
    top_1pct_cutoff = max(1, round(total_candidates * 0.01 + 0.499999))

    def p(row: Mapping[str, object]) -> Mapping[str, str]:
        return row["panel"] if isinstance(row["panel"], dict) else {}

    subset_predicates: dict[str, Callable[[Mapping[str, object]], bool]] = {
        "all_literature_positives": lambda row: True,
        "primary_four_anchors": lambda row: str(row["sequence"]) in PRIMARY_ANCHORS,
        "reported_antihypertensive": lambda row: "antihypertensive" in p(row).get("activity_class", ""),
        "reported_ic50_uM_le_30": lambda row: (reported_ic50_um(p(row)) or float("inf")) <= 30.0,
        "reported_ic50_uM_le_10": lambda row: (reported_ic50_um(p(row)) or float("inf")) <= 10.0,
        "short_len_le_5": lambda row: int(row["length"]) <= 5,
        "short_len_le_5_ic50_uM_le_30": lambda row: int(row["length"]) <= 5
        and (reported_ic50_um(p(row)) or float("inf")) <= 30.0,
        "naturally_released_digestome_positive": lambda row: row["label"] == "digestome_known_positive",
    }

    summary_rows = [
        summarize_subset(
            subset_id,
            [row for row in positives if predicate(row)],
            total_candidates,
            top_1pct_cutoff,
        )
        for subset_id, predicate in subset_predicates.items()
    ]
    write_csv_rows(
        Path(args.out),
        summary_rows,
        (
            "subset_id",
            "positive_count",
            "top_1pct_count",
            "top_1pct_fraction",
            "top_1pct_enrichment_over_random",
            "top_1pct_hypergeom_p",
            "top_100_count",
            "mean_rank",
            "median_rank",
            "best_rank",
        ),
    )

    table_rows = [
        [
            row["subset_id"],
            row["positive_count"],
            row["top_1pct_count"],
            row["top_1pct_fraction"],
            row["top_1pct_enrichment_over_random"],
            f"{float(row['top_1pct_hypergeom_p']):.3g}" if row["top_1pct_hypergeom_p"] != "" else "",
            row["top_100_count"],
            row["mean_rank"],
            row["median_rank"],
            row["best_rank"],
        ]
        for row in summary_rows
    ]
    report = f"""# ACE Rediscovery Stratified Report

Date: {date.today().isoformat()}

This report stratifies the unblinded known positives because the expanded literature panel contains heterogeneous evidence: short tripeptides, longer hydrolysate-derived peptides, reported antihypertensive peptides, and peptides with IC50 values reported in different units.

Total candidates: `{total_candidates}`  
Known positives: `{len(positives)}`  
Top 1 percent cutoff rank: `{top_1pct_cutoff}`

{markdown_table(["Subset", "N", "Top 1%", "Top 1% frac", "Enrichment", "Hypergeom p", "Top 100", "Mean rank", "Median rank", "Best rank"], table_rows)}

## Interpretation

The cleanest near-term claim should focus on recovery and enrichment of short ACE-like literature positives and the four primary anchors. The broader 81-sequence panel is valuable as a stress test, but it includes long and weak positives that the current sequence-only score does not yet model well.
"""
    Path(args.report_out).write_text(report, encoding="utf-8")
    print(f"wrote stratified report for {len(summary_rows)} subsets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

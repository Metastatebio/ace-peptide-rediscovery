#!/usr/bin/env python3
"""Summarize where an ACE holdout score-mode sweep succeeds or fails."""

from __future__ import annotations

import argparse
import math
import re
import statistics
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import PAPER_DIR, read_csv_rows, write_csv_rows  # noqa: E402

DEFAULT_SWEEP_DIR = PAPER_DIR / "experiments" / "2026-08-29-post-freeze-v3-primary-panel-001" / "score-mode-sweep"
DEFAULT_PANEL = PAPER_DIR / "known-positive-postfreeze-v3-primary-panel.csv"
DEFAULT_OUT = PAPER_DIR / "experiments" / "2026-08-29-post-freeze-v3-primary-panel-001" / "v3-failure-analysis.csv"
DEFAULT_REPORT = PAPER_DIR / "experiments" / "2026-08-29-post-freeze-v3-primary-panel-001" / "v3-failure-analysis.md"

OUTPUT_FIELDS = (
    "score_mode",
    "group_id",
    "positive_count",
    "top_1pct_count",
    "top_1pct_fraction",
    "top_10pct_count",
    "top_10pct_fraction",
    "mean_rank",
    "median_rank",
    "best_rank",
)

DISPLAY_GROUPS = (
    "all",
    "length_le_5",
    "length_le_7",
    "has_ic50_note",
    "ic50_explicit_uM",
    "ic50_uM_le_10",
    "ic50_uM_le_100",
    "percent_potency_rows",
    "potency_pct_ge_50",
    "potency_pct_ge_80",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep-dir", default=str(DEFAULT_SWEEP_DIR))
    parser.add_argument("--panel", default=str(DEFAULT_PANEL))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--report-out", default=str(DEFAULT_REPORT))
    return parser.parse_args()


def explicit_ic50_um(notes: str) -> float | None:
    match = re.search(r"IC50 (?:greater than )?([0-9.]+) uM", notes)
    if not match:
        return None
    return float(match.group(1))


def inhibitory_potency_percent(notes: str) -> float | None:
    match = re.search(r"inhibitory potency ([0-9.]+) percent", notes)
    if not match:
        return None
    return float(match.group(1))


def groups_for_positive(result_row: Mapping[str, str], panel_row: Mapping[str, str]) -> list[str]:
    sequence = result_row["sequence"]
    notes = panel_row.get("notes", "")
    groups = ["all", f"source:{panel_row.get('source_material', 'unknown')}"]
    if len(sequence) <= 5:
        groups.append("length_le_5")
    if len(sequence) <= 7:
        groups.append("length_le_7")
    if "IC50" in notes:
        groups.append("has_ic50_note")
    ic50 = explicit_ic50_um(notes)
    if ic50 is not None:
        groups.append("ic50_explicit_uM")
        if ic50 <= 10.0:
            groups.append("ic50_uM_le_10")
        if ic50 <= 100.0:
            groups.append("ic50_uM_le_100")
    potency = inhibitory_potency_percent(notes)
    if potency is not None:
        groups.append("percent_potency_rows")
        if potency >= 50.0:
            groups.append("potency_pct_ge_50")
        if potency >= 80.0:
            groups.append("potency_pct_ge_80")
    return groups


def summarize_group(score_mode: str, group_id: str, ranks: list[int], total_candidates: int) -> dict[str, object]:
    top_1pct_cutoff = max(1, math.ceil(total_candidates * 0.01))
    top_10pct_cutoff = max(1, math.ceil(total_candidates * 0.10))
    count = len(ranks)
    top_1pct = sum(rank <= top_1pct_cutoff for rank in ranks)
    top_10pct = sum(rank <= top_10pct_cutoff for rank in ranks)
    return {
        "score_mode": score_mode,
        "group_id": group_id,
        "positive_count": count,
        "top_1pct_count": top_1pct,
        "top_1pct_fraction": round(top_1pct / count, 6) if count else "",
        "top_10pct_count": top_10pct,
        "top_10pct_fraction": round(top_10pct / count, 6) if count else "",
        "mean_rank": round(sum(ranks) / count, 3) if count else "",
        "median_rank": statistics.median(ranks) if count else "",
        "best_rank": min(ranks) if count else "",
    }


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
    sweep_dir = Path(args.sweep_dir)
    panel_rows = read_csv_rows(Path(args.panel))
    panel_by_sequence = {row["sequence"]: row for row in panel_rows}
    mode_rows = read_csv_rows(sweep_dir / "score-mode-summary.csv")
    output_rows: list[dict[str, object]] = []

    for mode_row in mode_rows:
        mode = mode_row["score_mode"]
        result_path = sweep_dir / mode / "ace-rediscovery-results.csv"
        result_rows = read_csv_rows(result_path)
        total_candidates = len(result_rows)
        group_ranks: dict[str, list[int]] = defaultdict(list)
        for row in result_rows:
            if row.get("label") not in {"known_positive", "digestome_known_positive"}:
                continue
            panel_row = panel_by_sequence.get(row["sequence"], {})
            for group_id in groups_for_positive(row, panel_row):
                group_ranks[group_id].append(int(row["rank"]))
        for group_id, ranks in sorted(group_ranks.items()):
            output_rows.append(summarize_group(mode, group_id, ranks, total_candidates))

    write_csv_rows(Path(args.out), output_rows, OUTPUT_FIELDS)

    report_lines = [
        "# ACE Holdout Failure Analysis",
        "",
        f"Date: {date.today().isoformat()}",
        "",
        f"Sweep directory: `{sweep_dir}`",
        f"Panel: `{args.panel}`",
        "",
        "This report is diagnostic. It identifies which holdout strata the frozen scorer handled or missed; it is not a new validation claim.",
        "",
        "## Main Strata",
        "",
    ]
    for mode in [row["score_mode"] for row in mode_rows]:
        selected = [
            row
            for row in output_rows
            if row["score_mode"] == mode and row["group_id"] in DISPLAY_GROUPS
        ]
        table_rows = [
            [
                row["group_id"],
                row["positive_count"],
                row["top_1pct_count"],
                row["top_10pct_count"],
                row["mean_rank"],
                row["median_rank"],
                row["best_rank"],
            ]
            for row in selected
        ]
        report_lines.extend(
            [
                f"### {mode}",
                "",
                markdown_table(
                    ["Group", "N", "Top 1%", "Top 10%", "Mean rank", "Median rank", "Best rank"],
                    table_rows,
                ),
                "",
            ]
        )
    report_lines.extend(
        [
            "## Readout",
            "",
            "The holdout problem is top-rank resolution, not total absence of signal. Compare target-informed modes against target-agnostic controls to identify which strata cross the top-1-percent threshold and which remain broad-shortlist-only signals.",
            "",
            "For the next model iteration, use this file to choose pre-registered strata before curating a new validation panel.",
        ]
    )
    Path(args.report_out).write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(args.report_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

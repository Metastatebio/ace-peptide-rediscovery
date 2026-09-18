#!/usr/bin/env python3
"""Run multiple ACE scoring modes and compare rediscovery metrics."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import (  # noqa: E402
    BLINDING_KEY_PATH,
    BLINDED_UNIVERSE_PATH,
    PAPER_DIR,
    SCORING_MODE_DESCRIPTIONS,
    write_csv_rows,
)

DEFAULT_OUT_DIR = PAPER_DIR / "experiments" / "2026-08-29-score-mode-sweep-001"
DEFAULT_PANEL = PAPER_DIR / "known-positive-expanded-panel.csv"

SUMMARY_FIELDS = (
    "score_mode",
    "claim_tier",
    "total_candidates",
    "known_positive_count",
    "top_1pct_count",
    "top_1pct_fraction",
    "top_1pct_enrichment_over_random",
    "top_1pct_hypergeom_p",
    "top_10",
    "top_25",
    "top_100",
    "mean_positive_rank",
    "mean_matched_decoy_rank",
    "auroc_vs_all_negatives",
    "auroc_vs_matched_decoys",
    "average_precision_vs_all_negatives",
    "primary_success",
    "description",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universe", default=str(BLINDED_UNIVERSE_PATH))
    parser.add_argument("--key", default=str(BLINDING_KEY_PATH))
    parser.add_argument("--panel", default=str(DEFAULT_PANEL))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument(
        "--modes",
        nargs="+",
        default=[
            "generic_zero_shot",
            "oral_stability_zero_shot",
            "no_position_priors",
            "ace_no_exact_proline_motif",
            "ace_pharmacophore_zero_label",
            "literature_sar",
        ],
        choices=sorted(SCORING_MODE_DESCRIPTIONS),
    )
    parser.add_argument("--primary-min-top1pct-count", type=int, default=8)
    parser.add_argument("--primary-min-top1pct-fraction", type=float, default=0.10)
    parser.add_argument("--primary-min-auroc-matched", type=float, default=0.75)
    return parser.parse_args()


def claim_tier(mode: str) -> str:
    if mode == "generic_zero_shot":
        return "target_agnostic_zero_shot"
    if mode == "oral_stability_zero_shot":
        return "target_agnostic_developability_zero_shot"
    if mode == "no_position_priors":
        return "negative_control"
    if mode == "ace_pharmacophore_zero_label":
        return "target_informed_label_free"
    if mode == "ace_no_exact_proline_motif":
        return "target_informed_no_exact_terminal_proline_motif"
    if mode == "literature_sar":
        return "ace_sar_label_free_but_literature_informed"
    return "unknown"


def run_command(args: list[str]) -> None:
    subprocess.run(args, cwd=REPO_ROOT, check=True)


def load_summary(mode_dir: Path) -> dict[str, object]:
    path = mode_dir / "figures" / "ace-rediscovery-summary.json"
    return json.loads(path.read_text(encoding="utf-8"))


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
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []

    for mode in args.modes:
        mode_dir = out_dir / mode
        mode_dir.mkdir(parents=True, exist_ok=True)
        ranked = mode_dir / "ace-ranked-candidates.csv"
        scoring_manifest = mode_dir / "ace-scoring-manifest.json"
        results = mode_dir / "ace-rediscovery-results.csv"
        report = mode_dir / "ace-rediscovery-report.md"
        figures = mode_dir / "figures"
        stratified_results = mode_dir / "ace-stratified-results.csv"
        stratified_report = mode_dir / "ace-stratified-report.md"

        run_command(
            [
                sys.executable,
                "scripts/ace_score_blinded_universe.py",
                "--universe",
                args.universe,
                "--out",
                str(ranked),
                "--manifest-out",
                str(scoring_manifest),
                "--mode",
                mode,
            ]
        )
        run_command(
            [
                sys.executable,
                "scripts/ace_unblind_rediscovery_report.py",
                "--ranked",
                str(ranked),
                "--key",
                args.key,
                "--out",
                str(results),
                "--report-out",
                str(report),
                "--figure-dir",
                str(figures),
                "--primary-min-top1pct-count",
                str(args.primary_min_top1pct_count),
                "--primary-min-top1pct-fraction",
                str(args.primary_min_top1pct_fraction),
                "--primary-min-auroc-matched",
                str(args.primary_min_auroc_matched),
            ]
        )
        run_command(
            [
                sys.executable,
                "scripts/ace_stratified_rediscovery_report.py",
                "--results",
                str(results),
                "--panel",
                args.panel,
                "--out",
                str(stratified_results),
                "--report-out",
                str(stratified_report),
            ]
        )
        summary = load_summary(mode_dir)
        top_counts = summary.get("top_counts", {})
        if not isinstance(top_counts, dict):
            top_counts = {}
        summary_rows.append(
            {
                "score_mode": mode,
                "claim_tier": claim_tier(mode),
                "total_candidates": summary["total_candidates"],
                "known_positive_count": summary["known_positive_count"],
                "top_1pct_count": top_counts.get("top_1pct", ""),
                "top_1pct_fraction": summary["top_1pct_fraction"],
                "top_1pct_enrichment_over_random": summary["top_1pct_enrichment_over_random"],
                "top_1pct_hypergeom_p": summary["top_1pct_hypergeom_p"],
                "top_10": top_counts.get("top_10", ""),
                "top_25": top_counts.get("top_25", ""),
                "top_100": top_counts.get("top_100", ""),
                "mean_positive_rank": summary["mean_positive_rank"],
                "mean_matched_decoy_rank": summary["mean_matched_decoy_rank"],
                "auroc_vs_all_negatives": summary["auroc_vs_all_negatives"],
                "auroc_vs_matched_decoys": summary["auroc_vs_matched_decoys"],
                "average_precision_vs_all_negatives": summary["average_precision_vs_all_negatives"],
                "primary_success": summary["primary_success"],
                "description": SCORING_MODE_DESCRIPTIONS[mode],
            }
        )

    write_csv_rows(out_dir / "score-mode-summary.csv", summary_rows, SUMMARY_FIELDS)
    table_rows = [
        [
            row["score_mode"],
            row["claim_tier"],
            row["top_1pct_count"],
            f"{float(row['top_1pct_fraction']):.3f}",
            f"{float(row['top_1pct_enrichment_over_random']):.2f}x",
            f"{float(row['auroc_vs_matched_decoys']):.4f}",
            row["primary_success"],
        ]
        for row in summary_rows
    ]
    report = f"""# ACE Score Mode Sweep

Date: {date.today().isoformat()}

Universe: `{args.universe}`  
Key: `{args.key}`  
Modes: `{", ".join(args.modes)}`

## Summary

{markdown_table(["Mode", "Claim tier", "Top 1%", "Top 1% frac", "Enrichment", "AUROC matched", "Gate pass"], table_rows)}

## Interpretation

- `generic_zero_shot` and `oral_stability_zero_shot` are target-agnostic rediscovery settings.
- `ace_pharmacophore_zero_label` and `ace_no_exact_proline_motif` are cleaner than training on ACE literature labels, but they are still ACE-target-informed.
- `literature_sar` is a useful upper-bound heuristic, not a clean discovery claim.
- A publishable rediscovery claim should state exactly which tier passed.
"""
    (out_dir / "score-mode-sweep-report.md").write_text(report, encoding="utf-8")
    print(out_dir / "score-mode-sweep-report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

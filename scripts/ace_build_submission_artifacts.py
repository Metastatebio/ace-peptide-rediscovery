#!/usr/bin/env python3
"""Assemble submission-facing ACE tables and SVG figures from frozen artifacts."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import PAPER_DIR, read_csv_rows, write_csv_rows  # noqa: E402
from scripts.ace_unblind_rediscovery_report import hypergeom_sf  # noqa: E402

TABLE_DIR = PAPER_DIR / "tables"
FIGURE_DIR = PAPER_DIR / "figures" / "submission"

VALIDATION_RUNS = (
    {
        "run_id": "v4_early",
        "panel_id": "v4",
        "display_name": "v4 early",
        "summary_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v4-independent-validation-001"
        / "frozen-rerank"
        / "early-enrichment-v1"
        / "frozen-rerank-summary.csv",
        "positive_rank_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v4-independent-validation-001"
        / "frozen-rerank"
        / "early-enrichment-v1"
        / "frozen-rerank-positive-ranks.csv",
        "ranked_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v4-independent-validation-001"
        / "frozen-rerank"
        / "early-enrichment-v1"
        / "frozen-rerank-ranked-candidates.csv",
        "recipe_role": "aggressive early enrichment",
    },
    {
        "run_id": "v4_broad",
        "panel_id": "v4",
        "display_name": "v4 broad",
        "summary_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v4-independent-validation-001"
        / "frozen-rerank"
        / "broad-shortlist-v1"
        / "frozen-rerank-summary.csv",
        "positive_rank_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v4-independent-validation-001"
        / "frozen-rerank"
        / "broad-shortlist-v1"
        / "frozen-rerank-positive-ranks.csv",
        "ranked_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v4-independent-validation-001"
        / "frozen-rerank"
        / "broad-shortlist-v1"
        / "frozen-rerank-ranked-candidates.csv",
        "recipe_role": "broad shortlist",
    },
    {
        "run_id": "v5_broad",
        "panel_id": "v5",
        "display_name": "v5 broad",
        "summary_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v5-failure-mode-validation-001"
        / "frozen-rerank"
        / "broad-shortlist-v1"
        / "frozen-rerank-summary.csv",
        "positive_rank_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v5-failure-mode-validation-001"
        / "frozen-rerank"
        / "broad-shortlist-v1"
        / "frozen-rerank-positive-ranks.csv",
        "ranked_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v5-failure-mode-validation-001"
        / "frozen-rerank"
        / "broad-shortlist-v1"
        / "frozen-rerank-ranked-candidates.csv",
        "recipe_role": "v4 broad recipe on hard stress panel",
    },
    {
        "run_id": "v6_broad",
        "panel_id": "v6",
        "display_name": "v6 broad",
        "summary_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v6-independent-validation-001"
        / "frozen-rerank"
        / "broad-v6-candidate"
        / "frozen-rerank-summary.csv",
        "positive_rank_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v6-independent-validation-001"
        / "frozen-rerank"
        / "broad-v6-candidate"
        / "frozen-rerank-positive-ranks.csv",
        "ranked_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v6-independent-validation-001"
        / "frozen-rerank"
        / "broad-v6-candidate"
        / "frozen-rerank-ranked-candidates.csv",
        "recipe_role": "post-v5 broad recipe",
    },
    {
        "run_id": "v6_balanced",
        "panel_id": "v6",
        "display_name": "v6 balanced",
        "summary_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v6-independent-validation-001"
        / "frozen-rerank"
        / "balanced-v6-candidate"
        / "frozen-rerank-summary.csv",
        "positive_rank_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v6-independent-validation-001"
        / "frozen-rerank"
        / "balanced-v6-candidate"
        / "frozen-rerank-positive-ranks.csv",
        "ranked_path": PAPER_DIR
        / "experiments"
        / "2026-08-29-v6-independent-validation-001"
        / "frozen-rerank"
        / "balanced-v6-candidate"
        / "frozen-rerank-ranked-candidates.csv",
        "recipe_role": "post-v5 balanced recipe",
    },
)

VALIDATION_FIELDS = (
    "run_id",
    "panel_id",
    "display_name",
    "recipe_role",
    "recipe_id",
    "total_candidates",
    "known_positive_count",
    "top_1pct_count",
    "top_1pct_fraction",
    "top_1pct_hypergeom_p",
    "top_5pct_count",
    "top_5pct_fraction",
    "top_5pct_hypergeom_p",
    "top_10pct_count",
    "top_10pct_fraction",
    "top_10pct_hypergeom_p",
    "top_100_count",
    "best_positive_rank",
    "median_positive_rank",
    "mean_positive_rank",
    "auroc_vs_matched_decoys",
    "auroc_vs_all_negatives",
)

PANEL_FIELDS = (
    "panel_id",
    "panel_path",
    "peptide_count",
    "reference_count",
    "source_material_count",
    "length_le_5_count",
    "length_gt_5_count",
    "cterminal_proline_count",
    "non_cterminal_proline_count",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-dir", default=str(TABLE_DIR))
    parser.add_argument("--figure-dir", default=str(FIGURE_DIR))
    return parser.parse_args()


def read_single_row(path: Path) -> dict[str, str]:
    rows = read_csv_rows(path)
    if len(rows) != 1:
        raise ValueError(f"{path} expected one row, found {len(rows)}")
    return rows[0]


def count_data_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def positive_ranks(path: Path) -> list[int]:
    return [int(row["rank"]) for row in read_csv_rows(path)]


def validation_rows() -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for run in VALIDATION_RUNS:
        summary = read_single_row(run["summary_path"])  # type: ignore[arg-type]
        ranked_path = run["ranked_path"]  # type: ignore[assignment]
        ranks = positive_ranks(run["positive_rank_path"])  # type: ignore[arg-type]
        total_candidates = count_data_rows(ranked_path)  # type: ignore[arg-type]
        known_positive_count = len(ranks)
        top_5pct_cutoff = max(1, math.ceil(total_candidates * 0.05))
        top_10pct_cutoff = max(1, math.ceil(total_candidates * 0.10))
        top_5pct_count = sum(rank <= top_5pct_cutoff for rank in ranks)
        top_10pct_count = sum(rank <= top_10pct_cutoff for rank in ranks)
        output.append(
            {
                "run_id": run["run_id"],
                "panel_id": run["panel_id"],
                "display_name": run["display_name"],
                "recipe_role": run["recipe_role"],
                "recipe_id": summary["recipe_id"],
                "total_candidates": total_candidates,
                "known_positive_count": known_positive_count,
                "top_1pct_count": summary["top_1pct_count"],
                "top_1pct_fraction": summary["top_1pct_fraction"],
                "top_1pct_hypergeom_p": summary["top_1pct_hypergeom_p"],
                "top_5pct_count": top_5pct_count,
                "top_5pct_fraction": top_5pct_count / known_positive_count if known_positive_count else "",
                "top_5pct_hypergeom_p": hypergeom_sf(
                    top_5pct_count, total_candidates, known_positive_count, top_5pct_cutoff
                ),
                "top_10pct_count": top_10pct_count,
                "top_10pct_fraction": top_10pct_count / known_positive_count if known_positive_count else "",
                "top_10pct_hypergeom_p": hypergeom_sf(
                    top_10pct_count, total_candidates, known_positive_count, top_10pct_cutoff
                ),
                "top_100_count": summary["top_100_count"],
                "best_positive_rank": summary["best_positive_rank"],
                "median_positive_rank": summary["median_positive_rank"],
                "mean_positive_rank": summary["mean_positive_rank"],
                "auroc_vs_matched_decoys": summary["auroc_vs_matched_decoys"],
                "auroc_vs_all_negatives": summary["auroc_vs_all_negatives"],
            }
        )
    return output


def panel_rows() -> list[dict[str, object]]:
    panels = (
        ("v3", PAPER_DIR / "known-positive-postfreeze-v3-primary-panel.csv"),
        ("v4", PAPER_DIR / "known-positive-v4-independent-primary-panel.csv"),
        ("v5", PAPER_DIR / "known-positive-v5-failure-mode-primary-panel.csv"),
        ("v6", PAPER_DIR / "known-positive-v6-independent-primary-panel.csv"),
    )
    rows: list[dict[str, object]] = []
    for panel_id, path in panels:
        panel = read_csv_rows(path)
        sequences = [row["sequence"] for row in panel]
        rows.append(
            {
                "panel_id": panel_id,
                "panel_path": str(path),
                "peptide_count": len(sequences),
                "reference_count": len({row["reference_url"] for row in panel}),
                "source_material_count": len({row["source_material"] for row in panel}),
                "length_le_5_count": sum(len(sequence) <= 5 for sequence in sequences),
                "length_gt_5_count": sum(len(sequence) > 5 for sequence in sequences),
                "cterminal_proline_count": sum(sequence.endswith("P") for sequence in sequences),
                "non_cterminal_proline_count": sum(not sequence.endswith("P") for sequence in sequences),
            }
        )
    return rows


def external_baseline_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for panel_id, experiment_id in (
        ("v4", "2026-08-29-v4-independent-validation-001"),
        ("v5", "2026-08-29-v5-failure-mode-validation-001"),
        ("v6", "2026-08-29-v6-independent-validation-001"),
    ):
        path = PAPER_DIR / "experiments" / experiment_id / "external-baselines" / "peptideranker-style" / "peptideranker-style-summary.csv"
        if not path.exists():
            continue
        row = read_single_row(path)
        row["panel_id"] = panel_id
        rows.append(row)
    return rows


def svg_grouped_bars(
    *,
    title: str,
    labels: Sequence[str],
    series: Sequence[tuple[str, Sequence[float], str]],
    y_max: float,
    y_label: str,
    width: int = 1100,
    height: int = 560,
) -> str:
    margin_left = 82
    margin_right = 36
    margin_top = 76
    margin_bottom = 92
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    group_width = plot_width / max(1, len(labels))
    bar_gap = 5
    bar_width = max(10, (group_width - 28 - bar_gap * (len(series) - 1)) / max(1, len(series)))
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{margin_left}" y="38" font-family="Arial, sans-serif" font-size="24" font-weight="700" fill="#18212b">{title}</text>',
        f'<text x="24" y="{margin_top + plot_height / 2:.1f}" transform="rotate(-90 24 {margin_top + plot_height / 2:.1f})" font-family="Arial, sans-serif" font-size="14" fill="#394957">{y_label}</text>',
    ]
    for tick in range(0, 6):
        value = y_max * tick / 5
        y = margin_top + plot_height - plot_height * value / y_max
        lines.append(f'<line x1="{margin_left}" y1="{y:.1f}" x2="{width - margin_right}" y2="{y:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        lines.append(f'<text x="{margin_left - 10}" y="{y + 4:.1f}" text-anchor="end" font-family="Arial, sans-serif" font-size="12" fill="#526170">{value:.0f}</text>')
    for group_index, label in enumerate(labels):
        group_x = margin_left + group_index * group_width + 14
        for series_index, (_name, values, color) in enumerate(series):
            value = values[group_index]
            bar_height = plot_height * value / y_max
            x = group_x + series_index * (bar_width + bar_gap)
            y = margin_top + plot_height - bar_height
            lines.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" fill="{color}"/>')
            lines.append(f'<text x="{x + bar_width / 2:.1f}" y="{y - 6:.1f}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11" fill="#26323f">{value:.0f}</text>')
        lines.append(f'<text x="{group_x + (len(series) * (bar_width + bar_gap)) / 2 - bar_gap:.1f}" y="{height - 52}" text-anchor="middle" font-family="Arial, sans-serif" font-size="13" fill="#26323f">{label}</text>')
    legend_x = margin_left
    for index, (name, _values, color) in enumerate(series):
        x = legend_x + index * 180
        lines.append(f'<rect x="{x}" y="{height - 28}" width="14" height="14" fill="{color}"/>')
        lines.append(f'<text x="{x + 20}" y="{height - 16}" font-family="Arial, sans-serif" font-size="13" fill="#26323f">{name}</text>')
    lines.append("</svg>")
    return "\n".join(lines)


def svg_single_bars(
    *,
    title: str,
    labels: Sequence[str],
    values: Sequence[float],
    color: str,
    y_max: float,
    y_label: str,
    width: int = 980,
    height: int = 500,
) -> str:
    return svg_grouped_bars(
        title=title,
        labels=labels,
        series=((y_label, values, color),),
        y_max=y_max,
        y_label=y_label,
        width=width,
        height=height,
    )


def write_figures(rows: Sequence[Mapping[str, object]], panels: Sequence[Mapping[str, object]], figure_dir: Path) -> list[Path]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    selected = [row for row in rows if row["run_id"] in {"v4_broad", "v5_broad", "v6_broad", "v6_balanced"}]
    labels = [str(row["display_name"]) for row in selected]
    top1 = [float(row["top_1pct_count"]) for row in selected]
    top5 = [float(row["top_5pct_count"]) for row in selected]
    top10 = [float(row["top_10pct_count"]) for row in selected]
    auroc = [float(row["auroc_vs_matched_decoys"]) for row in selected]
    outputs = []
    top_recovery = figure_dir / "frozen-validation-recovery.svg"
    top_recovery.write_text(
        svg_grouped_bars(
            title="Frozen ACE holdout recovery",
            labels=labels,
            series=(
                ("Top 1%", top1, "#1f77b4"),
                ("Top 5%", top5, "#2ca02c"),
                ("Top 10%", top10, "#d62728"),
            ),
            y_max=max(top10) + 4,
            y_label="Recovered positives",
        ),
        encoding="utf-8",
    )
    outputs.append(top_recovery)
    auroc_chart = figure_dir / "frozen-validation-matched-auroc.svg"
    auroc_chart.write_text(
        svg_single_bars(
            title="Matched-decoy discrimination",
            labels=labels,
            values=auroc,
            color="#4b5563",
            y_max=1.0,
            y_label="AUROC",
        ),
        encoding="utf-8",
    )
    outputs.append(auroc_chart)
    panel_labels = [str(row["panel_id"]) for row in panels]
    short_counts = [float(row["length_le_5_count"]) for row in panels]
    long_counts = [float(row["length_gt_5_count"]) for row in panels]
    panel_chart = figure_dir / "panel-composition-length.svg"
    panel_chart.write_text(
        svg_grouped_bars(
            title="Primary-literature panel length mix",
            labels=panel_labels,
            series=(("Length <= 5", short_counts, "#007c89"), ("Length > 5", long_counts, "#f59e0b")),
            y_max=max(short_counts + long_counts) + 6,
            y_label="Peptides",
        ),
        encoding="utf-8",
    )
    outputs.append(panel_chart)
    return outputs


def main() -> int:
    args = parse_args()
    table_dir = Path(args.table_dir)
    figure_dir = Path(args.figure_dir)
    table_dir.mkdir(parents=True, exist_ok=True)
    rows = validation_rows()
    panels = panel_rows()
    write_csv_rows(table_dir / "ace-frozen-validation-summary.csv", rows, VALIDATION_FIELDS)
    write_csv_rows(table_dir / "ace-panel-composition-summary.csv", panels, PANEL_FIELDS)
    baseline = external_baseline_rows()
    if baseline:
        fields = ["panel_id"] + [field for field in baseline[0].keys() if field != "panel_id"]
        write_csv_rows(table_dir / "ace-external-baseline-summary.csv", baseline, fields)
    outputs = write_figures(rows, panels, figure_dir)
    print(table_dir / "ace-frozen-validation-summary.csv")
    print(table_dir / "ace-panel-composition-summary.csv")
    if baseline:
        print(table_dir / "ace-external-baseline-summary.csv")
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

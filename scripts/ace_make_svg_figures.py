#!/usr/bin/env python3
"""Generate lightweight SVG figures for ACE rediscovery reports."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import PAPER_DIR, read_csv_rows  # noqa: E402

DEFAULT_RUN_DIR = PAPER_DIR


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default=str(DEFAULT_RUN_DIR))
    return parser.parse_args()


def svg_bar_chart(title: str, labels: list[str], values: list[float], *, width: int = 980, height: int = 460) -> str:
    margin_left = 250
    margin_right = 40
    margin_top = 64
    margin_bottom = 54
    plot_width = width - margin_left - margin_right
    row_height = (height - margin_top - margin_bottom) / max(1, len(labels))
    max_value = max(values) if values else 1.0
    max_value = max(max_value, 1.0)
    rows: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{margin_left}" y="34" font-family="Arial, sans-serif" font-size="22" font-weight="700" fill="#1f2933">{title}</text>',
    ]
    for index, (label, value) in enumerate(zip(labels, values)):
        y = margin_top + index * row_height
        bar_width = plot_width * value / max_value
        rows.extend(
            [
                f'<text x="{margin_left - 12}" y="{y + row_height * 0.58:.1f}" text-anchor="end" font-family="Arial, sans-serif" font-size="14" fill="#26323f">{label}</text>',
                f'<rect x="{margin_left}" y="{y + row_height * 0.18:.1f}" width="{bar_width:.1f}" height="{row_height * 0.52:.1f}" fill="#2563eb"/>',
                f'<text x="{margin_left + bar_width + 8:.1f}" y="{y + row_height * 0.58:.1f}" font-family="Arial, sans-serif" font-size="14" fill="#26323f">{value:g}</text>',
            ]
        )
    axis_y = height - margin_bottom + 10
    rows.append(f'<line x1="{margin_left}" y1="{axis_y}" x2="{width - margin_right}" y2="{axis_y}" stroke="#9aa5b1" stroke-width="1"/>')
    rows.append("</svg>")
    return "\n".join(rows)


def write_stratified_figure(run_dir: Path) -> Path | None:
    path = run_dir / "ace-stratified-results.csv"
    if not path.exists():
        return None
    rows = read_csv_rows(path)
    labels = [row["subset_id"].replace("_", " ") for row in rows]
    values = [float(row["top_1pct_count"] or 0) for row in rows]
    out = run_dir / "figures" / "stratified-top1pct.svg"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg_bar_chart("Known ACE positives recovered in top 1%", labels, values), encoding="utf-8")
    return out


def write_ablation_figure(run_dir: Path) -> Path | None:
    path = run_dir / "ace-ablation-results.csv"
    if not path.exists():
        return None
    rows = read_csv_rows(path)
    total_results = len(read_csv_rows(run_dir / "ace-rediscovery-results.csv"))
    top_cutoff = max(1, round(total_results * 0.01 + 0.499999))
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        if int(row["ablation_rank"]) <= top_cutoff:
            counts[row["ablation_id"]] += 1
    order = [
        "full",
        "no_prodrug_triad",
        "no_cterm_proline_bonus",
        "no_position_priors",
        "composition_only",
        "random_hash_baseline",
    ]
    labels = [item.replace("_", " ") for item in order if item in counts]
    values = [float(counts[item]) for item in order if item in counts]
    out = run_dir / "figures" / "ablation-top1pct.svg"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg_bar_chart("Top-1% recovery by scoring ablation", labels, values), encoding="utf-8")
    return out


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir)
    outputs = [
        output
        for output in (
            write_stratified_figure(run_dir),
            write_ablation_figure(run_dir),
        )
        if output is not None
    ]
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

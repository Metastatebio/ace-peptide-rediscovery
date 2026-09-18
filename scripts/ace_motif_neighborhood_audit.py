#!/usr/bin/env python3
"""Audit ACE v6 positives for sequence-neighborhood leakage against prior panels."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from statistics import mean, median
from typing import Any, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO_ROOT / "docs" / "publication-engine" / "papers" / "01-ace-rediscovery"
EXPERIMENT_DIR = PAPER_DIR / "experiments" / "2026-08-29-v6-independent-validation-001"
DEFAULT_CANDIDATE_PANEL = PAPER_DIR / "known-positive-v6-independent-primary-panel.csv"
DEFAULT_RANKS = (
    EXPERIMENT_DIR
    / "frozen-rerank"
    / "balanced-v6-candidate"
    / "frozen-rerank-positive-ranks.csv"
)
DEFAULT_UNIVERSE_MANIFEST = EXPERIMENT_DIR / "universe-manifest.json"
DEFAULT_OUT_CSV = EXPERIMENT_DIR / "v6-motif-neighborhood-audit.csv"
DEFAULT_OUT_JSON = EXPERIMENT_DIR / "v6-motif-neighborhood-audit.json"
DEFAULT_OUT_REPORT = EXPERIMENT_DIR / "v6-motif-neighborhood-audit.md"
DEFAULT_EXCLUDE_PANELS = (
    PAPER_DIR / "known-positive-seed.csv",
    PAPER_DIR / "known-positive-expanded-panel.csv",
    PAPER_DIR / "known-positive-postfreeze-holdout.csv",
    PAPER_DIR / "known-positive-postfreeze-v3-primary-panel.csv",
    PAPER_DIR / "known-positive-v4-independent-primary-panel.csv",
    PAPER_DIR / "known-positive-v5-failure-mode-primary-panel.csv",
    PAPER_DIR / "known-positive-ahtpdb-ic50-training-panel.csv",
    PAPER_DIR / "known-positive-ahtpdb-ic50-le10-training-panel.csv",
    PAPER_DIR / "known-positive-ahtpdb-ic50-le30-training-panel.csv",
    PAPER_DIR / "known-positive-ahtpdb-ic50-le100-training-panel.csv",
    PAPER_DIR / "known-positive-ahtpdb-ic50-le100-min3-training-panel.csv",
)

ROW_FIELDS = (
    "sequence",
    "length",
    "rank",
    "top_1pct",
    "source_material",
    "reference_url",
    "nearest_prior_sequence",
    "nearest_prior_distance",
    "nearest_prior_identity",
    "nearest_prior_panels",
    "one_edit_neighbor_count",
    "two_edit_neighbor_count",
    "high_similarity_neighbor_count_identity_ge_0_8",
    "cterm2",
    "cterm2_prior_sequence_count",
    "cterm3",
    "cterm3_prior_sequence_count",
    "neighborhood_class",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-panel", default=str(DEFAULT_CANDIDATE_PANEL))
    parser.add_argument("--ranks", default=str(DEFAULT_RANKS))
    parser.add_argument("--universe-manifest", default=str(DEFAULT_UNIVERSE_MANIFEST))
    parser.add_argument("--exclude-panel", action="append", default=[])
    parser.add_argument("--out-csv", default=str(DEFAULT_OUT_CSV))
    parser.add_argument("--out-json", default=str(DEFAULT_OUT_JSON))
    parser.add_argument("--out-report", default=str(DEFAULT_OUT_REPORT))
    return parser.parse_args()


def clean_sequence(sequence: str) -> str:
    return re.sub(r"[^A-Z]", "", sequence.upper())


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


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def edit_distance(left: str, right: str) -> int:
    if left == right:
        return 0
    if not left:
        return len(right)
    if not right:
        return len(left)
    previous = list(range(len(right) + 1))
    for row_index, left_char in enumerate(left, start=1):
        current = [row_index]
        for col_index, right_char in enumerate(right, start=1):
            substitution_cost = 0 if left_char == right_char else 1
            current.append(
                min(
                    previous[col_index] + 1,
                    current[col_index - 1] + 1,
                    previous[col_index - 1] + substitution_cost,
                )
            )
        previous = current
    return previous[-1]


def normalized_identity(left: str, right: str, distance: int) -> float:
    denominator = max(len(left), len(right), 1)
    return max(0.0, 1.0 - (distance / denominator))


def load_prior_sequences(paths: Sequence[Path]) -> dict[str, set[str]]:
    panels_by_sequence: dict[str, set[str]] = defaultdict(set)
    for path in paths:
        for row in read_csv_rows(path):
            sequence = clean_sequence(row.get("sequence", ""))
            if sequence:
                panels_by_sequence[sequence].add(path.name)
    return panels_by_sequence


def load_rank_rows(path: Path) -> dict[str, dict[str, str]]:
    rows_by_sequence: dict[str, dict[str, str]] = {}
    for row in read_csv_rows(path):
        sequence = clean_sequence(row.get("sequence", ""))
        if sequence:
            rows_by_sequence[sequence] = row
    return rows_by_sequence


def cterminal(sequence: str, width: int) -> str:
    return sequence[-width:] if len(sequence) >= width else ""


def neighborhood_class(row: Mapping[str, Any]) -> str:
    if int(row["nearest_prior_distance"]) == 0:
        return "exact_overlap"
    if int(row["one_edit_neighbor_count"]) > 0:
        return "single_edit_neighbor"
    if int(row["high_similarity_neighbor_count_identity_ge_0_8"]) > 0:
        return "high_similarity_neighbor"
    if int(row["cterm3_prior_sequence_count"]) > 0:
        return "shared_cterm3_motif"
    if int(row["cterm2_prior_sequence_count"]) > 0:
        return "shared_cterm2_motif"
    return "distant_from_prior_positive_panels"


def top_threshold(universe_manifest: Path, rank_rows: Mapping[str, Mapping[str, str]]) -> int:
    if universe_manifest.exists():
        payload = json.loads(universe_manifest.read_text(encoding="utf-8"))
        total = int(payload.get("public_rows") or payload.get("key_rows") or 0)
        if total > 0:
            return max(1, math.ceil(total * 0.01))
    ranks = [int(float(row["rank"])) for row in rank_rows.values() if row.get("rank")]
    return max(1, math.ceil(max(ranks) * 0.01)) if ranks else 1


def build_rows(
    candidate_panel: Path,
    rank_path: Path,
    universe_manifest: Path,
    exclude_panels: Sequence[Path],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    prior_by_sequence = load_prior_sequences(exclude_panels)
    prior_sequences = sorted(prior_by_sequence)
    rank_rows = load_rank_rows(rank_path)
    threshold = top_threshold(universe_manifest, rank_rows)

    cterm2_map: dict[str, set[str]] = defaultdict(set)
    cterm3_map: dict[str, set[str]] = defaultdict(set)
    for sequence in prior_sequences:
        if cterminal(sequence, 2):
            cterm2_map[cterminal(sequence, 2)].add(sequence)
        if cterminal(sequence, 3):
            cterm3_map[cterminal(sequence, 3)].add(sequence)

    rows: list[dict[str, Any]] = []
    for row in read_csv_rows(candidate_panel):
        sequence = clean_sequence(row.get("sequence", ""))
        if not sequence:
            continue
        distances: list[tuple[int, float, str]] = []
        for prior_sequence in prior_sequences:
            distance = edit_distance(sequence, prior_sequence)
            distances.append((distance, normalized_identity(sequence, prior_sequence, distance), prior_sequence))
        nearest_distance, nearest_identity, nearest_sequence = min(
            distances,
            key=lambda item: (item[0], -item[1], abs(len(sequence) - len(item[2])), item[2]),
        )
        one_edit_count = sum(1 for distance, _, _ in distances if distance == 1)
        two_edit_count = sum(1 for distance, _, _ in distances if 0 < distance <= 2)
        high_similarity_count = sum(1 for distance, identity, _ in distances if distance > 0 and identity >= 0.8)
        rank_row = rank_rows.get(sequence, {})
        rank = int(float(rank_row.get("rank", "0"))) if rank_row.get("rank") else 0
        out_row = {
            "sequence": sequence,
            "length": len(sequence),
            "rank": rank or "",
            "top_1pct": "yes" if rank and rank <= threshold else "no",
            "source_material": row.get("source_material", ""),
            "reference_url": row.get("reference_url", ""),
            "nearest_prior_sequence": nearest_sequence,
            "nearest_prior_distance": nearest_distance,
            "nearest_prior_identity": round(nearest_identity, 4),
            "nearest_prior_panels": ";".join(sorted(prior_by_sequence[nearest_sequence])),
            "one_edit_neighbor_count": one_edit_count,
            "two_edit_neighbor_count": two_edit_count,
            "high_similarity_neighbor_count_identity_ge_0_8": high_similarity_count,
            "cterm2": cterminal(sequence, 2),
            "cterm2_prior_sequence_count": len(cterm2_map.get(cterminal(sequence, 2), set())),
            "cterm3": cterminal(sequence, 3),
            "cterm3_prior_sequence_count": len(cterm3_map.get(cterminal(sequence, 3), set())),
        }
        out_row["neighborhood_class"] = neighborhood_class(out_row)
        rows.append(out_row)

    top_rows = [row for row in rows if row["top_1pct"] == "yes"]
    classes = Counter(str(row["neighborhood_class"]) for row in rows)
    top_classes = Counter(str(row["neighborhood_class"]) for row in top_rows)
    identities = [float(row["nearest_prior_identity"]) for row in rows]
    summary = {
        "date": date.today().isoformat(),
        "candidate_panel": str(candidate_panel),
        "rank_file": str(rank_path),
        "universe_manifest": str(universe_manifest),
        "exclude_panels": [str(path) for path in exclude_panels],
        "positive_count": len(rows),
        "prior_sequence_count": len(prior_sequences),
        "top_1pct_rank_threshold": threshold,
        "top_1pct_count": len(top_rows),
        "exact_overlap_count": classes.get("exact_overlap", 0),
        "single_edit_neighbor_count": classes.get("single_edit_neighbor", 0),
        "high_similarity_neighbor_count_identity_ge_0_8": classes.get("high_similarity_neighbor", 0),
        "shared_cterm3_only_count": classes.get("shared_cterm3_motif", 0),
        "shared_cterm2_only_count": classes.get("shared_cterm2_motif", 0),
        "distant_count": classes.get("distant_from_prior_positive_panels", 0),
        "top_1pct_single_edit_neighbor_count": top_classes.get("single_edit_neighbor", 0),
        "top_1pct_high_similarity_neighbor_count_identity_ge_0_8": top_classes.get(
            "high_similarity_neighbor", 0
        ),
        "top_1pct_shared_cterm3_or_cterm2_count": top_classes.get("shared_cterm3_motif", 0)
        + top_classes.get("shared_cterm2_motif", 0),
        "top_1pct_distant_count": top_classes.get("distant_from_prior_positive_panels", 0),
        "mean_nearest_prior_identity": round(mean(identities), 4) if identities else None,
        "median_nearest_prior_identity": round(median(identities), 4) if identities else None,
        "neighborhood_class_counts": dict(sorted(classes.items())),
        "top_1pct_neighborhood_class_counts": dict(sorted(top_classes.items())),
    }
    return rows, summary


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def report_text(summary: Mapping[str, Any]) -> str:
    class_rows = [
        [key, value, summary.get("top_1pct_neighborhood_class_counts", {}).get(key, 0)]
        for key, value in summary.get("neighborhood_class_counts", {}).items()
    ]
    return f"""# ACE v6 Motif-Neighborhood Audit

Date: {summary['date']}

## Purpose

Check whether the funding-grade ACE v6 balanced rediscovery result is just exact reuse or near-neighbor reuse from prior positive panels.

## Inputs

- Candidate panel: `{summary['candidate_panel']}`
- Frozen positive ranks: `{summary['rank_file']}`
- Prior positive panels audited: `{len(summary['exclude_panels'])}`
- Prior unique positive sequences: `{summary['prior_sequence_count']}`

## Summary

- v6 positives audited: `{summary['positive_count']}`
- Exact overlaps with prior panels: `{summary['exact_overlap_count']}`
- Top-1% positives: `{summary['top_1pct_count']}`
- Top-1% positives with a one-edit prior neighbor: `{summary['top_1pct_single_edit_neighbor_count']}`
- Top-1% positives with a high-similarity prior neighbor, identity >= 0.8: `{summary['top_1pct_high_similarity_neighbor_count_identity_ge_0_8']}`
- Top-1% positives classified only by shared C-terminal motif: `{summary['top_1pct_shared_cterm3_or_cterm2_count']}`
- Top-1% positives distant from prior positive panels: `{summary['top_1pct_distant_count']}`
- Median nearest-prior identity: `{summary['median_nearest_prior_identity']}`

## Neighborhood Classes

{markdown_table(["Class", "All positives", "Top 1% positives"], class_rows)}

## Interpretation

The exact-sequence leakage result remains clean. The right claim boundary is narrower: this is a frozen ACE rediscovery and prioritization result under exact-sequence exclusion, not a proof of de novo molecular understanding. Sequence-neighborhood structure is now explicit and must travel with the manuscript.

For a funding-grade upgrade, the next holdout should be source-family and motif-neighborhood pre-registered before scoring, with success reported separately for close-neighbor and distant positives.
"""


def main() -> int:
    args = parse_args()
    exclude_panels = [Path(path) for path in args.exclude_panel] or list(DEFAULT_EXCLUDE_PANELS)
    rows, summary = build_rows(
        candidate_panel=Path(args.candidate_panel),
        rank_path=Path(args.ranks),
        universe_manifest=Path(args.universe_manifest),
        exclude_panels=exclude_panels,
    )
    write_csv(Path(args.out_csv), rows, ROW_FIELDS)
    write_json(Path(args.out_json), summary)
    Path(args.out_report).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_report).write_text(report_text(summary), encoding="utf-8")
    print(json.dumps({"rows": len(rows), "report": args.out_report}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

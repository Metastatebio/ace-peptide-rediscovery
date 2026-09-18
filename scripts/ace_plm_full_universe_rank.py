#!/usr/bin/env python3
"""Apply a frozen length-calibrated ACE/ProteinLM fusion to a full scored universe."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import date
from pathlib import Path
from statistics import mean, median
from typing import Any, Mapping, Sequence

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO_ROOT / "docs" / "publication-engine" / "papers" / "01-ace-rediscovery"
V6_EXPERIMENT_DIR = PAPER_DIR / "experiments" / "2026-08-29-v6-independent-validation-001"
DEFAULT_SCORED_PANEL = (
    PAPER_DIR
    / "experiments"
    / "2026-08-30-ace-protein-lm-fusion-004-esm8m-full-universe"
    / "ace_plm_scored_panel.csv"
)
DEFAULT_V6_SUMMARY = (
    V6_EXPERIMENT_DIR
    / "frozen-rerank"
    / "balanced-v6-candidate"
    / "frozen-rerank-summary.csv"
)
DEFAULT_OUT_DIR = PAPER_DIR / "experiments" / "2026-08-30-ace-plm-full-universe-rank-001"
POSITIVE_LABELS = {"known_positive", "digestome_known_positive"}
MATCHED_DECOY_LABEL = "matched_decoy"

RANKED_FIELDS = (
    "rank",
    "candidate_id",
    "sequence",
    "length",
    "label",
    "origin",
    "provenance",
    "ace_rank",
    "ace_score",
    "ace_norm",
    "plm_mean_log_likelihood",
    "plm_length_z",
    "plm_length_z_norm",
    "fused_score",
    "mol_wt",
    "mol_logp",
    "tpsa",
    "rotatable_bonds",
)
SUMMARY_FIELDS = (
    "fusion_id",
    "ace_weight",
    "plm_length_z_weight",
    "total_candidates",
    "positive_count",
    "top_0_5pct_count",
    "top_1pct_count",
    "top_1pct_fraction",
    "top_1pct_enrichment_over_random",
    "top_1pct_hypergeom_p",
    "top_2pct_count",
    "top_5pct_count",
    "top_10pct_count",
    "top_100_count",
    "best_positive_rank",
    "median_positive_rank",
    "mean_positive_rank",
    "auroc_vs_all_negatives",
    "auroc_vs_matched_decoys",
    "average_precision_vs_all_negatives",
    "average_precision_vs_matched_decoys",
)
POSITIVE_FIELDS = (
    "rank",
    "sequence",
    "label",
    "origin",
    "ace_rank",
    "ace_score",
    "plm_mean_log_likelihood",
    "plm_length_z",
    "fused_score",
    "provenance",
)
SHORTLIST_FIELDS = (
    "rank",
    "candidate_id",
    "sequence",
    "length",
    "label",
    "origin",
    "ace_rank",
    "ace_score",
    "plm_mean_log_likelihood",
    "plm_length_z",
    "fused_score",
    "mol_wt",
    "mol_logp",
    "tpsa",
    "rotatable_bonds",
    "provenance",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scored-panel", default=str(DEFAULT_SCORED_PANEL))
    parser.add_argument("--v6-summary", default=str(DEFAULT_V6_SUMMARY))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--ace-weight", type=float, default=0.65)
    parser.add_argument("--shortlist-n", type=int, default=200)
    return parser.parse_args()


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


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(payload), indent=2, sort_keys=True), encoding="utf-8")


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return default if math.isnan(number) else number


def minmax(values: Sequence[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    minimum = float(np.min(arr))
    maximum = float(np.max(arr))
    if math.isclose(minimum, maximum):
        return np.full(len(arr), 0.5, dtype=float)
    return (arr - minimum) / (maximum - minimum)


def length_z(rows: Sequence[Mapping[str, str]], field: str) -> np.ndarray:
    values = np.asarray([safe_float(row[field]) for row in rows], dtype=float)
    output = np.zeros(len(rows), dtype=float)
    global_sd = float(np.std(values)) or 1.0
    for length in sorted({int(row["length"]) for row in rows}):
        indices = [index for index, row in enumerate(rows) if int(row["length"]) == length]
        subset = values[indices]
        center = float(np.mean(subset))
        scale = float(np.std(subset)) or global_sd
        output[indices] = (subset - center) / scale
    return output


def hypergeom_tail(
    *,
    population_size: int,
    success_states: int,
    draws: int,
    observed_successes: int,
) -> float:
    if observed_successes <= 0:
        return 1.0
    denominator = (
        math.lgamma(population_size + 1)
        - math.lgamma(draws + 1)
        - math.lgamma(population_size - draws + 1)
    )
    terms = []
    for successes in range(observed_successes, min(success_states, draws) + 1):
        failures = draws - successes
        if failures > population_size - success_states:
            continue
        log_prob = (
            math.lgamma(success_states + 1)
            - math.lgamma(successes + 1)
            - math.lgamma(success_states - successes + 1)
            + math.lgamma(population_size - success_states + 1)
            - math.lgamma(failures + 1)
            - math.lgamma(population_size - success_states - failures + 1)
            - denominator
        )
        terms.append(log_prob)
    if not terms:
        return 0.0
    maximum = max(terms)
    return float(math.exp(maximum) * sum(math.exp(term - maximum) for term in terms))


def build_ranked_rows(rows: Sequence[Mapping[str, str]], ace_weight: float) -> list[dict[str, Any]]:
    ace_norm = minmax([safe_float(row["ace_score"]) for row in rows])
    plm_z = length_z(rows, "plm_mean_log_likelihood")
    plm_z_norm = minmax(plm_z)
    fused = (ace_weight * ace_norm) + ((1.0 - ace_weight) * plm_z_norm)
    order = sorted(range(len(rows)), key=lambda index: (-float(fused[index]), rows[index]["sequence"]))
    ranked_rows = []
    for rank, index in enumerate(order, start=1):
        row = rows[index]
        ranked_rows.append(
            {
                "rank": rank,
                "candidate_id": row.get("candidate_id", ""),
                "sequence": row.get("sequence", ""),
                "length": row.get("length", ""),
                "label": row.get("label", ""),
                "origin": row.get("origin", ""),
                "provenance": row.get("provenance", ""),
                "ace_rank": row.get("ace_rank", ""),
                "ace_score": row.get("ace_score", ""),
                "ace_norm": round(float(ace_norm[index]), 10),
                "plm_mean_log_likelihood": row.get("plm_mean_log_likelihood", ""),
                "plm_length_z": round(float(plm_z[index]), 10),
                "plm_length_z_norm": round(float(plm_z_norm[index]), 10),
                "fused_score": round(float(fused[index]), 10),
                "mol_wt": row.get("mol_wt", ""),
                "mol_logp": row.get("mol_logp", ""),
                "tpsa": row.get("tpsa", ""),
                "rotatable_bonds": row.get("rotatable_bonds", ""),
            }
        )
    return ranked_rows


def labels_for(rows: Sequence[Mapping[str, Any]]) -> np.ndarray:
    return np.asarray([1 if row.get("label") in POSITIVE_LABELS else 0 for row in rows], dtype=int)


def scores_for(rows: Sequence[Mapping[str, Any]]) -> np.ndarray:
    return np.asarray([safe_float(row["fused_score"]) for row in rows], dtype=float)


def summary_row(ranked_rows: Sequence[Mapping[str, Any]], ace_weight: float) -> dict[str, Any]:
    labels = labels_for(ranked_rows)
    scores = scores_for(ranked_rows)
    positive_count = int(np.sum(labels))
    positive_ranks = [int(row["rank"]) for row in ranked_rows if row.get("label") in POSITIVE_LABELS]
    top_1_cutoff = max(1, math.ceil(len(ranked_rows) * 0.01))
    matched_rows = [
        row
        for row in ranked_rows
        if row.get("label") in POSITIVE_LABELS or row.get("label") == MATCHED_DECOY_LABEL
    ]
    matched_labels = labels_for(matched_rows)
    matched_scores = scores_for(matched_rows)
    top_1_count = sum(1 for rank in positive_ranks if rank <= top_1_cutoff)
    random_fraction = top_1_cutoff / float(len(ranked_rows))
    return {
        "fusion_id": f"ace_{ace_weight:.2f}_plm_length_z_{1.0 - ace_weight:.2f}",
        "ace_weight": ace_weight,
        "plm_length_z_weight": round(1.0 - ace_weight, 4),
        "total_candidates": len(ranked_rows),
        "positive_count": positive_count,
        "top_0_5pct_count": sum(1 for rank in positive_ranks if rank <= max(1, math.ceil(len(ranked_rows) * 0.005))),
        "top_1pct_count": top_1_count,
        "top_1pct_fraction": top_1_count / float(max(positive_count, 1)),
        "top_1pct_enrichment_over_random": (
            (top_1_count / float(max(positive_count, 1))) / random_fraction if random_fraction else 0.0
        ),
        "top_1pct_hypergeom_p": hypergeom_tail(
            population_size=len(ranked_rows),
            success_states=positive_count,
            draws=top_1_cutoff,
            observed_successes=top_1_count,
        ),
        "top_2pct_count": sum(1 for rank in positive_ranks if rank <= max(1, math.ceil(len(ranked_rows) * 0.02))),
        "top_5pct_count": sum(1 for rank in positive_ranks if rank <= max(1, math.ceil(len(ranked_rows) * 0.05))),
        "top_10pct_count": sum(1 for rank in positive_ranks if rank <= max(1, math.ceil(len(ranked_rows) * 0.10))),
        "top_100_count": sum(1 for rank in positive_ranks if rank <= 100),
        "best_positive_rank": min(positive_ranks),
        "median_positive_rank": median(positive_ranks),
        "mean_positive_rank": mean(positive_ranks),
        "auroc_vs_all_negatives": float(roc_auc_score(labels, scores)),
        "auroc_vs_matched_decoys": float(roc_auc_score(matched_labels, matched_scores)),
        "average_precision_vs_all_negatives": float(average_precision_score(labels, scores)),
        "average_precision_vs_matched_decoys": float(average_precision_score(matched_labels, matched_scores)),
    }


def candidate_shortlist(
    ranked_rows: Sequence[Mapping[str, Any]],
    shortlist_n: int,
    *,
    digestome_only: bool = False,
) -> list[dict[str, Any]]:
    rows = [
        row
        for row in ranked_rows
        if row.get("label") not in POSITIVE_LABELS and row.get("label") != MATCHED_DECOY_LABEL
    ]
    if digestome_only:
        rows = [row for row in rows if row.get("label") == "digestome_background"]
    return [dict(row) for row in rows[:shortlist_n]]


def fmt(value: Any, digits: int = 4) -> str:
    number = safe_float(value, math.nan)
    if math.isnan(number):
        return ""
    return f"{number:.{digits}f}"


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def report_markdown(
    *,
    ranked_rows: Sequence[Mapping[str, Any]],
    summary: Mapping[str, Any],
    v6_summary: Mapping[str, str],
    shortlist: Sequence[Mapping[str, Any]],
    digestome_shortlist: Sequence[Mapping[str, Any]],
    out_dir: Path,
) -> str:
    positives = [row for row in ranked_rows if row.get("label") in POSITIVE_LABELS]
    top_positive_rows = [
        [
            row["rank"],
            row["sequence"],
            row["label"],
            row["ace_rank"],
            fmt(row["ace_score"], 4),
            fmt(row["plm_length_z"], 3),
            fmt(row["fused_score"], 4),
        ]
        for row in positives[:12]
    ]
    shortlist_rows = [
        [
            row["rank"],
            row["sequence"],
            row["label"],
            row["origin"],
            row["ace_rank"],
            fmt(row["plm_length_z"], 3),
            fmt(row["fused_score"], 4),
        ]
        for row in shortlist[:15]
    ]
    digestome_shortlist_rows = [
        [
            row["rank"],
            row["sequence"],
            row["origin"],
            row["ace_rank"],
            fmt(row["plm_length_z"], 3),
            fmt(row["fused_score"], 4),
        ]
        for row in digestome_shortlist[:15]
    ]
    return f"""# ACE PLM Full-Universe Rank

Date: {date.today().isoformat()}

## Frozen Fusion

- Fusion: `{summary["fusion_id"]}`
- Total candidates: `{summary["total_candidates"]}`
- Positive labels: `{summary["positive_count"]}`
- Input v6 top 1% positives: `{v6_summary.get("top_1pct_count", "")}`
- Input v6 top 100 positives: `{v6_summary.get("top_100_count", "")}`
- Input v6 best positive rank: `{v6_summary.get("best_positive_rank", "")}`
- Input v6 matched-decoy AUROC: `{fmt(v6_summary.get("auroc_vs_matched_decoys"), 4)}`
- Input v6 all-negative AUROC: `{fmt(v6_summary.get("auroc_vs_all_negatives"), 4)}`

## Full-Universe Result

| Metric | Value |
| --- | --- |
| Top 1% positives | {summary["top_1pct_count"]} |
| Top 1% enrichment | {fmt(summary["top_1pct_enrichment_over_random"], 2)}x |
| Top 1% hypergeometric p | {fmt(summary["top_1pct_hypergeom_p"], 6)} |
| Top 100 positives | {summary["top_100_count"]} |
| Best positive rank | {summary["best_positive_rank"]} |
| Matched-decoy AUROC | {fmt(summary["auroc_vs_matched_decoys"], 4)} |
| All-negative AUROC | {fmt(summary["auroc_vs_all_negatives"], 4)} |
| Median positive rank | {fmt(summary["median_positive_rank"], 1)} |

## Top Positive Recoveries

{markdown_table(["Rank", "Sequence", "Label", "v6 rank", "v6 score", "PLM z", "Fused"], top_positive_rows)}

## Computational Digestome Discovery Shortlist

This shortlist excludes known positives, matched decoys, and synthetic random peptides. It is the most relevant computational prioritization list for a food-peptide validation paper.

{markdown_table(["Rank", "Sequence", "Origin", "v6 rank", "PLM z", "Fused"], digestome_shortlist_rows)}

## Computational Discovery Shortlist

This shortlist excludes known positives and matched decoys. It is a computational prioritization list, not a binding or potency claim.

{markdown_table(["Rank", "Sequence", "Label", "Origin", "v6 rank", "PLM z", "Fused"], shortlist_rows)}

## Files

- Ranked full universe: `{out_dir / "ace_plm_lengthz_full_universe_ranked_candidates.csv"}`
- Positive ranks: `{out_dir / "ace_plm_lengthz_positive_ranks.csv"}`
- Summary: `{out_dir / "ace_plm_lengthz_full_universe_summary.csv"}`
- Digestome discovery shortlist: `{out_dir / "ace_plm_lengthz_digestome_discovery_shortlist.csv"}`
- Discovery shortlist: `{out_dir / "ace_plm_lengthz_discovery_shortlist.csv"}`
- Manifest: `{out_dir / "ace_plm_full_universe_rank_manifest.json"}`
"""


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = read_csv_rows(Path(args.scored_panel))
    ranked_rows = build_ranked_rows(rows, args.ace_weight)
    summary = summary_row(ranked_rows, args.ace_weight)
    v6_summary_rows = read_csv_rows(Path(args.v6_summary))
    v6_summary = v6_summary_rows[0] if v6_summary_rows else {}
    positive_rows = [row for row in ranked_rows if row.get("label") in POSITIVE_LABELS]
    shortlist = candidate_shortlist(ranked_rows, args.shortlist_n)
    digestome_shortlist = candidate_shortlist(ranked_rows, args.shortlist_n, digestome_only=True)

    write_csv(
        out_dir / "ace_plm_lengthz_full_universe_ranked_candidates.csv",
        ranked_rows,
        RANKED_FIELDS,
    )
    write_csv(out_dir / "ace_plm_lengthz_positive_ranks.csv", positive_rows, POSITIVE_FIELDS)
    write_csv(out_dir / "ace_plm_lengthz_full_universe_summary.csv", [summary], SUMMARY_FIELDS)
    write_csv(out_dir / "ace_plm_lengthz_discovery_shortlist.csv", shortlist, SHORTLIST_FIELDS)
    write_csv(
        out_dir / "ace_plm_lengthz_digestome_discovery_shortlist.csv",
        digestome_shortlist,
        SHORTLIST_FIELDS,
    )
    write_json(
        out_dir / "ace_plm_full_universe_rank_manifest.json",
        {
            "date": date.today().isoformat(),
            "scored_panel": str(Path(args.scored_panel)),
            "v6_summary": str(Path(args.v6_summary)),
            "ace_weight": args.ace_weight,
            "plm_length_z_weight": round(1.0 - args.ace_weight, 4),
            "label_counts": dict(Counter(row.get("label", "") for row in rows)),
            "summary": summary,
            "shortlist_n": args.shortlist_n,
            "digestome_shortlist_rows": len(digestome_shortlist),
            "claim_boundary": (
                "Uses the 0.65 ACE / 0.35 length-z PLM weight selected by parent-grouped CV. "
                "The full-universe rank is suitable for candidate prioritization and must be "
                "validated on a new frozen panel or wet lab before potency claims."
            ),
        },
    )
    (out_dir / "ace_plm_full_universe_rank_report.md").write_text(
        report_markdown(
            ranked_rows=ranked_rows,
            summary=summary,
            v6_summary=v6_summary,
            shortlist=shortlist,
            digestome_shortlist=digestome_shortlist,
            out_dir=out_dir,
        ),
        encoding="utf-8",
    )
    print(out_dir / "ace_plm_full_universe_rank_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

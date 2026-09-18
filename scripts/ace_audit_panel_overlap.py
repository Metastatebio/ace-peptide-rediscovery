#!/usr/bin/env python3
"""Audit exact sequence overlap between an ACE candidate panel and exclusion panels."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import read_csv_rows, save_json, write_csv_rows  # noqa: E402

DEFAULT_PAPER_DIR = REPO_ROOT / "docs" / "publication-engine" / "papers" / "01-ace-rediscovery"
DEFAULT_OUT_JSON = (
    DEFAULT_PAPER_DIR / "experiments" / "2026-08-29-v4-independent-validation-001" / "v4-exclusion-audit.json"
)
DEFAULT_OUT_CSV = (
    DEFAULT_PAPER_DIR / "experiments" / "2026-08-29-v4-independent-validation-001" / "v4-exclusion-overlaps.csv"
)
OVERLAP_FIELDS = ("sequence", "candidate_panel", "exclude_panel")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-panel", required=True)
    parser.add_argument("--exclude-panel", action="append", required=True)
    parser.add_argument("--out-json", default=str(DEFAULT_OUT_JSON))
    parser.add_argument("--out-csv", default=str(DEFAULT_OUT_CSV))
    return parser.parse_args()


def load_sequences(path: Path) -> set[str]:
    rows = read_csv_rows(path)
    if not rows or "sequence" not in rows[0]:
        raise ValueError(f"{path} does not contain a sequence column")
    return {row["sequence"].strip().upper() for row in rows if row.get("sequence", "").strip()}


def main() -> int:
    args = parse_args()
    candidate_panel = Path(args.candidate_panel)
    candidate_sequences = load_sequences(candidate_panel)
    sources_by_sequence: dict[str, list[str]] = defaultdict(list)
    for exclude_panel_text in args.exclude_panel:
        exclude_panel = Path(exclude_panel_text)
        for sequence in load_sequences(exclude_panel):
            if sequence in candidate_sequences:
                sources_by_sequence[sequence].append(str(exclude_panel))

    overlap_rows = [
        {
            "sequence": sequence,
            "candidate_panel": str(candidate_panel),
            "exclude_panel": exclude_panel,
        }
        for sequence, exclude_panels in sorted(sources_by_sequence.items())
        for exclude_panel in sorted(exclude_panels)
    ]
    write_csv_rows(Path(args.out_csv), overlap_rows, OVERLAP_FIELDS)
    save_json(
        Path(args.out_json),
        {
            "date": date.today().isoformat(),
            "candidate_panel": str(candidate_panel),
            "candidate_sequence_count": len(candidate_sequences),
            "exclude_panels": args.exclude_panel,
            "overlap_sequence_count": len(sources_by_sequence),
            "overlap_row_count": len(overlap_rows),
            "status": "pass" if not overlap_rows else "fail",
            "overlaps": {sequence: sorted(paths) for sequence, paths in sorted(sources_by_sequence.items())},
        },
    )
    print(f"{'PASS' if not overlap_rows else 'FAIL'}: {len(sources_by_sequence)} overlapping sequences")
    return 0 if not overlap_rows else 1


if __name__ == "__main__":
    raise SystemExit(main())

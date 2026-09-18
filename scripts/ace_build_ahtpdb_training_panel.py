#!/usr/bin/env python3
"""Build an ACE training-positive panel from the AHTPDB IC50 export."""

from __future__ import annotations

import argparse
import csv
import io
import re
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Mapping

import requests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import AA, PAPER_DIR, read_csv_rows, save_json, sha256_text, write_csv_rows  # noqa: E402

AHTPDB_IC50_URL = "http://crdd.osdd.net/raghava/ahtpdb/downloads/pepic50.txt"
DEFAULT_OUT = PAPER_DIR / "known-positive-ahtpdb-ic50-training-panel.csv"
DEFAULT_MANIFEST = PAPER_DIR / "known-positive-ahtpdb-ic50-training-panel.manifest.json"
DEFAULT_EXCLUDE = PAPER_DIR / "known-positive-postfreeze-v3-primary-panel.csv"

OUT_FIELDS = (
    "sequence",
    "target",
    "activity_class",
    "source_material",
    "generation_context",
    "evidence_type",
    "reference_url",
    "notes",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--manifest-out", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--exclude-panel", action="append", default=[str(DEFAULT_EXCLUDE)])
    parser.add_argument("--url", default=AHTPDB_IC50_URL)
    parser.add_argument("--min-length", type=int, default=2)
    parser.add_argument("--max-length", type=int, default=16)
    parser.add_argument(
        "--max-ic50-um",
        type=float,
        default=None,
        help="Keep only rows whose parsed IC50 is at or below this uM threshold.",
    )
    return parser.parse_args()


def normalized_row(row: Mapping[str, str]) -> dict[str, str]:
    return {str(key).strip(): str(value).strip() for key, value in row.items()}


def parse_ic50_to_um(ic50: str, molwt: str) -> float | None:
    value_match = re.search(r"(?<![A-Za-z0-9.])(-?[0-9]+(?:\.[0-9]+)?)", ic50)
    if not value_match:
        return None
    value = float(value_match.group(1))
    # A zero or negative reported IC50 is not a usable quantitative potency
    # measurement.  It must not enter thresholding or potency weighting as an
    # artificially maximal activity value.
    if value <= 0:
        return None
    unit = ic50.lower()
    unit = unit.replace("î¼", "µ").replace("âµ", "µ").replace("âî¼", "µ").replace("âµ", "µ")
    unit = unit.replace("μ", "µ").replace("Âµ".lower(), "µ")
    if "%" in unit:
        return None
    if "mm" in unit:
        return value * 1000.0
    if "µm" in unit or "um" in unit:
        return value
    mw_match = re.search(r"([0-9]+(?:\.[0-9]+)?)", molwt)
    if not mw_match:
        return None
    mw = float(mw_match.group(1))
    if mw <= 0:
        return None
    if "mg" in unit and "ml" in unit:
        return value / mw * 1_000_000.0
    if "µg" in unit or "ug" in unit:
        return (value * 0.001) / mw * 1_000_000.0
    return None


def load_exclusions(paths: list[str]) -> set[str]:
    excluded: set[str] = set()
    for path_text in paths:
        path = Path(path_text)
        if not path.exists():
            continue
        for row in read_csv_rows(path):
            excluded.add(row["sequence"].strip().upper())
    return excluded


def main() -> int:
    args = parse_args()
    response = requests.get(args.url, timeout=40, headers={"User-Agent": "PepLab AHTPDB training panel builder"})
    response.raise_for_status()
    raw_rows = [normalized_row(row) for row in csv.DictReader(io.StringIO(response.text), delimiter="\t")]
    excluded = load_exclusions(args.exclude_panel)
    by_sequence: dict[str, dict[str, object]] = {}
    skipped = defaultdict(int)

    for row in raw_rows:
        sequence = row.get("seq", "").upper()
        if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence):
            skipped["invalid_sequence"] += 1
            continue
        if set(sequence) - set(AA):
            skipped["noncanonical_sequence"] += 1
            continue
        if len(sequence) < args.min_length or len(sequence) > args.max_length:
            skipped["length_outside_range"] += 1
            continue
        if sequence in excluded:
            skipped["excluded_eval_sequence"] += 1
            continue
        ic50_um = parse_ic50_to_um(row.get("ic50", ""), row.get("molwt", ""))
        if ic50_um is None:
            skipped["unparsed_ic50"] += 1
            continue
        if args.max_ic50_um is not None and ic50_um > args.max_ic50_um:
            skipped["above_max_ic50_um"] += 1
            continue
        current = by_sequence.get(sequence)
        payload = {
            "sequence": sequence,
            "ic50_um": ic50_um,
            "source": row.get("source", "ND"),
            "method": row.get("method", "ND"),
            "assay": row.get("assay", "ND"),
            "mice": row.get("mice", "ND"),
            "ahtpdb_id": row.get("id", ""),
            "raw_ic50": row.get("ic50", ""),
            "molwt": row.get("molwt", ""),
        }
        if current is None or float(payload["ic50_um"]) < float(current["ic50_um"]):
            by_sequence[sequence] = payload

    output_rows: list[dict[str, object]] = []
    for sequence, row in sorted(by_sequence.items()):
        output_rows.append(
            {
                "sequence": sequence,
                "target": "ACE",
                "activity_class": "ACE inhibitory or antihypertensive peptide with IC50 in AHTPDB",
                "source_material": row["source"],
                "generation_context": row["method"],
                "evidence_type": "ahtpdb_ic50_training_positive",
                "reference_url": args.url,
                "notes": (
                    f"AHTPDB id {row['ahtpdb_id']}; retained best parsed IC50 {row['ic50_um']:.4g} uM "
                    f"from raw value {row['raw_ic50']}; assay {row['assay']}; mice {row['mice']}; molwt {row['molwt']}."
                ),
            }
        )

    write_csv_rows(Path(args.out), output_rows, OUT_FIELDS)
    save_json(
        Path(args.manifest_out),
        {
            "date": date.today().isoformat(),
            "source_url": args.url,
            "raw_rows": len(raw_rows),
            "output_rows": len(output_rows),
            "excluded_panels": args.exclude_panel,
            "excluded_sequence_count": len(excluded),
            "min_length": args.min_length,
            "max_length": args.max_length,
            "max_ic50_um": args.max_ic50_um,
            "skipped": dict(sorted(skipped.items())),
            "sequence_payload_sha256": sha256_text("\n".join(row["sequence"] for row in output_rows)),
            "claim_boundary": "AHTPDB rows are used for model development/training only, not as an independent validation panel.",
        },
    )
    print(f"wrote {len(output_rows)} AHTPDB training positives; skipped {dict(skipped)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

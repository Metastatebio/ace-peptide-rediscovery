#!/usr/bin/env python3
"""Build structured assay metadata tables from ACE positive-panel CSVs."""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import PAPER_DIR, read_csv_rows, save_json, write_csv_rows  # noqa: E402

DEFAULT_PANEL_PATHS = (
    PAPER_DIR / "known-positive-postfreeze-v3-primary-panel.csv",
    PAPER_DIR / "known-positive-v4-independent-primary-panel.csv",
    PAPER_DIR / "known-positive-v5-failure-mode-primary-panel.csv",
    PAPER_DIR / "known-positive-v6-independent-primary-panel.csv",
)
DEFAULT_OUT = PAPER_DIR / "tables" / "ace-positive-panel-assay-metadata.csv"
DEFAULT_SUMMARY_OUT = PAPER_DIR / "tables" / "ace-positive-panel-assay-metadata-summary.json"

OUTPUT_FIELDS = (
    "panel_id",
    "sequence",
    "length",
    "target",
    "activity_class",
    "source_material",
    "generation_context",
    "evidence_type",
    "reference_url",
    "assay_readout",
    "ic50_value",
    "ic50_unit",
    "ic50_um",
    "fixed_concentration_inhibition_percent",
    "digestive_context",
    "evidence_strength",
    "notes",
)

IC50_PATTERN = re.compile(
    r"IC50(?: values? (?:from|between|range from))?[^0-9]{0,40}"
    r"(?P<value>[0-9]+(?:\.[0-9]+)?)"
    r"(?:\s*(?:-|to)\s*(?P<upper>[0-9]+(?:\.[0-9]+)?))?"
    r"\s*(?P<unit>uM|µM|mM|ug/mL|ug per mL|ug/mL|mg/mL|mg per mL)",
    re.IGNORECASE,
)
INHIBITION_PATTERN = re.compile(
    r"(?P<value>[0-9]+(?:\.[0-9]+)?)\s*percent\s+ACE inhibition",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panel",
        action="append",
        default=None,
        help="Known-positive panel CSV. Defaults to v3/v4/v5/v6 claim-bearing panels.",
    )
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--summary-out", default=str(DEFAULT_SUMMARY_OUT))
    return parser.parse_args()


def panel_id_from_path(path: Path) -> str:
    stem = path.stem
    if "postfreeze-v3" in stem:
        return "v3_primary_development"
    if "v4" in stem:
        return "v4_frozen_validation"
    if "v5" in stem:
        return "v5_failure_mode"
    if "v6" in stem:
        return "v6_frozen_validation"
    return stem.replace("known-positive-", "").replace("-", "_")


def normalize_unit(unit: str) -> str:
    unit = unit.strip().replace("µ", "u")
    unit = unit.replace(" per ", "/")
    return unit


def ic50_to_um(value: float, unit: str) -> float | None:
    normalized = normalize_unit(unit).lower()
    if normalized == "um":
        return value
    if normalized == "mm":
        return value * 1000.0
    return None


def parse_ic50(notes: str) -> tuple[str, str, str]:
    match = IC50_PATTERN.search(notes)
    if not match:
        return "", "", ""
    value = float(match.group("value"))
    upper = match.group("upper")
    unit = normalize_unit(match.group("unit"))
    if upper:
        upper_value = float(upper)
        value_text = f"{value:g}-{upper_value:g}"
        um_lower = ic50_to_um(value, unit)
        um_upper = ic50_to_um(upper_value, unit)
        um_text = f"{um_lower:g}-{um_upper:g}" if um_lower is not None and um_upper is not None else ""
        return value_text, unit, um_text
    um = ic50_to_um(value, unit)
    return f"{value:g}", unit, f"{um:g}" if um is not None else ""


def parse_inhibition_percent(notes: str) -> str:
    match = INHIBITION_PATTERN.search(notes)
    return f"{float(match.group('value')):g}" if match else ""


def assay_readout(notes: str) -> str:
    if IC50_PATTERN.search(notes):
        return "IC50"
    if INHIBITION_PATTERN.search(notes):
        return "fixed_concentration_inhibition"
    return "activity_reported"


def digestive_context(generation_context: str, notes: str) -> str:
    text = f"{generation_context} {notes}".lower()
    if "before digestion" in text:
        return "pre_digestive_assay"
    if "after digestion" in text:
        return "post_digestive_assay"
    if "digestion-released" in text or "gi digestion release" in text:
        return "digestion_released_peptide"
    if "digestion" in text or "hydrolysate" in text:
        return "hydrolysate_or_digest_context"
    return "not_specified"


def evidence_strength(row: Mapping[str, str], readout: str) -> str:
    context = f"{row.get('generation_context', '')} {row.get('notes', '')}".lower()
    if readout == "IC50" and "synthesis" in context:
        return "synthesized_peptide_ic50"
    if readout == "IC50":
        return "reported_ic50"
    if readout == "fixed_concentration_inhibition":
        return "fixed_concentration_activity"
    return "reported_activity"


def structured_rows(panel_paths: Sequence[Path]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for path in panel_paths:
        panel_id = panel_id_from_path(path)
        for row in read_csv_rows(path):
            sequence = row["sequence"].strip().upper()
            notes = row.get("notes", "")
            readout = assay_readout(notes)
            ic50_value, ic50_unit, ic50_um = parse_ic50(notes)
            output.append(
                {
                    "panel_id": panel_id,
                    "sequence": sequence,
                    "length": len(sequence),
                    "target": row.get("target", ""),
                    "activity_class": row.get("activity_class", ""),
                    "source_material": row.get("source_material", ""),
                    "generation_context": row.get("generation_context", ""),
                    "evidence_type": row.get("evidence_type", ""),
                    "reference_url": row.get("reference_url", ""),
                    "assay_readout": readout,
                    "ic50_value": ic50_value,
                    "ic50_unit": ic50_unit,
                    "ic50_um": ic50_um,
                    "fixed_concentration_inhibition_percent": parse_inhibition_percent(notes),
                    "digestive_context": digestive_context(row.get("generation_context", ""), notes),
                    "evidence_strength": evidence_strength(row, readout),
                    "notes": notes,
                }
            )
    return output


def build_summary(rows: Sequence[Mapping[str, object]], panel_paths: Sequence[Path]) -> dict[str, object]:
    panel_counts = Counter(str(row["panel_id"]) for row in rows)
    readout_counts = Counter(str(row["assay_readout"]) for row in rows)
    strength_counts = Counter(str(row["evidence_strength"]) for row in rows)
    panel_source_counts: dict[str, int] = {}
    for panel_id in panel_counts:
        panel_source_counts[panel_id] = len(
            {str(row["source_material"]) for row in rows if row["panel_id"] == panel_id}
        )
    return {
        "date": date.today().isoformat(),
        "panel_paths": [str(path) for path in panel_paths],
        "row_count": len(rows),
        "panel_counts": dict(sorted(panel_counts.items())),
        "panel_source_material_counts": dict(sorted(panel_source_counts.items())),
        "assay_readout_counts": dict(sorted(readout_counts.items())),
        "evidence_strength_counts": dict(sorted(strength_counts.items())),
        "rows_with_numeric_ic50_um": sum(1 for row in rows if row.get("ic50_um")),
        "rows_with_mass_concentration_ic50": sum(
            1 for row in rows if str(row.get("ic50_unit", "")).lower() in {"ug/ml", "mg/ml"}
        ),
        "claim_boundary": "Derived metadata table. Parsed fields should be manually checked before final journal submission.",
    }


def main() -> int:
    args = parse_args()
    panel_paths = [Path(path) for path in args.panel] if args.panel else list(DEFAULT_PANEL_PATHS)
    rows = structured_rows(panel_paths)
    write_csv_rows(Path(args.out), rows, OUTPUT_FIELDS)
    save_json(Path(args.summary_out), build_summary(rows, panel_paths))
    print(args.out)
    print(args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Score the blinded ACE candidate universe without reading the blinding key."""

from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import (
    BLINDED_UNIVERSE_PATH,
    RANKED_CANDIDATES_PATH,
    SCORING_MANIFEST_PATH,
    SCORING_MODE_DESCRIPTIONS,
    ace_feature_row,
    rank_rows,
    read_csv_rows,
    save_json,
    score_features_for_mode,
    sha256_text,
    write_csv_rows,
)


FEATURE_COLUMNS = (
    "candidate_id",
    "rank",
    "sequence",
    "length",
    "score_mode",
    "ace_score",
    "ace_sar_score",
    "short_score",
    "cterm_anchor",
    "cterm_proline",
    "cterm_aromatic",
    "penultimate_proline",
    "antepenult_anchor",
    "nterm_hydrophobic",
    "second_basic_or_proline",
    "third_proline",
    "hydrophobic_fraction",
    "aromatic_fraction",
    "basic_fraction",
    "acidic_fraction",
    "proline_fraction",
    "hydrophobic_balance",
    "proline_balance",
    "net_charge",
    "charge_balance",
    "prodrug_triad",
    "acidic_penalty",
    "cysteine_penalty",
    "extreme_charge_penalty",
    "hydrophobicity",
    "frac_aromatic",
    "frac_proline",
    "frac_acidic",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universe", default=str(BLINDED_UNIVERSE_PATH))
    parser.add_argument("--out", default=str(RANKED_CANDIDATES_PATH))
    parser.add_argument("--manifest-out", default=str(SCORING_MANIFEST_PATH))
    parser.add_argument("--mode", choices=sorted(SCORING_MODE_DESCRIPTIONS), default="literature_sar")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    universe_path = Path(args.universe)
    rows = read_csv_rows(universe_path)
    scored: list[dict[str, object]] = []
    for row in rows:
        features = ace_feature_row(row["sequence"])
        features["candidate_id"] = row["candidate_id"]
        features["score_mode"] = args.mode
        features["ace_score"] = round(score_features_for_mode(features, args.mode), 6)
        scored.append(features)
    ranked = rank_rows(scored)
    write_csv_rows(Path(args.out), ranked, FEATURE_COLUMNS)

    universe_payload = universe_path.read_text(encoding="utf-8")
    save_json(
        Path(args.manifest_out),
        {
            "date": date.today().isoformat(),
            "universe": str(universe_path),
            "ranked_candidates": str(Path(args.out)),
            "candidate_count": len(ranked),
            "universe_sha256": sha256_text(universe_payload),
            "score_mode": args.mode,
            "score_mode_description": SCORING_MODE_DESCRIPTIONS[args.mode],
            "label_sources_read": [],
            "features": [
                "length preference favoring short ACE-like peptides",
                "C-terminal proline/hydrophobic/aromatic anchor",
                "penultimate proline",
                "antepenultimate hydrophobic/basic anchor",
                "N-terminal hydrophobic residue",
                "second-position basic/proline",
                "third-position proline",
                "hydrophobic, proline, and charge balance",
                "acidic, cysteine, and extreme-charge penalties",
            ],
            "claim_boundary": "Sequence-only score. No labels, PubMed evidence, or private blinding key were read during scoring.",
        },
    )
    print(f"wrote {len(ranked)} ranked ACE candidates with mode={args.mode}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

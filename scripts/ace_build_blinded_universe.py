#!/usr/bin/env python3
"""Build a blinded candidate universe for the ACE rediscovery pilot."""

from __future__ import annotations

import argparse
import random
import sys
from collections import Counter
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ace_rediscovery_common import (
    AA,
    BLINDED_UNIVERSE_PATH,
    BLINDING_KEY_PATH,
    DEFAULT_SOURCE_PROTEINS,
    KNOWN_POSITIVE_PATH,
    UNIVERSE_MANIFEST_PATH,
    candidate_id,
    clean_sequence,
    digest_sequence,
    fetch_uniprot_fasta,
    load_known_positives,
    random_decoys,
    random_peptide,
    save_json,
    sha256_text,
    shuffled_decoys,
    write_csv_rows,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(BLINDED_UNIVERSE_PATH))
    parser.add_argument("--key-out", default=str(BLINDING_KEY_PATH))
    parser.add_argument("--manifest-out", default=str(UNIVERSE_MANIFEST_PATH))
    parser.add_argument("--known-positive-csv", default=str(KNOWN_POSITIVE_PATH))
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--min-length", type=int, default=2)
    parser.add_argument("--max-length", type=int, default=12)
    parser.add_argument("--missed-cleavages", type=int, default=2)
    parser.add_argument("--random-background-count", type=int, default=12000)
    parser.add_argument("--decoys-per-positive", type=int, default=80)
    return parser.parse_args()


def add_candidate(
    *,
    sequence: str,
    label: str,
    origin: str,
    provenance: str,
    public_rows: list[dict[str, object]],
    key_rows: list[dict[str, object]],
    seen: set[str],
    label_counts: Counter[str] | None = None,
) -> bool:
    sequence = clean_sequence(sequence)
    if not sequence or sequence in seen:
        return False
    seen.add(sequence)
    cid = candidate_id(sequence)
    public_rows.append({"candidate_id": cid, "sequence": sequence, "length": len(sequence)})
    key_rows.append(
        {
            "candidate_id": cid,
            "sequence": sequence,
            "label": label,
            "origin": origin,
            "provenance": provenance,
        }
    )
    if label_counts is not None:
        label_counts[label] += 1
    return True


def main() -> int:
    args = parse_args()
    rng = random.Random(args.seed)
    known_positives = load_known_positives(Path(args.known_positive_csv))
    positive_sequences = {row["sequence"] for row in known_positives}
    public_rows: list[dict[str, object]] = []
    key_rows: list[dict[str, object]] = []
    seen: set[str] = set()
    live_label_counts: Counter[str] = Counter()
    source_status: list[dict[str, object]] = []
    digest_counts: Counter[str] = Counter()

    enzyme_recipes = ("trypsin", "chymotrypsin", "pepsin_like", "thermolysin_like")
    for source in DEFAULT_SOURCE_PROTEINS:
        try:
            header, protein_sequence = fetch_uniprot_fasta(source)
        except Exception as exc:  # pragma: no cover - network dependent
            source_status.append(
                {
                    "accession": source.accession,
                    "label": source.label,
                    "status": "fetch_failed",
                    "error": str(exc),
                }
            )
            continue
        before = len(public_rows)
        for enzyme in enzyme_recipes:
            for peptide in digest_sequence(
                protein_sequence,
                enzyme=enzyme,
                missed_cleavages=args.missed_cleavages,
                min_len=args.min_length,
                max_len=args.max_length,
            ):
                if set(peptide) - set(AA):
                    continue
                label = "digestome_known_positive" if peptide in positive_sequences else "digestome_background"
                add_candidate(
                    sequence=peptide,
                    label=label,
                    origin=f"{source.source_class}:{source.accession}:{enzyme}",
                    provenance=f"{source.label}; {header}",
                    public_rows=public_rows,
                    key_rows=key_rows,
                    seen=seen,
                    label_counts=live_label_counts,
                )
                digest_counts[f"{source.source_class}:{enzyme}"] += 1
        source_status.append(
            {
                "accession": source.accession,
                "label": source.label,
                "source_class": source.source_class,
                "organism": source.organism,
                "status": "ok",
                "sequence_length": len(protein_sequence),
                "unique_public_candidates_added": len(public_rows) - before,
                "uniprot_url": f"https://www.uniprot.org/uniprotkb/{source.accession}/entry",
            }
        )

    known_by_sequence = {row["sequence"]: row for row in known_positives}
    for sequence, row in known_by_sequence.items():
        add_candidate(
            sequence=sequence,
            label="known_positive",
            origin="blinded_known_positive_spike",
            provenance=f"{row.get('source_material', '')}; {row.get('reference_url', '')}",
            public_rows=public_rows,
            key_rows=key_rows,
            seen=seen,
            label_counts=live_label_counts,
        )

    for positive in sorted(positive_sequences):
        blocked = set(seen)
        for decoy in shuffled_decoys(positive, rng, max(12, args.decoys_per_positive // 4)):
            add_candidate(
                sequence=decoy,
                label="matched_decoy",
                origin=f"composition_shuffle:{positive}",
                provenance="private matched decoy",
                public_rows=public_rows,
                key_rows=key_rows,
                seen=seen,
                label_counts=live_label_counts,
            )
            blocked.add(decoy)
        remaining = args.decoys_per_positive - sum(
            1 for row in key_rows if row["label"] == "matched_decoy" and str(row["origin"]).endswith(f":{positive}")
        )
        for decoy in random_decoys(len(positive), rng, max(0, remaining), blocked):
            add_candidate(
                sequence=decoy,
                label="matched_decoy",
                origin=f"length_matched_random:{positive}",
                provenance="private matched decoy",
                public_rows=public_rows,
                key_rows=key_rows,
                seen=seen,
                label_counts=live_label_counts,
            )

    while live_label_counts["random_background"] < args.random_background_count:
        length = rng.randint(args.min_length, args.max_length)
        add_candidate(
            sequence=random_peptide(length, rng),
            label="random_background",
            origin="random_background",
            provenance="synthetic random peptide",
            public_rows=public_rows,
            key_rows=key_rows,
            seen=seen,
            label_counts=live_label_counts,
        )

    public_rows = sorted(public_rows, key=lambda row: str(row["candidate_id"]))
    key_rows = sorted(key_rows, key=lambda row: str(row["candidate_id"]))
    write_csv_rows(Path(args.out), public_rows, ("candidate_id", "sequence", "length"))
    write_csv_rows(Path(args.key_out), key_rows, ("candidate_id", "sequence", "label", "origin", "provenance"))

    label_counts = Counter(str(row["label"]) for row in key_rows)
    length_counts = Counter(int(row["length"]) for row in public_rows)
    public_payload = "\n".join(f"{row['candidate_id']},{row['sequence']},{row['length']}" for row in public_rows)
    key_payload = "\n".join(
        f"{row['candidate_id']},{row['sequence']},{row['label']},{row['origin']},{row['provenance']}" for row in key_rows
    )
    manifest = {
        "date": date.today().isoformat(),
        "seed": args.seed,
        "known_positive_source": str(Path(args.known_positive_csv)),
        "known_positive_sequences": sorted(positive_sequences),
        "public_rows": len(public_rows),
        "key_rows": len(key_rows),
        "label_counts_private": dict(sorted(label_counts.items())),
        "length_counts_public": {str(k): v for k, v in sorted(length_counts.items())},
        "digest_counts": dict(sorted(digest_counts.items())),
        "source_status": source_status,
        "blinding": {
            "public_columns": ["candidate_id", "sequence", "length"],
            "private_key_columns": ["candidate_id", "sequence", "label", "origin", "provenance"],
            "scoring_must_not_read": str(Path(args.key_out)),
        },
        "sha256": {
            "public_universe_payload": sha256_text(public_payload),
            "private_key_payload": sha256_text(key_payload),
        },
    }
    save_json(Path(args.manifest_out), manifest)
    print(
        f"wrote {len(public_rows)} blinded candidates, "
        f"{label_counts.get('known_positive', 0) + label_counts.get('digestome_known_positive', 0)} known positives"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

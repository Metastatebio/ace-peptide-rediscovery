#!/usr/bin/env python3
"""Verify the headline results from the public ACE reproducibility release."""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_LABEL_HASHES = {
    "v4-post-unblinding-labels.csv": "8d76847e95eebcb468063acf7f2f0b14d6307120103b9bac76015c391ed3fdc8",
    "v5-post-unblinding-labels.csv": "ab2c2ce9abfb0f4ad0ac70c6a22168c572169a8e642e5bb83538c5553a18b0c6",
    "v6-post-unblinding-labels.csv": "60baafe1bd5e477591afa1dcda304564cd6a322fe912e9d06b7336e56b53ecc8",
}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def average_ranks(values: list[float]) -> list[float]:
    ordered = sorted(range(len(values)), key=lambda index: values[index])
    output = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and values[ordered[end]] == values[ordered[start]]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for index in ordered[start:end]:
            output[index] = rank
        start = end
    return output


def auroc(labels: list[int], scores: list[float]) -> float:
    positives = sum(labels)
    negatives = len(labels) - positives
    if not positives or not negatives:
        raise ValueError("AUROC requires both positive and negative rows")
    ranks = average_ranks(scores)
    positive_rank_sum = sum(rank for rank, label in zip(ranks, labels) if label)
    return (positive_rank_sum - positives * (positives + 1) / 2.0) / (positives * negatives)


def assert_close(name: str, actual: float, expected: float, tolerance: float = 5e-5) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"{name}: expected {expected:.7f}, observed {actual:.7f}")


def main() -> int:
    labels_dir = ROOT / "data" / "labels-post-unblinding"
    for filename, expected in EXPECTED_LABEL_HASHES.items():
        observed = sha256(labels_dir / filename)
        if observed != expected:
            raise AssertionError(f"label hash mismatch for {filename}")

    balanced = rows(
        ROOT
        / "results"
        / "independent-validation"
        / "frozen-rerank"
        / "balanced-v6-candidate"
        / "frozen-rerank-ranked-candidates.csv"
    )
    positive_labels = {"known_positive", "digestome_known_positive"}
    labels = [int(row["label"] in positive_labels) for row in balanced]
    top_one = math.ceil(len(balanced) * 0.01)
    recovered = sum(label for label, row in zip(labels, balanced) if int(row["rank"]) <= top_one)
    matched = [row for row in balanced if row["label"] in positive_labels | {"matched_decoy"}]
    matched_labels = [int(row["label"] in positive_labels) for row in matched]
    matched_scores = [float(row["ace_score"]) for row in matched]
    assert len(balanced) == 58192
    assert sum(labels) == 52
    assert recovered == 6
    assert_close("v6 balanced matched-decoy AUROC", auroc(matched_labels, matched_scores), 0.7647235577)

    fused = rows(ROOT / "results" / "proteinlm" / "full-universe-rank" / "ace_plm_lengthz_full_universe_ranked_candidates.csv")
    fused_matched = [row for row in fused if row["label"] in positive_labels | {"matched_decoy"}]
    fused_labels = [int(row["label"] in positive_labels) for row in fused_matched]
    fused_scores = [float(row["fused_score"]) for row in fused_matched]
    assert_close("ESM2 fused matched-decoy AUROC", auroc(fused_labels, fused_scores), 0.7817492604)
    print("Release verification passed: corrected v6 and ProteinLM metrics match the public artifacts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Shared utilities for the ACE rediscovery pilot."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Iterable, Mapping, Sequence

import requests

try:
    from backend.app.analysis.features import charge, full_features, hydrophobicity
except ImportError:  # pragma: no cover - used when launched from scripts/
    import sys

    REPO_ROOT_FALLBACK = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(REPO_ROOT_FALLBACK))
    from backend.app.analysis.features import charge, full_features, hydrophobicity

REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO_ROOT / "docs" / "publication-engine" / "papers" / "01-ace-rediscovery"
KNOWN_POSITIVE_PATH = PAPER_DIR / "known-positive-seed.csv"
BLINDED_UNIVERSE_PATH = PAPER_DIR / "blinded-candidate-universe.csv"
BLINDING_KEY_PATH = PAPER_DIR / "blinding-key.private.csv"
UNIVERSE_MANIFEST_PATH = PAPER_DIR / "universe-manifest.json"
RANKED_CANDIDATES_PATH = PAPER_DIR / "ace-ranked-candidates.csv"
SCORING_MANIFEST_PATH = PAPER_DIR / "ace-scoring-manifest.json"
RESULTS_PATH = PAPER_DIR / "ace-rediscovery-results.csv"
REPORT_PATH = PAPER_DIR / "ace-rediscovery-report.md"
FIGURE_DIR = PAPER_DIR / "figures"
SOURCE_FASTA_DIR = PAPER_DIR / "source-fasta"

AA = "ACDEFGHIKLMNPQRSTVWY"
HYDROPHOBIC = set("AILMFWYV")
AROMATIC = set("FWY")
BASIC = set("KRH")
ACIDIC = set("DE")
CTERM_FAVORABLE = set("P" + "AILMFWYV")


@dataclass(frozen=True)
class SourceProtein:
    accession: str
    label: str
    source_class: str
    organism: str


DEFAULT_SOURCE_PROTEINS: tuple[SourceProtein, ...] = (
    SourceProtein("P02666", "bovine beta-casein", "milk", "Bos taurus"),
    SourceProtein("P02662", "bovine alpha-S1-casein", "milk", "Bos taurus"),
    SourceProtein("P02663", "bovine alpha-S2-casein", "milk", "Bos taurus"),
    SourceProtein("P02668", "bovine kappa-casein", "milk", "Bos taurus"),
    SourceProtein("P02754", "bovine beta-lactoglobulin", "milk", "Bos taurus"),
    SourceProtein("P02769", "bovine serum albumin", "milk", "Bos taurus"),
    SourceProtein("P01012", "chicken ovalbumin", "egg", "Gallus gallus"),
    SourceProtein("P00698", "chicken lysozyme C", "egg", "Gallus gallus"),
    SourceProtein("P02453", "bovine collagen alpha-1(I)", "collagen", "Bos taurus"),
    SourceProtein("Q9DGI8", "skipjack tuna myoglobin", "fish", "Katsuwonus pelamis"),
    SourceProtein("A9ZTF1", "skipjack tuna parvalbumin", "fish", "Katsuwonus pelamis"),
    SourceProtein("Q9IB34", "skipjack tuna myosin light chain 2", "fish", "Katsuwonus pelamis"),
    SourceProtein("Q9IB35", "skipjack tuna myosin light chain 3", "fish", "Katsuwonus pelamis"),
    SourceProtein("Q9IB36", "skipjack tuna myosin light chain 1", "fish", "Katsuwonus pelamis"),
)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def candidate_id(sequence: str, namespace: str = "ace-pilot-v1") -> str:
    digest = sha256_text(f"{namespace}:{sequence}")[:16]
    return f"ACE-{digest}"


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def save_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def load_known_positives(path: Path = KNOWN_POSITIVE_PATH) -> list[dict[str, str]]:
    rows = read_csv_rows(path)
    for row in rows:
        row["sequence"] = clean_sequence(row["sequence"])
    return rows


def clean_sequence(sequence: str) -> str:
    return re.sub(r"[^A-Z]", "", sequence.upper())


def parse_fasta(payload: str) -> tuple[str, str]:
    header = ""
    chunks: list[str] = []
    for line in payload.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if not header:
                header = line[1:]
            continue
        chunks.append(line)
    sequence = clean_sequence("".join(chunks))
    if not header or not sequence:
        raise ValueError("invalid FASTA payload")
    return header, sequence


def fetch_uniprot_fasta(source: SourceProtein, cache_dir: Path = SOURCE_FASTA_DIR) -> tuple[str, str]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{source.accession}.fasta"
    if cache_path.exists():
        return parse_fasta(cache_path.read_text(encoding="utf-8"))

    url = f"https://rest.uniprot.org/uniprotkb/{source.accession}.fasta"
    response = requests.get(url, timeout=25)
    response.raise_for_status()
    cache_path.write_text(response.text, encoding="utf-8")
    return parse_fasta(response.text)


def cleavage_sites(sequence: str, enzyme: str) -> list[int]:
    sequence = clean_sequence(sequence)
    sites = {0, len(sequence)}
    if enzyme == "none":
        return sorted(sites)
    for index, residue in enumerate(sequence):
        next_residue = sequence[index + 1] if index + 1 < len(sequence) else ""
        if enzyme == "trypsin" and residue in {"K", "R"} and next_residue != "P":
            sites.add(index + 1)
        elif enzyme == "chymotrypsin" and residue in {"F", "W", "Y", "L"} and next_residue != "P":
            sites.add(index + 1)
        elif enzyme == "pepsin_like" and residue in {"F", "L", "W", "Y"}:
            sites.add(index + 1)
        elif enzyme == "thermolysin_like" and index > 0 and residue in {"A", "F", "I", "L", "M", "V"}:
            sites.add(index)
    return sorted(sites)


def digest_sequence(
    sequence: str,
    *,
    enzyme: str,
    missed_cleavages: int,
    min_len: int,
    max_len: int,
) -> list[str]:
    cuts = cleavage_sites(sequence, enzyme)
    peptides: list[str] = []
    seen: set[str] = set()
    for start_index in range(len(cuts) - 1):
        for missed in range(missed_cleavages + 1):
            end_index = start_index + 1 + missed
            if end_index >= len(cuts):
                continue
            peptide = sequence[cuts[start_index] : cuts[end_index]]
            if min_len <= len(peptide) <= max_len and peptide not in seen:
                seen.add(peptide)
                peptides.append(peptide)
    return peptides


def shuffled_decoys(sequence: str, rng: random.Random, count: int) -> list[str]:
    sequence = clean_sequence(sequence)
    decoys: list[str] = []
    seen = {sequence}
    residues = list(sequence)
    max_attempts = max(200, count * 50)
    for _ in range(max_attempts):
        if len(decoys) >= count:
            break
        rng.shuffle(residues)
        decoy = "".join(residues)
        if decoy not in seen:
            seen.add(decoy)
            decoys.append(decoy)
    return decoys


def random_peptide(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(AA) for _ in range(length))


def random_decoys(length: int, rng: random.Random, count: int, blocked: set[str]) -> list[str]:
    decoys: list[str] = []
    max_attempts = max(200, count * 100)
    for _ in range(max_attempts):
        if len(decoys) >= count:
            break
        decoy = random_peptide(length, rng)
        if decoy in blocked:
            continue
        blocked.add(decoy)
        decoys.append(decoy)
    return decoys


def fraction(sequence: str, residues: set[str]) -> float:
    if not sequence:
        return 0.0
    return sum(1 for residue in sequence if residue in residues) / len(sequence)


def triangular_score(value: float, target: float, width: float) -> float:
    if width <= 0:
        return 0.0
    return clamp01(1.0 - abs(value - target) / width)


def clamp01(value: float) -> float:
    if math.isnan(value):
        return 0.0
    return min(1.0, max(0.0, value))


def length_preference(length: int) -> float:
    if length == 3:
        return 1.0
    if length == 2:
        return 0.78
    if length == 4:
        return 0.82
    if length == 5:
        return 0.72
    if 6 <= length <= 8:
        return 0.46
    if 9 <= length <= 12:
        return 0.22
    return 0.0


def bioactive_length_score(length: int) -> float:
    if length < 2 or length > 14:
        return 0.0
    return triangular_score(float(length), 5.0, 7.0)


def _feature_value(features: Mapping[str, object], key: str) -> float:
    try:
        return float(features.get(key, 0.0) or 0.0)
    except (TypeError, ValueError):
        return 0.0


def generic_zero_shot_score(features: Mapping[str, object]) -> float:
    """Target-agnostic peptide-likeness score with no ACE-position priors."""
    length = int(_feature_value(features, "length"))
    aromatic_fraction = _feature_value(features, "aromatic_fraction")
    acidic_penalty = _feature_value(features, "acidic_penalty")
    cysteine_penalty = _feature_value(features, "cysteine_penalty")
    extreme_charge_penalty = _feature_value(features, "extreme_charge_penalty")
    score = (
        0.20 * bioactive_length_score(length)
        + 0.20 * _feature_value(features, "hydrophobic_balance")
        + 0.18 * _feature_value(features, "charge_balance")
        + 0.16 * _feature_value(features, "proline_balance")
        + 0.12 * triangular_score(aromatic_fraction, 0.12, 0.22)
        + 0.08 * (1.0 - acidic_penalty)
        + 0.06 * (1.0 - cysteine_penalty)
        - 0.10 * extreme_charge_penalty
    )
    return clamp01(score)


def oral_stability_zero_shot_score(features: Mapping[str, object]) -> float:
    """Generic oral-peptide stability score without target labels."""
    length = int(_feature_value(features, "length"))
    score = (
        0.18 * bioactive_length_score(length)
        + 0.18 * _feature_value(features, "hydrophobic_balance")
        + 0.16 * _feature_value(features, "proline_balance")
        + 0.12 * _feature_value(features, "charge_balance")
        + 0.10 * _feature_value(features, "cterm_proline")
        + 0.07 * _feature_value(features, "penultimate_proline")
        + 0.07 * _feature_value(features, "nterm_hydrophobic")
        + 0.06 * (1.0 - _feature_value(features, "acidic_penalty"))
        + 0.06 * (1.0 - _feature_value(features, "cysteine_penalty"))
        - 0.08 * _feature_value(features, "extreme_charge_penalty")
    )
    return clamp01(score)


def ace_pharmacophore_zero_label_score(features: Mapping[str, object]) -> float:
    """ACE target-informed pharmacophore score without peptide activity labels."""
    score = (
        0.16 * _feature_value(features, "short_score")
        + 0.15 * _feature_value(features, "cterm_anchor")
        + 0.12 * _feature_value(features, "cterm_proline")
        + 0.07 * _feature_value(features, "cterm_aromatic")
        + 0.11 * _feature_value(features, "penultimate_proline")
        + 0.08 * _feature_value(features, "antepenult_anchor")
        + 0.08 * _feature_value(features, "nterm_hydrophobic")
        + 0.07 * _feature_value(features, "second_basic_or_proline")
        + 0.06 * _feature_value(features, "third_proline")
        + 0.08 * _feature_value(features, "hydrophobic_balance")
        + 0.07 * _feature_value(features, "proline_balance")
        + 0.05 * _feature_value(features, "charge_balance")
        + 0.10 * _feature_value(features, "prodrug_triad")
        - 0.10 * _feature_value(features, "acidic_penalty")
        - 0.16 * _feature_value(features, "cysteine_penalty")
        - 0.08 * _feature_value(features, "extreme_charge_penalty")
    )
    return clamp01(score)


def no_position_prior_score(features: Mapping[str, object]) -> float:
    """Control score that removes N/C-terminal position priors."""
    return clamp01(
        0.22 * bioactive_length_score(int(_feature_value(features, "length")))
        + 0.22 * _feature_value(features, "hydrophobic_balance")
        + 0.20 * _feature_value(features, "proline_balance")
        + 0.16 * _feature_value(features, "charge_balance")
        + 0.12 * triangular_score(_feature_value(features, "aromatic_fraction"), 0.12, 0.22)
        - 0.12 * _feature_value(features, "acidic_penalty")
        - 0.18 * _feature_value(features, "cysteine_penalty")
        - 0.10 * _feature_value(features, "extreme_charge_penalty")
    )


def ace_no_exact_proline_motif_score(features: Mapping[str, object]) -> float:
    """ACE target-informed score with direct terminal-proline motif priors removed."""
    sequence = str(features.get("sequence") or "")
    cterm = sequence[-1] if sequence else ""
    second = sequence[1] if len(sequence) > 1 else ""
    cterm_nonproline_hydrophobic = 1.0 if cterm in (HYDROPHOBIC - {"P"}) else 0.0
    second_basic = 1.0 if second in BASIC else 0.0
    score = (
        0.18 * _feature_value(features, "short_score")
        + 0.14 * cterm_nonproline_hydrophobic
        + 0.10 * _feature_value(features, "cterm_aromatic")
        + 0.12 * _feature_value(features, "antepenult_anchor")
        + 0.10 * _feature_value(features, "nterm_hydrophobic")
        + 0.08 * second_basic
        + 0.12 * _feature_value(features, "hydrophobic_balance")
        + 0.10 * _feature_value(features, "proline_balance")
        + 0.08 * _feature_value(features, "charge_balance")
        - 0.10 * _feature_value(features, "acidic_penalty")
        - 0.16 * _feature_value(features, "cysteine_penalty")
        - 0.08 * _feature_value(features, "extreme_charge_penalty")
    )
    return clamp01(score)


SCORING_MODE_DESCRIPTIONS = {
    "literature_sar": "ACE-informed SAR score used in the first pilot; no labels are read during scoring but priors are ACE-specific.",
    "ace_pharmacophore_zero_label": "ACE target-informed pharmacophore score; no ACE peptide activity labels are read during scoring.",
    "ace_no_exact_proline_motif": "ACE target-informed score with direct terminal-proline motif priors removed; no ACE peptide labels are read during scoring.",
    "oral_stability_zero_shot": "Target-agnostic oral peptide stability/developability score; no ACE target or activity labels.",
    "generic_zero_shot": "Target-agnostic generic peptide-likeness score; no ACE target, SAR, or activity labels.",
    "no_position_priors": "Negative-control score using global composition and no terminal-position priors.",
}


def score_features_for_mode(features: Mapping[str, object], mode: str) -> float:
    if mode == "literature_sar":
        return _feature_value(features, "ace_sar_score")
    if mode == "ace_pharmacophore_zero_label":
        return ace_pharmacophore_zero_label_score(features)
    if mode == "ace_no_exact_proline_motif":
        return ace_no_exact_proline_motif_score(features)
    if mode == "oral_stability_zero_shot":
        return oral_stability_zero_shot_score(features)
    if mode == "generic_zero_shot":
        return generic_zero_shot_score(features)
    if mode == "no_position_priors":
        return no_position_prior_score(features)
    raise ValueError(f"unknown ACE scoring mode: {mode}")


def ace_feature_row(sequence: str) -> dict[str, float | int | str]:
    sequence = clean_sequence(sequence)
    length = len(sequence)
    features = full_features(sequence)
    counts = Counter(sequence)
    hydrophobic_fraction = fraction(sequence, HYDROPHOBIC)
    aromatic_fraction = fraction(sequence, AROMATIC)
    basic_fraction = fraction(sequence, BASIC)
    acidic_fraction = fraction(sequence, ACIDIC)
    proline_fraction = counts["P"] / length if length else 0.0
    net_charge = charge(sequence)
    kd = hydrophobicity(sequence)
    nterm = sequence[0] if sequence else ""
    second = sequence[1] if len(sequence) > 1 else ""
    third = sequence[2] if len(sequence) > 2 else ""
    cterm = sequence[-1] if sequence else ""
    penultimate = sequence[-2] if len(sequence) > 1 else ""
    antepenult = sequence[-3] if len(sequence) > 2 else ""

    short_score = length_preference(length)
    cterm_anchor = 1.0 if cterm in CTERM_FAVORABLE else 0.0
    cterm_proline = 1.0 if cterm == "P" else 0.0
    cterm_aromatic = 1.0 if cterm in AROMATIC else 0.0
    penultimate_proline = 1.0 if penultimate == "P" else 0.0
    antepenult_anchor = 1.0 if antepenult in HYDROPHOBIC | BASIC else 0.0
    nterm_hydrophobic = 1.0 if nterm in HYDROPHOBIC else 0.0
    second_basic_or_proline = 1.0 if second in BASIC | {"P"} else 0.0
    third_proline = 1.0 if third == "P" else 0.0
    hydrophobic_balance = triangular_score(hydrophobic_fraction, 0.38, 0.38)
    proline_balance = triangular_score(proline_fraction, 0.32, 0.32)
    charge_balance = triangular_score(abs(net_charge), 0.45, 1.25)
    prodrug_triad = nterm_hydrophobic * second_basic_or_proline * third_proline * (1.0 if 3 <= length <= 6 else 0.0)
    acidic_penalty = acidic_fraction
    cysteine_penalty = 1.0 if counts["C"] else 0.0
    extreme_charge_penalty = max(0.0, abs(net_charge) - 2.0) / 3.0

    score = (
        0.18 * short_score
        + 0.16 * cterm_anchor
        + 0.12 * cterm_proline
        + 0.07 * cterm_aromatic
        + 0.12 * penultimate_proline
        + 0.08 * antepenult_anchor
        + 0.09 * nterm_hydrophobic
        + 0.08 * second_basic_or_proline
        + 0.07 * third_proline
        + 0.10 * hydrophobic_balance
        + 0.09 * proline_balance
        + 0.05 * charge_balance
        + 0.12 * prodrug_triad
        - 0.12 * acidic_penalty
        - 0.18 * cysteine_penalty
        - 0.10 * extreme_charge_penalty
    )

    output: dict[str, float | int | str] = {
        "sequence": sequence,
        "length": length,
        "ace_score": round(clamp01(score), 6),
        "ace_sar_score": round(clamp01(score), 6),
        "short_score": round(short_score, 6),
        "cterm_anchor": round(cterm_anchor, 6),
        "cterm_proline": round(cterm_proline, 6),
        "cterm_aromatic": round(cterm_aromatic, 6),
        "penultimate_proline": round(penultimate_proline, 6),
        "antepenult_anchor": round(antepenult_anchor, 6),
        "nterm_hydrophobic": round(nterm_hydrophobic, 6),
        "second_basic_or_proline": round(second_basic_or_proline, 6),
        "third_proline": round(third_proline, 6),
        "hydrophobic_fraction": round(hydrophobic_fraction, 6),
        "aromatic_fraction": round(aromatic_fraction, 6),
        "basic_fraction": round(basic_fraction, 6),
        "acidic_fraction": round(acidic_fraction, 6),
        "proline_fraction": round(proline_fraction, 6),
        "hydrophobic_balance": round(hydrophobic_balance, 6),
        "proline_balance": round(proline_balance, 6),
        "net_charge": round(net_charge, 6),
        "charge_balance": round(charge_balance, 6),
        "prodrug_triad": round(prodrug_triad, 6),
        "acidic_penalty": round(acidic_penalty, 6),
        "cysteine_penalty": round(cysteine_penalty, 6),
        "extreme_charge_penalty": round(extreme_charge_penalty, 6),
        "hydrophobicity": round(kd, 6),
        "frac_aromatic": round(float(features.get("frac_aromatic", 0.0)), 6),
        "frac_proline": round(float(features.get("frac_proline", 0.0)), 6),
        "frac_acidic": round(float(features.get("frac_acidic", 0.0)), 6),
    }
    for key, value in features.items():
        try:
            output[f"physchem_{key}"] = round(float(value), 6)
        except (TypeError, ValueError):
            continue
    return output


def rank_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    ranked = sorted(
        (dict(row) for row in rows),
        key=lambda row: (-float(row["ace_score"]), int(row["length"]), str(row["sequence"])),
    )
    for index, row in enumerate(ranked, start=1):
        row["rank"] = index
    return ranked


def mean_rank(values: Sequence[int]) -> float:
    return mean(values) if values else float("nan")


def auroc_from_scores(positive_scores: Sequence[float], negative_scores: Sequence[float]) -> float:
    if not positive_scores or not negative_scores:
        return float("nan")
    wins = 0.0
    total = 0
    for positive in positive_scores:
        for negative in negative_scores:
            total += 1
            if positive > negative:
                wins += 1.0
            elif positive == negative:
                wins += 0.5
    return wins / total if total else float("nan")


def average_precision(labels_and_scores: Sequence[tuple[int, float]]) -> float:
    ranked = sorted(labels_and_scores, key=lambda item: -item[1])
    positives = sum(label for label, _score in ranked)
    if positives == 0:
        return float("nan")
    hits = 0
    precision_sum = 0.0
    for index, (label, _score) in enumerate(ranked, start=1):
        if label:
            hits += 1
            precision_sum += hits / index
    return precision_sum / positives

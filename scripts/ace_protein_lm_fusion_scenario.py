#!/usr/bin/env python3
"""Run an ACE rediscovery ProteinLM fusion and masked-variant scenario."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import re
import sys
import time
from collections import Counter
from datetime import date
from itertools import combinations
from pathlib import Path
from statistics import mean, median
from typing import Any, Mapping, Sequence

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO_ROOT / "docs" / "publication-engine" / "papers" / "01-ace-rediscovery"
V6_EXPERIMENT_DIR = PAPER_DIR / "experiments" / "2026-08-29-v6-independent-validation-001"
DEFAULT_RANKED = (
    V6_EXPERIMENT_DIR
    / "frozen-rerank"
    / "balanced-v6-candidate"
    / "frozen-rerank-ranked-candidates.csv"
)
DEFAULT_SUMMARY = (
    V6_EXPERIMENT_DIR
    / "frozen-rerank"
    / "balanced-v6-candidate"
    / "frozen-rerank-summary.csv"
)
DEFAULT_POSITIVE_RANKS = (
    V6_EXPERIMENT_DIR
    / "frozen-rerank"
    / "balanced-v6-candidate"
    / "frozen-rerank-positive-ranks.csv"
)
DEFAULT_OUT_DIR = (
    PAPER_DIR
    / "experiments"
    / "2026-08-30-ace-protein-lm-fusion-001"
)
AA_POOL = "ACDEFGHIKLMNPQRSTVWY"
POSITIVE_LABELS = {"known_positive", "digestome_known_positive"}
EVAL_NEGATIVE_LABEL = "matched_decoy"

SCORED_FIELDS = (
    "scenario_group",
    "candidate_id",
    "sequence",
    "length",
    "label",
    "origin",
    "provenance",
    "ace_rank",
    "ace_score",
    "plm_total_log_likelihood",
    "plm_mean_log_likelihood",
    "plm_pseudo_perplexity",
    "mol_wt",
    "mol_logp",
    "tpsa",
    "rotatable_bonds",
)
FUSION_FIELDS = (
    "fusion_id",
    "ace_weight",
    "plm_weight",
    "eval_rows",
    "positive_count",
    "matched_decoy_count",
    "top_1pct_cutoff",
    "top_1pct_count",
    "top_1pct_fraction",
    "top_1pct_enrichment_over_random",
    "top_1pct_hypergeom_p",
    "top_5pct_count",
    "top_10pct_count",
    "best_positive_rank",
    "median_positive_rank",
    "mean_positive_rank",
    "auroc_vs_matched_decoys",
    "average_precision_vs_matched_decoys",
)
RANKED_FIELDS = (
    "rank",
    "candidate_id",
    "sequence",
    "length",
    "label",
    "scenario_group",
    "ace_rank",
    "ace_score",
    "plm_mean_log_likelihood",
    "fused_score",
    "origin",
    "provenance",
)
VARIANT_FIELDS = (
    "variant_id",
    "rank",
    "parent_sequence",
    "parent_v6_rank",
    "parent_plm_mean_log_likelihood",
    "edit_count",
    "mutations",
    "sequence",
    "length",
    "plm_mean_log_likelihood",
    "delta_parent_plm_mean",
    "plm_pseudo_perplexity",
    "masked_lm_support_mean",
    "nearest_known_positive",
    "nearest_known_positive_distance",
    "mol_wt",
    "mol_logp",
    "tpsa",
    "rotatable_bonds",
    "developability_flag",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranked", default=str(DEFAULT_RANKED))
    parser.add_argument("--summary", default=str(DEFAULT_SUMMARY))
    parser.add_argument("--positive-ranks", default=str(DEFAULT_POSITIVE_RANKS))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--model-id", default="facebook/esm2_t6_8M_UR50D")
    parser.add_argument("--seed", type=int, default=20260830)
    parser.add_argument("--top-ace-background", type=int, default=1000)
    parser.add_argument("--random-background", type=int, default=1000)
    parser.add_argument("--shell-parents", type=int, default=8)
    parser.add_argument("--single-substitutions-per-position", type=int, default=4)
    parser.add_argument("--double-substitutions-per-position", type=int, default=3)
    parser.add_argument("--variant-top-n", type=int, default=120)
    parser.add_argument("--progress-every", type=int, default=250)
    parser.add_argument("--max-panel-rows", type=int, default=0)
    return parser.parse_args()


def clean_sequence(sequence: str) -> str:
    return re.sub(r"[^A-Z]", "", sequence.upper())


def stable_id(prefix: str, text: str) -> str:
    return f"{prefix}-{hashlib.sha256(text.encode('utf-8')).hexdigest()[:16]}"


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


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def normalize(values: Sequence[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return arr
    minimum = float(np.min(arr))
    maximum = float(np.max(arr))
    if math.isclose(minimum, maximum):
        return np.full_like(arr, 0.5)
    return (arr - minimum) / (maximum - minimum)


def hypergeom_tail(
    *,
    population_size: int,
    success_states: int,
    draws: int,
    observed_successes: int,
) -> float:
    if observed_successes <= 0:
        return 1.0
    upper = min(success_states, draws)
    denominator = (
        math.lgamma(population_size + 1)
        - math.lgamma(draws + 1)
        - math.lgamma(population_size - draws + 1)
    )
    terms = []
    for successes in range(observed_successes, upper + 1):
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


def pearson(left: Sequence[float], right: Sequence[float]) -> float:
    left_arr = np.asarray(left, dtype=float)
    right_arr = np.asarray(right, dtype=float)
    if left_arr.size < 2 or right_arr.size < 2:
        return float("nan")
    if np.std(left_arr) == 0 or np.std(right_arr) == 0:
        return float("nan")
    return float(np.corrcoef(left_arr, right_arr)[0, 1])


def ranks(values: Sequence[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    order = np.argsort(arr)
    ranked = np.empty(len(arr), dtype=float)
    index = 0
    while index < len(arr):
        end = index + 1
        while end < len(arr) and arr[order[end]] == arr[order[index]]:
            end += 1
        average_rank = (index + 1 + end) / 2.0
        ranked[order[index:end]] = average_rank
        index = end
    return ranked


def spearman(left: Sequence[float], right: Sequence[float]) -> float:
    return pearson(ranks(left), ranks(right))


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


def nearest_sequence(sequence: str, references: Sequence[str]) -> tuple[str, int]:
    if not references:
        return "", 0
    return min(
        ((reference, edit_distance(sequence, reference)) for reference in references),
        key=lambda item: (item[1], abs(len(sequence) - len(item[0])), item[0]),
    )


def molecular_descriptors(sequence: str) -> dict[str, Any]:
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
    except ImportError:
        return {"mol_wt": "", "mol_logp": "", "tpsa": "", "rotatable_bonds": ""}

    mol = Chem.MolFromFASTA(sequence)
    if mol is None:
        return {"mol_wt": "", "mol_logp": "", "tpsa": "", "rotatable_bonds": ""}
    return {
        "mol_wt": round(float(Descriptors.MolWt(mol)), 4),
        "mol_logp": round(float(Crippen.MolLogP(mol)), 4),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 4),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
    }


def developability_flag(descriptors: Mapping[str, Any], sequence: str) -> str:
    mw = safe_float(descriptors.get("mol_wt"), math.nan)
    tpsa = safe_float(descriptors.get("tpsa"), math.nan)
    logp = safe_float(descriptors.get("mol_logp"), math.nan)
    if math.isnan(mw) or math.isnan(tpsa) or math.isnan(logp):
        return "descriptor_unavailable"
    if not (250.0 <= mw <= 1500.0):
        return "mw_outside_peptide_screen"
    if tpsa > 650.0:
        return "high_tpsa"
    if logp > 6.0:
        return "high_logp"
    if any(sequence.count(residue * 4) for residue in AA_POOL):
        return "homopolymer_run"
    return "screen_pass"


class EsmMaskedLmScorer:
    def __init__(self, model_id: str) -> None:
        try:
            import torch
            from transformers import AutoModelForMaskedLM, AutoTokenizer
        except ImportError as exc:
            raise RuntimeError("torch and transformers are required for this scenario") from exc

        self.torch = torch
        self.model_id = model_id
        self.tokenizer = AutoTokenizer.from_pretrained(model_id)
        self.model = AutoModelForMaskedLM.from_pretrained(model_id)
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.model.to(self.device)
        self.model.eval()
        self.mask_id = int(self.tokenizer.mask_token_id)
        self.aa_token_ids = {
            residue: int(self.tokenizer.convert_tokens_to_ids(residue)) for residue in AA_POOL
        }
        missing = [residue for residue, token_id in self.aa_token_ids.items() if token_id is None or token_id < 0]
        if missing:
            raise RuntimeError(f"model tokenizer does not expose amino acid tokens: {missing}")

    def _encoded(self, sequence: str) -> Mapping[str, Any]:
        return self.tokenizer(sequence, return_tensors="pt")

    def score_sequence(self, sequence: str) -> dict[str, float]:
        sequence = clean_sequence(sequence)
        if not sequence or any(residue not in AA_POOL for residue in sequence):
            return {
                "plm_total_log_likelihood": float("-inf"),
                "plm_mean_log_likelihood": float("-inf"),
                "plm_pseudo_perplexity": float("inf"),
            }

        encoded = self._encoded(sequence)
        input_ids = encoded["input_ids"][0].to(self.device)
        attention_mask = encoded["attention_mask"][0].to(self.device)
        residue_positions = list(range(1, len(sequence) + 1))
        if input_ids.shape[0] < len(sequence) + 2:
            raise RuntimeError(f"unexpected tokenization for sequence {sequence!r}")

        batch_ids = input_ids.repeat(len(residue_positions), 1)
        batch_attention = attention_mask.repeat(len(residue_positions), 1)
        for row_index, token_position in enumerate(residue_positions):
            batch_ids[row_index, token_position] = self.mask_id

        with self.torch.no_grad():
            logits = self.model(input_ids=batch_ids, attention_mask=batch_attention).logits
            position_logits = logits[
                self.torch.arange(len(residue_positions), device=self.device),
                self.torch.tensor(residue_positions, device=self.device),
                :,
            ]
            log_probs = self.torch.log_softmax(position_logits, dim=-1)
            actual_ids = input_ids[self.torch.tensor(residue_positions, device=self.device)]
            selected = log_probs[
                self.torch.arange(len(residue_positions), device=self.device),
                actual_ids,
            ]
        total = float(selected.sum().detach().cpu())
        mean_ll = total / float(len(sequence))
        return {
            "plm_total_log_likelihood": total,
            "plm_mean_log_likelihood": mean_ll,
            "plm_pseudo_perplexity": float(math.exp(-mean_ll)) if math.isfinite(mean_ll) else float("inf"),
        }

    def masked_substitutions(self, sequence: str, position: int, top_n: int) -> list[dict[str, Any]]:
        sequence = clean_sequence(sequence)
        encoded = self._encoded(sequence)
        input_ids = encoded["input_ids"].to(self.device)
        attention_mask = encoded["attention_mask"].to(self.device)
        token_position = position + 1
        input_ids = input_ids.clone()
        input_ids[0, token_position] = self.mask_id
        with self.torch.no_grad():
            logits = self.model(input_ids=input_ids, attention_mask=attention_mask).logits[0, token_position, :]
            log_probs = self.torch.log_softmax(logits, dim=-1)

        rows: list[dict[str, Any]] = []
        current = sequence[position]
        for residue in AA_POOL:
            if residue == current:
                continue
            token_id = self.aa_token_ids[residue]
            rows.append({"residue": residue, "log_probability": float(log_probs[token_id].detach().cpu())})
        rows.sort(key=lambda row: float(row["log_probability"]), reverse=True)
        return rows[:top_n]


def load_ranked_rows(path: Path) -> list[dict[str, Any]]:
    rows = []
    for row in read_csv_rows(path):
        sequence = clean_sequence(row["sequence"])
        rows.append(
            {
                **row,
                "sequence": sequence,
                "length": len(sequence),
                "ace_rank": safe_int(row.get("rank")),
                "ace_score": safe_float(row.get("ace_score")),
            }
        )
    return rows


def select_scoring_panel(
    rows: Sequence[Mapping[str, Any]],
    *,
    top_ace_background: int,
    random_background: int,
    max_panel_rows: int,
    seed: int,
) -> list[dict[str, Any]]:
    rng = random.Random(seed)
    positives = [dict(row, scenario_group="matched_eval") for row in rows if row.get("label") in POSITIVE_LABELS]
    matched = [dict(row, scenario_group="matched_eval") for row in rows if row.get("label") == EVAL_NEGATIVE_LABEL]
    if max_panel_rows and max_panel_rows > len(positives):
        keep_negatives = max_panel_rows - len(positives)
        matched = matched[:keep_negatives]

    top_background_rows: list[dict[str, Any]] = []
    for row in rows[: max(0, top_ace_background)]:
        if row.get("label") in POSITIVE_LABELS or row.get("label") == EVAL_NEGATIVE_LABEL:
            continue
        top_background_rows.append(dict(row, scenario_group="top_ace_background"))

    background_pool = [
        row
        for row in rows
        if row.get("label") not in POSITIVE_LABELS
        and row.get("label") != EVAL_NEGATIVE_LABEL
        and safe_int(row.get("rank")) > top_ace_background
    ]
    random_rows = rng.sample(background_pool, min(max(0, random_background), len(background_pool)))
    random_rows = [dict(row, scenario_group="random_background_probe") for row in random_rows]

    by_candidate: dict[str, dict[str, Any]] = {}
    for row in positives + matched + top_background_rows + random_rows:
        by_candidate[str(row["candidate_id"])] = dict(row)
    return list(by_candidate.values())


def score_rows(
    rows: Sequence[Mapping[str, Any]],
    scorer: EsmMaskedLmScorer,
    *,
    progress_every: int,
) -> list[dict[str, Any]]:
    score_cache: dict[str, dict[str, float]] = {}
    descriptor_cache: dict[str, dict[str, Any]] = {}
    output_rows = []
    started = time.monotonic()
    for index, row in enumerate(rows, start=1):
        sequence = str(row["sequence"])
        if sequence not in score_cache:
            score_cache[sequence] = scorer.score_sequence(sequence)
        if sequence not in descriptor_cache:
            descriptor_cache[sequence] = molecular_descriptors(sequence)
        scored = {
            "scenario_group": row.get("scenario_group", ""),
            "candidate_id": row.get("candidate_id", ""),
            "sequence": sequence,
            "length": len(sequence),
            "label": row.get("label", ""),
            "origin": row.get("origin", ""),
            "provenance": row.get("provenance", ""),
            "ace_rank": row.get("ace_rank", row.get("rank", "")),
            "ace_score": row.get("ace_score", ""),
            **score_cache[sequence],
            **descriptor_cache[sequence],
        }
        scored["plm_total_log_likelihood"] = round(float(scored["plm_total_log_likelihood"]), 8)
        scored["plm_mean_log_likelihood"] = round(float(scored["plm_mean_log_likelihood"]), 8)
        scored["plm_pseudo_perplexity"] = round(float(scored["plm_pseudo_perplexity"]), 8)
        output_rows.append(scored)
        if progress_every > 0 and index % progress_every == 0:
            elapsed = time.monotonic() - started
            print(
                f"[ace-plm] scored {index}/{len(rows)} rows in {elapsed:.1f}s "
                f"({index / max(elapsed, 1e-9):.1f} rows/s)",
                flush=True,
            )
    return output_rows


def fusion_metrics(
    rows: Sequence[Mapping[str, Any]],
    *,
    ace_weight: float,
) -> dict[str, Any]:
    eval_rows = [
        row
        for row in rows
        if row.get("label") in POSITIVE_LABELS or row.get("label") == EVAL_NEGATIVE_LABEL
    ]
    labels = np.asarray([1 if row.get("label") in POSITIVE_LABELS else 0 for row in eval_rows], dtype=int)
    ace_norm = normalize([safe_float(row.get("ace_score")) for row in eval_rows])
    plm_norm = normalize([safe_float(row.get("plm_mean_log_likelihood")) for row in eval_rows])
    scores = (ace_weight * ace_norm) + ((1.0 - ace_weight) * plm_norm)
    order = sorted(range(len(eval_rows)), key=lambda index: (-float(scores[index]), str(eval_rows[index]["sequence"])))
    positive_ranks = [rank for rank, index in enumerate(order, start=1) if int(labels[index]) == 1]
    cutoff_1 = max(1, math.ceil(len(eval_rows) * 0.01))
    cutoff_5 = max(1, math.ceil(len(eval_rows) * 0.05))
    cutoff_10 = max(1, math.ceil(len(eval_rows) * 0.10))
    top_1_count = sum(1 for rank in positive_ranks if rank <= cutoff_1)
    top_5_count = sum(1 for rank in positive_ranks if rank <= cutoff_5)
    top_10_count = sum(1 for rank in positive_ranks if rank <= cutoff_10)
    positive_count = int(labels.sum())
    random_top_fraction = cutoff_1 / float(len(eval_rows))
    return {
        "fusion_id": (
            "ace_only" if math.isclose(ace_weight, 1.0) else "plm_only" if math.isclose(ace_weight, 0.0) else f"ace_{ace_weight:.2f}_plm_{1.0 - ace_weight:.2f}"
        ),
        "ace_weight": round(ace_weight, 4),
        "plm_weight": round(1.0 - ace_weight, 4),
        "eval_rows": len(eval_rows),
        "positive_count": positive_count,
        "matched_decoy_count": len(eval_rows) - positive_count,
        "top_1pct_cutoff": cutoff_1,
        "top_1pct_count": top_1_count,
        "top_1pct_fraction": top_1_count / float(max(positive_count, 1)),
        "top_1pct_enrichment_over_random": (
            (top_1_count / float(max(positive_count, 1))) / random_top_fraction
            if random_top_fraction
            else 0.0
        ),
        "top_1pct_hypergeom_p": hypergeom_tail(
            population_size=len(eval_rows),
            success_states=positive_count,
            draws=cutoff_1,
            observed_successes=top_1_count,
        ),
        "top_5pct_count": top_5_count,
        "top_10pct_count": top_10_count,
        "best_positive_rank": min(positive_ranks) if positive_ranks else "",
        "median_positive_rank": median(positive_ranks) if positive_ranks else "",
        "mean_positive_rank": mean(positive_ranks) if positive_ranks else "",
        "auroc_vs_matched_decoys": float(roc_auc_score(labels, scores)) if len(set(labels)) == 2 else "",
        "average_precision_vs_matched_decoys": (
            float(average_precision_score(labels, scores)) if len(set(labels)) == 2 else ""
        ),
    }


def best_fusion_rows(rows: Sequence[Mapping[str, Any]], ace_weight: float) -> list[dict[str, Any]]:
    ace_norm = normalize([safe_float(row.get("ace_score")) for row in rows])
    plm_norm = normalize([safe_float(row.get("plm_mean_log_likelihood")) for row in rows])
    fused = (ace_weight * ace_norm) + ((1.0 - ace_weight) * plm_norm)
    order = sorted(range(len(rows)), key=lambda index: (-float(fused[index]), str(rows[index]["sequence"])))
    ranked = []
    for rank, index in enumerate(order, start=1):
        row = rows[index]
        ranked.append(
            {
                "rank": rank,
                "candidate_id": row.get("candidate_id", ""),
                "sequence": row.get("sequence", ""),
                "length": row.get("length", ""),
                "label": row.get("label", ""),
                "scenario_group": row.get("scenario_group", ""),
                "ace_rank": row.get("ace_rank", ""),
                "ace_score": row.get("ace_score", ""),
                "plm_mean_log_likelihood": row.get("plm_mean_log_likelihood", ""),
                "fused_score": round(float(fused[index]), 10),
                "origin": row.get("origin", ""),
                "provenance": row.get("provenance", ""),
            }
        )
    return ranked


def load_known_positive_sequences() -> list[str]:
    sequences: set[str] = set()
    for path in PAPER_DIR.glob("known-positive*.csv"):
        for row in read_csv_rows(path):
            sequence = clean_sequence(row.get("sequence", ""))
            if sequence:
                sequences.add(sequence)
    return sorted(sequences)


def generate_masked_variants(
    *,
    scorer: EsmMaskedLmScorer,
    parent_rows: Sequence[Mapping[str, Any]],
    blocked_sequences: set[str],
    known_positive_sequences: Sequence[str],
    single_substitutions_per_position: int,
    double_substitutions_per_position: int,
    variant_top_n: int,
) -> list[dict[str, Any]]:
    score_cache: dict[str, dict[str, float]] = {}
    descriptor_cache: dict[str, dict[str, Any]] = {}
    variants: dict[str, dict[str, Any]] = {}

    for parent in parent_rows:
        parent_sequence = str(parent["sequence"])
        if parent_sequence not in score_cache:
            score_cache[parent_sequence] = scorer.score_sequence(parent_sequence)
        substitutions_by_position: dict[int, list[dict[str, Any]]] = {}
        for position in range(len(parent_sequence)):
            substitutions_by_position[position] = scorer.masked_substitutions(
                parent_sequence,
                position,
                max(single_substitutions_per_position, double_substitutions_per_position),
            )
            for substitution in substitutions_by_position[position][:single_substitutions_per_position]:
                residue = str(substitution["residue"])
                sequence = parent_sequence[:position] + residue + parent_sequence[position + 1 :]
                if sequence in blocked_sequences:
                    continue
                variants[sequence] = {
                    "parent_sequence": parent_sequence,
                    "parent_v6_rank": parent.get("ace_rank", parent.get("rank", "")),
                    "parent_plm_mean_log_likelihood": score_cache[parent_sequence]["plm_mean_log_likelihood"],
                    "edit_count": 1,
                    "mutations": f"{parent_sequence[position]}{position + 1}{residue}",
                    "masked_lm_support_mean": substitution["log_probability"],
                }

        for left_position, right_position in combinations(range(len(parent_sequence)), 2):
            left_subs = substitutions_by_position[left_position][:double_substitutions_per_position]
            right_subs = substitutions_by_position[right_position][:double_substitutions_per_position]
            for left_sub in left_subs:
                for right_sub in right_subs:
                    sequence_list = list(parent_sequence)
                    sequence_list[left_position] = str(left_sub["residue"])
                    sequence_list[right_position] = str(right_sub["residue"])
                    sequence = "".join(sequence_list)
                    if sequence in blocked_sequences:
                        continue
                    support = mean([float(left_sub["log_probability"]), float(right_sub["log_probability"])])
                    existing = variants.get(sequence)
                    if existing and safe_int(existing.get("edit_count")) < 2:
                        continue
                    variants[sequence] = {
                        "parent_sequence": parent_sequence,
                        "parent_v6_rank": parent.get("ace_rank", parent.get("rank", "")),
                        "parent_plm_mean_log_likelihood": score_cache[parent_sequence][
                            "plm_mean_log_likelihood"
                        ],
                        "edit_count": 2,
                        "mutations": (
                            f"{parent_sequence[left_position]}{left_position + 1}{left_sub['residue']};"
                            f"{parent_sequence[right_position]}{right_position + 1}{right_sub['residue']}"
                        ),
                        "masked_lm_support_mean": support,
                    }

    rows = []
    for sequence, base in variants.items():
        if sequence not in score_cache:
            score_cache[sequence] = scorer.score_sequence(sequence)
        if sequence not in descriptor_cache:
            descriptor_cache[sequence] = molecular_descriptors(sequence)
        nearest, distance = nearest_sequence(sequence, known_positive_sequences)
        descriptors = descriptor_cache[sequence]
        row = {
            "variant_id": stable_id("ACEVAR", f"{base['parent_sequence']}:{sequence}:{base['mutations']}"),
            **base,
            "sequence": sequence,
            "length": len(sequence),
            **score_cache[sequence],
            "delta_parent_plm_mean": (
                score_cache[sequence]["plm_mean_log_likelihood"]
                - safe_float(base.get("parent_plm_mean_log_likelihood"))
            ),
            "nearest_known_positive": nearest,
            "nearest_known_positive_distance": distance,
            **descriptors,
            "developability_flag": developability_flag(descriptors, sequence),
        }
        for key in (
            "plm_mean_log_likelihood",
            "parent_plm_mean_log_likelihood",
            "delta_parent_plm_mean",
            "plm_pseudo_perplexity",
            "masked_lm_support_mean",
        ):
            row[key] = round(float(row[key]), 8)
        rows.append(row)

    rows.sort(
        key=lambda row: (
            safe_int(row["edit_count"]) < 2,
            str(row["developability_flag"]) != "screen_pass",
            -safe_float(row["plm_mean_log_likelihood"]),
            -safe_float(row["delta_parent_plm_mean"]),
            str(row["sequence"]),
        )
    )
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows[:variant_top_n]


def report_markdown(
    *,
    args: argparse.Namespace,
    label_counts: Mapping[str, int],
    scored_rows: Sequence[Mapping[str, Any]],
    fusion_rows: Sequence[Mapping[str, Any]],
    ranked_rows: Sequence[Mapping[str, Any]],
    variant_rows: Sequence[Mapping[str, Any]],
    v6_summary: Mapping[str, str],
    correlation_summary: Mapping[str, Any],
    elapsed_seconds: float,
) -> str:
    best = max(fusion_rows, key=lambda row: safe_float(row["auroc_vs_matched_decoys"], -1.0))
    ace_only = next(row for row in fusion_rows if row["fusion_id"] == "ace_only")
    plm_only = next(row for row in fusion_rows if row["fusion_id"] == "plm_only")
    improvement = safe_float(best["auroc_vs_matched_decoys"]) - safe_float(ace_only["auroc_vs_matched_decoys"])
    v6_reference_auroc = safe_float(v6_summary.get("auroc_vs_matched_decoys"), math.nan)
    v6_reference_top1 = safe_int(v6_summary.get("top_1pct_count"))
    if improvement >= 0.02 and safe_int(best["top_1pct_count"]) >= safe_int(ace_only["top_1pct_count"]):
        verdict = "strengthens"
        interpretation = (
            "The small ESM2 pilot adds independent sequence-language signal on the matched panel. "
            "This is worth scaling to a larger ESM2 checkpoint and a preregistered wider universe."
        )
    elif safe_float(plm_only["auroc_vs_matched_decoys"]) >= 0.70:
        verdict = "supports"
        interpretation = (
            "ProteinLM signal is meaningful but not yet dominant. Use it as an orthogonal reranker and scale the "
            "model/checkpoint before making a discovery claim."
        )
    else:
        verdict = "does_not_strengthen"
        interpretation = (
            "This ESM2 checkpoint does not materially improve the frozen v6 matched-panel result. "
            "This is a concrete negative scenario: scaling or receptor-conditioned structure is required."
        )

    def fmt(value: Any, digits: int = 4) -> str:
        number = safe_float(value, math.nan)
        if math.isnan(number):
            return ""
        return f"{number:.{digits}f}"

    top_variant_lines = []
    for row in variant_rows[:12]:
        top_variant_lines.append(
            "| "
            + " | ".join(
                [
                    str(row["rank"]),
                    str(row["sequence"]),
                    str(row["parent_sequence"]),
                    str(row["edit_count"]),
                    str(row["mutations"]),
                    fmt(row["plm_mean_log_likelihood"], 4),
                    fmt(row["delta_parent_plm_mean"], 4),
                    str(row["developability_flag"]),
                ]
            )
            + " |"
        )
    top_variants = "\n".join(
        [
            "| Rank | Variant | Parent | Edits | Mutations | PLM mean LL | Delta parent | Developability |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
            *top_variant_lines,
        ]
    )

    top_discovery_lines = []
    for row in ranked_rows[:12]:
        top_discovery_lines.append(
            "| "
            + " | ".join(
                [
                    str(row["rank"]),
                    str(row["sequence"]),
                    str(row["label"]),
                    str(row["scenario_group"]),
                    str(row["ace_rank"]),
                    fmt(row["plm_mean_log_likelihood"], 4),
                    fmt(row["fused_score"], 4),
                ]
            )
            + " |"
        )
    top_discovery = "\n".join(
        [
            "| Rank | Sequence | Label | Group | v6 rank | PLM mean LL | Fused score |",
            "| --- | --- | --- | --- | --- | --- | --- |",
            *top_discovery_lines,
        ]
    )

    return f"""# ACE ProteinLM Fusion Scenario

Date: {date.today().isoformat()}

## Scenario

- Model: `{args.model_id}`
- Device: `{correlation_summary.get("device", "")}`
- Scored rows: `{len(scored_rows)}`
- Label counts in scored panel: `{dict(sorted(label_counts.items()))}`
- Runtime seconds: `{elapsed_seconds:.1f}`
- Reference v6 matched-decoy AUROC from frozen full run: `{fmt(v6_reference_auroc, 4)}`
- Reference v6 full-universe top 1% positives: `{v6_reference_top1}`

## Matched-Panel Result

| Scenario | ACE weight | PLM weight | Top 1% positives | AUROC matched | AP matched | Median positive rank |
| --- | --- | --- | --- | --- | --- | --- |
| Frozen v6 scores on scored panel | {ace_only["ace_weight"]} | {ace_only["plm_weight"]} | {ace_only["top_1pct_count"]} | {fmt(ace_only["auroc_vs_matched_decoys"], 4)} | {fmt(ace_only["average_precision_vs_matched_decoys"], 4)} | {fmt(ace_only["median_positive_rank"], 1)} |
| ESM2 only | {plm_only["ace_weight"]} | {plm_only["plm_weight"]} | {plm_only["top_1pct_count"]} | {fmt(plm_only["auroc_vs_matched_decoys"], 4)} | {fmt(plm_only["average_precision_vs_matched_decoys"], 4)} | {fmt(plm_only["median_positive_rank"], 1)} |
| Best fusion `{best["fusion_id"]}` | {best["ace_weight"]} | {best["plm_weight"]} | {best["top_1pct_count"]} | {fmt(best["auroc_vs_matched_decoys"], 4)} | {fmt(best["average_precision_vs_matched_decoys"], 4)} | {fmt(best["median_positive_rank"], 1)} |

Verdict: `{verdict}`.

{interpretation}

Correlation between frozen ACE score and ESM2 mean pseudo-log-likelihood on the matched panel:
Pearson `{fmt(correlation_summary.get("pearson_ace_plm"), 4)}`, Spearman `{fmt(correlation_summary.get("spearman_ace_plm"), 4)}`.

## Best-Fusion Top Rows

{top_discovery}

## Masked-Variant Candidates

The variant scenario uses ESM2 masked substitutions from the best v6 positives, generates one-edit and two-edit candidates, excludes exact sequences already present in the ACE universe or known-positive panels, then ranks by ProteinLM plausibility with RDKit descriptor flags.

{top_variants}

## Files

- Scored panel: `{Path(args.out_dir) / "ace_plm_scored_panel.csv"}`
- Fusion summary: `{Path(args.out_dir) / "ace_plm_fusion_summary.csv"}`
- Best-fusion ranked subset: `{Path(args.out_dir) / "ace_plm_best_fusion_ranked_subset.csv"}`
- Masked variants: `{Path(args.out_dir) / "ace_plm_masked_variant_candidates.csv"}`
- Manifest: `{Path(args.out_dir) / "ace_plm_scenario_manifest.json"}`
"""


def main() -> int:
    args = parse_args()
    started = time.monotonic()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ranked_input = load_ranked_rows(Path(args.ranked))
    v6_summary_rows = read_csv_rows(Path(args.summary))
    v6_summary = v6_summary_rows[0] if v6_summary_rows else {}
    panel = select_scoring_panel(
        ranked_input,
        top_ace_background=args.top_ace_background,
        random_background=args.random_background,
        max_panel_rows=args.max_panel_rows,
        seed=args.seed,
    )
    print(f"[ace-plm] selected {len(panel)} rows for PLM scoring", flush=True)

    model_started = time.monotonic()
    scorer = EsmMaskedLmScorer(args.model_id)
    model_load_seconds = time.monotonic() - model_started
    print(
        f"[ace-plm] loaded {args.model_id} on {scorer.device} in {model_load_seconds:.1f}s",
        flush=True,
    )

    scored_rows = score_rows(panel, scorer, progress_every=args.progress_every)
    write_csv(out_dir / "ace_plm_scored_panel.csv", scored_rows, SCORED_FIELDS)

    eval_rows = [
        row
        for row in scored_rows
        if row.get("label") in POSITIVE_LABELS or row.get("label") == EVAL_NEGATIVE_LABEL
    ]
    correlation_summary = {
        "device": scorer.device,
        "pearson_ace_plm": pearson(
            [safe_float(row.get("ace_score")) for row in eval_rows],
            [safe_float(row.get("plm_mean_log_likelihood")) for row in eval_rows],
        ),
        "spearman_ace_plm": spearman(
            [safe_float(row.get("ace_score")) for row in eval_rows],
            [safe_float(row.get("plm_mean_log_likelihood")) for row in eval_rows],
        ),
    }

    weights = [round(index / 20.0, 2) for index in range(21)]
    fusion_rows = [fusion_metrics(scored_rows, ace_weight=weight) for weight in weights]
    write_csv(out_dir / "ace_plm_fusion_summary.csv", fusion_rows, FUSION_FIELDS)

    best = max(fusion_rows, key=lambda row: safe_float(row["auroc_vs_matched_decoys"], -1.0))
    ranked_subset = best_fusion_rows(scored_rows, safe_float(best["ace_weight"]))
    write_csv(out_dir / "ace_plm_best_fusion_ranked_subset.csv", ranked_subset, RANKED_FIELDS)

    positive_parent_sequences = {row.get("sequence") for row in read_csv_rows(Path(args.positive_ranks))}
    parent_rows = [
        row
        for row in ranked_input
        if row.get("sequence") in positive_parent_sequences and row.get("label") in POSITIVE_LABELS
    ]
    parent_rows.sort(key=lambda row: safe_int(row.get("ace_rank")))
    parent_rows = parent_rows[: args.shell_parents]
    known_positive_sequences = load_known_positive_sequences()
    blocked_sequences = {str(row["sequence"]) for row in ranked_input}
    blocked_sequences.update(known_positive_sequences)
    variant_rows = generate_masked_variants(
        scorer=scorer,
        parent_rows=parent_rows,
        blocked_sequences=blocked_sequences,
        known_positive_sequences=known_positive_sequences,
        single_substitutions_per_position=args.single_substitutions_per_position,
        double_substitutions_per_position=args.double_substitutions_per_position,
        variant_top_n=args.variant_top_n,
    )
    write_csv(out_dir / "ace_plm_masked_variant_candidates.csv", variant_rows, VARIANT_FIELDS)

    label_counts = Counter(str(row.get("label", "")) for row in scored_rows)
    elapsed_seconds = time.monotonic() - started
    manifest = {
        "date": date.today().isoformat(),
        "model_id": args.model_id,
        "device": scorer.device,
        "model_load_seconds": round(model_load_seconds, 3),
        "elapsed_seconds": round(elapsed_seconds, 3),
        "ranked_input": str(Path(args.ranked)),
        "v6_summary": str(Path(args.summary)),
        "positive_ranks": str(Path(args.positive_ranks)),
        "selected_rows": len(panel),
        "label_counts": dict(sorted(label_counts.items())),
        "top_ace_background": args.top_ace_background,
        "random_background": args.random_background,
        "shell_parents": args.shell_parents,
        "single_substitutions_per_position": args.single_substitutions_per_position,
        "double_substitutions_per_position": args.double_substitutions_per_position,
        "variant_top_n": args.variant_top_n,
        "best_fusion": best,
        "correlation_summary": correlation_summary,
        "python": sys.version,
        "torch_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
    }
    write_json(out_dir / "ace_plm_scenario_manifest.json", manifest)
    (out_dir / "ace_plm_fusion_report.md").write_text(
        report_markdown(
            args=args,
            label_counts=label_counts,
            scored_rows=scored_rows,
            fusion_rows=fusion_rows,
            ranked_rows=ranked_subset,
            variant_rows=variant_rows,
            v6_summary=v6_summary,
            correlation_summary=correlation_summary,
            elapsed_seconds=elapsed_seconds,
        ),
        encoding="utf-8",
    )
    print(out_dir / "ace_plm_fusion_report.md", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

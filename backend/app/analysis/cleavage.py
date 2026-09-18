from __future__ import annotations

import math
from typing import Dict, Iterable, List, Tuple

# Logistic preference matrices distilled from the PROSPERous cleavage models.
# We encode position-specific log-odds contributions (P4 -> P4') for the
# supported proteases. Positions are indexed relative to the scissile bond such
# that 0 represents the residue immediately N-terminal of the cleavage site (P1)
# and +1 represents the residue immediately C-terminal (P1'). Negative offsets
# move further toward the N-terminus (P2 = -1, P3 = -2, …).
_PROTEASE_MODELS: Dict[str, Dict[str, object]] = {
    "trypsin": {
        "intercept": -6.5,
        "positions": {
            -3: {"K": 0.3, "R": 0.2, "P": -0.4, "*": -0.1, "X": -0.2},
            -2: {"K": 0.4, "R": 0.3, "P": -0.6, "D": -0.4, "E": -0.4, "*": -0.1, "X": -0.3},
            -1: {"K": 1.0, "R": 0.9, "H": 0.6, "P": -1.2, "D": -0.6, "E": -0.6, "*": -0.2, "X": -0.5},
            0: {
                "K": 6.0,
                "R": 5.8,
                "H": 1.5,
                "P": -3.5,
                "D": -2.5,
                "E": -2.5,
                "*": -1.0,
                "X": -1.4,
            },
            1: {
                "P": -5.5,
                "D": -1.0,
                "E": -1.0,
                "K": 0.6,
                "R": 0.6,
                "G": 0.3,
                "S": 0.3,
                "_": -1.2,
                "*": -0.3,
                "X": -0.8,
            },
            2: {"P": -0.8, "G": 0.5, "S": 0.4, "T": 0.3, "K": 0.2, "R": 0.2, "*": 0.0, "X": -0.2},
            3: {"P": -0.6, "G": 0.2, "S": 0.2, "T": 0.2, "*": 0.0, "X": -0.1},
        },
        "min_probability": 0.01,
    },
    "chymotrypsin": {
        "intercept": -5.5,
        "positions": {
            -3: {"F": 0.3, "W": 0.3, "Y": 0.3, "L": 0.2, "I": 0.2, "M": 0.2, "*": 0.0, "X": -0.2},
            -2: {
                "F": 0.6,
                "W": 0.6,
                "Y": 0.6,
                "L": 0.5,
                "I": 0.4,
                "M": 0.4,
                "P": -0.6,
                "*": -0.1,
                "X": -0.4,
            },
            -1: {
                "F": 1.4,
                "W": 1.3,
                "Y": 1.2,
                "L": 1.1,
                "I": 1.0,
                "M": 1.1,
                "V": 0.8,
                "A": 0.7,
                "P": -1.2,
                "D": -0.8,
                "E": -0.8,
                "*": -0.2,
                "X": -0.6,
            },
            0: {
                "F": 5.2,
                "W": 5.4,
                "Y": 5.0,
                "L": 3.0,
                "I": 2.7,
                "M": 2.8,
                "H": 1.7,
                "V": 1.5,
                "A": 1.2,
                "P": -3.0,
                "D": -2.2,
                "E": -2.2,
                "*": -1.2,
                "X": -1.6,
            },
            1: {
                "P": -4.2,
                "G": 0.7,
                "A": 0.7,
                "S": 0.7,
                "T": 0.6,
                "V": 0.5,
                "Y": -0.5,
                "_": -1.0,
                "*": -0.2,
                "X": -0.5,
            },
            2: {"P": -1.4, "G": 0.3, "A": 0.3, "S": 0.3, "T": 0.2, "*": 0.0, "X": -0.2},
            3: {"P": -0.8, "G": 0.1, "A": 0.1, "S": 0.1, "*": 0.0},
        },
        "min_probability": 0.01,
    },
    "pepsin": {
        "intercept": -4.8,
        "positions": {
            -2: {
                "F": 0.6,
                "W": 0.6,
                "Y": 0.6,
                "L": 0.5,
                "I": 0.5,
                "M": 0.5,
                "V": 0.4,
                "P": -0.7,
                "*": -0.2,
                "X": -0.4,
            },
            -1: {
                "F": 1.0,
                "W": 1.0,
                "Y": 1.0,
                "L": 0.8,
                "I": 0.7,
                "M": 0.7,
                "V": 0.6,
                "A": 0.4,
                "P": -1.0,
                "*": -0.2,
                "X": -0.5,
            },
            0: {
                "F": 4.2,
                "W": 4.3,
                "Y": 4.0,
                "L": 3.5,
                "I": 3.3,
                "M": 3.2,
                "V": 3.0,
                "A": 2.0,
                "P": -2.5,
                "D": -2.0,
                "E": -2.0,
                "*": -1.2,
                "X": -1.6,
            },
            1: {
                "F": 3.0,
                "W": 3.2,
                "Y": 3.0,
                "L": 2.8,
                "I": 2.5,
                "M": 2.4,
                "V": 2.2,
                "A": 1.5,
                "P": -2.5,
                "D": -2.0,
                "E": -2.0,
                "_": -1.0,
                "*": -0.5,
                "X": -0.8,
            },
            2: {"F": 0.6, "W": 0.6, "Y": 0.6, "L": 0.5, "I": 0.5, "M": 0.5, "V": 0.4, "P": -0.8, "*": 0.0},
        },
        "min_probability": 0.01,
    },
    "elastase": {
        "intercept": -6.0,
        "positions": {
            -2: {"A": 0.6, "G": 0.6, "S": 0.5, "V": 0.5, "T": 0.4, "P": -0.8, "*": -0.1, "X": -0.3},
            -1: {
                "A": 1.4,
                "G": 1.3,
                "S": 1.1,
                "V": 1.0,
                "T": 0.9,
                "L": 0.8,
                "I": 0.8,
                "M": 0.7,
                "P": -1.2,
                "Y": -0.8,
                "W": -0.8,
                "*": -0.2,
                "X": -0.5,
            },
            0: {
                "A": 5.0,
                "G": 4.8,
                "S": 4.5,
                "V": 4.0,
                "T": 3.8,
                "L": 3.6,
                "I": 3.5,
                "M": 3.2,
                "P": -2.8,
                "Y": -1.5,
                "W": -1.5,
                "F": -1.2,
                "*": -1.0,
                "X": -1.4,
            },
            1: {"G": 1.2, "A": 1.0, "S": 0.9, "V": 0.7, "T": 0.7, "P": -3.5, "_": -1.0, "*": -0.2, "X": -0.5},
            2: {"G": 0.6, "A": 0.5, "S": 0.5, "V": 0.4, "T": 0.4, "P": -1.0, "*": 0.0},
        },
        "min_probability": 0.01,
    },
    "dpp4": {
        "intercept": -9.0,
        "positions": {
            -1: {
                "P": 6.2,
                "A": 5.0,
                "S": 4.0,
                "G": 3.6,
                "V": 2.4,
                "L": 2.0,
                "I": 2.0,
                "T": 2.0,
                "D": -1.8,
                "E": -1.8,
                "*": -0.8,
                "X": -1.2,
            },
            0: {
                "P": 6.8,
                "A": 5.8,
                "S": 5.0,
                "G": 4.8,
                "V": 3.0,
                "L": 2.8,
                "I": 2.8,
                "T": 2.6,
                "D": -2.5,
                "E": -2.5,
                "*": -1.5,
                "X": -1.8,
            },
            1: {
                "P": -2.5,
                "A": 1.2,
                "S": 1.2,
                "G": 1.2,
                "V": 0.8,
                "T": 0.8,
                "_": -1.0,
                "*": -0.2,
                "X": -0.4,
            },
            2: {"P": -1.5, "A": 0.6, "S": 0.5, "G": 0.5, "T": 0.4, "*": 0.0},
        },
        "min_probability": 0.005,
        "min_index": 1,
        "max_index": 6,
        "distance_penalty": -1.2,
    },
}

_DEFAULT_PROTEASES = ("trypsin", "chymotrypsin", "elastase", "dpp4", "pepsin")


def _sigmoid(x: float) -> float:
    if x >= 60:
        return 1.0
    if x <= -60:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def _weight_for(position_weights: Dict[str, float], residue: str) -> float:
    if residue in position_weights:
        return position_weights[residue]
    if residue == "" or residue is None:
        return position_weights.get("_", 0.0)
    return position_weights.get("*", 0.0)


def _score_site(peptide: str, idx: int, model: Dict[str, object]) -> float:
    score = float(model.get("intercept", 0.0))
    for offset, weights in model["positions"].items():
        pos = idx + offset
        residue = peptide[pos] if 0 <= pos < len(peptide) else "_"
        score += _weight_for(weights, residue)
    distance_penalty = model.get("distance_penalty")
    if isinstance(distance_penalty, (int, float)):
        score += float(distance_penalty) * idx
    return _sigmoid(score)


def _scan_sites(peptide: str, protease: str, model: Dict[str, object]) -> Tuple[float, List[Dict[str, object]]]:
    if not peptide:
        return 0.0, []
    length = len(peptide)
    if length < 2:
        return 0.0, []

    min_idx = int(model.get("min_index", 0))
    max_idx = int(model.get("max_index", length - 2))
    max_idx = min(max_idx, length - 2)
    min_idx = max(min_idx, 0)

    threshold = float(model.get("min_probability", 0.0))
    sites: List[Tuple[int, float]] = []

    for idx in range(min_idx, max_idx + 1):
        prob = _score_site(peptide, idx, model)
        if prob >= threshold:
            sites.append((idx, prob))

    sites.sort(key=lambda x: x[0])

    overall = 1.0
    for _, prob in sites:
        overall *= (1.0 - prob)
    overall = 1.0 - overall
    if overall < 0.0:
        overall = 0.0
    elif overall > 1.0:
        overall = 1.0

    site_payload: List[Dict[str, object]] = []
    for idx, prob in sites:
        residue = peptide[idx] if idx < length else ""
        next_residue = peptide[idx + 1] if idx + 1 < length else ""
        site_payload.append(
            {
                "position": idx + 1,  # 1-based index for clarity
                "probability": prob,
                "residue": residue,
                "nextResidue": next_residue,
            }
        )

    return overall, site_payload


def protease_cleavage_profile(
    peptide: str,
    proteases: Iterable[str] | None = None,
) -> Dict[str, Dict[str, object]]:
    clean = "".join(ch for ch in peptide.upper() if ch.isalpha())
    if not clean:
        return {}
    targets = list(proteases) if proteases is not None else list(_DEFAULT_PROTEASES)
    profile: Dict[str, Dict[str, object]] = {}

    for protease in targets:
        if protease == "none":
            profile[protease] = {"probability": 0.0, "sites": []}
            continue
        model = _PROTEASE_MODELS.get(protease)
        if not model:
            continue
        overall, site_payload = _scan_sites(clean, protease, model)
        profile[protease] = {
            "probability": overall,
            "sites": site_payload,
        }

    return profile


def generic_cleavage_risk(peptide: str, enzyme: str) -> float:
    profile = protease_cleavage_profile(peptide, [enzyme])
    data = profile.get(enzyme)
    if not data:
        return 0.0
    probability = data.get("probability")
    try:
        value = float(probability)
    except (TypeError, ValueError):
        return 0.0
    return max(0.0, min(1.0, value))


__all__ = ["protease_cleavage_profile", "generic_cleavage_risk"]

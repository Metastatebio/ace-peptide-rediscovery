from __future__ import annotations

import hashlib
import math
import random
from collections import Counter
from functools import lru_cache
from pathlib import Path
from statistics import mean, pstdev
from typing import Dict, List, Optional, Tuple

try:
    import yaml
except ImportError:  # pragma: no cover - PyYAML is optional but listed in requirements
    yaml = None

from .cleavage import protease_cleavage_profile

_AA = "ACDEFGHIKLMNPQRSTVWY"
_KD = {
    'I': 4.5,
    'V': 4.2,
    'L': 3.8,
    'F': 2.8,
    'C': 2.5,
    'M': 1.9,
    'A': 1.8,
    'G': -0.4,
    'T': -0.7,
    'S': -0.8,
    'W': -0.9,
    'Y': -1.3,
    'P': -1.6,
    'H': -3.2,
    'E': -3.5,
    'Q': -3.5,
    'D': -3.5,
    'N': -3.5,
    'K': -3.9,
    'R': -4.5,
}
_pKa = {'Cterm': 3.1, 'Nterm': 8.0, 'C': 8.5, 'D': 3.9, 'E': 4.1, 'H': 6.5, 'K': 10.5, 'R': 12.5, 'Y': 10.1}
_HYDROPHOBIC_AA = set("AVILMFWY")
_BASIC_AA = set("KRH")
_AROMATIC_AA = set("FWY")

def _clamp01(value: float) -> float:
    if not isinstance(value, (int, float)):
        return 0.0
    if math.isnan(value):
        return 0.0
    if value < 0.0:
        return 0.0
    if value > 1.0:
        return 1.0
    return float(value)

def aac(peptide: str) -> Dict[str, float]:
    counts = {aa:0 for aa in _AA}
    for ch in peptide:
        if ch in counts: counts[ch]+=1
    n = len(peptide) or 1
    return {aa: counts[aa]/n for aa in _AA}

def hydrophobicity(peptide: str) -> float:
    if not peptide: return 0.0
    return sum(_KD.get(a,0.0) for a in peptide) / len(peptide)

def charge(peptide: str, pH: float=7.4) -> float:
    if not peptide: return 0.0
    nterm = 1.0/(1.0 + 10**(pH-_pKa['Nterm'])); cterm = -1.0/(1.0 + 10**(_pKa['Cterm']-pH))
    pos = neg = 0.0
    for a in peptide:
        if a == 'K': pos += 1.0/(1.0 + 10**(pH-_pKa['K']))
        elif a == 'R': pos += 1.0/(1.0 + 10**(pH-_pKa['R']))
        elif a == 'H': pos += 1.0/(1.0 + 10**(pH-_pKa['H']))
        elif a == 'D': neg += 1.0/(1.0 + 10**(_pKa['D']-pH))
        elif a == 'E': neg += 1.0/(1.0 + 10**(_pKa['E']-pH))
        elif a == 'C': neg += 1.0/(1.0 + 10**(_pKa['C']-pH))
        elif a == 'Y': neg += 1.0/(1.0 + 10**(_pKa['Y']-pH))
    return nterm + pos + cterm - neg

try:
    from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
    _HAVE_MODLAMP = True
except Exception:
    _HAVE_MODLAMP = False

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    _HAVE_PROTPARAM = True
except Exception:
    _HAVE_PROTPARAM = False


def _first_float(value: object) -> Optional[float]:
    """Extract the first numerical value from possibly nested descriptor output."""
    if value is None:
        return None
    try:
        if hasattr(value, "size") and getattr(value, "size") == 1 and hasattr(value, "item"):
            return float(value.item())  # numpy scalar
    except Exception:
        pass
    try:
        if hasattr(value, "flat"):
            iterator = iter(value.flat)
            first = next(iterator)
            return float(first)
    except Exception:
        pass
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        return _first_float(value[0])
    try:
        return float(value)
    except Exception:
        return None

def modlamp_features(peptide: str) -> Dict[str, float]:
    if not _HAVE_MODLAMP or not peptide:
        return {}
    try:
        result: Dict[str, float] = {}
        pd = PeptideDescriptor([peptide], "eisenberg")
        pd.calculate_moment()
        hm = _first_float(pd.descriptor)
        if hm is not None:
            result["modlamp_hydrophobic_moment"] = hm

        gd = GlobalDescriptor([peptide])
        modlamp_methods = [
            ("isoelectric_point", "modlamp_pI"),
            ("charge_density", "modlamp_charge_density"),
            ("aromaticity", "modlamp_aromaticity"),
            ("aliphatic_index", "modlamp_aliphatic_index"),
            ("hydrophobic_ratio", "modlamp_hydrophobic_ratio"),
            ("instability_index", "modlamp_instability_index"),
            ("boman_index", "modlamp_boman_index"),
            ("calculate_MW", "modlamp_molecular_weight"),
        ]
        for method_name, feature_name in modlamp_methods:
            try:
                getattr(gd, method_name)()
                value = _first_float(gd.descriptor)
                if value is not None:
                    result[feature_name] = value
            except Exception:
                continue
        return result
    except Exception:
        return {}


def protparam_features(peptide: str) -> Dict[str, float]:
    if not _HAVE_PROTPARAM or not peptide:
        return {}
    try:
        analysis = ProteinAnalysis(peptide.upper())
    except Exception:
        return {}
    try:
        features: Dict[str, float] = {
            "protparam_molecular_weight": float(analysis.molecular_weight()),
            "protparam_aromaticity": float(analysis.aromaticity()),
            "protparam_instability_index": float(analysis.instability_index()),
            "protparam_isoelectric_point": float(analysis.isoelectric_point()),
            "protparam_gravy": float(analysis.gravy()),
        }
    except Exception:
        return {}
    try:
        helix, turn, sheet = analysis.secondary_structure_fraction()
        features.update({
            "protparam_secstruct_helix": float(helix),
            "protparam_secstruct_turn": float(turn),
            "protparam_secstruct_sheet": float(sheet),
        })
    except Exception:
        pass
    try:
        flexibility = analysis.flexibility()
        if flexibility:
            features["protparam_flexibility_mean"] = float(sum(flexibility) / len(flexibility))
    except Exception:
        pass
    return features

def full_features(peptide: str) -> Dict[str, float]:
    a = aac(peptide)
    d = hydrophobicity(peptide)
    q = charge(peptide)
    arom = sum(a[x] for x in ('F', 'W', 'Y'))
    pro = a['P']
    acidic = a['D'] + a['E']
    base = {
        "len": len(peptide),
        "hydrophobicity": d,
        "charge": q,
        "frac_aromatic": arom,
        "frac_proline": pro,
        "frac_acidic": acidic,
        **{f"aac_{k}": v for k, v in a.items()},
    }
    base.update(modlamp_features(peptide))
    base.update(protparam_features(peptide))
    return base


def _neutral_positive_score(window: str, pH: float = 6.8) -> float:
    value = abs(charge(window, pH=pH))
    return max(0.0, 1.0 - min(value / 2.0, 1.0))


def _hydrophobic_nterminus_score(residue: str) -> float:
    kd = _KD.get(residue, 0.0)
    if kd <= 0.0:
        return 0.0
    if kd >= 2.5:
        return 1.0
    return max(0.0, kd / 2.5)


def _moderate_hydropathy_score(window: str) -> float:
    avg = hydrophobicity(window)
    return max(0.0, 1.0 - min(abs(avg - 0.5) / 3.0, 1.0))


def pept1_uptake_index(peptide: str) -> float:
    clean = ''.join(res for res in peptide.upper() if res.isalpha())
    if len(clean) < 2:
        return 0.0

    scores: List[float] = []
    for window_size in (2, 3):
        if len(clean) < window_size:
            continue
        for idx in range(len(clean) - window_size + 1):
            window = clean[idx: idx + window_size]
            charge_score = _neutral_positive_score(window)
            nterm_score = _hydrophobic_nterminus_score(window[0])
            hydro_score = _moderate_hydropathy_score(window)
            scores.append((charge_score + nterm_score + hydro_score) / 3.0)

    if not scores:
        return 0.0
    return _clamp01(sum(scores) / len(scores))


def chemosensor_engagement(peptide: str) -> float:
    clean = ''.join(res for res in peptide.upper() if res.isalpha())
    length = len(clean)
    if length == 0:
        return 0.0
    basic_fraction = sum(1 for res in clean if res in _BASIC_AA) / length
    aromatic_fraction = sum(1 for res in clean if res in _AROMATIC_AA) / length
    hydrophobic_fraction = sum(1 for res in clean if res in _HYDROPHOBIC_AA) / length

    basic_term = min(basic_fraction / 0.35, 1.0)
    aromatic_term = min(aromatic_fraction / 0.25, 1.0)
    hydrophobic_term = min(hydrophobic_fraction / 0.6, 1.0)
    score = 0.4 * basic_term + 0.3 * aromatic_term + 0.3 * hydrophobic_term
    return _clamp01(score)


def _count_gi_cleavage_sites(profile: Dict[str, Dict[str, object]]) -> Tuple[int, Dict[str, int]]:
    total = 0
    per_protease: Dict[str, int] = {}
    for protease, payload in profile.items():
        sites = payload.get('sites') if isinstance(payload, dict) else None
        if not isinstance(sites, list):
            continue
        count = 0
        for site in sites:
            if not isinstance(site, dict):
                continue
            probability = site.get('probability')
            next_residue = site.get('nextResidue')
            try:
                prob_value = float(probability)
            except (TypeError, ValueError):
                continue
            if prob_value < 0.2:
                continue
            if isinstance(next_residue, str) and next_residue.upper() == 'P':
                continue
            count += 1
        if count:
            per_protease[protease] = count
            total += count
    return total, per_protease


def gi_protease_stability_score(peptide: str, profile: Dict[str, Dict[str, object]] | None = None) -> Tuple[float, Dict[str, int]]:
    clean = ''.join(res for res in peptide.upper() if res.isalpha())
    length = len(clean)
    if length <= 1:
        return 1.0, {}

    if profile is None:
        profile = protease_cleavage_profile(clean, ['trypsin', 'chymotrypsin', 'pepsin'])
    else:
        profile = {key: value for key, value in profile.items() if key in {'trypsin', 'chymotrypsin', 'pepsin'}}

    total_sites, per_protease = _count_gi_cleavage_sites(profile)
    max_sites = max(1, length - 1)
    score = 1.0 - min(total_sites / max_sites, 1.0)
    return _clamp01(score), per_protease


def hydrophobic_moment_alpha(peptide: str) -> Tuple[float, float]:
    clean = ''.join(res for res in peptide.upper() if res.isalpha())
    length = len(clean)
    if length == 0:
        return 0.0, 0.0

    angle_step = math.radians(100.0)
    cos_sum = 0.0
    sin_sum = 0.0
    for index, residue in enumerate(clean):
        value = _KD.get(residue, 0.0)
        angle = angle_step * index
        cos_sum += value * math.cos(angle)
        sin_sum += value * math.sin(angle)

    moment = math.sqrt(cos_sum ** 2 + sin_sum ** 2) / length
    normalized = _clamp01(moment / 4.5)
    return normalized, moment


def solubility_exposure_score(peptide: str) -> float:
    gravy = hydrophobicity(peptide)
    score = (4.5 - gravy) / 9.0
    return _clamp01(score)


def _residue_charge(residue: str, pH: float = 6.8) -> float:
    if residue == 'K':
        return 1.0 / (1.0 + 10 ** (pH - _pKa['K']))
    if residue == 'R':
        return 1.0 / (1.0 + 10 ** (pH - _pKa['R']))
    if residue == 'H':
        return 1.0 / (1.0 + 10 ** (pH - _pKa['H']))
    if residue == 'D':
        return -1.0 / (1.0 + 10 ** (_pKa['D'] - pH))
    if residue == 'E':
        return -1.0 / (1.0 + 10 ** (_pKa['E'] - pH))
    if residue == 'C':
        return -1.0 / (1.0 + 10 ** (_pKa['C'] - pH))
    if residue == 'Y':
        return -1.0 / (1.0 + 10 ** (_pKa['Y'] - pH))
    return 0.0


def sequence_property_trace(peptide: str, window: int = 5) -> List[Dict[str, float]]:
    clean = ''.join(res for res in peptide.upper() if res.isalpha())
    if not clean:
        return []
    length = len(clean)
    hydro_values = [_KD.get(res, 0.0) for res in clean]
    charge_values = [_residue_charge(res) for res in clean]

    def _moving_average(values: List[float]) -> List[float]:
        if length == 1:
            return values[:]
        radius = max(1, window)
        half = radius // 2
        smoothed: List[float] = []
        for index in range(length):
            start = max(0, index - half)
            end = min(length, index + half + 1)
            segment = values[start:end]
            smoothed.append(sum(segment) / len(segment))
        return smoothed

    hydro_smooth = _moving_average(hydro_values)
    charge_smooth = _moving_average(charge_values)

    trace: List[Dict[str, float]] = []
    for index, residue in enumerate(clean):
        trace.append(
            {
                'position': index + 1,
                'hydropathy': hydro_smooth[index],
                'charge': charge_smooth[index],
                'residue': residue,
            }
        )
    return trace


_GPCR_DEFAULTS: Dict[str, Dict[str, object]] = {
    'logistic': {
        'aromatic': {'center': 0.22, 'scale': 0.06},
        'basic': {'center': 0.18, 'scale': 0.06},
        'solubility': {'center': 0.0, 'scale': 0.2},
        'toxicity': {'center': 0.15, 'scale': 0.05},
        'bitterness': {'center': 0.20, 'scale': 0.07},
    },
    'weights': {
        'aromatic': 0.12,
        'basic': 0.10,
        'turn': 0.10,
        'length': 0.06,
        'charge': 0.04,
        'solubility': 0.04,
        'toxicity': -0.03,
        'bitterness': 0.03,
    },
    'turn': {
        'proline_bonus': 0.6,
        'glycine_adjacent_bonus': 0.2,
        'aromatic_i_i3_bonus': 0.2,
        'cap': 1.0,
    },
    'basic_fraction': {'histidine_cap': 0.5},
}


def _gpcr_config_path() -> Path:
    return Path(__file__).with_name('weights_gpcr.yaml')


def _merge_dict(base: Dict[str, object], update: Dict[str, object]) -> Dict[str, object]:
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            base[key] = _merge_dict(dict(base[key]), value)  # type: ignore[arg-type]
        else:
            base[key] = value
    return base


@lru_cache(maxsize=1)
def _load_gpcr_config() -> Dict[str, Dict[str, object]]:
    config: Dict[str, Dict[str, object]] = {
        key: dict(value) for key, value in _GPCR_DEFAULTS.items()
    }

    path = _gpcr_config_path()
    if path.exists() and yaml is not None:
        try:
            with path.open('r', encoding='utf-8') as handle:
                loaded = yaml.safe_load(handle)
            if isinstance(loaded, dict):
                config = _merge_dict(config, loaded)
        except Exception:
            pass
    return config


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def _apply_logistic(value: float, key: str, config: Dict[str, Dict[str, object]]) -> float:
    defaults = _GPCR_DEFAULTS['logistic'][key]
    params = config.get('logistic', {}).get(key, {})
    center = float(params.get('center', defaults['center']))
    scale = float(params.get('scale', defaults['scale'])) or 1.0
    adjusted = (value - center) / scale
    return _sigmoid(adjusted)


def _fraction(count: float, length: int) -> float:
    if length <= 0:
        return 0.0
    return count / float(length)


def _gpcr_feature_values(peptide: str) -> Dict[str, float]:
    clean = ''.join(res for res in peptide.upper() if res.isalpha())
    length = len(clean)
    counts = Counter(clean)

    aromatic_set = {'Y', 'F', 'W'}
    polar_set = {'S', 'T', 'N', 'Q', 'Y'}
    hydrophobic_set = {'A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y'}
    strong_hydrophobic = {'F', 'W'}

    histidine_cap = float(_load_gpcr_config()['basic_fraction'].get('histidine_cap', 0.5))
    aromatic_fraction = _fraction(sum(counts[res] for res in aromatic_set), length)
    histidine_effective = min(counts['H'], histidine_cap * length)
    basic_fraction = _fraction(counts['R'] + counts['K'] + histidine_effective, length)

    turn_config = _load_gpcr_config()['turn']
    turn_score = 0.0
    if counts['P']:
        turn_score += float(turn_config.get('proline_bonus', 0.6))

    glycine_bonus = float(turn_config.get('glycine_adjacent_bonus', 0.2))
    if glycine_bonus > 0 and length:
        for index, residue in enumerate(clean):
            if residue != 'P':
                continue
            if (index > 0 and clean[index - 1] == 'G') or (index < length - 1 and clean[index + 1] == 'G'):
                turn_score += glycine_bonus
                break

    aromatic_bonus = float(turn_config.get('aromatic_i_i3_bonus', 0.2))
    if aromatic_bonus > 0 and length:
        for index in range(length - 3):
            if clean[index] in aromatic_set and clean[index + 3] in aromatic_set:
                turn_score += aromatic_bonus
                break

    turn_score = min(turn_score, float(turn_config.get('cap', 1.0)))

    length_score = max(0.0, 1.0 - abs(length - 8.0) / 4.0) if length else 0.0

    net_charge = counts['R'] + counts['K'] + 0.1 * counts['H'] - counts['D'] - counts['E']
    charge_score = math.exp(-abs(net_charge) / 3.0)

    polar_fraction = _fraction(sum(counts[res] for res in polar_set), length)
    hydrophobic_fraction = _fraction(sum(counts[res] for res in hydrophobic_set), length)
    solv_input = polar_fraction - hydrophobic_fraction

    strong_run = 0.0
    streak = 0
    strong_streak = 0
    for residue in clean:
        if residue in hydrophobic_set:
            streak += 1
        else:
            if streak >= 3:
                strong_run += streak - 2
            streak = 0

        if residue in strong_hydrophobic:
            strong_streak += 1
        else:
            if strong_streak >= 2:
                strong_run += strong_streak - 1
            strong_streak = 0

    if streak >= 3:
        strong_run += streak - 2
    if strong_streak >= 2:
        strong_run += strong_streak - 1

    bitterness_raw = _fraction(strong_run, length)

    basic_cluster_score = 0.0
    cluster = 0
    for residue in clean:
        if residue in {'K', 'R'}:
            cluster += 1
        else:
            if cluster >= 3:
                basic_cluster_score += cluster - 2
            cluster = 0
    if cluster >= 3:
        basic_cluster_score += cluster - 2

    rare_motifs = 0
    for index in range(length - 1):
        if clean[index:index + 2] in {'CW', 'HW'}:
            rare_motifs += 1

    toxicity_raw = _fraction(basic_cluster_score + 0.5 * rare_motifs, length)

    return {
        'clean_sequence': clean,
        'length': float(length),
        'aromatic_fraction': aromatic_fraction,
        'basic_fraction': basic_fraction,
        'turn_score': turn_score,
        'length_score': length_score,
        'net_charge': net_charge,
        'charge_score': charge_score,
        'solv_input': solv_input,
        'polar_fraction': polar_fraction,
        'hydrophobic_fraction': hydrophobic_fraction,
        'toxicity_raw': toxicity_raw,
        'bitterness_raw': bitterness_raw,
    }


def _gpcr_transforms(values: Dict[str, float]) -> Dict[str, float]:
    config = _load_gpcr_config()

    s_aromatic = _apply_logistic(values['aromatic_fraction'], 'aromatic', config)
    s_basic = _apply_logistic(values['basic_fraction'], 'basic', config)
    s_turn = _clamp01(values['turn_score'])
    s_length = _clamp01(values['length_score'])
    s_charge = _clamp01(values['charge_score'])
    s_solv = _apply_logistic(values['solv_input'], 'solubility', config)
    s_tox = _apply_logistic(values['toxicity_raw'], 'toxicity', config)
    s_bitt = _apply_logistic(values['bitterness_raw'], 'bitterness', config)

    return {
        'aromatic': s_aromatic,
        'basic': s_basic,
        'turn': s_turn,
        'length': s_length,
        'charge': s_charge,
        'solubility': s_solv,
        'toxicity': s_tox,
        'bitterness': s_bitt,
    }


def _gpcr_score(transforms: Dict[str, float]) -> float:
    config = _load_gpcr_config()
    weights = config.get('weights', {})
    total_weight = 0.0
    weighted_sum = 0.0

    for key, value in transforms.items():
        weight = float(weights.get(key, _GPCR_DEFAULTS['weights'].get(key, 0.0)))
        magnitude = abs(weight)
        if magnitude == 0:
            continue

        contribution = value if weight >= 0 else 1.0 - value
        weighted_sum += magnitude * contribution
        total_weight += magnitude

    if total_weight == 0:
        return 0.0

    normalized = weighted_sum / total_weight
    return _clamp01(normalized)


def compute_glico_metrics(peptide: str, profile: Dict[str, Dict[str, object]] | None = None) -> Dict[str, object]:
    del profile  # unused in GPCR mode

    values = _gpcr_feature_values(peptide)
    transforms = _gpcr_transforms(values)
    glico_score = _gpcr_score(transforms)

    clean = values['clean_sequence']
    scramble_scores: List[float] = []
    if len(clean) > 1:
        seed_bytes = hashlib.sha256(clean.encode('utf-8')).digest()
        seed = int.from_bytes(seed_bytes[:8], 'big')
        rng = random.Random(seed)
        chars = list(clean)
        for _ in range(200):
            rng.shuffle(chars)
            scrambled = ''.join(chars)
            scrambled_values = _gpcr_feature_values(scrambled)
            scramble_transforms = _gpcr_transforms(scrambled_values)
            scramble_scores.append(_gpcr_score(scramble_transforms))

    z_score = None
    if scramble_scores:
        m = mean(scramble_scores)
        try:
            sd = pstdev(scramble_scores)
        except Exception:
            sd = 0.0
        if sd and sd > 1e-6:
            z_score = (glico_score - m) / sd
        elif sd == 0:
            z_score = 0.0

    return {
        'aromatic': transforms['aromatic'],
        'basic': transforms['basic'],
        'turn': transforms['turn'],
        'length': transforms['length'],
        'charge': transforms['charge'],
        'solubility': transforms['solubility'],
        'toxicity': transforms['toxicity'],
        'bitterness': transforms['bitterness'],
        'score': glico_score,
        'zScore': z_score,
        'scrambleScores': scramble_scores,
    }

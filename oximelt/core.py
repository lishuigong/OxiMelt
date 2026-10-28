from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from itertools import combinations, product
import csv
import json
import math
import re
import sysconfig
from statistics import median
from typing import Dict, List, Tuple, Any

PACKAGE_ROOT = Path(__file__).resolve().parent.parent
_SOURCE_DATA_DIR = PACKAGE_ROOT / "data"
_INSTALLED_DATA_DIR = Path(sysconfig.get_path("data")) / "data"
DATA_DIR = _SOURCE_DATA_DIR if _SOURCE_DATA_DIR.exists() else _INSTALLED_DATA_DIR

SUBSCRIPT_TRANS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
TOKEN_RE = re.compile(r"[A-Z][a-z]?|\(|\)|(?:\d+(?:\.\d*)?|\.\d+)")

class OxiMeltError(ValueError):
    pass

@dataclass
class PairContribution:
    pair: str
    x_i: float
    x_j: float
    E_strain: float
    E_chem: float
    Omega_K: float
    weighted_contribution_K: float
    seen_in_two_cation_training: bool

@dataclass
class LocalResidualSupport:
    dopant: str
    dopant_share: float
    supported: bool
    mode: str
    residual_K: float
    distance: float | None
    support_n: int
    support_formula: str
    support_x_Hf: float | None = None
    support_x_Hf_low: float | None = None
    support_x_Hf_high: float | None = None
    target_vacancy_fraction: float | None = None
    support_vacancy_fraction: float | None = None
    vacancy_scaling_factor: float | None = None
    unscaled_residual_K: float | None = None
    vacancy_fraction_distance: float | None = None
    concentration_extrapolation: bool = False

@dataclass
class PredictionResult:
    formula: str
    # `predicted_Tm_*` is the selected estimate, not a claim of superior accuracy.
    predicted_Tm_K: float
    predicted_Tm_C: float
    global_Tm_K: float
    global_Tm_C: float
    local_correction_K: float
    corrected_Tm_K: float
    corrected_Tm_C: float
    prediction_mode: str
    hf_local_domain_eligible: bool
    hf_local_correction_applied: bool
    hf_local_support_coverage: float
    T0_K: float
    pair_correction_K: float
    n_cations: int
    cation_fractions: Dict[str, float]
    formal_oxidation_states: Dict[str, float]
    applicability: str
    applicability_detail: str
    warnings: List[str]
    pair_contributions: List[PairContribution]
    local_support_details: List[LocalResidualSupport]
    nominal_oxygen_vacancy_site_fraction: float
    hf_local_concentration_extrapolation: bool

    def to_flat_dict(self) -> Dict[str, Any]:
        return {
            "formula": self.formula,
            "predicted_Tm_K": self.predicted_Tm_K,
            "predicted_Tm_C": self.predicted_Tm_C,
            "global_Tm_K": self.global_Tm_K,
            "global_Tm_C": self.global_Tm_C,
            "local_correction_K": self.local_correction_K,
            "corrected_Tm_K": self.corrected_Tm_K,
            "corrected_Tm_C": self.corrected_Tm_C,
            "prediction_mode": self.prediction_mode,
            "hf_local_domain_eligible": int(self.hf_local_domain_eligible),
            "hf_local_correction_applied": int(self.hf_local_correction_applied),
            "hf_local_support_coverage": self.hf_local_support_coverage,
            "nominal_oxygen_vacancy_site_fraction": self.nominal_oxygen_vacancy_site_fraction,
            "hf_local_concentration_extrapolation": int(self.hf_local_concentration_extrapolation),
            "hf_local_support_details_json": json.dumps(
                [asdict(d) for d in self.local_support_details], ensure_ascii=False,
                allow_nan=False,
            ),
            "hf_local_support_detail": "; ".join(
                f"{d.dopant}:{d.mode}:{d.residual_K:+.3f}K:share={d.dopant_share:.3f}"
                for d in self.local_support_details
            ),
            "T0_K": self.T0_K,
            "pair_correction_K": self.pair_correction_K,
            "n_cations": self.n_cations,
            "formal_oxidation_states": "; ".join(
                f"{k}({v:+g})" for k, v in sorted(self.formal_oxidation_states.items())
            ),
            "applicability": self.applicability,
            "applicability_detail": self.applicability_detail,
            "warnings": " | ".join(self.warnings),
        }

def _load_csv_dict(path: Path, key: str) -> Dict[str, Dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return {row[key]: row for row in csv.DictReader(f)}

def _load_resources():
    props = _load_csv_dict(DATA_DIR / "element_properties.csv", "element")
    for e, row in props.items():
        for k in (
            "ionization_energy_eV", "pauling_en", "ionic_radius_A",
            "cohesive_energy_eV", "bulk_modulus_GPa"
        ):
            row[k] = float(row[k])

    refs: Dict[str, List[Tuple[float, float]]] = {}
    with (DATA_DIR / "simple_oxide_reference_states.csv").open(
        newline="", encoding="utf-8-sig"
    ) as f:
        for row in csv.DictReader(f):
            refs.setdefault(row["element"], []).append(
                (float(row["oxidation_state"]), float(row["Tm_anchor_K"]))
            )

    seen_pairs = set()
    with (DATA_DIR / "training_pair_landscape.csv").open(
        newline="", encoding="utf-8-sig"
    ) as f:
        for row in csv.DictReader(f):
            a, b = row["pair"].replace("–", "-").split("-")
            seen_pairs.add(tuple(sorted((a, b))))

    model = json.loads((DATA_DIR / "model.json").read_text(encoding="utf-8"))

    local_support: Dict[str, List[Dict[str, Any]]] = {}
    local_path = DATA_DIR / "hf_local_residual_support.csv"
    if local_path.exists():
        with local_path.open(newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                item = {
                    "dopant": row["dopant"],
                    "x_Hf": float(row["x_Hf"]),
                    "residual_K": float(row["residual_K"]),
                    "n_records_at_x": int(row["n_records_at_x"]),
                    "support_formula": row["support_formula"],
                    "dopant_formal_valence": float(row["dopant_formal_valence"]),
                    "nominal_vacancy_site_fraction": float(row["nominal_vacancy_site_fraction"]),
                }
                fv = item["nominal_vacancy_site_fraction"]
                expected = (1.0 - item["x_Hf"]) * (4.0 - item["dopant_formal_valence"]) / 4.0
                if (not all(math.isfinite(item[k]) for k in (
                        "x_Hf", "residual_K", "dopant_formal_valence", "nominal_vacancy_site_fraction"))
                        or not 0 < fv < 1 or not 0 < item["x_Hf"] < 1
                        or not math.isclose(fv, expected, rel_tol=1e-8, abs_tol=1e-10)
                        or item["n_records_at_x"] < 1):
                    raise OxiMeltError(f"Invalid vacancy/valence metadata for support: {item['support_formula']}")
                local_support.setdefault(item["dopant"], []).append(item)
        for dopant in local_support:
            local_support[dopant].sort(key=lambda r: r["x_Hf"])

    return props, refs, seen_pairs, model, local_support

PROPS, REFS, SEEN_PAIRS, MODEL, HF_LOCAL_SUPPORT = _load_resources()

def _merge(dst: Dict[str, float], src: Dict[str, float], mult: float = 1.0):
    for e, v in src.items():
        dst[e] = dst.get(e, 0.0) + v * mult

def parse_formula(formula: str) -> Dict[str, float]:
    """Parse ordinary oxide formulas with optional decimal stoichiometry and parentheses."""
    if not formula or not formula.strip():
        raise OxiMeltError("Please enter a chemical formula.")

    s = formula.strip().translate(SUBSCRIPT_TRANS).replace(" ", "")
    if any(x in s for x in ("·", "•", "[", "]", "{", "}", "+", "-")):
        raise OxiMeltError(
            "Hydrates, brackets, explicit charges, and dot-separated formulas are not supported."
        )

    tokens = TOKEN_RE.findall(s)
    if "".join(tokens) != s:
        raise OxiMeltError(f"Unsupported formula syntax: {formula}")

    pos = 0

    def parse_group(stop_at_paren=False):
        nonlocal pos
        comp: Dict[str, float] = {}
        while pos < len(tokens):
            tok = tokens[pos]
            if tok == ")":
                if not stop_at_paren:
                    raise OxiMeltError("Unmatched ')' in formula.")
                pos += 1
                return comp
            if tok == "(":
                pos += 1
                inner = parse_group(True)
                mult = 1.0
                if pos < len(tokens) and re.fullmatch(r"(?:\d+(?:\.\d*)?|\.\d+)", tokens[pos]):
                    mult = float(tokens[pos]); pos += 1
                _merge(comp, inner, mult)
                continue
            if re.fullmatch(r"[A-Z][a-z]?", tok):
                element = tok
                pos += 1
                amount = 1.0
                if pos < len(tokens) and re.fullmatch(r"(?:\d+(?:\.\d*)?|\.\d+)", tokens[pos]):
                    amount = float(tokens[pos]); pos += 1
                if amount <= 0:
                    raise OxiMeltError("Stoichiometric coefficients must be positive.")
                comp[element] = comp.get(element, 0.0) + amount
                continue
            raise OxiMeltError(f"Unexpected token '{tok}'.")
        if stop_at_paren:
            raise OxiMeltError("Unmatched '(' in formula.")
        return comp

    comp = parse_group(False)
    if pos != len(tokens):
        raise OxiMeltError("Could not fully parse formula.")
    return comp

def assign_formal_oxidation_states(comp: Dict[str, float]) -> Tuple[Dict[str, float], Dict[str, float]]:
    if "O" not in comp:
        raise OxiMeltError("The current model is oxide-specific: oxygen must be present.")

    non_oxygen = [e for e in comp if e != "O"]
    unsupported = sorted(e for e in non_oxygen if e not in REFS or e not in PROPS)
    if unsupported:
        raise OxiMeltError(
            "Unsupported cation(s) or missing frozen reference/property data: "
            + ", ".join(unsupported)
        )

    nO = comp["O"]
    choices = [REFS[e] for e in non_oxygen]
    solutions = []
    scale = max(1.0, 2.0 * nO)

    for combo in product(*choices):
        charge = sum(comp[e] * combo[i][0] for i, e in enumerate(non_oxygen))
        if abs(charge - 2.0 * nO) <= 1e-6 * scale:
            states = {e: combo[i][0] for i, e in enumerate(non_oxygen)}
            anchors = {e: combo[i][1] for i, e in enumerate(non_oxygen)}
            solutions.append((states, anchors))

    if not solutions:
        raise OxiMeltError(
            "No unique conventional O(-II) formal-valence assignment matches this formula. "
            "The composition may be mixed-valent, peroxide/superoxide-like, non-stoichiometric "
            "outside the supported rule set, or outside the frozen reference-state domain."
        )

    unique = {}
    for states, anchors in solutions:
        k = tuple(sorted(states.items()))
        unique[k] = (states, anchors)
    solutions = list(unique.values())

    if len(solutions) != 1:
        opts = "; ".join(
            ", ".join(f"{e}{z:+g}" for e, z in sorted(states.items()))
            for states, _ in solutions[:6]
        )
        raise OxiMeltError(
            "Formal oxidation-state assignment is ambiguous under the frozen reference rules. "
            f"Candidate assignments: {opts}"
        )
    return solutions[0]

def nominal_vacancy_fraction(comp: Dict[str, float]) -> float:
    """Stoichiometric deficiency relative to MO2; not a measured defect content."""
    ncat = sum(v for e, v in comp.items() if e != "O")
    return max(0.0, (2.0 - comp["O"] / ncat) / 2.0) if ncat > 0 else 0.0


def _hf_defect_like(comp, states, cation_fractions) -> bool:
    """Composition-only Hf defect-family rule.

    The same vacancy-scaled binary residual rule applies to two-/three-cation
    targets. Other inputs use the global equation without a family correction.
    """
    cfg = MODEL.get("hf_local_transfer", {})
    if not cfg.get("enabled", False):
        return False
    if "Hf" not in cation_fractions or len(cation_fractions) not in (2, 3):
        return False
    ncat = sum(v for e, v in comp.items() if e != "O")
    if ncat <= 0:
        return False
    vacancy_site_fraction_proxy = nominal_vacancy_fraction(comp)
    tol = float(cfg.get("hf_valence_tol", 1.0e-6))
    dopant_states = [z for e, z in states.items() if e != "Hf"]
    return (
        cation_fractions["Hf"] >= float(cfg.get("host_min_fraction", 0.55))
        and abs(states.get("Hf", math.nan) - 4.0) <= tol
        and vacancy_site_fraction_proxy > float(cfg.get("vacancy_min", 0.0))
        and bool(dopant_states)
        and min(dopant_states) < 4.0 - tol
    )

def _estimate_hf_binary_residual(
    dopant: str, target_x_hf: float, target_vacancy_fraction: float,
    dopant_valence: float,
) -> LocalResidualSupport:
    """Nearest vacancy-fraction support, then signed residual scaling.

    No experimental target temperature is used. Equal-distance ties select the
    lower vacancy fraction. Duplicate support fractions use the median residual.
    """
    cfg = MODEL.get("hf_local_transfer", {})
    tol = float(cfg.get("hf_valence_tol", 1e-6))
    numerical_tol = 1e-12
    profile = HF_LOCAL_SUPPORT.get(str(dopant), [])
    prof = [r for r in profile
            if r["x_Hf"] >= float(cfg.get("host_min_fraction", 0.55))
            and r["nominal_vacancy_site_fraction"] > float(cfg.get("vacancy_min", 0.0))
            and r["dopant_formal_valence"] < 4.0 - tol
            and abs(r["dopant_formal_valence"] - dopant_valence) <= tol]
    if not prof:
        return LocalResidualSupport(
            dopant=dopant, dopant_share=0.0, supported=False,
            mode="no_binary_support" if not profile else "no_eligible_binary_support",
            residual_K=0.0, distance=None, support_n=0, support_formula="",
            target_vacancy_fraction=target_vacancy_fraction,
        )

    best_distance = min(abs(r["nominal_vacancy_site_fraction"] - target_vacancy_fraction)
                        for r in prof)
    nearest = [r for r in prof
               if abs(abs(r["nominal_vacancy_site_fraction"] - target_vacancy_fraction)
                      - best_distance) <= numerical_tol]
    reference_fv = min(r["nominal_vacancy_site_fraction"] for r in nearest)
    selected = [r for r in prof
                if abs(r["nominal_vacancy_site_fraction"] - reference_fv) <= numerical_tol]
    raw_residual = median(r["residual_K"] for r in selected)
    ratio = target_vacancy_fraction / reference_fv
    reference_x = median(r["x_Hf"] for r in selected)
    minimum = min(r["nominal_vacancy_site_fraction"] for r in prof)
    maximum = max(r["nominal_vacancy_site_fraction"] for r in prof)
    extrapolation = (target_vacancy_fraction < minimum - numerical_tol
                     or target_vacancy_fraction > maximum + numerical_tol)
    if abs(target_vacancy_fraction - reference_fv) <= numerical_tol:
        mode = "vacancy_scaled_exact"
    elif extrapolation:
        mode = "vacancy_scaled_extrapolation"
    else:
        mode = "vacancy_scaled_nearest"
    return LocalResidualSupport(
        dopant=dopant, dopant_share=0.0, supported=True, mode=mode,
        residual_K=raw_residual * ratio,
        distance=abs(target_x_hf - reference_x),
        support_n=sum(r["n_records_at_x"] for r in selected),
        support_formula=" || ".join(sorted(set(r["support_formula"] for r in selected))),
        support_x_Hf=reference_x,
        target_vacancy_fraction=target_vacancy_fraction,
        support_vacancy_fraction=reference_fv,
        vacancy_scaling_factor=ratio, unscaled_residual_K=raw_residual,
        vacancy_fraction_distance=abs(target_vacancy_fraction - reference_fv),
        concentration_extrapolation=extrapolation,
    )

def _hf_local_transfer(comp, states, cation_fractions):
    eligible = _hf_defect_like(comp, states, cation_fractions)
    if not eligible:
        return False, 0.0, False, 0.0, []

    xhf = float(cation_fractions["Hf"])
    dopants = [e for e in sorted(cation_fractions) if e != "Hf"]
    denom = sum(cation_fractions[e] for e in dopants)
    if denom <= 0:
        return True, 0.0, False, 0.0, []

    details: List[LocalResidualSupport] = []
    coverage = 0.0
    correction = 0.0
    for dopant in dopants:
        alpha = float(cation_fractions[dopant] / denom)
        est = _estimate_hf_binary_residual(
            dopant, xhf, nominal_vacancy_fraction(comp), states[dopant]
        )
        est.dopant_share = alpha
        details.append(est)
        if est.supported:
            coverage += alpha
            correction += alpha * est.residual_K

    cfg = MODEL.get("hf_local_transfer", {})
    applied = coverage + 1e-12 >= float(cfg.get("min_support_coverage", 0.80))
    if not applied:
        correction = 0.0
    return True, coverage, applied, correction, details

def predict_formula(formula: str) -> PredictionResult:
    comp = parse_formula(formula)

    extra_anions = [e for e in comp if e != "O" and e not in PROPS]
    # Unsupported species are also caught by the frozen reference check below.

    states, anchors = assign_formal_oxidation_states(comp)
    cations = sorted(states)
    ncat = sum(comp[e] for e in cations)
    if ncat <= 0:
        raise OxiMeltError("No cations found.")

    x = {e: comp[e] / ncat for e in cations}
    T0 = sum(x[e] * anchors[e] for e in cations)

    b0 = float(MODEL["b0_K"])
    b1 = float(MODEL["b1_K_per_descriptor"])
    b2 = float(MODEL["b2_K_per_descriptor"])
    lam = float(MODEL["lambda"])
    chi_O = float(MODEL["oxygen_pauling_en"])

    pair_rows: List[PairContribution] = []
    correction = 0.0
    seen_count = 0

    for i, j in combinations(cations, 2):
        pi, pj = PROPS[i], PROPS[j]
        rbar = 0.5 * (pi["ionic_radius_A"] + pj["ionic_radius_A"])
        kbar = 0.5 * (pi["bulk_modulus_GPa"] + pj["bulk_modulus_GPa"])
        E_strain = (
            kbar * rbar
            * (abs(pi["ionic_radius_A"] - pj["ionic_radius_A"]) / rbar) ** 2
        )
        chibar = 0.5 * (pi["pauling_en"] + pj["pauling_en"])
        denom = chi_O - chibar
        if denom <= 0:
            raise OxiMeltError(
                f"Invalid O-referenced electronegativity denominator for {i}-{j}."
            )
        E_chem = (
            0.5 * (pi["ionization_energy_eV"] + pj["ionization_energy_eV"])
            - lam * 0.5 * (pi["cohesive_energy_eV"] + pj["cohesive_energy_eV"])
        ) * (abs(pi["pauling_en"] - pj["pauling_en"]) / denom) ** 2

        omega = b0 + b1 * E_strain + b2 * E_chem
        weighted = x[i] * x[j] * omega
        correction += weighted
        seen = tuple(sorted((i, j))) in SEEN_PAIRS
        seen_count += int(seen)
        pair_rows.append(
            PairContribution(
                pair=f"{i}-{j}",
                x_i=x[i],
                x_j=x[j],
                E_strain=E_strain,
                E_chem=E_chem,
                Omega_K=omega,
                weighted_contribution_K=weighted,
                seen_in_two_cation_training=seen,
            )
        )

    pred = T0 + correction
    n_cations = len(cations)
    total_pairs = n_cations * (n_cations - 1) // 2

    warnings: List[str] = []
    if n_cations == 1:
        applicability = "REFERENCE ONLY"
        detail = (
            "Single-cation oxide: the output equals the frozen simple-oxide reference. "
            "No pair correction is applied."
        )
    elif n_cations == 2:
        if seen_count == 1:
            applicability = "TWO-CATION DOMAIN"
            detail = "The cation pair occurs in the two-cation training chemistry."
        else:
            applicability = "PAIR EXTRAPOLATION"
            detail = (
                "The cation pair was not present in the two-cation training chemistry. "
                "The frozen analytical relation is being extrapolated to this pair."
            )
    elif n_cations == 3:
        applicability = "THREE-CATION ZERO-SHOT DOMAIN"
        detail = (
            f"Three-cation composition; {seen_count}/{total_pairs} constituent pair(s) "
            "occur in the two-cation training chemistry. The model was validated on "
            "three-cation compositions by frozen zero-shot transfer."
        )
    else:
        applicability = "HIGHER-ORDER EXTRAPOLATION"
        detail = (
            f"{n_cations}-cation composition. The pairwise equation can be evaluated, "
            "but direct validation in the paper was limited to two- and three-cation oxides."
        )
        warnings.append(
            "Higher-order compositions are mathematical extrapolations beyond the directly "
            "validated two-/three-cation composition orders."
        )

    global_pred = pred
    local_eligible, local_coverage, local_applied, local_correction, local_details = \
        _hf_local_transfer(comp, states, x)
    corrected_pred = global_pred + local_correction if local_applied else global_pred

    concentration_extrapolation = any(
        d.supported and d.concentration_extrapolation for d in local_details
    )
    if local_applied:
        prediction_mode = "GLOBAL 3P + Hf-FAMILY VACANCY SCALED"
        applicability = f"{n_cations}-CATION Hf DEFECT FAMILY — VACANCY SCALED"
        detail += (
            f" Vacancy-scaled binary OOF residuals were combined with "
            f"{local_coverage:.0%} dopant-support coverage."
        )
        warnings.append(
            "The family estimate scales empirical OOF residuals using nominal "
            "stoichiometric vacancy fractions, not measured vacancy concentrations "
            "or an isolated thermodynamic vacancy contribution."
        )
        if concentration_extrapolation:
            prediction_mode += " (CONCENTRATION EXTRAPOLATION)"
            detail += " The target vacancy fraction is outside at least one dopant support range."
            warnings.append(
                "Concentration extrapolation: proportional residual scaling outside "
                "the support range has not been independently validated."
            )
        if n_cations == 2:
            warnings.append(
                "This binary family estimate uses the binary OOF support library; "
                "it is not an independent external-validation prediction."
            )
    elif local_eligible:
        prediction_mode = "GLOBAL 3P (LOCAL SUPPORT INSUFFICIENT)"
        warnings.append(
            "HfO2-rich oxygen-vacancy family detected, but frozen binary residual support "
            f"coverage is only {local_coverage:.0%}; the local correction was not applied."
        )
    else:
        prediction_mode = "GLOBAL 3P"

    return PredictionResult(
        formula=formula,
        predicted_Tm_K=corrected_pred,
        predicted_Tm_C=corrected_pred - 273.15,
        global_Tm_K=global_pred,
        global_Tm_C=global_pred - 273.15,
        local_correction_K=local_correction if local_applied else 0.0,
        corrected_Tm_K=corrected_pred,
        corrected_Tm_C=corrected_pred - 273.15,
        prediction_mode=prediction_mode,
        hf_local_domain_eligible=local_eligible,
        hf_local_correction_applied=local_applied,
        hf_local_support_coverage=local_coverage,
        T0_K=T0,
        pair_correction_K=correction,
        n_cations=n_cations,
        cation_fractions=x,
        formal_oxidation_states=states,
        applicability=applicability,
        applicability_detail=detail,
        warnings=warnings,
        pair_contributions=pair_rows,
        local_support_details=local_details,
        nominal_oxygen_vacancy_site_fraction=nominal_vacancy_fraction(comp),
        hf_local_concentration_extrapolation=concentration_extrapolation,
    )

def predict_csv(input_csv: str | Path, output_csv: str | Path) -> int:
    input_csv, output_csv = Path(input_csv), Path(output_csv)
    with input_csv.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows or "formula" not in rows[0]:
        raise OxiMeltError("Batch CSV must contain a 'formula' column.")

    out = []
    for row in rows:
        formula = row.get("formula", "").strip()
        try:
            r = predict_formula(formula)
            flat = r.to_flat_dict()
            flat["error"] = ""
        except Exception as exc:
            flat = {
                "formula": formula,
                "predicted_Tm_K": "",
                "predicted_Tm_C": "",
                "global_Tm_K": "",
                "global_Tm_C": "",
                "local_correction_K": "",
                "corrected_Tm_K": "",
                "corrected_Tm_C": "",
                "prediction_mode": "UNAVAILABLE",
                "hf_local_domain_eligible": "",
                "hf_local_correction_applied": "",
                "hf_local_support_coverage": "",
                "hf_local_support_detail": "",
                "nominal_oxygen_vacancy_site_fraction": "",
                "hf_local_concentration_extrapolation": "",
                "hf_local_support_details_json": "[]",
                "T0_K": "",
                "pair_correction_K": "",
                "n_cations": "",
                "formal_oxidation_states": "",
                "applicability": "UNAVAILABLE",
                "applicability_detail": "",
                "warnings": "",
                "error": str(exc),
            }
        out.append(flat)

    fields = list(out[0].keys())
    with output_csv.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out)
    return len(out)

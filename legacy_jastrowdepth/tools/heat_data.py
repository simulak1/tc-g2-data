#!/usr/bin/env python3
"""
Shared data model + I/O utilities for HEAT / depth-sweep analyses.

This module is intended to be reused by multiple plotting / analysis scripts.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

# Public constant (importable by plotting scripts)
KJM_PER_HARTREE = 2625.5  # kJ/mol per Hartree


# ---------------------------------------------------------------------------
# 1) HEAT & formation reference tables (copied from your script)
# ---------------------------------------------------------------------------

HEAT_TABLE = """
MOLECULE METHOD EHF     dECCSD(T)   dECCSDT dEHLC   dErel   dEZPE   dEDBOC  dESO    E       REACTION
N2       HEAT   482.94  472.26      -3.14   3.61    -0.56   -14.12  0.08    0.00    941.07  1*N2>2*N
H2       HEAT   350.81  107.41      0.00    0.00    -0.01   -25.97  0.22    0.00    432.46  1*H2>2*H
F2       HEAT   -154.98 316.00      -1.31   3.42    -0.20   -5.50   0.02    -3.01   154.44  1*F2>2*F
CO       HEAT   729.86  356.91      -2.37   2.21    -0.69   -12.98  0.07    -1.20   1071.82 1*CO>1*C+1*O
O2       HEAT   140.36  362.56      -2.07   4.46    -0.82   -9.56   0.06    -1.61   493.39  1*O2>2*O
C2H2     HEAT   1228.97 468.20      -3.03   2.23    -1.14   -68.93  0.50    -0.76   1626.06 1*C2H2>2*H+2*C
CCH      HEAT   777.34  333.02      0.56    2.28    -1.12   -36.34  0.23    -0.76   1075.23 1*CCH>1*H+2*C
CF       HEAT   338.46  216.63      -0.60   1.00    -0.72   -7.88   0.03    -1.42   545.50  1*CF>1*C+1*F
CH2      HEAT   649.11  149.13      -0.02   0.14    -0.63   -45.07  0.11    -0.38   752.39  1*CH2>2*H+1*C
CH       HEAT   238.28  113.57      0.43    0.12    -0.18   -16.97  -0.34   -0.20   334.70  1*CH>1*H+1*C
CH3      HEAT   1017.30 271.11      -0.15   0.24    -0.72   -77.64  0.19    -0.38   1209.93 1*CH3>3*H+1*C
CN       HEAT   379.38  371.64      3.14    3.43    -0.63   -12.75  0.10    -0.38   743.92  1*CN>1*C+1*N
CO2      HEAT   1033.08 599.66      -4.49   4.11    -2.05   -30.40  0.20    -2.02   1598.10 1*CO2>1*C+2*O
H2O2     HEAT   562.96  562.77      -2.37   2.95    -1.63   -68.77  0.52    -1.64   1054.81 1*H2O2>2*H+2*O
H2O      HEAT   652.40  323.02      -0.97   0.99    -1.14   -55.73  0.53    -0.82   918.26  1*H2O>2*H+1*O
HCN      HEAT   833.50  477.62      -3.21   3.04    -0.94   -41.74  0.32    -0.38   1268.20 1*HCN>1*H+1*C+1*N
HCO      HEAT   764.53  403.66      -1.78   2.19    -1.16   -34.03  -0.10   -1.20   1132.11 1*HCO>1*H+1*C+1*O
HF       HEAT   405.62  187.00      -0.67   0.72    -0.84   -24.66  0.34    -1.51   566.01  1*HF>1*H+1*F
HNO      HEAT   331.69  528.35      -2.45   3.54    -1.19   -35.92  -0.23   -0.82   822.97  1*HNO>1*H+1*N+1*O
HO2      HEAT   300.93  431.57      -0.60   2.52    -1.20   -37.05  0.02    -1.64   694.56  1*HO2>1*H+2*O
NH2      HEAT   492.92  270.76      0.06    0.48    -0.67   -49.57  -0.03   0.00    713.95  1*NH2>2*H+1*N
NH3      HEAT   841.38  406.03      -0.63   0.72    -1.08   -89.45  0.56    0.00    1157.53 1*NH3>3*H+1*N
NH       HEAT   214.88  132.40      0.29    0.22    -0.31   -19.46  -0.19   0.00    327.83  1*NH>1*H+1*N
NO       HEAT   225.52  411.50      -1.54   3.07    -0.83   -11.46  0.06    -0.10   626.22  1*NO>1*N+1*O
OF       HEAT   -64.04  282.44      1.05    2.10    -0.41   -6.37   0.02    -1.24   213.56  1*OF>1*O+1*F
OH       HEAT   286.41  162.12      -0.13   0.48    -0.54   -22.21  0.05    -0.04   426.14  1*OH>1*H+1*O
""".strip()

FORMATION_TABLE = """
C2H2 356.157
CCH 491.011
CH2 284.324
CH 366.917
CH3 184.906
CO2 -15.4
H2O2 -61.999
H2O -100.761
HCO 152.272
HF -107.4065
HO2 0
NO 35.0185
OH 12.526
HNO 37.5095
CN 307.3225
HCN 183.8215
CF 233.3675
NH2 65.4465
NH3 -31.3565
NH 136.7945
OF 42.9115
""".strip()


# ---------------------------------------------------------------------------
# 2) Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class HeatRef:
    molecule: str         # e.g. CO2 (with CCH mapped to C2H)
    reaction: str         # e.g. "1*CO2>1*C+2*O"
    elec_au: float        # EHF + dECCSD(T) + dECCSDT + dEHLC [Hartree]


@dataclass
class FormRef:
    molecule: str     # e.g. C2H2 (with CCH mapped to C2H)
    form_au: float    # formation energy [Hartree] (note: currently kJ→/1000)


@dataclass
class CCRecord:
    project: str          # from 'projectname'
    core: str             # 'eCEPP', 'ccECP', 'AE', ...
    molecule: str
    basis: str
    method: str
    energy_au: float
    depth: float | None   # from 'ee_depth' (numeric) or None for 'free'


# ---------------------------------------------------------------------------
# 3) Small helpers for parsing
# ---------------------------------------------------------------------------

_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> Dict[str, int]:
    """
    Parse e.g. 'C2H2' -> {'C': 2, 'H': 2}, 'HO2' -> {'H': 1, 'O': 2}.
    """
    counts: Dict[str, int] = {}
    for elem, num_str in _FORMULA_RE.findall(formula):
        n = int(num_str) if num_str else 1
        counts[elem] = counts.get(elem, 0) + n
    return counts


def parse_formation_table(raw: str = FORMATION_TABLE) -> Dict[str, FormRef]:
    """
    Parse FORMATION_TABLE text.

    Each line: MOL  E_form[kJ/mol]

    Convert to "Hartree-like" units here by dividing by 1000 as in your original
    script. (You can refine that later if you want strict kJ/mol→Ha via KJM_PER_HARTREE.)
    """
    refs: Dict[str, FormRef] = {}
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(f"Unexpected formation line:\n{line}")
        mol = parts[0]
        val_kj = float(parts[1])
        mol_key = "C2H" if mol == "CCH" else mol
        refs[mol_key] = FormRef(molecule=mol_key, form_au=val_kj / 1000.0)
    return refs


def parse_heat_table(raw: str = HEAT_TABLE) -> Dict[str, HeatRef]:
    """
    Parse HEAT_TABLE text.

    Columns (kJ/mol):
      EHF, dECCSD(T), dECCSDT, dEHLC, dErel, dEZPE, dEDBOC, dESO, E
    We define the electronic reference as:
      E_elec = EHF + dECCSD(T) + dECCSDT + dEHLC
    and convert to Hartree.
    """
    refs: Dict[str, HeatRef] = {}

    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith("MOLECULE"):
            continue

        parts = line.split()
        mol = parts[0]
        method = parts[1]
        if method != "HEAT":
            continue

        numeric_strs = parts[2:-1]
        if len(numeric_strs) != 9:
            raise ValueError(f"Unexpected numeric column count in line:\n{line}")

        nums = list(map(float, numeric_strs))
        EHF = nums[0]
        dECCSDT_corr = nums[1]
        dECCSDT_triples = nums[2]
        dEHLC = nums[3]

        elec_kjmol = EHF + dECCSDT_corr + dECCSDT_triples + dEHLC
        elec_au = elec_kjmol / KJM_PER_HARTREE

        reaction = parts[-1]
        mol_key = "C2H" if mol == "CCH" else mol

        refs[mol_key] = HeatRef(molecule=mol_key, reaction=reaction, elec_au=elec_au)

    return refs


# ---------------------------------------------------------------------------
# 4) CSV reading
# ---------------------------------------------------------------------------

def read_ccsd_csv(path: Path) -> List[CCRecord]:
    """
    Read the CSV format with columns like:

        projectname, ee_depth, en_depth, een_depth,
        molecule, core, basis, method, energy_au,
        spin_dependent_ee_cusp, spin_dependent_ee,
        spin_dependent_en, spin_dependent_een

    We currently ignore:
        en_depth, een_depth, and all spin_dependent_* columns.
    """
    records: List[CCRecord] = []
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            depth: float | None = None
            if row["core"] == "AE":
                continue
            if "ee_depth" in row and row["ee_depth"] is not None:
                s = row["ee_depth"].strip()
                if s and s.lower() != "free":
                    depth = float(s)

            core = row.get("core", "unknown").strip() if row.get("core") is not None else "unknown"

            rec = CCRecord(
                project=row["projectname"],
                core=core,
                molecule=row["molecule"],
                basis=row["basis"],
                method=row["method"],
                energy_au=float(row["energy_au"]),
                depth=depth,
            )
            records.append(rec)
    return records


# ---------------------------------------------------------------------------
# 5) Helpers for reactions & energy construction
# ---------------------------------------------------------------------------

def parse_side(side: str) -> Dict[str, int]:
    """Parse e.g. '2*H+1*C' → {'H': 2, 'C': 1}."""
    species: Dict[str, int] = {}
    side = side.strip()
    if not side:
        return species

    for term in side.split("+"):
        term = term.strip()
        if not term:
            continue
        coef_str, mol = term.split("*")
        coef = int(coef_str)
        species[mol] = species.get(mol, 0) + coef
    return species


def parse_reaction(reaction: str) -> Tuple[Dict[str, int], Dict[str, int]]:
    """'1*CO2>1*C+2*O' → (lhs, rhs) dicts."""
    lhs_str, rhs_str = reaction.split(">")
    lhs = parse_side(lhs_str)
    rhs = parse_side(rhs_str)
    return lhs, rhs


def build_energy_index(
    records: List[CCRecord],
    basis: str | None = None,
) -> tuple[Dict[Tuple[str, str], float], Dict[str, float | None]]:
    """
    Build:
      energies[(project, molecule)] = energy_au
      project_depths[project] = depth (float or None)

    If basis is not None, filter to that basis.
    """
    energies: Dict[Tuple[str, str], float] = {}
    project_depths: Dict[str, float | None] = {}

    for rec in records:
        if basis is not None and rec.basis != basis:
            continue

        key = (rec.project, rec.molecule)
        energies[key] = rec.energy_au

        if rec.project not in project_depths:
            project_depths[rec.project] = rec.depth
        else:
            if project_depths[rec.project] is None and rec.depth is not None:
                project_depths[rec.project] = rec.depth

    return energies, project_depths


def compute_our_atomization(
    energies: Dict[Tuple[str, str], float],
    heat_refs: Dict[str, HeatRef],
) -> Dict[Tuple[str, str], float]:
    """
    Return our atomization energies in Hartree for each (project, molecule),
    using the HEAT reaction definition:

        E_form(project, mol) = sum_nu E(prod) - sum_nu E(react)   [Ha]
    """
    result: Dict[Tuple[str, str], float] = {}

    projects = sorted({p for (p, _) in energies.keys()})

    for mol, href in heat_refs.items():
        lhs, rhs = parse_reaction(href.reaction)

        for proj in projects:
            E_reac_au = 0.0
            E_prod_au = 0.0
            ok = True

            for sp, coef in lhs.items():
                key = (proj, sp)
                if key not in energies:
                    ok = False
                    break
                E_reac_au += coef * energies[key]
            if not ok:
                continue

            for sp, coef in rhs.items():
                key = (proj, sp)
                if key not in energies:
                    if key[1] == "H":
                        E_prod_au += coef * -0.5
                        continue
                    else:
                        ok = False
                        break
                E_prod_au += coef * energies[key]
            if not ok:
                continue

            E_form_au = E_prod_au - E_reac_au
            result[(proj, mol)] = E_form_au

    return result


def compute_our_formation(
    energies: Dict[Tuple[str, str], float],
    form_refs: Dict[str, FormRef],
) -> Dict[Tuple[str, str], float]:
    """
    Compute our formation energies in Hartree for each (project, molecule) in
    form_refs using elemental references:

        μ_H = 0.5 * E(H2)
        μ_N = 0.5 * E(N2)
        μ_O = 0.5 * E(O2)
        μ_F = 0.5 * E(F2)
        μ_C = E(CO) - E(O)

      E_form(project, mol) = E(mol) - sum_elem n_elem(mol) * μ_elem
    """
    result: Dict[Tuple[str, str], float] = {}

    projects = sorted({p for (p, _) in energies.keys()})
    stoich: Dict[str, Dict[str, int]] = {
        mol: parse_formula(mol) for mol in form_refs.keys()
    }

    for proj in projects:
        def get_E(sp: str) -> float | None:
            return energies.get((proj, sp))

        mu: Dict[str, float] = {}

        E_H2 = get_E("H2")
        if E_H2 is not None:
            mu["H"] = 0.5 * E_H2

        E_N2 = get_E("N2")
        if E_N2 is not None:
            mu["N"] = 0.5 * E_N2

        E_O2 = get_E("O2")
        if E_O2 is not None:
            mu["O"] = 0.5 * E_O2

        E_F2 = get_E("F2")
        if E_F2 is not None:
            mu["F"] = 0.5 * E_F2

        E_CO = get_E("CO")
        E_O_atom = get_E("O")
        if E_CO is not None and E_O_atom is not None:
            mu["C"] = E_CO - E_O_atom

        if not mu:
            continue

        for mol, fref in form_refs.items():
            E_mol = energies.get((proj, mol))
            if E_mol is None:
                continue

            elems = stoich[mol]
            if any(elem not in mu for elem in elems):
                continue

            mu_tot = 0.0
            for elem, n_e in elems.items():
                mu_tot += n_e * mu[elem]

            E_form_au = E_mol - mu_tot
            result[(proj, mol)] = E_form_au

    return result

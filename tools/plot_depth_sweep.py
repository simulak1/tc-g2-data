#!/usr/bin/env python3
"""
Compare HEAT reference atomization energies to CCSD(T) results from ccsd_energies.csv,
and also compare formation energies to given benchmarks.

- Atomization:
    E_form_atom(project, mol) from HEAT reaction, as before.
    ΔE_atom = E_form_atom(ours) - E_form_atom(HEAT).

- Formation (new):
    Our formation energies are defined with elemental references:
        μ_H = 0.5 * E(H2)
        μ_N = 0.5 * E(N2)
        μ_O = 0.5 * E(O2)
        μ_F = 0.5 * E(F2)
        μ_C = E(CO) - E(O)

      E_form(project, mol) = E(mol) - sum_elem n_elem(mol) * μ_elem

    These are compared to the given benchmark formation energies (kJ/mol).

Layout
------
If --basis-short is given:
    3×2 grid:
      row 1: no-H molecules            (ΔE_atom | ΔE_form)
      row 2: H-containing molecules    (ΔE_atom | ΔE_form)
      row 3: MAE(ΔE_atom) vs project   (left only)

If --basis-short omitted:
    5×2 grid:
      rows 1–2: avdz (no-H / with-H)
      rows 3–4: avtz (no-H / with-H)
      row 5: MAE(ΔE_atom) vs project (avdz | avtz)

Usage
-----
    # Single basis
    python plot_heat_delta_vs_project.py ccsd_energies.csv \
        --basis-short avdz \
        --unit mHa \
        --out-csv heat_delta_vs_project.csv

    # avdz + avtz in one figure
    python plot_heat_delta_vs_project.py ccsd_energies.csv \
        --unit mHa \
        --out-csv heat_delta_vs_project.csv
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np


KJM_PER_HARTREE = 2625.5  # kJ/mol per Hartree


# ---------------------------------------------------------------------------
# 1) HEAT reference table (for atomization)
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


@dataclass
class HeatRef:
    molecule: str
    reaction: str
    elec_au: float  # EHF + dECCSD(T) + dECCSDT + dEHLC [Ha]


def parse_heat_table(raw: str) -> Dict[str, HeatRef]:
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
# 2) Formation energy benchmarks (kJ/mol)
# ---------------------------------------------------------------------------

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


def parse_formation_table(raw: str) -> Dict[str, float]:
    """
    Return formation reference energies in Hartree:
        form_refs[mol] = E_form_ref [Ha]
    """
    refs: Dict[str, float] = {}
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(f"Bad formation line: {line}")
        mol = parts[0]
        val_kj = float(parts[1])
        mol_key = "C2H" if mol == "CCH" else mol
        refs[mol_key] = -1.*val_kj / KJM_PER_HARTREE
    return refs


# ---------------------------------------------------------------------------
# 3) CCSD(T) energies from CSV
# ---------------------------------------------------------------------------

@dataclass
class CCRecord:
    project: str        # from 'projectname'
    core: str           # 'AE', 'eCEPP', 'ccECP', ...
    molecule: str
    basis_short: str
    method: str
    energy_au: float
    # ee_depth, en_depth, een_depth, spin_dependent_* are ignored for the moment

def read_ccsd_csv(path: Path) -> List[CCRecord]:
    records: List[CCRecord] = []
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rec = CCRecord(
                project=row["projectname"],
                core=row["core"],
                molecule=row["molecule"],
                basis_short=row["basis_short"],
                method=row["method"],
                energy_au=float(row["energy_au"]),
            )
            records.append(rec)
    return records

# ---------------------------------------------------------------------------
# 4) HEAT reaction parsing
# ---------------------------------------------------------------------------

def parse_side(side: str) -> Dict[str, int]:
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
    lhs_str, rhs_str = reaction.split(">")
    return parse_side(lhs_str), parse_side(rhs_str)


# ---------------------------------------------------------------------------
# 5) Energy index
# ---------------------------------------------------------------------------

def build_energy_index(records: List[CCRecord], basis_short: str | None = None):
    energies: Dict[Tuple[str, str], float] = {}
    for rec in records:
        if basis_short is not None and rec.basis_short != basis_short:
            continue
        key = (rec.project, rec.molecule)
        energies[key] = rec.energy_au
    return energies

# ---------------------------------------------------------------------------
# 6) Our atomization energies (HEAT reaction)
# ---------------------------------------------------------------------------

def compute_our_atomization(
    energies: Dict[Tuple[str, str], float],
    heat_refs: Dict[str, HeatRef],
) -> Dict[Tuple[str, str], float]:
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
                    ok = False
                    break
                E_prod_au += coef * energies[key]
            if not ok:
                continue

            result[(proj, mol)] = E_prod_au - E_reac_au
    return result


# ---------------------------------------------------------------------------
# 7) Our formation energies from elemental references
# ---------------------------------------------------------------------------

FORMULA_COUNTS: Dict[str, Dict[str, int]] = {
    "N2": {"N": 2},
    "H2": {"H": 2},
    "F2": {"F": 2},
    "CO": {"C": 1, "O": 1},
    "O2": {"O": 2},
    "C2H2": {"C": 2, "H": 2},
    "C2H": {"C": 2, "H": 1},
    "CF": {"C": 1, "F": 1},
    "CH2": {"C": 1, "H": 2},
    "CH": {"C": 1, "H": 1},
    "CH3": {"C": 1, "H": 3},
    "CN": {"C": 1, "N": 1},
    "CO2": {"C": 1, "O": 2},
    "H2O2": {"H": 2, "O": 2},
    "H2O": {"H": 2, "O": 1},
    "HCN": {"H": 1, "C": 1, "N": 1},
    "HCO": {"H": 1, "C": 1, "O": 1},
    "HF": {"H": 1, "F": 1},
    "HNO": {"H": 1, "N": 1, "O": 1},
    "HO2": {"H": 1, "O": 2},
    "NH2": {"N": 1, "H": 2},
    "NH3": {"N": 1, "H": 3},
    "NH": {"N": 1, "H": 1},
    "NO": {"N": 1, "O": 1},
    "OF": {"O": 1, "F": 1},
    "OH": {"O": 1, "H": 1},
}


def compute_our_formation(
    energies: Dict[Tuple[str, str], float]
) -> Dict[Tuple[str, str], float]:
    """
    Our formation energies (Hartree) for each (project, molecule) where possible,
    using the elemental reference scheme specified in the user’s message.
    """
    result: Dict[Tuple[str, str], float] = {}
    projects = sorted({p for (p, _) in energies.keys()})

    for proj in projects:
        def E(m: str):
            return energies.get((proj, m), None)

        required = ["H2", "N2", "O2", "F2", "CO", "O"]
        if any(E(m) is None for m in required):
            continue

        mu_H = 0.5 * E("H2")
        mu_N = 0.5 * E("N2")
        mu_O = 0.5 * E("O2")
        mu_F = 0.5 * E("F2")
        mu_C = E("CO") - E("O")

        mus = {"H": mu_H, "N": mu_N, "O": mu_O, "F": mu_F, "C": mu_C}

        for (p, mol), E_mol in energies.items():
            if p != proj:
                continue
            comp = FORMULA_COUNTS.get(mol)
            if comp is None:
                continue

            ref_E = 0.0
            ok = True
            for elem, count in comp.items():
                if elem not in mus:
                    ok = False
                    break
                ref_E += count * mus[elem]
            if not ok:
                continue

            result[(proj, mol)] = E_mol - ref_E
    return result


# ---------------------------------------------------------------------------
# 8) Plotting helpers
# ---------------------------------------------------------------------------

def _get_scale_and_labels(unit: str):
    if unit == "Ha":
        scale = 1.0
        delta_label = r"$\Delta E_\mathrm{atom}$ (ours − HEAT) [Ha]"
        form_label = r"$\Delta E_\mathrm{form}$ (ours − benchmark) [Ha]"
        mae_label = r"MAE(ΔE_atom) [Ha]"
    elif unit == "mHa":
        scale = 1000.0
        delta_label = r"$\Delta E_\mathrm{atom}$ (ours − HEAT) [mHa]"
        form_label = r"$\Delta E_\mathrm{form}$ (ours − benchmark) [mHa]"
        mae_label = r"MAE(ΔE_atom) [mHa]"
    else:
        raise ValueError("unit must be 'Ha' or 'mHa'")
    return scale, delta_label, form_label, mae_label


def _split_molecules(heat_refs: Dict[str, HeatRef]):
    molecules = sorted(heat_refs.keys())
    mol_no_H = [m for m in molecules if "H" not in m]
    mol_with_H = [m for m in molecules if "H" in m]
    return molecules, mol_no_H, mol_with_H


def _plot_group_delta_and_form(
    ax_delta,
    ax_form,
    mol_list: List[str],
    projects: List[str],
    project_index: Dict[str, int],
    our_atom_E: Dict[Tuple[str, str], float],
    heat_refs: Dict[str, HeatRef],
    our_form_E: Dict[Tuple[str, str], float],
    form_refs: Dict[str, float],
    scale: float,
    basis_label: str | None,
    csv_rows: List[Dict[str, str]] | None,
    deltas_by_project: Dict[str, List[float]] | None,
):
    """
    For each molecule in mol_list:
      left axis:  ΔE_atom(project, mol) = our_atom - HEAT_elec  [scaled]
      right axis: ΔE_form(project, mol) = our_form - benchmark [scaled]
    """
    for mol in mol_list:
        href = heat_refs[mol]
        E_ref_atom = href.elec_au

        xs_delta, ys_delta = [], []
        xs_form, ys_form = [], []

        form_ref_au = form_refs.get(mol)

        for proj in projects:
            key = (proj, mol)
            if key not in our_atom_E:
                continue

            E_our_atom = our_atom_E[key]
            delta_atom_au = E_our_atom - E_ref_atom

            xs_delta.append(project_index[proj])
            ys_delta.append(delta_atom_au * scale)

            if deltas_by_project is not None:
                deltas_by_project[proj].append(delta_atom_au)

            if csv_rows is not None:
                row = {
                    "molecule": mol,
                    "project": proj,
                    "E_our_Ha": f"{E_our_atom:.10f}",
                    "E_HEAT_Ha": f"{E_ref_atom:.10f}",
                    "delta_Ha": f"{delta_atom_au:.10f}",
                }
                if basis_label is not None:
                    row["basis_short"] = basis_label
                csv_rows.append(row)

            if form_ref_au is not None and key in our_form_E:
                delta_form_au = our_form_E[key] - form_ref_au
                xs_form.append(project_index[proj])
                ys_form.append(delta_form_au * scale)

        if xs_delta:
            xs_d = np.array(xs_delta)
            ys_d = np.array(ys_delta)
            order = np.argsort(xs_d)
            ax_delta.plot(xs_d[order], ys_d[order], marker="o", label=mol)

        if xs_form:
            xs_f = np.array(xs_form)
            ys_f = np.array(ys_form)
            order = np.argsort(xs_f)
            ax_form.plot(xs_f[order], ys_f[order], marker="o", label=mol)


# ---------------------------------------------------------------------------
# 9) Single-basis plotting (3×2 grid)
# ---------------------------------------------------------------------------

def plot_deltas_single_basis(
    our_atom_E: Dict[Tuple[str, str], float],
    our_form_E: Dict[Tuple[str, str], float],
    heat_refs: Dict[str, HeatRef],
    form_refs: Dict[str, float],
    out_csv: Path | None,
    unit: str,
    basis_short: str,
    core_label: str | None = None,
) -> None:

    core_tag = f", core={core_label}" if core_label is not None else ""

    scale, delta_label, form_label, mae_label = _get_scale_and_labels(unit)
    molecules, mol_no_H, mol_with_H = _split_molecules(heat_refs)
    projects = sorted({p for (p, _) in our_atom_E.keys()})
    project_index = {p: i for i, p in enumerate(projects)}

    fig, axes = plt.subplots(3, 2, sharex=True, figsize=(12, 9))
    ax_noH_delta, ax_noH_form = axes[0]
    ax_H_delta, ax_H_form = axes[1]
    ax_stats_delta, ax_stats_dummy = axes[2]

    csv_rows: List[Dict[str, str]] = []
    deltas_by_project: Dict[str, List[float]] = {p: [] for p in projects}

    # Row 1: non-H
    if mol_no_H:
        _plot_group_delta_and_form(
            ax_noH_delta,
            ax_noH_form,
            mol_no_H,
            projects,
            project_index,
            our_atom_E,
            heat_refs,
            our_form_E,
            form_refs,
            scale,
            basis_short,
            csv_rows,
            deltas_by_project,
        )
        ax_noH_delta.axhline(0.0, linestyle="--")
        ax_noH_delta.set_ylabel(delta_label)
        ax_noH_delta.set_title(    f"ΔE_atom vs HEAT (no H), basis={basis_short}{core_tag}"
)
        ax_noH_delta.legend(fontsize="x-small", ncol=3)

        ax_noH_form.axhline(0.0, linestyle="--")
        ax_noH_form.set_ylabel(form_label)
        ax_noH_form.set_title(f"ΔE_form vs benchmark (no H), basis={basis_short}{core_tag}"
)
        ax_noH_form.legend(fontsize="x-small", ncol=3)

    # Row 2: H-containing
    if mol_with_H:
        _plot_group_delta_and_form(
            ax_H_delta,
            ax_H_form,
            mol_with_H,
            projects,
            project_index,
            our_atom_E,
            heat_refs,
            our_form_E,
            form_refs,
            scale,
            basis_short,
            csv_rows,
            deltas_by_project,
        )
        ax_H_delta.axhline(0.0, linestyle="--")
        ax_H_delta.set_ylabel(delta_label)
        ax_H_delta.set_title(
            f"ΔE_atom vs HEAT (with H), basis={basis_short}{core_tag}"
        )        
        ax_H_delta.legend(fontsize="x-small", ncol=3)

        ax_H_form.axhline(0.0, linestyle="--")
        ax_H_form.set_ylabel(form_label)
        ax_H_form.set_title(
            f"ΔE_form vs benchmark (with H), basis={basis_short}{core_tag}"
        )
        ax_H_form.legend(fontsize="x-small", ncol=3)

    # Row 3: MAE(ΔE_atom)
    xs_stats, mae_vals = [], []
    for proj in projects:
        deltas = deltas_by_project[proj]
        if not deltas:
            continue
        d_arr = np.array(deltas)
        mae = np.mean(np.abs(d_arr)) * scale
        xs_stats.append(project_index[proj])
        mae_vals.append(mae)

    if xs_stats:
        xs_stats = np.array(xs_stats)
        order = np.argsort(xs_stats)
        ax_stats_delta.plot(xs_stats[order], np.array(mae_vals)[order],
                            marker="o", label=f"MAE({unit})")
        ax_stats_delta.axhline(0.0, linestyle="--")
        ax_stats_delta.set_ylabel(mae_label)
        ax_stats_delta.set_title(
            f"MAE(ΔE_atom) vs project, basis={basis_short}{core_tag}"
        )
        ax_stats_delta.legend(fontsize="small")

    if core_label is not None:
        fig.suptitle(f"core = {core_label}")

    ax_stats_delta.set_xlabel("Project")
    ax_stats_delta.set_xticks(np.arange(len(projects)))
    ax_stats_delta.set_xticklabels(projects, rotation=45, ha="right")
    ax_stats_dummy.axis("off")

    fig.tight_layout()
    plt.show()

    if out_csv is not None:
        with out_csv.open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "molecule",
                    "project",
                    "basis_short",
                    "E_our_Ha",
                    "E_HEAT_Ha",
                    "delta_Ha",
                ],
            )
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"Wrote delta CSV to {out_csv}")


# ---------------------------------------------------------------------------
# 10) Two-basis plotting (5×2 grid for avdz & avtz)
# ---------------------------------------------------------------------------

def plot_deltas_two_bases(
    our_atom_avdz: Dict[Tuple[str, str], float],
    our_form_avdz: Dict[Tuple[str, str], float],
    our_atom_avtz: Dict[Tuple[str, str], float],
    our_form_avtz: Dict[Tuple[str, str], float],
    heat_refs: Dict[str, HeatRef],
    form_refs: Dict[str, float],
    out_csv: Path | None,
    unit: str,
    core_label: str | None = None,
) -> None:

    core_tag = f", core={core_label}" if core_label is not None else ""

    scale, delta_label, form_label, mae_label = _get_scale_and_labels(unit)
    molecules, mol_no_H, mol_with_H = _split_molecules(heat_refs)

    fig, axes = plt.subplots(5, 2, figsize=(12, 14))
    csv_rows: List[Dict[str, str]] = []

    # --- AVDZ block (rows 0–1)
    if our_atom_avdz:
        projects_avdz = sorted({p for (p, _) in our_atom_avdz.keys()})
        proj_idx_avdz = {p: i for i, p in enumerate(projects_avdz)}
        deltas_by_project_avdz: Dict[str, List[float]] = {p: [] for p in projects_avdz}

        ax_noH_delta, ax_noH_form = axes[0]
        ax_H_delta, ax_H_form = axes[1]

        if mol_no_H:
            _plot_group_delta_and_form(
                ax_noH_delta,
                ax_noH_form,
                mol_no_H,
                projects_avdz,
                proj_idx_avdz,
                our_atom_avdz,
                heat_refs,
                our_form_avdz,
                form_refs,
                scale,
                "avdz",
                csv_rows,
                deltas_by_project_avdz,
            )
            ax_noH_delta.axhline(0.0, linestyle="--")
            ax_noH_delta.set_ylabel(delta_label)
            ax_noH_delta.set_title(f"ΔE_atom vs HEAT (no H), avdz{core_tag}")
            ax_noH_delta.legend(fontsize="x-small", ncol=3)

            ax_noH_form.axhline(0.0, linestyle="--")
            ax_noH_form.set_ylabel(form_label)
            ax_noH_form.set_title(f"ΔE_form vs benchmark (no H), avdz{core_tag}")
            ax_noH_form.legend(fontsize="x-small", ncol=3)

        if mol_with_H:
            _plot_group_delta_and_form(
                ax_H_delta,
                ax_H_form,
                mol_with_H,
                projects_avdz,
                proj_idx_avdz,
                our_atom_avdz,
                heat_refs,
                our_form_avdz,
                form_refs,
                scale,
                "avdz",
                csv_rows,
                deltas_by_project_avdz,
            )
            ax_H_delta.axhline(0.0, linestyle="--")
            ax_H_delta.set_ylabel(delta_label)
            ax_H_delta.set_title(f"ΔE_atom vs HEAT (with H), avdz{core_tag}")
            ax_H_delta.legend(fontsize="x-small", ncol=3)

            ax_H_form.axhline(0.0, linestyle="--")
            ax_H_form.set_ylabel(form_label)
            ax_H_form.set_title(f"ΔE_form vs benchmark (with H), avdz{core_tag}")
            ax_H_form.legend(fontsize="x-small", ncol=3)
    else:
        axes[0, 0].axis("off")
        axes[0, 1].axis("off")
        axes[1, 0].axis("off")
        axes[1, 1].axis("off")
        projects_avdz = []
        proj_idx_avdz = {}
        deltas_by_project_avdz = {}

    # --- AVTZ block (rows 2–3)
    if our_atom_avtz:
        projects_avtz = sorted({p for (p, _) in our_atom_avtz.keys()})
        proj_idx_avtz = {p: i for i, p in enumerate(projects_avtz)}
        deltas_by_project_avtz: Dict[str, List[float]] = {p: [] for p in projects_avtz}

        ax_noH_delta_tz, ax_noH_form_tz = axes[2]
        ax_H_delta_tz, ax_H_form_tz = axes[3]

        if mol_no_H:
            _plot_group_delta_and_form(
                ax_noH_delta_tz,
                ax_noH_form_tz,
                mol_no_H,
                projects_avtz,
                proj_idx_avtz,
                our_atom_avtz,
                heat_refs,
                our_form_avtz,
                form_refs,
                scale,
                "avtz",
                csv_rows,
                deltas_by_project_avtz,
            )
            ax_noH_delta_tz.axhline(0.0, linestyle="--")
            ax_noH_delta_tz.set_ylabel(delta_label)
            ax_noH_delta_tz.set_title(f"ΔE_atom vs HEAT (no H), avtz{core_tag}")
            ax_noH_delta_tz.legend(fontsize="x-small", ncol=3)

            ax_noH_form_tz.axhline(0.0, linestyle="--")
            ax_noH_form_tz.set_ylabel(form_label)
            ax_noH_form_tz.set_title(f"ΔE_form vs benchmark (no H), avtz{core_tag}")
            ax_noH_form_tz.legend(fontsize="x-small", ncol=3)

        if mol_with_H:
            _plot_group_delta_and_form(
                ax_H_delta_tz,
                ax_H_form_tz,
                mol_with_H,
                projects_avtz,
                proj_idx_avtz,
                our_atom_avtz,
                heat_refs,
                our_form_avtz,
                form_refs,
                scale,
                "avtz",
                csv_rows,
                deltas_by_project_avtz,
            )
            ax_H_delta_tz.axhline(0.0, linestyle="--")
            ax_H_delta_tz.set_ylabel(delta_label)
            ax_H_delta_tz.set_title(f"ΔE_atom vs HEAT (with H), avtz{core_tag}")
            ax_H_delta_tz.legend(fontsize="x-small", ncol=3)

            ax_H_form_tz.axhline(0.0, linestyle="--")
            ax_H_form_tz.set_ylabel(form_label)
            ax_H_form_tz.set_title(f"ΔE_form vs benchmark (with H), avtz{core_tag}")
            ax_H_form_tz.legend(fontsize="x-small", ncol=3)
    else:
        axes[2, 0].axis("off")
        axes[2, 1].axis("off")
        axes[3, 0].axis("off")
        axes[3, 1].axis("off")
        projects_avtz = []
        proj_idx_avtz = {}
        deltas_by_project_avtz = {}

    # --- Row 5: MAE(ΔE_atom) vs project for avdz & avtz
    ax_mae_avdz, ax_mae_avtz = axes[4]

    if our_atom_avdz:
        xs, maes = [], []
        for proj in projects_avdz:
            d = deltas_by_project_avdz[proj]
            if not d:
                continue
            arr = np.array(d)
            maes.append(np.mean(np.abs(arr)) * scale)
            xs.append(proj_idx_avdz[proj])
        if xs:
            xs = np.array(xs)
            order = np.argsort(xs)
            ax_mae_avdz.plot(xs[order], np.array(maes)[order], marker="o")
        ax_mae_avdz.axhline(0.0, linestyle="--")
        ax_mae_avdz.set_ylabel(mae_label)
        ax_mae_avdz.set_title(f"MAE(ΔE_atom) vs project (avdz{core_tag})")
        ax_mae_avdz.set_xlabel("Project")
        ax_mae_avdz.set_xticks(np.arange(len(projects_avdz)))
        ax_mae_avdz.set_xticklabels(projects_avdz, rotation=45, ha="right")
    else:
        ax_mae_avdz.axis("off")

    if our_atom_avtz:
        xs, maes = [], []
        for proj in projects_avtz:
            d = deltas_by_project_avtz[proj]
            if not d:
                continue
            arr = np.array(d)
            maes.append(np.mean(np.abs(arr)) * scale)
            xs.append(proj_idx_avtz[proj])
        if xs:
            xs = np.array(xs)
            order = np.argsort(xs)
            ax_mae_avtz.plot(xs[order], np.array(maes)[order], marker="o")
        ax_mae_avtz.axhline(0.0, linestyle="--")
        ax_mae_avtz.set_ylabel(mae_label)
        ax_mae_avtz.set_title(f"MAE(ΔE_atom) vs project (avtz{core_tag})")
        ax_mae_avtz.set_xlabel("Project")
        ax_mae_avtz.set_xticks(np.arange(len(projects_avtz)))
        ax_mae_avtz.set_xticklabels(projects_avtz, rotation=45, ha="right")
    else:
        ax_mae_avtz.axis("off")

    if core_label is not None:
        fig.suptitle(f"core = {core_label}")

    fig.tight_layout()
    plt.show()

    if out_csv is not None and csv_rows:
        with out_csv.open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "molecule",
                    "project",
                    "basis_short",
                    "E_our_Ha",
                    "E_HEAT_Ha",
                    "delta_Ha",
                ],
            )
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"Wrote delta CSV to {out_csv}")


# ---------------------------------------------------------------------------
# 11) CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compare CCSD(T) atomization energies vs HEAT references and "
            "formation energies vs given benchmarks (Ha / mHa)."
        )
    )
    parser.add_argument(
        "ccsd_csv",
        type=Path,
        help="Path to ccsd_energies.csv (from collect_ccsd_csv.py)",
    )
    parser.add_argument(
        "--basis-short",
        default=None,
        help=(
            "If given, plot only this basis in a 3×2 grid. "
            "If omitted, plot avdz+avtz in a 5×2 grid."
        ),
    )
    parser.add_argument(
        "--unit",
        choices=["Ha", "mHa"],
        default="mHa",
        help="Unit for y-axis (default: mHa).",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=None,
        help="Output CSV filename for ΔE_atom values (Hartree).",
    )

    args = parser.parse_args()

    heat_refs = parse_heat_table(HEAT_TABLE)
    form_refs = parse_formation_table(FORMATION_TABLE)
    records = read_ccsd_csv(args.ccsd_csv)

    # Unique cores in the CSV (e.g. 'eCEPP', 'ccECP')
    cores = sorted({rec.core for rec in records})

    for core in cores:
        core_records = [rec for rec in records if rec.core == core]
        print(f"Processing core = {core} with {len(core_records)} records")

        if args.basis_short is not None:
            energies = build_energy_index(core_records, basis_short=args.basis_short)
            our_atom = compute_our_atomization(energies, heat_refs)
            our_form = compute_our_formation(energies)
            plot_deltas_single_basis(
                our_atom_E=our_atom,
                our_form_E=our_form,
                heat_refs=heat_refs,
                form_refs=form_refs,
                out_csv=args.out_csv,
                unit=args.unit,
                basis_short=args.basis_short,
                core_label=core,
            )
        else:
            energies_avdz = build_energy_index(core_records, basis_short="avdz")
            energies_avtz = build_energy_index(core_records, basis_short="avtz")

            our_atom_avdz = compute_our_atomization(energies_avdz, heat_refs)
            our_atom_avtz = compute_our_atomization(energies_avtz, heat_refs)

            our_form_avdz = compute_our_formation(energies_avdz)
            our_form_avtz = compute_our_formation(energies_avtz)

            plot_deltas_two_bases(
                our_atom_avdz=our_atom_avdz,
                our_form_avdz=our_form_avdz,
                our_atom_avtz=our_atom_avtz,
                our_form_avtz=our_form_avtz,
                heat_refs=heat_refs,
                form_refs=form_refs,
                out_csv=args.out_csv,
                unit=args.unit,
                core_label=core,
            )

if __name__ == "__main__":
    main()

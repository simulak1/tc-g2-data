#!/usr/bin/env python3
"""
Compare HEAT reference atomization energies to CCSD(T) results from ccsd_energies.csv.

- Reads ccsd_energies.csv (produced by collect_ccsd_csv.py).
- Embeds the HEAT reference table and parses:
    * the reaction (e.g. 1*CO2>1*C+2*O)
    * E_ref_elec = EHF + dECCSD(T) + dECCSDT + dEHLC  [kJ/mol] → [Ha]
- For each *project* and molecule, computes our CCSD(T) atomization energy:
    E_form(project, mol) = sum_nu E(prod) - sum_nu E(react)  [Ha]
  using the same reaction definition as in the HEAT table.
- Computes ΔE = E_form_ours - E_form_HEAT  [Ha].
- Plots ΔE vs project for each molecule (one line per molecule) in Ha or mHa.
- Optionally writes heat_delta_vs_project.csv with the numerical values in Hartree.

Usage
-----
    python plot_heat_delta_vs_depth.py ccsd_energies.csv \
        --basis-short avdz \
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


# More precise if you like, but 2625.5 is fine for plotting
KJM_PER_HARTREE = 2625.5  # kJ/mol per Hartree


# ---------------------------------------------------------------------------
# 1) HEAT reference table (your pasted table)
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


@dataclass
class HeatRef:
    molecule: str         # e.g. CO2 (with CCH mapped to C2H)
    reaction: str         # e.g. "1*CO2>1*C+2*O"
    elec_au: float        # EHF + dECCSD(T) + dECCSDT + dEHLC [Hartree]

@dataclass
class FormRef:
    molecule: str     # e.g. C2H2 (with CCH mapped to C2H)
    form_au: float    # formation energy [Hartree]


def parse_formation_table(raw: str) -> Dict[str, FormRef]:
    """
    Parse FORMATION_TABLE text.

    Each line: MOL  E_form[kJ/mol]

    Convert to Hartree and map HEAT naming to CSV naming (CCH -> C2H).
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
        refs[mol_key] = FormRef(molecule=mol_key, form_au=val_kj /1000.0)
    return refs

import re

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


def parse_heat_table(raw: str) -> Dict[str, HeatRef]:
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

        # All numeric columns except the final REACTION string
        numeric_strs = parts[2:-1]
        if len(numeric_strs) != 9:
            raise ValueError(f"Unexpected numeric column count in line:\n{line}")

        nums = list(map(float, numeric_strs))
        # Indices:
        # 0: EHF
        # 1: dECCSD(T)
        # 2: dECCSDT
        # 3: dEHLC
        # 4: dErel
        # 5: dEZPE
        # 6: dEDBOC
        # 7: dESO
        # 8: E (total formation energy including thermal corrections)

        EHF = nums[0]
        dECCSDT_corr = nums[1]
        dECCSDT_triples = nums[2]
        dEHLC = nums[3]

        elec_kjmol = EHF + dECCSDT_corr + dECCSDT_triples + dEHLC
        elec_au = elec_kjmol / KJM_PER_HARTREE

        reaction = parts[-1]

        # Map HEAT naming to your CSV naming (CCH -> C2H)
        mol_key = "C2H" if mol == "CCH" else mol

        refs[mol_key] = HeatRef(molecule=mol_key, reaction=reaction, elec_au=elec_au)

    return refs

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

    Only projects that have the required reference species are used.
    Molecules with elements lacking μ for a given project are skipped.
    """
    result: Dict[Tuple[str, str], float] = {}

    projects = sorted({p for (p, _) in energies.keys()})

    # Precompute stoichiometries
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
            continue  # cannot form anything for this project

        for mol, fref in form_refs.items():
            # Need total energy of the molecule itself
            E_mol = energies.get((proj, mol))
            if E_mol is None:
                continue

            elems = stoich[mol]
            # Require μ for all elements in this molecule
            if any(elem not in mu for elem in elems):
                continue

            mu_tot = 0.0
            for elem, n_e in elems.items():
                mu_tot += n_e * mu[elem]

            E_form_au = E_mol - mu_tot
            result[(proj, mol)] = E_form_au

    return result


# ---------------------------------------------------------------------------
# 2) CCSD(T) energies from CSV
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 2) CCSD(T) energies from CSV (new format)
# ---------------------------------------------------------------------------

@dataclass
class CCRecord:
    project: str          # from 'projectname'
    core: str             # 'eCEPP', 'ccECP', 'AE', ...
    molecule: str
    basis_short: str
    method: str
    energy_au: float
    depth: float | None   # from 'ee_depth' (numeric) or None for 'free'


def read_ccsd_csv(path: Path) -> List[CCRecord]:
    """
    Read the new CSV format with columns like:

        projectname, ee_depth, en_depth, een_depth,
        molecule, core, basis_short, method, energy_au,
        spin_dependent_ee_cusp, spin_dependent_ee,
        spin_dependent_en, spin_dependent_een

    We currently ignore:
        en_depth, een_depth, and all spin_dependent_* columns.
    """
    records: List[CCRecord] = []
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            # ee_depth: "free" or numeric; treat "free" as None
            depth: float | None = None
            if "ee_depth" in row and row["ee_depth"] is not None:
                s = row["ee_depth"].strip()
                if s and s.lower() != "free":
                    depth = float(s)

            core = row.get("core", "unknown").strip() if row.get("core") is not None else "unknown"

            rec = CCRecord(
                project=row["projectname"],
                core=core,
                molecule=row["molecule"],
                basis_short=row["basis_short"],
                method=row["method"],
                energy_au=float(row["energy_au"]),
                depth=depth,
            )
            records.append(rec)
    return records


# ---------------------------------------------------------------------------
# 3) Parse reactions like "1*CO2>1*C+2*O"
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


# ---------------------------------------------------------------------------
# 4) Build index of our CCSD(T) energies
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 4) Build index of our CCSD(T) energies
# ---------------------------------------------------------------------------

def build_energy_index(
    records: List[CCRecord],
    basis_short: str | None = None,
) -> tuple[Dict[Tuple[str, str], float], Dict[str, float | None]]:
    """
    Build:
      energies[(project, molecule)] = energy_au
      project_depths[project] = depth (float or None)

    If basis_short is not None, filter to that basis.
    """
    energies: Dict[Tuple[str, str], float] = {}
    project_depths: Dict[str, float | None] = {}

    for rec in records:
        if basis_short is not None and rec.basis_short != basis_short:
            continue

        key = (rec.project, rec.molecule)
        energies[key] = rec.energy_au

        if rec.project not in project_depths:
            project_depths[rec.project] = rec.depth
        else:
            # If we previously had None but now see a real depth, upgrade it
            if project_depths[rec.project] is None and rec.depth is not None:
                project_depths[rec.project] = rec.depth

    return energies, project_depths

# ---------------------------------------------------------------------------
# 5) Compute our atomization energies (Hartree)
# ---------------------------------------------------------------------------

def compute_our_atomization(
    energies: Dict[Tuple[str, str], float],                         # <<< CHANGED
    heat_refs: Dict[str, HeatRef],
) -> Dict[Tuple[str, str], float]:                                  # <<< CHANGED
    """
    Return our atomization energies in Hartree for each (project, molecule),
    using the HEAT reaction definition:

        E_form(project, mol) = sum_nu E(prod) - sum_nu E(react)   [Ha]  # <<< CHANGED
    """
    result: Dict[Tuple[str, str], float] = {}                       # <<< CHANGED

    projects = sorted({p for (p, _) in energies.keys()})            # <<< CHANGED

    for mol, href in heat_refs.items():
        lhs, rhs = parse_reaction(href.reaction)

        for proj in projects:
            E_reac_au = 0.0
            E_prod_au = 0.0

            ok = True

            # Reactants
            for sp, coef in lhs.items():
                key = (proj, sp)
                if key not in energies:
                    ok = False
                    break
                E_reac_au += coef * energies[key]
            if not ok:
                continue

            # Products
            for sp, coef in rhs.items():
                key = (proj, sp)
                if key not in energies:
                    if key[1] == 'H':
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


# ---------------------------------------------------------------------------
# 6) ΔE vs project and plotting (Ha / mHa)
# ---------------------------------------------------------------------------
def plot_deltas(
    atom_E_by_basis: Dict[str, Dict[Tuple[str, str], float]],
    form_E_by_basis: Dict[str, Dict[Tuple[str, str], float]],
    project_depths: Dict[str, Dict[str, float | None]],
    heat_refs: Dict[str, HeatRef],
    form_refs: Dict[str, FormRef],
    out_csv: Path | None = None,
    unit: str = "mHa",
    save_path: Path | None = None,  # <--- NEW
) -> None:
    """
    atom_E_by_basis:
        basis_short -> {(project, molecule): atomization E_form_our [Ha]}
    form_E_by_basis:
        basis_short -> {(project, molecule): formation E_form_our [Ha]}
    project_depths:
        basis_short -> {project: depth (float or None)}

    heat_refs:
        molecule -> HeatRef (electronic atomization reference [Ha])
    form_refs:
        molecule -> FormRef (benchmark formation energies [Ha])

    unit: "Ha" or "mHa" for plotting.

    Layout:
        rows: 0 = no-H molecules, 1 = H-containing molecules, 2 = MAE/RMSE
        columns: 0..(B-1) = atomization per basis
                 B..(2B-1) = formation per basis

    x-axis:
        - depth in [-1, 1] → label depth
        - otherwise → label "free"
    """
    if unit == "Ha":
        scale = 1.0
        y_label = r"$\Delta E$ (ours − ref) [Ha]"
        stats_label = r"error metric [Ha]"
    elif unit == "mHa":
        scale = 1000.0
        y_label = r"$\Delta E$ (ours − ref) [mHa]"
        stats_label = r"error metric [mHa]"
    else:
        raise ValueError("unit must be 'Ha' or 'mHa'")

    # Union of all bases that have either atomization or formation data
    basis_list = sorted(set(atom_E_by_basis.keys()) | set(form_E_by_basis.keys()))
    if not basis_list:
        print("No data to plot.")
        return

    # Union of all molecules for a global colour map
    all_molecules = set(heat_refs.keys()) | set(form_refs.keys())
    molecules_all = sorted(all_molecules)

    import matplotlib.pyplot as plt

    if len(molecules_all) <= 10:
        cmap = plt.get_cmap("tab10")
    else:
        cmap = plt.get_cmap("tab20")
    denom = max(len(molecules_all) - 1, 1)
    mol_colors = {mol: cmap(i / denom) for i, mol in enumerate(molecules_all)}

    n_basis = len(basis_list)
    n_cols = 2 * n_basis  # atom block + formation block

    fig, axes = plt.subplots(
        3,
        n_cols,
        sharex="col",
        sharey="row",
        figsize=(4.0 * n_cols, 9.0),
        squeeze=False,
    )

    csv_rows: List[Dict[str, str]] = []

    def plot_block(
        kind: str,                      # "atom" or "form"
        basis: str,
        col_index: int,
        our_E: Dict[Tuple[str, str], float],
        refs: Dict[str, object],
        get_ref_val,
    ) -> None:
        """
        Plot one vertical block (3 rows) for a given basis and kind.
        """
        ax_noH = axes[0, col_index]
        ax_H = axes[1, col_index]
        ax_stats = axes[2, col_index]

        depths = project_depths.get(basis, {})
        if not our_E:
            ax_noH.text(0.5, 0.5, f"No {kind} data", ha="center", va="center")
            ax_H.text(0.5, 0.5, f"No {kind} data", ha="center", va="center")
            ax_stats.text(0.5, 0.5, f"No {kind} data", ha="center", va="center")
            for ax in (ax_noH, ax_H, ax_stats):
                ax.set_axis_off()
            return

        molecules = sorted(refs.keys())
        mol_no_H = [m for m in molecules if "H" not in m]
        mol_with_H = [m for m in molecules if "H" in m]

        projects = sorted(
            {p for (p, _) in our_E.keys()},
            key=lambda p: (
                1
                if (
                    depths.get(p) is None
                    or depths[p] < -1.0
                    or depths[p] > 1.0
                )
                else 0,
                depths.get(p) if depths.get(p) is not None else 0.0,
                p,
            ),
        )

        project_to_x = {p: i for i, p in enumerate(projects)}
        depth_labels: List[str] = []
        for proj in projects:
            d = depths.get(proj)
            if d is None or d < -1.0 or d > 1.0:
                depth_labels.append("free")
            else:
                depth_labels.append(f"{d:g}")

        deltas_by_project: Dict[str, List[float]] = {p: [] for p in projects}

        def plot_group(ax, mol_list, group_suffix: str) -> None:
            for mol in mol_list:
                href = refs[mol]
                E_ref_au = get_ref_val(href)

                xs: List[float] = []
                ys: List[float] = []

                for proj in projects:
                    key = (proj, mol)
                    if key not in our_E:
                        continue
                    E_our_au = our_E[key]
                    delta_au = E_our_au - E_ref_au

                    xs.append(project_to_x[proj])
                    ys.append(delta_au * scale)

                    d = depths.get(proj)
                    if d is None or d < -1.0 or d > 1.0:
                        depth_label = "free"
                        depth_str = ""
                    else:
                        depth_label = f"{d:g}"
                        depth_str = f"{d:.8f}"

                    csv_rows.append(
                        {
                            "basis_short": basis,
                            "kind": kind,
                            "molecule": mol,
                            "project": proj,
                            "depth": depth_str,
                            "depth_label": depth_label,
                            "E_our_Ha": f"{E_our_au:.10f}",
                            "E_HEAT_Ha": f"{E_ref_au:.10f}",  # ref energy (HEAT or formation)
                            "delta_Ha": f"{delta_au:.10f}",
                        }
                    )

                    deltas_by_project[proj].append(delta_au)

                if xs:
                    xs_arr = np.array(xs)
                    ys_arr = np.array(ys)
                    order = np.argsort(xs_arr)
                    xs_arr = xs_arr[order]
                    ys_arr = ys_arr[order]

                    color = mol_colors.get(mol, None)
                    ax.plot(
                        xs_arr,
                        ys_arr,
                        marker="o",
                        label=mol,
                        color=color,
                    )

            ax.axhline(0.0, linestyle="--")
            ax.set_ylabel(y_label)
            ax.set_title(f"{basis} {kind}: {group_suffix}")
            ax.legend(fontsize="x-small", ncol=3)
            ax.grid(True, linestyle=":", alpha=0.5)

        # Row 0: non-H molecules
        if mol_no_H:
            plot_group(ax_noH, mol_no_H, "ΔE (no H)")
        else:
            ax_noH.text(0.5, 0.5, "No H-free molecules", ha="center", va="center")
            ax_noH.set_axis_off()

        # Row 1: H-containing molecules
        if mol_with_H:
            plot_group(ax_H, mol_with_H, "ΔE (with H)")
        else:
            ax_H.text(0.5, 0.5, "No H-containing molecules", ha="center", va="center")
            ax_H.set_axis_off()

        # Row 2: MAE & RMSE over all molecules in this ref set
        xs_stats: List[int] = []
        mae_vals: List[float] = []
        rmse_vals: List[float] = []

        for proj in projects:
            deltas = deltas_by_project[proj]
            if not deltas:
                continue
            deltas_arr = np.array(deltas)

            mae = np.mean(np.abs(deltas_arr)) * scale
            rmse = np.sqrt(np.mean(deltas_arr**2)) * scale

            xs_stats.append(project_to_x[proj])
            mae_vals.append(mae)
            rmse_vals.append(rmse)

        if xs_stats:
            xs_stats_arr = np.array(xs_stats)
            order = np.argsort(xs_stats_arr)
            xs_stats_arr = xs_stats_arr[order]
            mae_vals_arr = np.array(mae_vals)[order]
            rmse_vals_arr = np.array(rmse_vals)[order]

            ax_stats.plot(xs_stats_arr, mae_vals_arr, marker="o", label=f"MAE ({unit})")
            ax_stats.plot(xs_stats_arr, rmse_vals_arr, marker="s", label=f"RMSE ({unit})")
            ax_stats.set_ylabel(stats_label)
            ax_stats.set_title(f"{basis} {kind}: MAE & RMSE")
            ax_stats.axhline(0.0, linestyle="--")
            ax_stats.legend(fontsize="small")
        else:
            ax_stats.text(0.5, 0.5, "No data for stats", ha="center", va="center")
            ax_stats.set_axis_off()

        ax_stats.set_xlabel("Jastrow depth")
        ax_stats.set_xticks(np.arange(len(projects)))
        ax_stats.set_xticklabels(depth_labels, rotation=45, ha="right")
        ax_stats.grid(True, linestyle=":", alpha=0.5)

        # Hide x tick labels on upper rows; they still share the same x-axis
        plt.setp(ax_noH.get_xticklabels(), visible=False)
        plt.setp(ax_H.get_xticklabels(), visible=False)

    # --- Loop over bases: left block = atomization, right block = formation ---
    for j, basis in enumerate(basis_list):
        # Atomization columns
        atom_E = atom_E_by_basis.get(basis, {})
        plot_block(
            kind="atom",
            basis=basis,
            col_index=j,
            our_E=atom_E,
            refs=heat_refs,
            get_ref_val=lambda href: href.elec_au,
        )

        # Formation columns
        form_E = form_E_by_basis.get(basis, {})
        plot_block(
            kind="form",
            basis=basis,
            col_index=j + n_basis,
            our_E=form_E,
            refs=form_refs,
            get_ref_val=lambda fref: fref.form_au,
        )

    fig.tight_layout()
    # Save if requested
    if save_path is not None:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300)
        print(f"Saved figure to {save_path}")

    # Still show interactively
    plt.show()

    if out_csv is not None:
        with out_csv.open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "basis_short",
                    "kind",         # "atom" or "form"
                    "molecule",
                    "project",
                    "depth",
                    "depth_label",
                    "E_our_Ha",
                    "E_HEAT_Ha",    # ref energy (HEAT electronic or formation benchmark)
                    "delta_Ha",
                ],
            )
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"Wrote delta CSV to {out_csv}")


# ---------------------------------------------------------------------------
# 7) CLI
# --------------------------------------------------------------------------- 
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare CCSD(T) atomization and formation energies vs references (Ha / mHa)."
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="If set, save each core figure into the figures/ directory.",
    )
    parser.add_argument(
        "ccsd_csv",
        type=Path,
        help="Path to heat_tc_energies.csv (from collect_ccsd_csv.py)",
    )
    parser.add_argument(
        "--basis-short",
        default=None,
        help=(
            "Filter by basis_short. "
            "Use a single value (e.g. avdz) or a comma-separated list "
            "(e.g. avdz,avtz). Each basis gets two columns: atom and form."
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
        help="Output CSV filename for ΔE values (Hartree). If omitted, no CSV is written.",
    )

    args = parser.parse_args()

    heat_refs = parse_heat_table(HEAT_TABLE)
    form_refs = parse_formation_table(FORMATION_TABLE)
    records = read_ccsd_csv(args.ccsd_csv)

    if not records:
        print("No records found in CSV.")
        return

    # Unique cores: e.g. ['AE', 'ccECP', 'eCEPP']
    cores = sorted({rec.core for rec in records})

    for core in cores:
        core_records = [rec for rec in records if rec.core == core]
        if not core_records:
            continue

        print(f"Processing core = {core} with {len(core_records)} records")

        # Determine which bases to use for *this* core
        if args.basis_short is None:
            basis_list = sorted({rec.basis_short for rec in core_records})
        else:
            basis_list = [b.strip() for b in args.basis_short.split(",") if b.strip()]

        atom_E_by_basis: Dict[str, Dict[Tuple[str, str], float]] = {}
        form_E_by_basis: Dict[str, Dict[Tuple[str, str], float]] = {}
        project_depths: Dict[str, Dict[str, float | None]] = {}

        for basis in basis_list:
            energies, depths = build_energy_index(core_records, basis_short=basis)
            if not energies:
                continue

            project_depths[basis] = depths
            atom_E_by_basis[basis] = compute_our_atomization(energies, heat_refs)
            form_E_by_basis[basis] = compute_our_formation(energies, form_refs)

        if not atom_E_by_basis and not form_E_by_basis:
            print(f"No matching data for core={core} and requested basis_short selection.")
            continue

        print(f"Plotting for core = {core}, bases = {', '.join(atom_E_by_basis.keys())}")

        save_path = None
        if args.save:
            # figures/depth_sweep_<core>.png
            safe_core = core.replace(" ", "_")
            save_dir = Path("figures")
            save_path = save_dir / f"depth_sweep_{safe_core}.png"

        plot_deltas(
            atom_E_by_basis=atom_E_by_basis,
            form_E_by_basis=form_E_by_basis,
            project_depths=project_depths,
            heat_refs=heat_refs,
            form_refs=form_refs,
            out_csv=args.out_csv,
            unit=args.unit,
            save_path=save_path,
        )

if __name__ == "__main__":
    main()
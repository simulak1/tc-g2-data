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
from heat_data import (
    KJM_PER_HARTREE,
    HEAT_TABLE,
    FORMATION_TABLE,
    HeatRef,
    FormRef,
    CCRecord,
    parse_heat_table,
    parse_formation_table,
    read_ccsd_csv,
    build_energy_index,
    compute_our_atomization,
    compute_our_formation,
)
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
    else:
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
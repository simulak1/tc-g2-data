#!/usr/bin/env python3
"""
Plot depth parameters for AE vs ccECP from a HEAT-style CSV.

- Filters rows with projectname in {"heat_ae_default", "heat_ccecp_default"}
  AND basis_short == "avdz".
- Parses:
    * ee_depth
    * e{atom}_depth  for atoms C, N, O, F, H
    * ee{atom}_depth for atoms C, N, O, F, H
- Produces a 2x2 subplot figure:
    (a) ee_depth per molecule (molecule on x-axis, avdz only)
    (b) e-n depths per atom (C, N, O, F, H on x-axis)
    (c) "een-depths" per atom, using ee{atom}_depth columns
    (d) sum_{atom}( e{atom}_depth + ee{atom}_depth ) per atom

AE and ccECP points are horizontally offset and use different markers.
"""

import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


PROJECTS = {
    "heat_ae_default": {
        "label": "AE",
        "offset": -0.15,
        "marker": "o",
    },
    "heat_ccecp_default": {
        "label": "ccECP",
        "offset": +0.15,
        "marker": "s",
    },
}

ELEMENTS = ["C", "N", "O", "F", "H"]


def load_and_filter(csv_file: str) -> pd.DataFrame:
    df = pd.read_csv(csv_file)

    # Keep only AE / ccECP + avdz
    df = df[
        (df["projectname"].isin(PROJECTS.keys())) &
        (df["basis_short"] == "avdz")
    ].copy()

    # With avdz-only, molecule label can just be the molecule name
    df["mol_label"] = df["molecule"].astype(str)
    return df


def panel_ee_depth(ax, df: pd.DataFrame):
    """Subplot (a): ee_depth vs molecule."""

    mol_labels = sorted(df["mol_label"].unique())
    x_base = np.arange(len(mol_labels))
    mol_index = {m: i for i, m in enumerate(mol_labels)}

    # To avoid duplicate legend entries
    seen_proj = set()

    for proj_name, info in PROJECTS.items():
        df_p = df[df["projectname"] == proj_name]
        x_vals = []
        y_vals = []

        for _, row in df_p.iterrows():
            ml = row["mol_label"]
            if ml not in mol_index:
                continue
            x = mol_index[ml] + info["offset"]
            y = row["ee_depth"]
            x_vals.append(x)
            y_vals.append(y)

        if not x_vals:
            continue

        label = info["label"] if proj_name not in seen_proj else None
        seen_proj.add(proj_name)

        ax.scatter(x_vals, y_vals, marker=info["marker"], label=label)

    ax.set_xticks(x_base)
    ax.set_xticklabels(mol_labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("ee_depth (a.u.)")
    ax.set_title("(a) e–e depths per system (avdz)")
    ax.legend(frameon=False)
    ax.grid(True, linestyle=":", alpha=0.4)


def _scatter_atom_depths(ax, df: pd.DataFrame, value_type: str, add_legend: bool = False):
    """
    Helper for panels (b)-(d).

    value_type:
      "e"   -> e{X}_depth
      "ee"  -> ee{X}_depth
      "sum" -> e{X}_depth + ee{X}_depth + ee_depth

    add_legend:
      If True, add a legend mapping colors to molecule names.
    """

    x_base = np.arange(len(ELEMENTS))

    # Color map for molecules (consistent across subplots)
    mol_labels = sorted(df["mol_label"].unique())
    cmap = plt.get_cmap("tab20")
    mol_color = {m: cmap(i % 20) for i, m in enumerate(mol_labels)}

    for proj_name, info in PROJECTS.items():
        df_p = df[df["projectname"] == proj_name]

        for _, row in df_p.iterrows():
            mol = row["molecule"]
            mol_label = row["mol_label"]

            for i, el in enumerate(ELEMENTS):
                col_e = f"e{el}_depth"
                col_ee = f"ee{el}_depth"
                col_2e = f"ee_depth"

                val_e = row[col_e] if col_e in df.columns else np.nan
                val_ee = row[col_ee] if col_ee in df.columns else np.nan
                val_2e = row[col_2e] if col_2e in df.columns else np.nan

                if value_type == "e":
                    val = val_e
                elif value_type == "ee":
                    val = val_ee
                elif value_type == "sum":
                    # Sum of contributions for this atom
                    val = val_2e if pd.notna(val_2e) else 0.0
                    
                    if pd.notna(val_e):
                        val += val_e
                    if pd.notna(val_ee):
                        val += val_ee
                    if pd.isna(val_e) and pd.isna(val_ee):
                        val = np.nan
                else:
                    raise ValueError(f"Unknown value_type {value_type!r}")

                if pd.isna(val):
                    continue

                x = x_base[i] + info["offset"]
                ax.scatter(
                    x,
                    val,
                    marker=info["marker"],
                    color=mol_color[mol_label],
                    alpha=0.8,
                )

    ax.set_xticks(x_base)
    ax.set_xticklabels(ELEMENTS)
    ax.grid(True, linestyle=":", alpha=0.4)

    if add_legend:
        # Legend: colors → molecules (marker shape still distinguishes AE/ccECP)
        from matplotlib.lines import Line2D

        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                color=mol_color[m],
                label=m,
                markersize=5,
            )
            for m in mol_labels
        ]

        # Put it outside to keep plot area clear
        ax.legend(
            handles=handles,
            title="Molecules",
            bbox_to_anchor=(1.02, 1.0),
            loc="upper left",
            borderaxespad=0.0,
            fontsize=6,
            title_fontsize=7,
        )


def panel_en_depth(ax, df: pd.DataFrame):
    """Subplot (b): e–n depths per atom (eX_depth)."""
    _scatter_atom_depths(ax, df, value_type="e", add_legend=True)
    ax.set_ylabel("e–n depth (a.u.)")
    ax.set_title("(b) e–n depths per atom (avdz)")


def panel_een_depth(ax, df: pd.DataFrame):
    """
    Subplot (c): 'een-depths' per atom.

    Here we use the eeX_depth columns as the atom-resolved 'een' contributions.
    """
    _scatter_atom_depths(ax, df, value_type="ee", add_legend=False)
    ax.set_ylabel("een-depth (a.u.)")
    ax.set_title("(c) een-depths per atom (avdz)")


def panel_sum_depth(ax, df: pd.DataFrame):
    """Subplot (d): sum depths per atom: eX_depth + eeX_depth."""
    _scatter_atom_depths(ax, df, value_type="sum", add_legend=False)
    ax.set_ylabel("sum depth (a.u.)")
    ax.set_title("(d) sum of depths per atom (avdz)")


def main():
    parser = argparse.ArgumentParser(
        description="Plot AE vs ccECP depth parameters from HEAT-style CSV."
    )
    parser.add_argument("csv_file", help="Input CSV file")
    parser.add_argument(
        "-o",
        "--output",
        default="none",
        help="Output image filename (default: show interactively)",
    )
    args = parser.parse_args()

    df = load_and_filter(args.csv_file)
    if df.empty:
        raise SystemExit(
            "No avdz rows with projectname in {heat_ae_default, heat_ccecp_default}."
        )

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_a, ax_b, ax_c, ax_d = axes.flatten()

    panel_ee_depth(ax_a, df)
    panel_en_depth(ax_b, df)
    panel_een_depth(ax_c, df)
    panel_sum_depth(ax_d, df)

    plt.tight_layout()
    if args.output != "none":
        fig.savefig(args.output, dpi=300, bbox_inches="tight")
        print(f"Saved figure to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()

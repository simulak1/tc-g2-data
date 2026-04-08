#!/usr/bin/env python3
"""
tools/analyze_cutoff_sweep.py

Goal: find best (en_cutoff, een_cutoff, ee_order, en_order, een_order) parameters
from data/energies_vs_order_cutoff.csv by analyzing error distributions.

Generates EXACTLY 5 figures:

  1) Heatmap: best achievable score vs (en_cutoff, een_cutoff) after optimizing orders.
  2) Bar: best achievable score vs order-tuple after optimizing cutoffs.
  3) Sensitivity lines: for best order-tuple, score vs en_cutoff (one line per een_cutoff).
  4) Histogram: error distribution overall vs best parameter combo.
  5) Pareto scatter: mean error vs worst-case error across species for every parameter combo,
     with Pareto front highlighted and best few annotated.

No seaborn; matplotlib only.

Usage
-----
python tools/analyze_cutoff_sweep.py data/energies_vs_order_cutoff.csv
python tools/analyze_cutoff_sweep.py data/energies_vs_order_cutoff.csv --outdir figures
python tools/analyze_cutoff_sweep.py data/energies_vs_order_cutoff.csv --metric median --topk 20
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


REQUIRED = [
    "en_cutoff", "een_cutoff",
    "ee_order", "en_order", "een_order",
    "species", "error_mHa",
]


def metric_fn(metric: str):
    metric = metric.lower()
    if metric == "mean":
        return np.mean
    if metric == "median":
        return np.median
    if metric == "rms":
        return lambda x: np.sqrt(np.mean(np.square(x)))
    raise SystemExit(f"Unknown --metric {metric}. Use mean|median|rms.")


def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in ["en_cutoff", "een_cutoff", "error_mHa", "ee_order", "en_order", "een_order"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["species"] = df["species"].astype(str)
    df = df.dropna(subset=["en_cutoff", "een_cutoff", "error_mHa", "ee_order", "en_order", "een_order"])
    return df



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", help="Input CSV (e.g. data/energies_vs_order_cutoff.csv)")
    ap.add_argument("--outdir", default=None, help="If given, save 5 figures into this dir; else show interactively.")
    ap.add_argument("--metric", default="mean", help="Score metric over errors: mean|median|rms (default: mean)")
    ap.add_argument("--topk", type=int, default=15, help="How many best combos to print / annotate (default: 15)")
    ap.add_argument("--worstq", type=float, default=1.0,
                    help="Robustness metric across species: 1.0=max, 0.9=90th percentile, etc. (default: 1.0)")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    df = coerce_numeric(df)

    score = metric_fn(args.metric)
    worstq = float(args.worstq)
    if not (0 < worstq <= 1.0):
        raise SystemExit("--worstq must be in (0,1].")

    outdir = Path(args.outdir) if args.outdir else None
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)

    # Parameter keys
    CUTOFFS = ["en_cutoff", "een_cutoff"]
    ORDERS = ["ee_order", "en_order", "een_order"]
    ALL = CUTOFFS + ORDERS

    # ------------------------------------------------------------------
    # Build an aggregated table for each full parameter combo:
    #   central = metric(errors over all rows/species)  [mean/median/rms]
    #   worst   = max (or quantile) of per-species metric (robustness)
    # ------------------------------------------------------------------
    def per_species_metric(g: pd.DataFrame) -> pd.Series:
        # g is one full param combo, multiple species
        # compute metric per species then take worst quantile
        by_sp = g.groupby("species")["error_mHa"].apply(lambda x: score(x.to_numpy()))
        if worstq == 1.0:
            worst = float(by_sp.max())
        else:
            worst = float(np.quantile(by_sp.to_numpy(), worstq))
        return pd.Series({
            "central": float(score(g["error_mHa"].to_numpy())),
            "worst": worst,
            "n_species": int(by_sp.shape[0]),
            "n_rows": int(g.shape[0]),
        })

    agg = df.groupby(ALL, dropna=False).apply(per_species_metric).reset_index()

    # Rank combos by central, then worst
    agg_sorted = agg.sort_values(["central", "worst"], ascending=[True, True]).reset_index(drop=True)

    # Pick "best" as the minimum central (tie-break worst)
    best = agg_sorted.iloc[0].to_dict()

    # Print best few
    topk = min(args.topk, len(agg_sorted))
    print("\nTop parameter combos (sorted by central, then worst):")
    cols_print = ALL + ["central", "worst", "n_species", "n_rows"]
    print(agg_sorted.loc[: topk - 1, cols_print].to_string(index=False))

    best_cut = (best["en_cutoff"], best["een_cutoff"])
    best_ord = (best["ee_order"], best["en_order"], best["een_order"])
    print("\nSelected BEST combo:")
    print(f"  en_cutoff={best_cut[0]}, een_cutoff={best_cut[1]}, "
          f"ee_order={best_ord[0]}, en_order={best_ord[1]}, een_order={best_ord[2]}")
    print(f"  central({args.metric})={best['central']:.3f} mHa, "
          f"worst(q={worstq})={best['worst']:.3f} mHa")

    # ==============================================================
    # FIGURE 1: Heatmap of best achievable central score vs cutoffs,
    #           after optimizing orders for each cutoff pair.
    # ==============================================================
    cut_opt = agg.groupby(CUTOFFS, dropna=False)[["central", "worst"]].min().reset_index()
    en_vals = sorted(cut_opt["en_cutoff"].unique())
    een_vals = sorted(cut_opt["een_cutoff"].unique())

    mat = np.full((len(een_vals), len(en_vals)), np.nan, dtype=float)
    lookup = {(r["en_cutoff"], r["een_cutoff"]): r["central"] for _, r in cut_opt.iterrows()}
    for iy, een in enumerate(een_vals):
        for ix, en in enumerate(en_vals):
            mat[iy, ix] = lookup.get((en, een), np.nan)

    fig1, ax1 = plt.subplots(figsize=(1.2 * max(6, len(en_vals)), 1.0 * max(5, len(een_vals))))
    im1 = ax1.imshow(np.ma.masked_invalid(mat), origin="lower", aspect="auto")
    ax1.set_xticks(np.arange(len(en_vals)))
    ax1.set_xticklabels([f"{x:g}" for x in en_vals], rotation=45, ha="right")
    ax1.set_yticks(np.arange(len(een_vals)))
    ax1.set_yticklabels([f"{y:g}" for y in een_vals])
    ax1.set_xlabel("en_cutoff")
    ax1.set_ylabel("een_cutoff")
    ax1.set_title(f"Best achievable {args.metric} error (mHa) vs cutoffs (orders optimized)")

    # annotate cells
    for iy in range(len(een_vals)):
        for ix in range(len(en_vals)):
            v = mat[iy, ix]
            if np.isfinite(v):
                ax1.text(ix, iy, f"{v:.1f}", ha="center", va="center", fontsize=8)

    # mark global best cutoff point
    best_en, best_een = best_cut
    if best_en in en_vals and best_een in een_vals:
        ix = en_vals.index(best_en)
        iy = een_vals.index(best_een)
        ax1.scatter([ix], [iy], marker="o")  # default color

    cbar1 = fig1.colorbar(im1, ax=ax1, shrink=0.9)
    cbar1.set_label("central error (mHa)")
    fig1.tight_layout()

    # ==============================================================
    # FIGURE 2: Bar plot — best achievable central score vs order tuple,
    #           after optimizing cutoffs for each order tuple.
    # ==============================================================
    ord_opt = agg.groupby(ORDERS, dropna=False)[["central", "worst"]].min().reset_index()
    ord_opt = ord_opt.sort_values(["central", "worst"], ascending=[True, True]).reset_index(drop=True)

    top_ord = ord_opt.head(min(20, len(ord_opt))).copy()
    labels = [f"({int(r.ee_order)},{int(r.en_order)},{int(r.een_order)})" for r in top_ord.itertuples()]
    vals = top_ord["central"].to_numpy()

    fig2, ax2 = plt.subplots(figsize=(max(10, 0.45 * len(labels)), 5.5))
    ax2.bar(np.arange(len(labels)), vals)
    ax2.set_xticks(np.arange(len(labels)))
    ax2.set_xticklabels(labels, rotation=45, ha="right")
    ax2.set_ylabel(f"best achievable {args.metric} error (mHa)")
    ax2.set_title(f"Top order-tuples by best achievable {args.metric} error")
    fig2.tight_layout()

    # ==============================================================
    # FIGURE 3: Sensitivity (best order fixed) — score vs en_cutoff,
    #           with one line per een_cutoff.
    # ==============================================================
    ee_o, en_o, een_o = best_ord
    sub_best_ord = agg[(agg["ee_order"] == ee_o) & (agg["en_order"] == en_o) & (agg["een_order"] == een_o)].copy()

    # Pivot: rows een_cutoff, cols en_cutoff -> central
    piv = sub_best_ord.pivot_table(index="een_cutoff", columns="en_cutoff", values="central", aggfunc="min")
    piv = piv.sort_index().sort_index(axis=1)

    fig3, ax3 = plt.subplots(figsize=(9.5, 6.0))
    x = piv.columns.to_numpy(dtype=float)
    for een in piv.index.to_numpy(dtype=float):
        y = piv.loc[een].to_numpy(dtype=float)
        ax3.plot(x, y, marker="o", label=f"een_cutoff={een:g}")
    ax3.set_xlabel("en_cutoff")
    ax3.set_ylabel(f"{args.metric} error (mHa)")
    ax3.set_title(f"Sensitivity at best orders (ee,en,een)=({int(ee_o)},{int(en_o)},{int(een_o)})")
    ax3.legend(fontsize=8, ncol=2)
    fig3.tight_layout()

    # ==============================================================
    # FIGURE 4: Histogram — overall error distribution vs best combo
    # ==============================================================
    overall = df["error_mHa"].to_numpy(dtype=float)
    best_rows = df[
        (df["en_cutoff"] == best_cut[0]) &
        (df["een_cutoff"] == best_cut[1]) &
        (df["ee_order"] == best_ord[0]) &
        (df["en_order"] == best_ord[1]) &
        (df["een_order"] == best_ord[2])
    ]["error_mHa"].to_numpy(dtype=float)

    fig4, ax4 = plt.subplots(figsize=(9.0, 5.5))
    bins = 30
    ax4.hist(overall, bins=bins, alpha=0.6, label="All samples")
    ax4.hist(best_rows, bins=bins, alpha=0.6, label="Best combo samples")
    ax4.set_xlabel("error (mHa)")
    ax4.set_ylabel("count")
    ax4.set_title("Error distribution: all vs best parameter combo")
    ax4.legend()
    fig4.tight_layout()
    # ==============================================================
    # FIGURE 5: central vs worst (top-K only, labeled)
    # ==============================================================
    K = 9  # number of best candidates to consider

    sub = agg_sorted.head(K).copy()
    pts = sub[["central", "worst"]].to_numpy(dtype=float)

    fig5, ax5 = plt.subplots(figsize=(8.5, 6.5))

    ax5.scatter(
        pts[:, 0],
        pts[:, 1],
        s=60,
        alpha=0.8,
        label=f"Top {K} parameter combos"
    )

    # Compact labels: index + orders + cutoffs
    for i, r in enumerate(sub.itertuples(index=False), start=1):
        label = (
            f"{i}: "
            f"en={r.en_cutoff:g},een={r.een_cutoff:g}\n"
            f"ord=({int(r.ee_order)},{int(r.en_order)},{int(r.een_order)})"
        )
        ax5.annotate(
            label,
            (r.central, r.worst),
            xytext=(5, 5),           # small offset in points
            textcoords="offset points",
            fontsize=8
        )

    ax5.set_xlabel(f"central {args.metric} error (mHa)")
    ax5.set_ylabel(f"worst per-species error (q={worstq}) (mHa)")
    ax5.set_title(f"Accuracy vs robustness (top {K} combos)")
    ax5.legend()

    fig5.tight_layout()


    # ------------------------------------------------------------------
    # Save or show: exactly 5 figures
    # ------------------------------------------------------------------
    if outdir:
        fig1.savefig(outdir / "01_heatmap_best_vs_cutoffs.png", dpi=200, bbox_inches="tight")
        fig2.savefig(outdir / "02_top_orders_best_possible.png", dpi=200, bbox_inches="tight")
        fig3.savefig(outdir / "03_sensitivity_best_orders.png", dpi=200, bbox_inches="tight")
        fig4.savefig(outdir / "04_hist_all_vs_best.png", dpi=200, bbox_inches="tight")
        fig5.savefig(outdir / "05_pareto_central_vs_worst.png", dpi=200, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()

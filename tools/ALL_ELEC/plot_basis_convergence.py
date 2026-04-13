#!/usr/bin/env python3
"""Plot MAE, RMSE, MaxE vs basis set size (D→T→Q→5) for ANO and CC families."""

import csv
import re
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Configuration ──────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# (csv_path_relative_to_SCRIPT_DIR, legend_label)
SOURCES = [
    ("F12/Diag/basis_error_stats.csv", "F12"),
    ("Mix/basis_error_stats.csv",      "STANDARD"),
    ("Ref-corr/basis_error_stats.csv", "Ref-corr"),
    ("SVD/basis_error_stats.csv",      "SVD"),
]

ZETA_ORDER = ["D", "T", "Q", "5"]
ZETA_LABELS = ["DZ", "TZ", "QZ", "5Z"]

METRICS = [
    ("MAE_kcal_mol",     "MAE (kcal/mol)"),
    ("RMSE_kcal_mol",    "RMSE (kcal/mol)"),
    ("abs_MaxE_kcal_mol", "Max|E| (kcal/mol)"),
]

# ── Helpers ────────────────────────────────────────────────────────────────
ZETA_RE = re.compile(r"pV\(([DTQ5])\+d\)Z", re.IGNORECASE)

def detect_family(basis_name):
    """Return 'ano' or 'cc' depending on the basis set family."""
    low = basis_name.lower()
    if "ano" in low:
        return "ano"
    return "cc"

def extract_zeta(basis_name):
    m = ZETA_RE.search(basis_name)
    if m:
        return m.group(1).upper()
    return None

def load_source(csv_path, label):
    """Return list of dicts: {family, zeta, label, MAE, RMSE, MaxE}."""
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            z = extract_zeta(r["basis"])
            if z is None:
                continue
            rows.append({
                "family": detect_family(r["basis"]),
                "zeta":   z,
                "label":  label,
                "MAE_kcal_mol":      float(r["MAE_kcal_mol"]),
                "RMSE_kcal_mol":     float(r["RMSE_kcal_mol"]),
                "abs_MaxE_kcal_mol": float(r["abs_MaxE_kcal_mol"]),
                "n_species":         int(r["n_species"]),
            })
    return rows

# ── Load data ──────────────────────────────────────────────────────────────
all_rows = []
for rel_path, legend in SOURCES:
    full = os.path.join(SCRIPT_DIR, rel_path)
    if not os.path.isfile(full):
        print(f"WARNING: {full} not found, skipping")
        continue
    all_rows.extend(load_source(full, legend))

# Group: family → label → {zeta: row}
from collections import defaultdict
grouped = defaultdict(lambda: defaultdict(dict))
for r in all_rows:
    grouped[r["family"]][r["label"]][r["zeta"]] = r

# ── Plot ───────────────────────────────────────────────────────────────────
STYLE_MAP = {
    "F12/Diag": dict(marker="s", linestyle="-"),
    "STANDARD": dict(marker="o", linestyle="-"),
    "Ref-corr": dict(marker="^", linestyle="--"),
    "SVD":      dict(marker="D", linestyle="-."),
}

def plot_family(family, title_suffix, filename):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False)
    fig.suptitle(f"Basis convergence — {title_suffix}", fontsize=14, y=1.02)

    data = grouped[family]
    if not data:
        print(f"No data for family '{family}', skipping plot")
        return

    x_pos = list(range(len(ZETA_ORDER)))

    for ax, (metric_key, metric_label) in zip(axes, METRICS):
        for label in [s[1] for s in SOURCES]:
            if label not in data:
                continue
            zeta_dict = data[label]
            xs, ys, ns = [], [], []
            for i, z in enumerate(ZETA_ORDER):
                if z in zeta_dict:
                    xs.append(i)
                    ys.append(zeta_dict[z][metric_key])
                    ns.append(zeta_dict[z]["n_species"])
            if not xs:
                continue
            style = STYLE_MAP.get(label, {})
            line, = ax.plot(xs, ys, label=label, markersize=7, linewidth=1.8, **style)
            # annotate n_species
            for xi, yi, ni in zip(xs, ys, ns):
                ax.annotate(f"n={ni}", (xi, yi), textcoords="offset points",
                            xytext=(4, 6), fontsize=7, color=line.get_color(), alpha=0.7)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(ZETA_LABELS)
        ax.set_xlabel("Basis set size")
        ax.set_ylabel(metric_label)
        ax.set_title(metric_label)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out = os.path.join(SCRIPT_DIR, filename)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved: {out}")
    plt.close(fig)

plot_family("cc",  "CC basis sets",  "basis_convergence_cc.png")
plot_family("ano", "ANO basis sets", "basis_convergence_ano.png")

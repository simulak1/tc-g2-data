#!/usr/bin/env python3
"""Plot pipeline waterfall: stacked step times per system.

Data sources (from svd_ecp_timings.csv):
  - 00 RHF + 02a+02b VMC: g2_largebasis_rhf, av5z (shared across all bases)
  - 00 projection + 03 + 04: g2_dfunc_projection_5z (per basis)
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from g2_structures import G2_GEOMETRIES

CSV_PATH = Path(__file__).resolve().parents[2] / "data" / "ECP" / "svd_ecp_timings.csv"

PIPELINE_STEPS = [
    ("00_rhf",  "#4C78A8", "00 RHF"),
    ("02a+02b", "#F58518", "02a+02b VMC opt"),
    ("00_proj", "#9D755D", "00 projection"),
    ("03",      "#E45756", "03 TCHInt"),
    ("04",      "#54A24B", "04 CCSD(T)"),
]

BASES = ["pvdz", "avdz", "pvtz", "avtz", "pvqz", "avqz"]


def load_csv():
    with open(str(CSV_PATH)) as f:
        return list(csv.DictReader(f))


def nelec_for(system):
    g = G2_GEOMETRIES.get(system)
    return g.nelec if g is not None else None


def parse_float(s):
    try:
        v = float(s)
        return v if v > 0 else None
    except (ValueError, TypeError):
        return None


def collect_data(rows):
    """Return (step_times, rhf_times, vmc_times)."""
    step_times = defaultdict(lambda: defaultdict(lambda: None))
    rhf_times = {}
    vmc_times = {}

    for row in rows:
        project = row["project"]
        system, basis = row["system"], row["basis"]

        if project == "g2_new_projection_5z":
            t00 = parse_float(row["00"])
            if t00 is not None:
                step_times[(system, basis)]["00_proj"] = t00
            t03 = parse_float(row["03_no_fci"])
            if t03 is not None:
                step_times[(system, basis)]["03"] = t03
            t04 = parse_float(row["04"])
            if t04 is not None:
                step_times[(system, basis)]["04"] = t04

        if project == "g2_largebasis_rhf" and basis == "av5z":
            t00 = parse_float(row["00"])
            if t00 is not None:
                rhf_times[system] = t00
            t02a = parse_float(row["02a"])
            t02b = parse_float(row["02b"])
            if t02a is not None and t02b is not None:
                vmc_times[system] = t02a + t02b

    return step_times, rhf_times, vmc_times


def build_systems(step_times, rhf_times, vmc_times, basis):
    """Build sorted list of (system, total, combined_dict) for a given basis."""
    systems = []
    for (sys, bas), steps in step_times.items():
        if bas != basis:
            continue
        combined = dict(steps)
        if sys in rhf_times:
            combined["00_rhf"] = rhf_times[sys]
        if sys in vmc_times:
            combined["02a+02b"] = vmc_times[sys]
        total = sum(v for v in combined.values() if v is not None)
        if total > 0 and sum(1 for v in combined.values() if v is not None) >= 2:
            systems.append((sys, total, combined))
    systems.sort(key=lambda x: x[1], reverse=True)
    return systems


def save_latex_table(step_times, rhf_times, vmc_times, output_path: Path):
    step_labels = {sk: lbl for sk, _, lbl in PIPELINE_STEPS}
    lines = []
    lines.append(r'\begin{table}[h]')
    lines.append(r'\centering')
    lines.append(
        r'\caption{Pipeline wall-time breakdown per basis set for the G2/55 set '
        r'(ECP calculations). RHF and VMC optimisation are performed once at '
        r'aug-cc-pV5Z and shared across all projection bases.}')
    lines.append(r'\label{tab:timing}')

    ncols = len(PIPELINE_STEPS)
    lines.append(r'\begin{tabular}{@{}l r ' + 'r ' * ncols + r'@{}}')
    lines.append(r'\toprule')
    header_steps = ' & '.join(
        r'\multicolumn{1}{c}{' + step_labels[sk].replace(' ', r'\ ') + r'}'
        for sk, _, _ in PIPELINE_STEPS)
    lines.append(r'Basis & Total (nh) & ' + header_steps + r' \\')
    lines.append(r' & & ' + ' & '.join([r'\multicolumn{1}{c}{(\%)}'] * ncols) + r' \\')
    lines.append(r'\midrule')

    for basis in BASES:
        systems = build_systems(step_times, rhf_times, vmc_times, basis)
        if not systems:
            continue
        step_totals = {sk: sum(s[2].get(sk) or 0 for s in systems)
                       for sk, _, _ in PIPELINE_STEPS}
        grand_total = sum(step_totals.values())
        node_hours = grand_total / 3600
        pcts = []
        for sk, _, _ in PIPELINE_STEPS:
            pct = 100 * step_totals[sk] / grand_total if grand_total else 0
            pcts.append(f'{pct:.1f}')
        lines.append(f'{basis} & {node_hours:.1f} & ' + ' & '.join(pcts) + r' \\')

    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular}')
    lines.append(r'\end{table}')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text('\n'.join(lines) + '\n')
    print(f'Saved {output_path}')


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path,
                        default=Path(__file__).resolve().parents[2] / 'figures' / 'ECP')
    parser.add_argument('--html-output-dir', type=Path, default=None)
    parser.add_argument('--tables-dir', type=Path,
                        default=Path(__file__).resolve().parents[2] / 'tables' / 'ECP')
    args = parser.parse_args()

    rows = load_csv()
    step_times, rhf_times, vmc_times = collect_data(rows)

    outdir = args.output_dir
    html_outdir = args.html_output_dir or outdir
    outdir.mkdir(parents=True, exist_ok=True)
    html_outdir.mkdir(parents=True, exist_ok=True)

    nbases = len(BASES)
    titles = ["{} (RHF+VMC from av5z)".format(b) for b in BASES]

    fig = make_subplots(
        rows=nbases, cols=1,
        subplot_titles=titles,
        vertical_spacing=0.06,
    )

    for row_idx, basis in enumerate(BASES, start=1):
        systems = build_systems(step_times, rhf_times, vmc_times, basis)
        sys_names = [s[0] for s in systems]

        for step_key, color, label in PIPELINE_STEPS:
            is_av5z = step_key in ("00_rhf", "02a+02b")
            note = " (av5z)" if is_av5z else ""
            vals = [s[2].get(step_key) or 0 for s in systems]
            hovers = [
                "{}<br>{}{}={:.1f}s<br>total={:.1f}s<br>nelec={}".format(
                    s[0], label, note, v, s[1], nelec_for(s[0]))
                for s, v in zip(systems, vals)
            ]
            fig.add_trace(go.Bar(
                x=sys_names, y=vals,
                name=label + note,
                marker_color=color,
                hovertext=hovers, hoverinfo="text",
                legendgroup=step_key,
                showlegend=(row_idx == 1),
            ), row=row_idx, col=1)

        grand_total_s = sum(s[1] for s in systems)
        node_hours = grand_total_s / 3600
        fig.add_annotation(
            text="Total: {:.1f} node-hours".format(node_hours),
            xref="paper", yref="paper",
            x=1.0, y=1.0 - (row_idx - 1) / nbases + 0.005,
            xanchor="right", yanchor="bottom",
            showarrow=False,
            font=dict(size=12),
            bordercolor="#666", borderwidth=1, borderpad=4,
            bgcolor="rgba(255,255,255,0.85)",
        )

        fig.update_yaxes(title_text="Time (s)", row=row_idx, col=1)
        fig.update_xaxes(tickangle=60, tickfont_size=9, row=row_idx, col=1)

    max_systems = max(
        len(build_systems(step_times, rhf_times, vmc_times, b)) for b in BASES)

    fig.update_layout(
        barmode="stack",
        title_text="Pipeline wall time per system by basis (ECP)",
        height=500 * nbases,
        width=max(900, 22 * max_systems),
        legend=dict(orientation="h", yanchor="bottom", y=1.01),
    )

    html_path = html_outdir / "pipeline_waterfall.html"
    fig.write_html(str(html_path))
    if html_outdir != outdir:
        (outdir / "pipeline_waterfall.html").write_text(html_path.read_text())
    print("Saved:", html_path)

    save_latex_table(step_times, rhf_times, vmc_times,
                     args.tables_dir / 'timing_summary.tex')

    # Print per-basis step breakdown
    print()
    col_w = 18
    step_labels = {sk: lbl for sk, _, lbl in PIPELINE_STEPS}
    print("{:<8}  {}  {}".format(
        "Basis",
        "  ".join(f"{step_labels[sk]:>{col_w}}" for sk, _, _ in PIPELINE_STEPS),
        f"{'Total (node-h)':>16}",
    ))
    print("-" * (8 + 2 + (col_w + 2) * len(PIPELINE_STEPS) + 18))
    for basis in BASES:
        systems = build_systems(step_times, rhf_times, vmc_times, basis)
        if not systems:
            continue
        step_totals = {sk: sum(s[2].get(sk) or 0 for s in systems) for sk, _, _ in PIPELINE_STEPS}
        grand_total = sum(step_totals.values())
        node_hours = grand_total / 3600
        parts = []
        for sk, _, _ in PIPELINE_STEPS:
            pct = 100 * step_totals[sk] / grand_total if grand_total else 0
            parts.append(f"{pct:>{col_w}.1f}%")
        print("{:<8}  {}  {:>16.2f}".format(basis, "  ".join(parts), node_hours))


if __name__ == "__main__":
    main()

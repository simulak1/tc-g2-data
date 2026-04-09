#!/usr/bin/env python3
"""
Plot MAE, RMSE, and MaxE as a function of basis set for each method family.

Reads summary_vs_shci_pbe.csv and summary_vs_experiment.csv produced by
atomization_plots.py and creates convergence figures.
"""

import argparse
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path


# Unified x-axis: interleave cc-pVxZ and aug-cc-pVxZ by cardinal number.
UNIFIED_BASIS_ORDER = ['pvdz', 'avdz', 'pvtz', 'avtz', 'pvqz', 'avqz']

# CCSD(T) uses 2z/3z/4z naming for the same cardinal number as pVxZ.
BASIS_REMAP = {'2z': 'pvdz', '3z': 'pvtz', '4z': 'pvqz'}

# Display labels
BASIS_LABELS = {
    'avdz': 'avDZ', 'avtz': 'avTZ', 'avqz': 'avQZ',
    'pvdz': 'pvDZ', 'pvtz': 'pvTZ', 'pvqz': 'pvQZ',
    'CBS': 'CBS',
}

# Method family styling
FAMILY_STYLES = {
    'xTC-CCSD(T)':      {'color': 'darkorange',  'symbol': 'diamond',       'dash': 'solid'},
    'xTC-DC-CCSDT':     {'color': 'firebrick',   'symbol': 'star',          'dash': 'solid'},
    'md-xTC-CCSD(T)':   {'color': 'goldenrod',   'symbol': 'square',        'dash': 'dash'},
    'md-xTC-DC-CCSDT':  {'color': 'darkred',     'symbol': 'star-square',   'dash': 'dash'},
    'xTC-FCIQMC':       {'color': 'purple',       'symbol': 'circle',        'dash': 'solid'},
    'md-xTC-FCIQMC':    {'color': 'mediumpurple', 'symbol': 'circle-open',   'dash': 'dash'},
    'CCSD(T)':          {'color': 'dodgerblue',   'symbol': 'triangle-down', 'dash': 'dot'},
    'SHCI':             {'color': 'forestgreen',  'symbol': 'cross',         'dash': 'solid'},
    'SHCI+PBE':         {'color': 'gold',         'symbol': 'cross',         'dash': 'dash'},
    'F12':              {'color': 'teal',         'symbol': 'hexagram',      'dash': 'solid'},
}


def _get_style(family):
    return FAMILY_STYLES.get(family, {'color': 'gray', 'symbol': 'circle', 'dash': 'solid'})


def plot_convergence(df: pd.DataFrame, reference_label: str, output_path: Path,
                     png_path: Path = None,
                     metrics=('MAE', 'RMSE', 'MaxE'), qz_corr: bool = False):
    """
    Plot convergence of error metrics across basis sets.

    One subplot per metric (MAE, RMSE, MaxE), with lines for each method family.
    aug-cc-pV*Z and cc-pV*Z families on the same x-axis (cardinal number).
    """
    # Filter by QZ correction flag
    df = df[df['qz_corr'] == qz_corr].copy()

    if len(df) == 0:
        print(f"No data for qz_corr={qz_corr}, skipping.")
        return None

    n_metrics = len(metrics)
    fig = make_subplots(
        rows=1, cols=n_metrics,
        subplot_titles=[m + f' vs {reference_label}' for m in metrics],
        shared_yaxes=False,
        horizontal_spacing=0.08,
    )

    # Remap CCSD(T)/{2,3,4}z onto the pVxZ positions
    df = df.copy()
    df['x_basis'] = df['basis'].replace(BASIS_REMAP)
    df = df[df['x_basis'].isin(UNIFIED_BASIS_ORDER)]

    order_map = {b: i for i, b in enumerate(UNIFIED_BASIS_ORDER)}
    x_tick_labels = [BASIS_LABELS.get(b, b) for b in UNIFIED_BASIS_ORDER]

    families = df['family'].unique()

    for family in sorted(families):
        fam_df = df[df['family'] == family].copy()
        if len(fam_df) == 0:
            continue
        fam_df['_order'] = fam_df['x_basis'].map(order_map)
        fam_df = fam_df.sort_values('_order')

        style = _get_style(family)
        x_labels = [BASIS_LABELS.get(b, b) for b in fam_df['x_basis']]

        for col_idx, metric in enumerate(metrics):
            fig.add_trace(
                go.Scatter(
                    x=x_labels,
                    y=fam_df[metric],
                    mode='lines+markers',
                    name=family,
                    legendgroup=family,
                    showlegend=(col_idx == 0),
                    marker=dict(
                        color=style['color'],
                        symbol=style['symbol'],
                        size=10,
                        line=dict(width=1, color='black'),
                    ),
                    line=dict(
                        color=style['color'],
                        dash=style['dash'],
                        width=2,
                    ),
                    hovertemplate=f'{family}<br>{metric}=%{{y:.3f}}<extra></extra>',
                ),
                row=1, col=col_idx + 1,
            )

    # Add horizontal reference lines at chemical accuracy (1 kcal/mol)
    for col_idx in range(n_metrics):
        fig.add_hline(y=1.0, line_dash='dash', line_color='gray', line_width=1,
                      row=1, col=col_idx + 1,
                      annotation_text='1 kcal/mol' if col_idx == 0 else None,
                      annotation_position='top left')

    qz_tag = ' (QZ-HF corr)' if qz_corr else ''
    fig.update_layout(
        title=f'Basis Set Convergence of Atomization Energy Errors vs {reference_label}{qz_tag}',
        height=500,
        width=400 * n_metrics,
        legend=dict(orientation='h', yanchor='bottom', y=1.05, xanchor='center', x=0.5),
        hovermode='x unified',
    )

    for col_idx, metric in enumerate(metrics):
        fig.update_yaxes(title_text=f'{metric} (kcal/mol)', row=1, col=col_idx + 1)
        fig.update_xaxes(
            title_text='Basis set',
            categoryorder='array',
            categoryarray=x_tick_labels,
            row=1, col=col_idx + 1,
        )

    if output_path:
        fig.write_html(output_path)
        png_out = png_path or output_path.with_suffix('.png')
        fig.write_image(png_out, width=400 * n_metrics, height=500, scale=2)
        print(f"Saved {output_path} / {png_out}")

    return fig


def main():
    parser = argparse.ArgumentParser(description='Plot basis set convergence of atomization energy errors')
    parser.add_argument('--data-dir', type=Path,
                        default=Path(__file__).parent.parent / 'data',
                        help='Directory containing summary CSV files')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help='Output directory for PNG plots (default: figures/)')
    parser.add_argument('--html-output-dir', type=Path, default=None,
                        help='Output directory for HTML plots (default: same as --output-dir)')
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display plots interactively')
    args = parser.parse_args()

    data_dir = args.data_dir
    output_dir = args.output_dir or Path(__file__).parent.parent / 'figures'
    html_output_dir = args.html_output_dir or output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    html_output_dir.mkdir(parents=True, exist_ok=True)

    figs = []

    # Plot vs SHCI+PBE+CV
    pbe_csv = data_dir / 'summary_vs_shci_pbe.csv'
    if pbe_csv.exists():
        df_pbe = pd.read_csv(pbe_csv)
        fig = plot_convergence(df_pbe, 'SHCI+PBE+CV',
                               html_output_dir / 'convergence_vs_shci_pbe.html',
                               png_path=output_dir / 'convergence_vs_shci_pbe.png')
        if fig:
            figs.append(fig)
        fig_qz = plot_convergence(df_pbe, 'SHCI+PBE+CV',
                                  html_output_dir / 'convergence_vs_shci_pbe_QZ_corr.html',
                                  png_path=output_dir / 'convergence_vs_shci_pbe_QZ_corr.png',
                                  qz_corr=True)
        if fig_qz:
            figs.append(fig_qz)
    else:
        print(f"Not found: {pbe_csv}")

    if not args.no_show:
        for fig in figs:
            fig.show()


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Plotting library for atomization energy comparison data.

Reads atomization_comparison.csv and generates various plots:
- Bar chart comparison of atomization energies
- Error distribution histograms
- Overlaid error distributions
- Summary statistics
"""

import argparse
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path


def load_data(csv_path: Path) -> pd.DataFrame:
    """Load atomization comparison data from CSV."""
    return pd.read_csv(csv_path)


# xTC method slugs → (display label, marker symbol, marker size)
XTC_METHOD_DEFS = [
    ('xtc_ccsdt', 'xTC-CCSD(T)', 'diamond', 10),
    ('md_xtc_ccsdt', 'md-xTC-CCSD(T)', 'square', 10),
    ('xtc_dc_ccsdt', 'xTC-DC-CCSDT', 'star', 12),
    ('md_xtc_dc_ccsdt', 'md-xTC-DC-CCSDT', 'star-square', 10),
    ('xtc_fciqmc', 'xTC-FCIQMC', 'circle', 10),
    ('md_xtc_fciqmc', 'md-xTC-FCIQMC', 'circle-open', 12),
]

# Colors for xTC bar charts (one per method, in order of XTC_METHOD_DEFS)
_XTC_BAR_COLORS = ['darkorange', 'orangered', 'tomato', 'firebrick', 'darkred', 'maroon']


def compute_best_xtc_de(df: pd.DataFrame, basis: str = 'avtz') -> pd.Series:
    """
    For each molecule, pick the best available xTC De value.

    Hierarchy (highest priority first):
      md-xTC-FCIQMC > xTC-FCIQMC > md-xTC-DC-CCSDT > xTC-DC-CCSDT >
      md-xTC-CCSD(T) > xTC-CCSD(T)
    """
    hierarchy = [
        f'De_md_xtc_fciqmc_{basis}',
        f'De_xtc_fciqmc_{basis}',
        f'De_md_xtc_dc_ccsdt_{basis}',
        f'De_xtc_dc_ccsdt_{basis}',
        f'De_md_xtc_ccsdt_{basis}',
        f'De_xtc_ccsdt_{basis}',
    ]
    best = pd.Series(np.nan, index=df.index)
    for col in reversed(hierarchy):  # lowest priority first, overwritten by higher
        if col in df.columns:
            mask = df[col].notna()
            best[mask] = df.loc[mask, col]
    return best


def compute_best_xtc_d0(df: pd.DataFrame, basis: str = 'avtz') -> pd.Series:
    """
    For each molecule, pick the best available xTC D0 (corrected) value.

    Same hierarchy as compute_best_xtc_de but uses D0_ columns.
    """
    hierarchy = [
        f'D0_md_xtc_fciqmc_{basis}',
        f'D0_xtc_fciqmc_{basis}',
        f'D0_md_xtc_dc_ccsdt_{basis}',
        f'D0_xtc_dc_ccsdt_{basis}',
        f'D0_md_xtc_ccsdt_{basis}',
        f'D0_xtc_ccsdt_{basis}',
    ]
    best = pd.Series(np.nan, index=df.index)
    for col in reversed(hierarchy):  # lowest priority first, overwritten by higher
        if col in df.columns:
            mask = df[col].notna()
            best[mask] = df.loc[mask, col]
    return best


def _xtc_col_suffix(basis: str, qz_corr: bool) -> str:
    """Build the column suffix for xTC columns, e.g. 'avtz' or 'avdz_QZ_corr'."""
    return f'{basis}_QZ_corr' if qz_corr else basis


def plot_atomization_comparison(merged_df: pd.DataFrame, output_path: Path = None,
                                png_path: Path = None,
                                basis: str = 'avtz', qz_corr: bool = False,
                                show_f12: bool = False):
    """
    Plot atomization energies: raw De, corrected D0, and experimental values.
    """
    merged_df = merged_df.sort_values('experiment').reset_index(drop=True)
    sfx = _xtc_col_suffix(basis, qz_corr)
    title_tag = f' (xTC basis: {basis}, QZ-HF corr)' if qz_corr else f' (xTC basis: {basis})'

    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=(f'G2 Set Atomization Energies: SHCI (CBS) vs Experiment{title_tag}',
                       'Deviation from Experiment'),
        row_heights=[0.65, 0.35],
        vertical_spacing=0.12
    )

    # Top panel: Absolute energies
    fig.add_trace(
        go.Bar(name='D0_SHCI (ref)', x=merged_df['formula'], y=merged_df['D0_shci_ref'],
               marker_color='forestgreen', opacity=0.8),
        row=1, col=1
    )
    fig.add_trace(
        go.Bar(name='Experiment', x=merged_df['formula'], y=merged_df['experiment'],
               marker_color='coral', opacity=0.8),
        row=1, col=1
    )

    # Add F12 data if available
    if show_f12 and 'D0_f12' in merged_df.columns:
        fig.add_trace(
            go.Bar(name='D0_F12', x=merged_df['formula'], y=merged_df['D0_f12'],
                   marker_color='purple', opacity=0.8),
            row=1, col=1
        )
    if 'D0_shci_pbe' in merged_df.columns:
        fig.add_trace(
            go.Bar(name='D0_SHCI+PBE', x=merged_df['formula'], y=merged_df['D0_shci_pbe'],
                   marker_color='gold', opacity=0.8),
            row=1, col=1
        )

    # Add xTC method data if available
    for (slug, label, _, _), color in zip(XTC_METHOD_DEFS, _XTC_BAR_COLORS):
        d0_col = f'D0_{slug}_{sfx}'
        if d0_col in merged_df.columns and merged_df[d0_col].notna().any():
            fig.add_trace(
                go.Bar(name=f'D0 {label}/{sfx}', x=merged_df['formula'], y=merged_df[d0_col],
                       marker_color=color, opacity=0.8),
                row=1, col=1
            )

    # Bottom panel: Deviations from experiment
    deviation_corr = merged_df['D0_shci_ref'] - merged_df['experiment']
    mad_corr = np.abs(deviation_corr).mean()

    fig.add_trace(
        go.Bar(name=f'D0_SHCI - Exp (MAD={mad_corr:.2f})', x=merged_df['formula'], y=deviation_corr,
               marker_color='forestgreen', opacity=0.8, showlegend=False),
        row=2, col=1
    )

    if show_f12 and 'D0_f12' in merged_df.columns:
        mask_f12 = merged_df['D0_f12'].notna()
        deviation_f12 = merged_df.loc[mask_f12, 'D0_f12'] - merged_df.loc[mask_f12, 'experiment']
        mad_f12 = np.abs(deviation_f12).mean()
        fig.add_trace(
            go.Bar(name=f'D0_F12 - Exp (MAD={mad_f12:.2f})',
                   x=merged_df.loc[mask_f12, 'formula'], y=deviation_f12,
                   marker_color='purple', opacity=0.8, showlegend=False),
            row=2, col=1
        )
    if 'D0_shci_pbe' in merged_df.columns:
        mask_pbe = merged_df['D0_shci_pbe'].notna()
        deviation_exp_f12 = merged_df.loc[mask_pbe, 'D0_shci_pbe'] - merged_df.loc[mask_pbe, 'experiment']
        mad_pbe = np.abs(deviation_exp_f12).mean()
        fig.add_trace(
            go.Bar(name=f'D0_SHCI+PBE - Exp (MAD={mad_pbe:.2f})',
                   x=merged_df.loc[mask_pbe, 'formula'], y=deviation_exp_f12,
                   marker_color='gold', opacity=0.8, showlegend=False),
            row=2, col=1
        )

    # Add xTC method deviations
    for (slug, label, _, _), color in zip(XTC_METHOD_DEFS, _XTC_BAR_COLORS):
        d0_col = f'D0_{slug}_{sfx}'
        if d0_col in merged_df.columns and merged_df[d0_col].notna().any():
            mask = merged_df[d0_col].notna()
            deviation = merged_df.loc[mask, d0_col] - merged_df.loc[mask, 'experiment']
            fig.add_trace(
                go.Bar(name=f'D0 {label}/{sfx} - Exp (n={mask.sum()})',
                       x=merged_df.loc[mask, 'formula'],
                       y=deviation,
                       marker_color=color, opacity=0.8, showlegend=False),
                row=2, col=1
            )

    # Aggregate D0 MAD annotations (vs Experiment)
    _ref_col = 'experiment'
    _mad_lines = []

    # Reference methods: SHCI CBS, CC-F12, SHCI+PBE
    _ref_methods = [('D0_shci_ref', 'SHCI CBS'), ('D0_shci_pbe', 'SHCI+PBE')]
    if show_f12:
        _ref_methods.insert(1, ('D0_f12', 'CC-F12'))
    for _col, _lbl in _ref_methods:
        if _col in merged_df.columns:
            _m = merged_df[_col].notna() & merged_df[_ref_col].notna()
            _v = np.abs(merged_df.loc[_m, _col] - merged_df.loc[_m, _ref_col]).mean() if _m.sum() else float('nan')
        else:
            _v = float('nan')
        _mad_lines.append(f'{_lbl} MAD={_v:.2f}')

    # xTC methods
    _ccsdt_d0 = f'D0_xtc_ccsdt_{sfx}'
    if _ccsdt_d0 in merged_df.columns:
        _m = merged_df[_ccsdt_d0].notna() & merged_df[_ref_col].notna()
        _mad1 = np.abs(merged_df.loc[_m, _ccsdt_d0] - merged_df.loc[_m, _ref_col]).mean() if _m.sum() else float('nan')
    else:
        _mad1 = float('nan')
    _best_d0 = compute_best_xtc_d0(merged_df, sfx)
    _mb = _best_d0.notna() & merged_df[_ref_col].notna()
    _mad2 = np.abs(_best_d0[_mb] - merged_df.loc[_mb, _ref_col]).mean() if _mb.sum() else float('nan')
    _mad_lines.append(f'xTC-CCSD(T)/{sfx} MAD={_mad1:.2f}')
    _mad_lines.append(f'xTC-best/{sfx} MAD={_mad2:.2f}')

    fig.add_annotation(
        text='<br>'.join(_mad_lines),
        xref='paper', yref='paper', x=0.98, y=0.02,
        showarrow=False, font=dict(size=11), align='right',
        bgcolor='rgba(255,255,255,0.85)', bordercolor='darkorange', borderwidth=1,
        xanchor='right', yanchor='bottom'
    )

    fig.add_hline(y=0, line_dash="solid", line_color="black", line_width=1, row=2, col=1)

    fig.update_layout(
        height=800,
        barmode='group',
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        hovermode='x unified'
    )

    fig.update_yaxes(title_text="Atomization Energy (kcal/mol)", row=1, col=1)
    fig.update_yaxes(title_text="Deviation (kcal/mol)", row=2, col=1)
    fig.update_xaxes(title_text="Molecule", row=2, col=1, tickangle=45)

    if output_path:
        fig.write_html(output_path)
        png_out = png_path or output_path.with_suffix('.png')
        fig.write_image(png_out, width=1600, height=800, scale=2)
        print(f"Saved plot to {output_path} / {png_out}")

    return fig


def plot_method_discrepancy_shci_pbe(merged_df: pd.DataFrame, output_path: Path = None,
                                     png_path: Path = None,
                                     basis: str = 'avtz', qz_corr: bool = False,
                                     show_f12: bool = False):
    """
    Plot discrepancies between raw atomization energies and SHCI+PBE as reference.

    Y-axis: Discrepancy from SHCI+PBE (De_shci_pbe = 0 reference line)
    X-axis: Molecule names
    Methods compared: SHCI CBS, F12 CBS, xTC-PP-CCSD(T)/<basis>, CCSDT in various basis sets
    """
    if 'De_shci_pbe' not in merged_df.columns:
        print("WARNING: De_shci_pbe column not found, skipping SHCI+PBE discrepancy plot.")
        return None

    sfx = _xtc_col_suffix(basis, qz_corr)
    title_tag = f' (xTC basis: {basis}, QZ-HF corr)' if qz_corr else f' (xTC basis: {basis})'

    # Restrict to molecules where SHCI+PBE data is available
    merged_df = merged_df[merged_df['De_shci_pbe'].notna()].copy()
    # Add CV correction to SHCI+PBE reference for fair comparison
    merged_df['De_shci_pbe_cv'] = merged_df['De_shci_pbe'] + merged_df['CV']
    merged_df = merged_df.sort_values('De_shci_pbe_cv').reset_index(drop=True)

    fig = go.Figure()

    # Reference line at y=0 (SHCI+PBE+CV)
    fig.add_hline(y=0, line_dash="solid", line_color="gold", line_width=2,
                  annotation_text="SHCI+PBE+CV (reference)", annotation_position="top left")

    methods = []

    # SHCI CBS (add CV correction to raw De for fair comparison with SHCI+PBE+CV)
    mask = merged_df['De_shci_ref'].notna()
    if mask.sum() > 0:
        discrepancy = (merged_df.loc[mask, 'De_shci_ref'] + merged_df.loc[mask, 'CV']) - merged_df.loc[mask, 'De_shci_pbe_cv']
        methods.append({
            'name': 'SHCI CBS+CV',
            'x': merged_df.loc[mask, 'formula'],
            'y': discrepancy,
            'color': 'forestgreen',
            'symbol': 'circle',
            'size': 10
        })

    # F12 CBS
    if show_f12 and 'De_f12' in merged_df.columns:
        mask = merged_df['De_f12'].notna()
        if mask.sum() > 0:
            discrepancy = (merged_df.loc[mask, 'De_f12'] + merged_df.loc[mask, 'CV']) - merged_df.loc[mask, 'De_shci_pbe_cv']
            methods.append({
                'name': 'F12 CC/CBS+CV',
                'x': merged_df.loc[mask, 'formula'],
                'y': discrepancy,
                'color': 'purple',
                'symbol': 'star',
                'size': 12
            })

    # xTC methods
    for slug, label, symbol, sz in XTC_METHOD_DEFS:
        de_col = f'De_{slug}_{sfx}'
        if de_col in merged_df.columns:
            mask = merged_df[de_col].notna()
            if mask.sum() > 0:
                discrepancy = merged_df.loc[mask, de_col] - merged_df.loc[mask, 'De_shci_pbe_cv']
                methods.append({
                    'name': f'{label}/{sfx}',
                    'x': merged_df.loc[mask, 'formula'],
                    'y': discrepancy,
                    'color': 'darkorange',
                    'symbol': symbol,
                    'size': sz,
                    'is_xtc': True,
                })

    # CCSD(T) in various basis sets (skip 2z)
    ccsdt_colors = {'3z': 'dodgerblue', '4z': 'blue', '5z': 'darkblue'}
    ccsdt_symbols = {'3z': 'triangle-down', '4z': 'square', '5z': 'pentagon'}
    for ccsdt_basis in ['3z', '4z', '5z']:
        de_col = f'De_ccsdt_{ccsdt_basis}'
        if de_col in merged_df.columns:
            mask = merged_df[de_col].notna()
            if mask.sum() > 0:
                discrepancy = (merged_df.loc[mask, de_col] + merged_df.loc[mask, 'CV']) - merged_df.loc[mask, 'De_shci_pbe_cv']
                methods.append({
                    'name': f'CCSD(T)+CV/{ccsdt_basis}',
                    'x': merged_df.loc[mask, 'formula'],
                    'y': discrepancy,
                    'color': ccsdt_colors[ccsdt_basis],
                    'symbol': ccsdt_symbols[ccsdt_basis],
                    'size': 8
                })

    # Add traces
    for method in methods:
        n_pts = len(method['y'])
        is_xtc = method.get('is_xtc', False)
        if is_xtc:
            label = f"{method['name']} (n={n_pts})"
        else:
            mad = np.abs(method['y']).mean()
            label = f"{method['name']} (MAD={mad:.2f})"
        fig.add_trace(go.Scatter(
            x=method['x'],
            y=method['y'],
            mode='markers',
            name=label,
            marker=dict(
                color=method['color'],
                symbol=method['symbol'],
                size=method['size'],
                line=dict(width=1, color='black')
            ),
            hovertemplate='%{x}<br>Δ = %{y:.2f} kcal/mol<extra></extra>'
        ))

    # Aggregate xTC MAD annotations (vs SHCI+PBE+CV)
    _ref = 'De_shci_pbe_cv'
    _ccsdt_col = f'De_xtc_ccsdt_{sfx}'
    if _ccsdt_col in merged_df.columns:
        _m = merged_df[_ccsdt_col].notna() & merged_df[_ref].notna()
        _mad1 = np.abs(merged_df.loc[_m, _ccsdt_col] - merged_df.loc[_m, _ref]).mean() if _m.sum() else float('nan')
    else:
        _mad1 = float('nan')
    _best = compute_best_xtc_de(merged_df, sfx)
    _mb = _best.notna() & merged_df[_ref].notna()
    _mad2 = np.abs(_best[_mb] - merged_df.loc[_mb, _ref]).mean() if _mb.sum() else float('nan')
    fig.add_annotation(
        text=f'xTC-CCSD(T)/{sfx} MAD={_mad1:.2f}<br>xTC-best/{sfx} MAD={_mad2:.2f}',
        xref='paper', yref='paper', x=0.98, y=0.02,
        showarrow=False, font=dict(size=11), align='right',
        bgcolor='rgba(255,255,255,0.85)', bordercolor='darkorange', borderwidth=1,
        xanchor='right', yanchor='bottom'
    )

    # Chemical accuracy band
    fig.add_hrect(y0=-1, y1=1, fillcolor="lightyellow", opacity=0.3,
                  line_width=0, annotation_text="±1 kcal/mol",
                  annotation_position="top right")

    fig.update_layout(
        title=f'Atomization Energy Discrepancy from SHCI+PBE+CV (Raw De){title_tag}',
        xaxis_title='Molecule',
        yaxis_title='Discrepancy (kcal/mol): Method - SHCI+PBE+CV',
        height=600,
        hovermode='x unified',
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        xaxis=dict(tickangle=45)
    )

    if output_path:
        fig.write_html(output_path)
        png_out = png_path or output_path.with_suffix('.png')
        fig.write_image(png_out, width=1600, height=600, scale=2)
        print(f"Saved plot to {output_path} / {png_out}")

    return fig


def print_statistics(merged_df: pd.DataFrame, show_f12: bool = False):
    """Print summary statistics."""
    deviation_raw = merged_df['De_shci_ref'] - merged_df['experiment']
    deviation_corr = merged_df['D0_shci_ref'] - merged_df['experiment']
    
    print("\n" + "="*60)
    print("ATOMIZATION ENERGY STATISTICS (kcal/mol)")
    print("="*60)
    
    print(f"\nNumber of molecules: {len(merged_df)}")
    
    print("\nRaw SHCI CBS (De_shci_ref) vs Experiment:")
    print(f"  Mean Absolute Deviation (MAD): {np.abs(deviation_raw).mean():.3f}")
    print(f"  Root Mean Square Deviation:    {np.sqrt((deviation_raw**2).mean()):.3f}")
    print(f"  Max Deviation:                 {deviation_raw.abs().max():.3f}")
    print(f"  Mean Signed Error:             {deviation_raw.mean():.3f}")
    
    print("\nCorrected SHCI CBS (D0_shci_ref) vs Experiment:")
    print(f"  Mean Absolute Deviation (MAD): {np.abs(deviation_corr).mean():.3f}")
    print(f"  Root Mean Square Deviation:    {np.sqrt((deviation_corr**2).mean()):.3f}")
    print(f"  Max Deviation:                 {deviation_corr.abs().max():.3f}")
    print(f"  Mean Signed Error:             {deviation_corr.mean():.3f}")
    
    # F12 statistics
    if show_f12 and 'D0_f12' in merged_df.columns:
        deviation_f12 = merged_df['D0_f12'].dropna() - merged_df.loc[merged_df['D0_f12'].notna(), 'experiment']
        n_f12 = len(deviation_f12)
        print(f"\nCorrected CC F12 vs Experiment ({n_f12} molecules):")
        print(f"  Mean Absolute Deviation (MAD): {np.abs(deviation_f12).mean():.3f}")
        print(f"  Root Mean Square Deviation:    {np.sqrt((deviation_f12**2).mean()):.3f}")
        print(f"  Max Deviation:                 {deviation_f12.abs().max():.3f}")
        print(f"  Mean Signed Error:             {deviation_f12.mean():.3f}")
    
    if 'D0_shci_pbe' in merged_df.columns:
        deviation_exp_f12 = merged_df['D0_shci_pbe'].dropna() - merged_df.loc[merged_df['D0_shci_pbe'].notna(), 'experiment']
        n_exp_f12 = len(deviation_exp_f12)
        print(f"\nCorrected SHCI+PBE vs Experiment ({n_exp_f12} molecules):")
        print(f"  Mean Absolute Deviation (MAD): {np.abs(deviation_exp_f12).mean():.3f}")
        print(f"  Root Mean Square Deviation:    {np.sqrt((deviation_exp_f12**2).mean()):.3f}")
        print(f"  Max Deviation:                 {deviation_exp_f12.abs().max():.3f}")
        print(f"  Mean Signed Error:             {deviation_exp_f12.mean():.3f}")
    
    # xTC method statistics (raw + QZ-corrected)
    all_sfx = ['avdz', 'avtz', 'avqz', 'pvdz', 'pvtz', 'pvqz', 'avdz_QZ_corr', 'avtz_QZ_corr', 'pvdz_QZ_corr', 'pvtz_QZ_corr']
    for sfx in all_sfx:
        for slug, label, _, _ in XTC_METHOD_DEFS:
            de_col = f'De_{slug}_{sfx}'
            d0_col = f'D0_{slug}_{sfx}'

            if de_col in merged_df.columns and merged_df[de_col].notna().any():
                mask_de = merged_df[de_col].notna()
                deviation_de = merged_df.loc[mask_de, de_col] - merged_df.loc[mask_de, 'experiment']
                n_de = len(deviation_de)
                print(f"\nRaw {label}/{sfx} vs Experiment ({n_de} molecules):")
                print(f"  Mean Absolute Deviation (MAD): {np.abs(deviation_de).mean():.3f}")
                print(f"  Root Mean Square Deviation:    {np.sqrt((deviation_de**2).mean()):.3f}")
                print(f"  Max Deviation:                 {deviation_de.abs().max():.3f}")
                print(f"  Mean Signed Error:             {deviation_de.mean():.3f}")

            if d0_col in merged_df.columns and merged_df[d0_col].notna().any():
                mask_d0 = merged_df[d0_col].notna()
                deviation_d0 = merged_df.loc[mask_d0, d0_col] - merged_df.loc[mask_d0, 'experiment']
                n_d0 = len(deviation_d0)
                print(f"\nCorrected D0 {label}/{sfx} vs Experiment ({n_d0} molecules):")
                print(f"  Mean Absolute Deviation (MAD): {np.abs(deviation_d0).mean():.3f}")
                print(f"  Root Mean Square Deviation:    {np.sqrt((deviation_d0**2).mean()):.3f}")
                print(f"  Max Deviation:                 {deviation_d0.abs().max():.3f}")
                print(f"  Mean Signed Error:             {deviation_d0.mean():.3f}")
    
    print("\n" + "="*60)


def print_summary_table(merged_df: pd.DataFrame, show_f12: bool = False):
    """Print a clean summary table of MAE, RMSE, MaxE for all methods."""
    print("\n" + "="*80)
    print("SUMMARY: Corrected D0 vs Experiment (kcal/mol)")
    print("="*80)
    print(f"{'Method':<35} {'n':>5} {'MAE':>10} {'RMSE':>10} {'MaxE':>10}")
    print("-"*80)
    
    methods = []
    
    # SHCI CBS
    dev = merged_df['D0_shci_ref'] - merged_df['experiment']
    methods.append(('SHCI CBS', len(dev), np.abs(dev).mean(), np.sqrt((dev**2).mean()), dev.abs().max()))
    
    # F12 CBS
    if show_f12 and 'D0_f12' in merged_df.columns:
        mask = merged_df['D0_f12'].notna()
        dev = merged_df.loc[mask, 'D0_f12'] - merged_df.loc[mask, 'experiment']
        methods.append(('F12 CBS', len(dev), np.abs(dev).mean(), np.sqrt((dev**2).mean()), dev.abs().max()))
    
    # SHCI+PBE CBS
    if 'D0_shci_pbe' in merged_df.columns:
        mask = merged_df['D0_shci_pbe'].notna()
        dev = merged_df.loc[mask, 'D0_shci_pbe'] - merged_df.loc[mask, 'experiment']
        methods.append(('SHCI+PBE CBS', len(dev), np.abs(dev).mean(), np.sqrt((dev**2).mean()), dev.abs().max()))
    
    # xTC methods (raw + QZ-corrected)
    for sfx in ['avdz', 'avtz', 'avqz', 'pvdz', 'pvtz', 'pvqz', 'avdz_QZ_corr', 'avtz_QZ_corr', 'pvdz_QZ_corr', 'pvtz_QZ_corr']:
        for slug, label, _, _ in XTC_METHOD_DEFS:
            d0_col = f'D0_{slug}_{sfx}'
            if d0_col in merged_df.columns and merged_df[d0_col].notna().any():
                mask = merged_df[d0_col].notna()
                dev = merged_df.loc[mask, d0_col] - merged_df.loc[mask, 'experiment']
                methods.append((f'{label}/{sfx}', len(dev), np.abs(dev).mean(), np.sqrt((dev**2).mean()), dev.abs().max()))
    
    for name, n, mae, rmse, max_e in methods:
        print(f"{name:<35} {n:>5} {mae:>10.3f} {rmse:>10.3f} {max_e:>10.3f}")

    print("="*80)


def print_detailed_table(merged_df: pd.DataFrame):
    """Print detailed results table."""
    print("\nDetailed Results:")
    print("-"*100)
    print(f"{'Molecule':<12} {'De(ref)':<10} {'D0(ref)':<10} {'D0(PBE)':<10} {'Exp':<10} {'D0-Exp':>8}")
    print("-"*100)
    for _, row in merged_df.iterrows():
        d0_pbe = row.get('D0_shci_pbe', float('nan'))
        print(f"{row['formula']:<12} {row['De_shci_ref']:<10.2f} {row['D0_shci_ref']:<10.2f} {d0_pbe:<10.2f} {row['experiment']:<10.2f} {row['D0_shci_ref'] - row['experiment']:>8.2f}")
    print("-"*100)


def _collect_shci_pbe_summary(merged_df: pd.DataFrame, show_f12: bool = False):
    """Collect summary statistics of De errors vs SHCI+PBE+CV reference.

    Returns list of dicts with keys: method, basis, family, qz_corr, n, MAE, RMSE, MaxE, MSE.
    """
    if 'De_shci_pbe' not in merged_df.columns:
        return []

    df = merged_df[merged_df['De_shci_pbe'].notna()].copy()
    ref = df['De_shci_pbe'] + df['CV']  # SHCI+PBE+CV

    rows = []

    def _add(name, family, basis, qz_corr, dev):
        rows.append({
            'method': name, 'family': family, 'basis': basis, 'qz_corr': qz_corr,
            'n': len(dev), 'MAE': np.abs(dev).mean(), 'RMSE': np.sqrt((dev**2).mean()),
            'MaxE': dev.abs().max(), 'MSE': dev.mean(),
        })

    # SHCI CBS + CV
    dev = (df['De_shci_ref'] + df['CV']) - ref
    _add('SHCI CBS+CV', 'SHCI', 'CBS', False, dev)

    # F12 CBS + CV
    if show_f12 and 'De_f12' in df.columns:
        mask = df['De_f12'].notna()
        dev = (df.loc[mask, 'De_f12'] + df.loc[mask, 'CV']) - ref[mask]
        _add('F12 CBS+CV', 'F12', 'CBS', False, dev)

    # xTC methods
    all_sfx = ['avdz', 'avtz', 'avqz', 'pvdz', 'pvtz', 'pvqz',
               'avdz_QZ_corr', 'avtz_QZ_corr', 'pvdz_QZ_corr', 'pvtz_QZ_corr']
    for sfx in all_sfx:
        qz = sfx.endswith('_QZ_corr')
        raw_basis = sfx.replace('_QZ_corr', '') if qz else sfx
        for slug, label, _, _ in XTC_METHOD_DEFS:
            de_col = f'De_{slug}_{sfx}'
            if de_col in df.columns and df[de_col].notna().any():
                mask = df[de_col].notna()
                dev = df.loc[mask, de_col] - ref[mask]
                _add(f'{label}/{sfx}', label, raw_basis, qz, dev)

    # CCSDT (Molpro) + CV
    for basis in ['2z', '3z', '4z', '5z']:
        de_col = f'De_ccsdt_{basis}'
        if de_col in df.columns and df[de_col].notna().any():
            mask = df[de_col].notna()
            dev = (df.loc[mask, de_col] + df.loc[mask, 'CV']) - ref[mask]
            _add(f'CCSD(T)+CV/{basis}', 'CCSD(T)', basis, False, dev)

    return rows


def _collect_experiment_summary(merged_df: pd.DataFrame, show_f12: bool = False):
    """Collect summary statistics of D0 errors vs Experiment.

    Returns list of dicts with keys: method, basis, family, qz_corr, n, MAE, RMSE, MaxE, MSE.
    """
    rows = []

    def _add(name, family, basis, qz_corr, dev):
        rows.append({
            'method': name, 'family': family, 'basis': basis, 'qz_corr': qz_corr,
            'n': len(dev), 'MAE': np.abs(dev).mean(), 'RMSE': np.sqrt((dev**2).mean()),
            'MaxE': dev.abs().max(), 'MSE': dev.mean(),
        })

    # SHCI CBS
    dev = merged_df['D0_shci_ref'] - merged_df['experiment']
    _add('SHCI CBS', 'SHCI', 'CBS', False, dev)

    # F12 CBS
    if show_f12 and 'D0_f12' in merged_df.columns:
        mask = merged_df['D0_f12'].notna()
        dev = merged_df.loc[mask, 'D0_f12'] - merged_df.loc[mask, 'experiment']
        _add('F12 CBS', 'F12', 'CBS', False, dev)

    # SHCI+PBE CBS
    if 'D0_shci_pbe' in merged_df.columns:
        mask = merged_df['D0_shci_pbe'].notna()
        dev = merged_df.loc[mask, 'D0_shci_pbe'] - merged_df.loc[mask, 'experiment']
        _add('SHCI+PBE CBS', 'SHCI+PBE', 'CBS', False, dev)

    # xTC methods
    for sfx in ['avdz', 'avtz', 'avqz', 'pvdz', 'pvtz', 'pvqz',
                'avdz_QZ_corr', 'avtz_QZ_corr', 'pvdz_QZ_corr', 'pvtz_QZ_corr']:
        qz = sfx.endswith('_QZ_corr')
        raw_basis = sfx.replace('_QZ_corr', '') if qz else sfx
        for slug, label, _, _ in XTC_METHOD_DEFS:
            d0_col = f'D0_{slug}_{sfx}'
            if d0_col in merged_df.columns and merged_df[d0_col].notna().any():
                mask = merged_df[d0_col].notna()
                dev = merged_df.loc[mask, d0_col] - merged_df.loc[mask, 'experiment']
                _add(f'{label}/{sfx}', label, raw_basis, qz, dev)

    # CCSDT (Molpro)
    for basis in ['2z', '3z', '4z', '5z']:
        d0_col = f'D0_ccsdt_{basis}'
        if d0_col in merged_df.columns and merged_df[d0_col].notna().any():
            mask = merged_df[d0_col].notna()
            dev = merged_df.loc[mask, d0_col] - merged_df.loc[mask, 'experiment']
            _add(f'CCSD(T)/{basis}', 'CCSD(T)', basis, False, dev)

    return rows


def save_summary_csv(merged_df: pd.DataFrame, output_dir: Path, show_f12: bool = False):
    """Save summary statistics to CSV files for downstream convergence plots."""
    rows_pbe = _collect_shci_pbe_summary(merged_df, show_f12)
    if rows_pbe:
        df_pbe = pd.DataFrame(rows_pbe)
        path = output_dir / 'summary_vs_shci_pbe.csv'
        df_pbe.to_csv(path, index=False)
        print(f"Saved {path}")


def print_shci_pbe_summary_table(merged_df: pd.DataFrame, show_f12: bool = False):
    """Print summary table of De errors vs SHCI+PBE+CV reference (MAE, RMSE, MaxE)."""
    rows = _collect_shci_pbe_summary(merged_df, show_f12)
    if not rows:
        return

    print("\n" + "="*80)
    print("SUMMARY: Raw De vs SHCI+PBE+CV (kcal/mol)")
    print("="*80)
    print(f"{'Method':<35} {'n':>5} {'MAE':>10} {'RMSE':>10} {'MaxE':>10} {'MSE':>10}")
    print("-"*80)

    for r in rows:
        print(f"{r['method']:<35} {r['n']:>5} {r['MAE']:>10.3f} {r['RMSE']:>10.3f} {r['MaxE']:>10.3f} {r['MSE']:>10.3f}")

    print("="*80)


def print_shci_pbe_detailed_table(merged_df: pd.DataFrame):
    """Print per-molecule De discrepancy from SHCI+PBE+CV for each xTC method/basis."""
    if 'De_shci_pbe' not in merged_df.columns:
        return

    df = merged_df[merged_df['De_shci_pbe'].notna()].copy()
    ref = df['De_shci_pbe'] + df['CV']  # SHCI+PBE+CV

    all_sfx = ['avdz', 'avtz', 'avqz', 'pvdz', 'pvtz', 'pvqz',
               'avdz_QZ_corr', 'avtz_QZ_corr', 'pvdz_QZ_corr', 'pvtz_QZ_corr']

    for sfx in all_sfx:
        # Collect available method columns for this basis suffix
        available = []
        for slug, label, _, _ in XTC_METHOD_DEFS:
            de_col = f'De_{slug}_{sfx}'
            if de_col in df.columns and df[de_col].notna().any():
                available.append((slug, label, de_col))

        if not available:
            continue

        # Build header
        method_labels = [label for _, label, _ in available]
        header_methods = ''.join(f'{lbl:>14}' for lbl in method_labels)
        disc_headers = ''.join(f'{"disc":>10}' for _ in method_labels)

        print(f"\n{'='*120}")
        print(f"Detailed De vs SHCI+PBE+CV — basis: {sfx}")
        print(f"{'='*120}")
        print(f"{'Molecule':<14}{'SHCI+PBE+CV':>14}{header_methods}{disc_headers}")
        print("-"*120)

        for idx, row in df.iterrows():
            line = f"{row['formula']:<14}{ref[idx]:>14.2f}"
            disc_line = ""
            for slug, label, de_col in available:
                val = row[de_col]
                if pd.notna(val):
                    line += f"{val:>14.2f}"
                    disc_line += f"{val - ref[idx]:>10.2f}"
                else:
                    line += f"{'---':>14}"
                    disc_line += f"{'---':>10}"
            print(line + disc_line)

        # Print MAE/RMSE/MaxE footer
        print("-"*120)
        footer_ref = f"{'MAE/RMSE/MaxE':<14}{' ':>14}"
        footer_disc = ""
        for slug, label, de_col in available:
            mask = df[de_col].notna()
            dev = df.loc[mask, de_col] - ref[mask]
            mae = np.abs(dev).mean()
            rmse = np.sqrt((dev**2).mean())
            maxe = dev.abs().max()
            footer_ref += f"{' ':>14}"
            footer_disc += f"{mae:>10.2f}"
        print(footer_ref + footer_disc)

        footer_ref2 = f"{' ':<14}{' ':>14}"
        footer_disc2 = ""
        for slug, label, de_col in available:
            mask = df[de_col].notna()
            dev = df.loc[mask, de_col] - ref[mask]
            rmse = np.sqrt((dev**2).mean())
            footer_ref2 += f"{' ':>14}"
            footer_disc2 += f"{rmse:>10.2f}"
        print(footer_ref2 + footer_disc2 + "  (RMSE)")

        footer_ref3 = f"{' ':<14}{' ':>14}"
        footer_disc3 = ""
        for slug, label, de_col in available:
            mask = df[de_col].notna()
            dev = df.loc[mask, de_col] - ref[mask]
            maxe = dev.abs().max()
            footer_ref3 += f"{' ':>14}"
            footer_disc3 += f"{maxe:>10.2f}"
        print(footer_ref3 + footer_disc3 + "  (MaxE)")

        print("="*120)


def main():
    parser = argparse.ArgumentParser(description='Plot atomization energy comparison data')
    parser.add_argument('csv_file', type=Path, nargs='?', 
                        default=Path(__file__).parent / 'atomization_comparison.csv',
                        help='Path to atomization_comparison.csv')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help='Output directory for PNG plots (default: same as CSV file)')
    parser.add_argument('--html-output-dir', type=Path, default=None,
                        help='Output directory for HTML plots (default: same as --output-dir)')
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display plots interactively')
    parser.add_argument('--stats-only', action='store_true',
                        help='Only print statistics, no plots')
    parser.add_argument('--basis', type=str, default='avtz',
                        choices=['avdz', 'avtz', 'avqz', 'pvdz', 'pvtz', 'pvqz'],
                        help='Basis set for xTC columns (default: avtz)')
    parser.add_argument('--f12', action='store_true',
                        help='Include F12 CBS data in plots and statistics (disabled by default)')
    args = parser.parse_args()

    csv_path = args.csv_file
    output_dir = args.output_dir or csv_path.parent
    html_output_dir = args.html_output_dir or output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    html_output_dir.mkdir(parents=True, exist_ok=True)
    basis = args.basis

    print(f"Loading data from {csv_path}...")
    print(f"xTC basis: {basis}")
    merged_df = load_data(csv_path)
    print(f"Loaded {len(merged_df)} molecules")

    # Print statistics (SHCI+PBE reference only)
    print_shci_pbe_summary_table(merged_df, show_f12=args.f12)
    print_shci_pbe_detailed_table(merged_df)
    save_summary_csv(merged_df, csv_path.parent, show_f12=args.f12)

    if args.stats_only:
        return

    # Generate plots for the raw basis data
    print("\nGenerating plots...")

    fig1 = plot_atomization_comparison(
        merged_df,
        output_path=html_output_dir / f'atomization_comparison_{basis}.html',
        png_path=output_dir / f'atomization_comparison_{basis}.png',
        basis=basis, show_f12=args.f12)
    fig5 = plot_method_discrepancy_shci_pbe(
        merged_df,
        output_path=html_output_dir / f'atomization_method_discrepancy_shci_pbe_{basis}.html',
        png_path=output_dir / f'atomization_method_discrepancy_shci_pbe_{basis}.png',
        basis=basis, show_f12=args.f12)

    figs = [fig1, fig5]

    # For avdz/avtz, also generate QZ-HF-corrected plots
    if basis in ('avdz', 'avtz'):
        fig1q = plot_atomization_comparison(
            merged_df,
            output_path=html_output_dir / f'atomization_comparison_{basis}_QZ_corr.html',
            png_path=output_dir / f'atomization_comparison_{basis}_QZ_corr.png',
            basis=basis, qz_corr=True, show_f12=args.f12)
        fig5q = plot_method_discrepancy_shci_pbe(
            merged_df,
            output_path=html_output_dir / f'atomization_method_discrepancy_shci_pbe_{basis}_QZ_corr.html',
            png_path=output_dir / f'atomization_method_discrepancy_shci_pbe_{basis}_QZ_corr.png',
            basis=basis, qz_corr=True, show_f12=args.f12)
        figs.extend([fig1q, fig5q])

    if not args.no_show:
        for fig in figs:
            if fig:
                fig.show()


if __name__ == '__main__':
    main()

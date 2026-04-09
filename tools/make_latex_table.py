#!/usr/bin/env python3
"""
Generate a booktabs LaTeX table of MAE/RMSE/MaxE vs SHCI+PBE+CV
from data/summary_vs_shci_pbe.csv.
"""

import argparse
import pandas as pd
from pathlib import Path

# Method families in preferred display order
FAMILY_ORDER = [
    'SHCI',
    'xTC-CCSD(T)',
    'xTC-DC-CCSDT',
    'md-xTC-CCSD(T)',
    'md-xTC-DC-CCSDT',
    'xTC-FCIQMC',
    'md-xTC-FCIQMC',
    'CCSD(T)',
]

# Unified basis order (interleaved pV/apV by cardinal number)
BASIS_ORDER = ['pvdz', 'avdz', 'pvtz', 'avtz', 'pvqz', 'avqz', 'CBS']

# CCSD(T) uses Xz notation; map to canonical pVxZ names for ordering
CCSDT_BASIS_REMAP = {'2z': 'pvdz', '3z': 'pvtz', '4z': 'pvqz', '5z': 'pv5z'}

# Display labels
BASIS_LABELS = {
    'pvdz': 'pvDZ', 'avdz': 'avDZ',
    'pvtz': 'pvTZ', 'avtz': 'avTZ',
    'pvqz': 'pvQZ', 'avqz': 'avQZ',
    'pv5z': 'pv5Z', 'CBS':  'CBS',
}


def _basis_sort_key(basis):
    full_order = BASIS_ORDER + ['pv5z']
    try:
        return full_order.index(basis)
    except ValueError:
        return 99


def build_table(csv_path: Path, output_path: Path, caption: str, label: str):
    df = pd.read_csv(csv_path)

    # Remap CCSD(T) basis names to canonical forms for sorting/display
    df['x_basis'] = df['basis'].replace(CCSDT_BASIS_REMAP)

    lines = []
    lines.append(r'\begin{table}[h]')
    lines.append(r'\centering')
    lines.append(rf'\caption{{{caption}}}')
    lines.append(rf'\label{{{label}}}')
    lines.append(r'\begin{tabular}{@{}llrrr@{}}')
    lines.append(r'\toprule')
    lines.append(r'Method & Basis & MAE & RMSE & MaxE \\')
    lines.append(r' & & \multicolumn{3}{c}{(kcal/mol)} \\')
    lines.append(r'\midrule')

    families_present = [f for f in FAMILY_ORDER if f in df['family'].values]
    for i, family in enumerate(families_present):
        fam_df = df[df['family'] == family].copy()
        fam_df['_order'] = fam_df['x_basis'].apply(_basis_sort_key)
        fam_df = fam_df.sort_values('_order')

        lines.append(rf'\multicolumn{{5}}{{l}}{{\textit{{{family}}}}} \\')
        for _, row in fam_df.iterrows():
            blabel = BASIS_LABELS.get(row['x_basis'], row['x_basis'])
            lines.append(
                rf' & {blabel} & {row["MAE"]:.2f} & {row["RMSE"]:.2f} & {row["MaxE"]:.2f} \\'
            )

        if i < len(families_present) - 1:
            lines.append(r'\midrule')

    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular}')
    lines.append(r'\end{table}')

    output_path.write_text('\n'.join(lines) + '\n')
    print(f'Saved {output_path}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=Path,
                        default=Path(__file__).parent.parent / 'data')
    parser.add_argument('--output-dir', type=Path,
                        default=Path(__file__).parent.parent / 'tables')
    parser.add_argument('--caption',
                        default='MAE, RMSE, and MaxE of atomization energies '
                                r'vs SHCI+PBE+CV reference (kcal/mol), '
                                'for each method and basis set.')
    parser.add_argument('--label', default='tab:errors_shci_pbe')
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv = args.data_dir / 'summary_vs_shci_pbe.csv'
    out = args.output_dir / 'shci_pbe_summary.tex'
    build_table(csv, out, args.caption, args.label)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Plot atomization-energy discrepancies against G2 Umrigar reference values.

Usage:
    python plot_atomization_energy_discrepancy.py [csv_file]

    Without arguments, the default atomization_energies.csv is used.
    The CSV may contain multiple basis-set datasets, distinguished by the basis_folder column.
"""

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go


SPECIES_TO_SHCI = {
    'li_h': 'LiH',
    'be_h': 'BeH',
    'c_h': 'CH',
    'c_h2_triplet': 'CH2(3B1)',
    'c_h2_singlet': 'CH2(1A1)',
    'c_h3': 'CH3',
    'c_h4': 'CH4',
    'n_h': 'NH',
    'n_h2': 'NH2',
    'n_h3': 'NH3',
    'o_h': 'OH',
    'h2_o': 'H2O',
    'h_f': 'HF',
    'si_h2_singlet': 'SiH2(1A1)',
    'si_h2_triplet': 'SiH2(3B1)',
    'si_h3': 'SiH3',
    'si_h4': 'SiH4',
    'p_h2': 'PH2',
    'p_h3': 'PH3',
    'h2_s': 'H2S',
    'h_cl': 'HCl',
    'li2': 'Li2',
    'li_f': 'LiF',
    'c2_h2': 'C2H2',
    'c2_h4': 'C2H4',
    'c2_h6': 'C2H6',
    'c_n': 'CN',
    'h_c_n': 'HCN',
    'c_o': 'CO',
    'h_c_o': 'HCO',
    'h2_c_o': 'H2CO',
    'h3_c_o_h': 'H3COH',
    'n2': 'N2',
    'n2_h4': 'N2H4',
    'n_o': 'NO',
    'o2': 'O2',
    'h2_o2': 'H2O2',
    'f2': 'F2',
    'c_o2': 'CO2',
    'na2': 'Na2',
    'si2': 'Si2',
    'p2': 'P2',
    's2': 'S2',
    'cl2': 'Cl2',
    'na_cl': 'NaCl',
    'si_o': 'SiO',
    'c_s': 'CS',
    's_o': 'SO',
    'cl_o': 'ClO',
    'cl_f': 'ClF',
    'si2_h6': 'Si2H6',
    'c_h3_cl': 'CH3Cl',
    'h3_c_s_h': 'H3CSH',
    'h_o_cl': 'HOCl',
    's_o2': 'SO2',
}


def _get_fallback_basis(basis_name):
    """For a +d basis, return candidate non-+d fallback basis names."""
    pattern = r'\((\w)\+d\)'
    if not re.search(pattern, basis_name, re.IGNORECASE):
        return []

    fallback = re.sub(pattern, r'\1', basis_name)
    candidates = []

    def add_candidate(name):
        if name not in candidates:
            candidates.append(name)
        lower = name.lower()
        upper = name.upper()
        if lower not in candidates:
            candidates.append(lower)
        if upper not in candidates:
            candidates.append(upper)

    add_candidate(fallback)
    if fallback.lower().startswith('nacv-'):
        add_candidate(fallback[5:])
    return candidates


def resolve_output_file(base_file, species_filter=None, exclude_list=None):
    """Resolve a final output path, optionally nesting by species/exclusion filters."""
    final_file = Path(base_file)
    if species_filter:
        species_tag = '_'.join(species_filter)
        output_dir = final_file.parent / f'{species_tag}_'
        output_dir.mkdir(parents=True, exist_ok=True)
        final_file = output_dir / final_file.name
    if exclude_list:
        exclude_tag = 'exclude_' + '_'.join(exclude_list)
        output_dir = final_file.parent / exclude_tag
        output_dir.mkdir(parents=True, exist_ok=True)
        final_file = output_dir / final_file.name
    return final_file


def load_data(analysis_dir, csv_files=None, basis_filter=None, basis_order=None, umrigar_file=None):
    """Load one or more data files and group them by basis set."""
    datasets = []

    if not csv_files:
        csv_paths = [analysis_dir / 'atomization_energies.csv']
    else:
        csv_paths = []
        for csv_file in csv_files:
            csv_path = Path(csv_file)
            if not csv_path.is_absolute():
                csv_path = analysis_dir / csv_path
            csv_paths.append(csv_path)

    fallback_needed = set()
    if basis_filter:
        for basis_name in basis_filter:
            fallback_needed.update(_get_fallback_basis(basis_name))

    grouped = {}
    for csv_path in csv_paths:
        if not csv_path.exists():
            print(f'Warning: file does not exist: {csv_path}')
            continue

        with open(csv_path, 'r') as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                species = row.get('species')
                if species is None:
                    continue
                try:
                    ae = float(row['atomization_energy_kcal_mol'])
                except (KeyError, ValueError):
                    continue

                label = row.get('basis_folder') or csv_path.stem
                if basis_filter and label not in basis_filter and label not in fallback_needed:
                    continue
                grouped.setdefault(label, []).append({
                    'species': species,
                    'atomization_energy_kcal_mol': ae,
                })

    for basis_name in list(grouped.keys()):
        fallback_candidates = _get_fallback_basis(basis_name)
        if not fallback_candidates:
            continue
        fallback_key = None
        for candidate in fallback_candidates:
            if candidate in grouped:
                fallback_key = candidate
                break
        if fallback_key is None:
            continue

        existing_species = {row['species'] for row in grouped[basis_name]}
        filled = []
        for row in grouped[fallback_key]:
            if row['species'] not in existing_species:
                filled.append(row)
        if filled:
            species_names = [row['species'] for row in filled]
            print(f'  [{basis_name}] filled {len(filled)} molecules from {fallback_key}: {species_names}')
            grouped[basis_name].extend(filled)

    if basis_filter:
        for key in list(grouped.keys()):
            if key not in basis_filter:
                del grouped[key]

    keys = sorted(grouped.keys())
    if basis_order:
        ordered_keys = [key for key in basis_order if key in grouped]
        remaining_keys = [key for key in keys if key not in ordered_keys]
        keys = ordered_keys + remaining_keys

    for label in keys:
        datasets.append((label, grouped[label]))

    if umrigar_file is None:
        umrigar_file = analysis_dir / 'G2_Umrigar.csv'
    umrigar_data = {}
    with open(umrigar_file, 'r') as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            molecule = row.get('Molecule')
            if not molecule:
                continue
            umrigar_data[molecule] = {
                'expt': float(row['Expt']),
                'zpe': float(row['ZPE']),
                'sr_so': float(row['SR_SO']),
                'cv': float(row['CV']),
                'shci_d0': float(row['SHCI_D0']),
                'shci_pbe_d0': float(row['SHCI_PBE_D0']),
            }

    return datasets, umrigar_data


def _ref_key_from_theo_ref(theo_ref):
    if theo_ref == 'extrap':
        return 'shci_d0'
    if theo_ref == 'PBE':
        return 'shci_pbe_d0'
    return 'expt'


def _ref_label_from_theo_ref(theo_ref):
    if theo_ref == 'extrap':
        return 'SHCI_D0'
    if theo_ref == 'PBE':
        return 'SHCI_PBE_D0'
    return 'Expt'


def calculate_discrepancy(data, umrigar_dict, theo_ref):
    """Compute discrepancy = corrected method value - theoretical reference."""
    results = []
    ref_key = _ref_key_from_theo_ref(theo_ref)
    ref_label = _ref_label_from_theo_ref(theo_ref)

    for row in data:
        species = row['species']
        ae = row['atomization_energy_kcal_mol']
        if species not in SPECIES_TO_SHCI:
            continue

        shci_name = SPECIES_TO_SHCI[species]
        if shci_name not in umrigar_dict:
            continue

        ref = umrigar_dict[shci_name]
        correction = ref['zpe'] + ref['sr_so']
        corrected_ae = ae + correction
        ref_value = ref[ref_key]
        discrepancy = corrected_ae - ref_value
        results.append({
            'species': species,
            'shci_name': shci_name,
            'ae_raw': ae,
            'correction': correction,
            'ae_corrected': corrected_ae,
            'ref_label': ref_label,
            'ref_value': ref_value,
            'expt': ref['expt'],
            'discrepancy': discrepancy,
        })

    return results


def calculate_all_discrepancies(datasets, umrigar_dict, theo_ref):
    all_results = []
    for label, data in datasets:
        all_results.append((label, calculate_discrepancy(data, umrigar_dict, theo_ref)))
    return all_results


def _get_union_shci_names(all_data, umrigar_dict, theo_ref, species_filter=None, exclude_list=None):
    if len(all_data) == 0:
        return []

    union_names = set()
    for _, data in all_data:
        union_names.update({entry['shci_name'] for entry in data})

    if species_filter:
        union_names = {name for name in union_names if name in species_filter}

    if exclude_list:
        print(f'Excluded molecules: {exclude_list}')
        print(union_names)
        union_names = {name for name in union_names if name not in exclude_list}
        print(f'Remaining molecules: {union_names}')

    if not union_names:
        return []

    ref_key = _ref_key_from_theo_ref(theo_ref)
    return sorted(union_names, key=lambda name: umrigar_dict.get(name, {}).get(ref_key, float('inf')))


def build_shci_reference_series(shci_names, umrigar_dict, theo_ref):
    selected_results = []
    ref_key = _ref_key_from_theo_ref(theo_ref)
    ref_label = _ref_label_from_theo_ref(theo_ref)

    if theo_ref is None:
        shci_d0_results = []
        shci_pbe_d0_results = []
        for name in shci_names:
            if name not in umrigar_dict:
                continue
            ref = umrigar_dict[name]
            expt = ref['expt']
            shci_d0 = ref['shci_d0']
            shci_pbe_d0 = ref['shci_pbe_d0']

            shci_d0_results.append({
                'species': None,
                'shci_name': name,
                'ae_raw': shci_d0,
                'correction': 0.0,
                'ae_corrected': shci_d0,
                'ref_label': 'Expt',
                'ref_value': expt,
                'expt': expt,
                'discrepancy': shci_d0 - expt,
            })
            shci_pbe_d0_results.append({
                'species': None,
                'shci_name': name,
                'ae_raw': shci_pbe_d0,
                'correction': 0.0,
                'ae_corrected': shci_pbe_d0,
                'ref_label': 'Expt',
                'ref_value': expt,
                'expt': expt,
                'discrepancy': shci_pbe_d0 - expt,
            })

        return [('SHCI_D0', shci_d0_results), ('SHCI_PBE_D0', shci_pbe_d0_results)]

    for name in shci_names:
        if name not in umrigar_dict:
            continue
        ref = umrigar_dict[name]
        ref_value = ref[ref_key]
        selected_results.append({
            'species': None,
            'shci_name': name,
            'ae_raw': ref_value,
            'correction': 0.0,
            'ae_corrected': ref_value,
            'ref_label': ref_label,
            'ref_value': ref_value,
            'expt': ref['expt'],
            'discrepancy': 0.0,
        })

    return [(ref_label, selected_results)]


def plot_discrepancy_interactive(all_data, output_file, umrigar_dict, theo_ref, species_filter=None, exclude_list=None):
    colors = ['#8B008B', 'orange', '#00CED1', 'blue', '#4B0082', 'red', 'green', 'brown', 'pink', 'gray']
    symbols = ['star', 'diamond', 'triangle-down', 'square', 'hexagon', 'circle', 'cross', 'x', 'pentagon', 'octagon']

    if len(all_data) == 0:
        print('No data available for plotting')
        return []

    shci_names = _get_union_shci_names(all_data, umrigar_dict, theo_ref, species_filter, exclude_list)
    if not shci_names:
        print('No molecules available for plotting')
        return []

    all_data_to_plot = list(all_data)
    if theo_ref is None:
        all_data_to_plot.extend(build_shci_reference_series(shci_names, umrigar_dict, theo_ref))

    fig = go.Figure()
    fig.add_hrect(
        y0=-1, y1=1,
        fillcolor='green', opacity=0.2,
        layer='below', line_width=0,
        annotation_text='±1 kcal/mol', annotation_position='top right'
    )
    fig.add_hline(
        y=0, line_color='green', line_width=2,
        annotation_text='SHCI CBS (reference)', annotation_position='top left'
    )

    all_discrepancies = []
    for i, (label, data) in enumerate(all_data_to_plot):
        color = colors[i % len(colors)]
        symbol = symbols[i % len(symbols)]
        data_dict = {entry['shci_name']: entry for entry in data}

        x_vals = []
        y_vals = []
        hover_texts = []
        for name in shci_names:
            if name not in data_dict:
                continue
            entry = data_dict[name]
            x_vals.append(name)
            y_vals.append(entry['discrepancy'])
            hover_texts.append(
                f'<b>{name}</b><br>'
                f'Method AE(raw): {entry["ae_raw"]:.2f} kcal/mol<br>'
                f'Correction(ZPE+SR_SO): {entry["correction"]:.2f} kcal/mol<br>'
                f'Method AE(corrected): {entry["ae_corrected"]:.2f} kcal/mol<br>'
                f'{entry["ref_label"]}: {entry["ref_value"]:.2f} kcal/mol<br>'
                f'Discrepancy: {entry["discrepancy"]:.2f} kcal/mol'
            )

        all_discrepancies.extend(y_vals)
        mad = np.mean(np.abs(y_vals)) if y_vals else 0
        fig.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            mode='markers',
            name=f'{label} (MAD={mad:.2f})',
            marker=dict(symbol=symbol, size=12, color=color, line=dict(width=1, color='black')),
            hovertemplate='%{text}<extra></extra>',
            text=hover_texts,
        ))

    y_limit = 15
    if all_discrepancies:
        y_max = max(abs(min(all_discrepancies)), abs(max(all_discrepancies)))
        y_limit = max(15, y_max * 1.1)

    fig.update_layout(
        title=dict(
            text=f'Atomization Energy Discrepancy from {_ref_label_from_theo_ref(theo_ref)} (with method ZPE+SR_SO correction)',
            font=dict(size=16),
        ),
        xaxis=dict(
            title='Molecule',
            tickangle=45,
            tickfont=dict(size=10),
            categoryorder='array',
            categoryarray=shci_names,
        ),
        yaxis=dict(
            title=f'Discrepancy (kcal/mol): Method - {_ref_label_from_theo_ref(theo_ref)}',
            range=[-y_limit, y_limit],
            gridcolor='white',
            gridwidth=1,
        ),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='left', x=0, font=dict(size=11)),
        plot_bgcolor='#e8f0f8',
        paper_bgcolor='white',
        hovermode='closest',
        width=1400,
        height=700,
    )

    html_file = resolve_output_file(Path(output_file).with_suffix('.html'), species_filter, exclude_list)
    fig.write_html(str(html_file), include_plotlyjs=True)
    print(f'Interactive chart saved to: {html_file}')
    return fig


def plot_discrepancy(all_data, output_file, umrigar_dict, theo_ref, species_filter=None, exclude_list=None):
    colors = ['orange', 'blue', 'red', 'purple', 'green', 'brown', 'pink', 'gray', 'olive', 'cyan']
    markers = ['D', 'o', 's', '^', 'v', '<', '>', 'p', 'h', '*']

    if len(all_data) == 0:
        print('No data available for plotting')
        return []

    shci_names = _get_union_shci_names(all_data, umrigar_dict, theo_ref, species_filter, exclude_list)
    if not shci_names:
        print('No molecules available for plotting')
        return []

    all_data_to_plot = list(all_data)
    if theo_ref is None:
        all_data_to_plot.extend(build_shci_reference_series(shci_names, umrigar_dict, theo_ref))

    name_to_idx = {name: idx for idx, name in enumerate(shci_names)}
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.set_facecolor('#e8f0f8')
    ax.axhline(y=0, color='green', linewidth=2, label='SHCI CBS (reference)')
    ax.axhspan(-1, 1, alpha=0.3, color='green', label='±1 kcal/mol')

    all_discrepancies = []
    length = 3
    for i, (label, data) in enumerate(all_data_to_plot):
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]
        data_dict = {entry['shci_name']: entry for entry in data}

        x_vals = []
        y_vals = []
        for name in shci_names:
            if name not in data_dict:
                continue
            x_vals.append(name_to_idx[name] * length)
            y_vals.append(data_dict[name]['discrepancy'])

        all_discrepancies.extend(y_vals)
        mad = np.mean(np.abs(y_vals)) if y_vals else 0
        offset = (i - (len(all_data) - 1) / 2) * 0.15
        x_offset = [x + offset for x in x_vals]
        ax.scatter(
            x_offset,
            y_vals,
            marker=marker,
            s=80,
            c=color,
            edgecolors='black',
            linewidths=0.5,
            label=f'{label} (MAD={mad:.2f})',
            zorder=5 + i,
        )

    x_positions = np.arange(len(shci_names)) * length
    ax.set_xticks(x_positions)
    ax.set_xticklabels(shci_names, rotation=45, ha='right', fontsize=9)
    ax.set_xlabel('Molecule', fontsize=12)
    ax.set_ylabel(f'Discrepancy (kcal/mol): Method - {_ref_label_from_theo_ref(theo_ref)}', fontsize=12)
    ax.set_title(f'Atomization Energy Discrepancy from {_ref_label_from_theo_ref(theo_ref)} (with method ZPE+SR_SO correction)', fontsize=14)
    ax.legend(loc='upper left', fontsize=10)

    if all_discrepancies:
        y_max = max(abs(min(all_discrepancies)), abs(max(all_discrepancies)))
        y_limit = max(15, y_max * 1.1)
        ax.set_ylim(-y_limit, y_limit)

    for xpos in x_positions:
        ax.axvline(x=xpos, color='gray', linestyle='--', linewidth=0.5, alpha=0.5, zorder=1)

    ax.grid(True, axis='y', alpha=0.3, linestyle='-', color='white')
    plt.tight_layout()

    png_file = resolve_output_file(output_file, species_filter, exclude_list)
    plt.savefig(png_file, dpi=150, bbox_inches='tight', facecolor='white')
    print(f'Image saved to: {png_file}')
    plt.show()
    return all_data


def print_statistics(all_data):
    print('\nStatistics:')
    for label, data in all_data:
        discrepancies = [entry['discrepancy'] for entry in data]
        mad = np.mean(np.abs(discrepancies))
        print(f'\n  [{label}]')
        print(f'    Molecule count: {len(data)}')
        print(f'    MAD (Mean Absolute Deviation): {mad:.2f} kcal/mol')
        print(f'    Maximum positive discrepancy: {max(discrepancies):.2f} kcal/mol')
        print(f'    Maximum negative discrepancy: {min(discrepancies):.2f} kcal/mol')
        print(f'    Mean signed discrepancy: {np.mean(discrepancies):.2f} kcal/mol')


def save_statistics_csv(all_data, output_dir):
    rows = []
    for label, data in all_data:
        discrepancies = np.array([entry['discrepancy'] for entry in data], dtype=float)
        if discrepancies.size == 0:
            print(f'  [WARN] {label}: no discrepancy data available, skipping CSV write')
            continue
        rows.append({
            'basis': label,
            'n_species': int(discrepancies.size),
            'MAE_kcal_mol': float(np.mean(np.abs(discrepancies))),
            'RMSE_kcal_mol': float(np.sqrt(np.mean(discrepancies ** 2))),
            'abs_MaxE_kcal_mol': float(np.max(np.abs(discrepancies))),
        })

    stats_csv = output_dir / 'basis_error_stats.csv'
    with open(stats_csv, 'w', newline='') as handle:
        fieldnames = ['basis', 'n_species', 'MAE_kcal_mol', 'RMSE_kcal_mol', 'abs_MaxE_kcal_mol']
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f'Statistics summary saved to: {stats_csv}')
    return stats_csv


def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot atomization-energy discrepancies against G2 Umrigar Expt reference values',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
    python plot_atomization_energy_discrepancy.py
    python plot_atomization_energy_discrepancy.py xTC_CCSDpT_aug-ano-pvdz_AE.csv
    python plot_atomization_energy_discrepancy.py xTC_CCSDpT_aug-ano-pvdz_AE.csv xTC_CCSDpT_aug-ano-pvtz_AE.csv
        ''',
    )
    parser.add_argument('csv_files', nargs='*', help='CSV files containing atomization-energy data')
    parser.add_argument('-o', '--output', default='output', help='Relative output subdirectory under the PNG/HTML/CSV output roots')
    parser.add_argument('--png-output-root', type=Path, default=None, help='Root directory for PNG output (default: <repo>/figures/ALL_ELEC)')
    parser.add_argument('--html-output-root', type=Path, default=None, help='Root directory for HTML output (default: <repo>/docs/figures/ALL_ELEC)')
    parser.add_argument('--csv-output-root', type=Path, default=None, help='Root directory for CSV output (default: <repo>/docs/figures/ALL_ELEC)')
    parser.add_argument('--no-html', action='store_true', help='Do not generate the interactive HTML chart (Plotly)')
    parser.add_argument('--species', default=None, help='Comma-separated list of molecule names to plot')
    parser.add_argument('--exclude', default=None, help='Comma-separated list of molecule names to exclude')
    parser.add_argument('--basis-folder', nargs='*', default=None, help='Plot only the specified basis_folder names')
    parser.add_argument('--basis-order', nargs='*', default=None, help='Display basis sets in the specified order')
    parser.add_argument('--no-png', action='store_true', help='Do not generate the static PNG image')
    parser.add_argument('--theo-ref', choices=['extrap', 'PBE'], default=None, help='Theoretical reference type')
    parser.add_argument('--exclude-elements', default=None, help='Exclude molecules containing the specified elements, as comma-separated symbols')
    parser.add_argument('--stats', action='store_true', help='Write MAE, RMSE, and abs(MaxE) to CSV')
    return parser.parse_args()


def main():
    args = parse_args()

    def parse_csv_list(value):
        if not value:
            return None
        return [item.strip() for item in value.split(',') if item.strip()]

    def parse_element_list(value):
        if not value:
            return None
        return {item.strip() for item in value.split(',') if item.strip()}

    script_dir = Path(__file__).resolve().parent
    analysis_dir = script_dir
    repo_root = script_dir.parents[1]

    output_subdir = Path(args.output)
    png_output_dir = (args.png_output_root or (repo_root / 'figures' / 'ALL_ELEC')) / output_subdir
    html_output_dir = (args.html_output_root or (repo_root / 'docs' / 'figures' / 'ALL_ELEC')) / output_subdir
    csv_output_dir = (args.csv_output_root or (repo_root / 'docs' / 'figures' / 'ALL_ELEC')) / output_subdir
    png_output_dir.mkdir(parents=True, exist_ok=True)
    html_output_dir.mkdir(parents=True, exist_ok=True)
    csv_output_dir.mkdir(parents=True, exist_ok=True)

    png_output_file = png_output_dir / 'atomization_energy_discrepancy.png'
    html_output_file = html_output_dir / 'atomization_energy_discrepancy.html'

    print('Loading data...')
    umrigar_file = repo_root / 'reference_data' / 'G2_Umrigar.csv'
    datasets, umrigar_dict = load_data(
        analysis_dir,
        args.csv_files,
        args.basis_folder,
        args.basis_order,
        umrigar_file=umrigar_file,
    )
    if len(datasets) == 0:
        print('Error: no data was loaded')
        return

    for label, data in datasets:
        print(f'  {label}: {len(data)} molecules')
    print(f'G2 Umrigar reference data: {len(umrigar_dict)} molecules')

    print('\nComputing discrepancies...')
    exclude_elements = parse_element_list(args.exclude_elements)
    if exclude_elements:
        atom_names = {
            'h': 'H', 'c': 'C', 'n': 'N', 'o': 'O', 'f': 'F',
            'li': 'Li', 'be': 'Be', 'na': 'Na', 'si': 'Si',
            'p': 'P', 's': 'S', 'cl': 'Cl',
        }

        def get_elements_from_species(species):
            name = species.replace('_singlet', '').replace('_triplet', '')
            elements = set()
            for part in name.split('_'):
                token = part.rstrip('0123456789')
                if token in atom_names:
                    elements.add(atom_names[token])
            return elements

        filtered_datasets = []
        for label, data in datasets:
            filtered = []
            for row in data:
                elems = get_elements_from_species(row['species'])
                if elems & exclude_elements:
                    continue
                filtered.append(row)
            removed = len(data) - len(filtered)
            if removed:
                print(f'  [{label}] excluded {removed} molecules containing elements {sorted(exclude_elements)}')
            filtered_datasets.append((label, filtered))
        datasets = filtered_datasets

    all_data = calculate_all_discrepancies(datasets, umrigar_dict, args.theo_ref)
    for label, data in all_data:
        print(f'  {label}: matched {len(data)} molecules')

    species_filter = parse_csv_list(args.species)
    exclude_list = parse_csv_list(args.exclude)

    print('\nGenerating plots...')
    if not args.no_html:
        plot_discrepancy_interactive(all_data, html_output_file, umrigar_dict, args.theo_ref, species_filter, exclude_list)
    if not args.no_png:
        plot_discrepancy(all_data, png_output_file, umrigar_dict, args.theo_ref, species_filter, exclude_list)

    shci_names = _get_union_shci_names(all_data, umrigar_dict, args.theo_ref, species_filter, exclude_list)
    all_stats_data = list(all_data)
    all_stats_data.extend(build_shci_reference_series(shci_names, umrigar_dict, args.theo_ref))
    print_statistics(all_stats_data)

    if args.stats:
        save_statistics_csv(all_data, csv_output_dir)

    for label, data in all_data:
        output_csv = csv_output_dir / f'{label}_discrepancy.csv'
        with open(output_csv, 'w', newline='') as handle:
            fieldnames = ['species', 'shci_name', 'ae_raw', 'correction', 'ae_corrected', 'ref_label', 'ref_value', 'expt', 'discrepancy']
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(sorted(data, key=lambda entry: entry['ref_value']))
        print(f'Discrepancy data saved to: {output_csv}')

    for label, data in build_shci_reference_series(shci_names, umrigar_dict, args.theo_ref):
        output_csv = csv_output_dir / f'{label}_discrepancy.csv'
        with open(output_csv, 'w', newline='') as handle:
            fieldnames = ['species', 'shci_name', 'ae_raw', 'correction', 'ae_corrected', 'ref_label', 'ref_value', 'expt', 'discrepancy']
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(sorted(data, key=lambda entry: entry['ref_value']))
        print(f'Discrepancy data saved to: {output_csv}')


if __name__ == '__main__':
    main()

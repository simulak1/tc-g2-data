#!/usr/bin/env python3
"""
Build merged atomization energy data from multiple sources:
- SHCI CBS extrapolated values from combined_reference_data.csv
- F12 reference data from g2_reference.py
- CCSD(T) energies from ccsd_energies.csv
- Corrections (ZPE, SR_SO, CV) from g2_atomization_reference.csv

Output: atomization_comparison.csv with all atomization energies and corrections.
"""

import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

# Add reference_data directory to path for f12_reference import
sys.path.insert(0, str(Path(__file__).parent.parent / 'reference_data'))
from f12_reference import G2_97_ATOMIZATION

# Conversion factor: Hartree to kcal/mol
HA_TO_KCAL = 627.5094740631

# Conversion: kJ/mol to kcal/mol
KJ_TO_KCAL = 1.0 / 4.184

# Mapping from formula to constituent atoms
# Format: {formula: {atom: count}}
FORMULA_TO_ATOMS = {
    'LiH': {'Li': 1, 'H': 1},
    'BeH': {'Be': 1, 'H': 1},
    'CH': {'C': 1, 'H': 1},
    'CH2': {'C': 1, 'H': 2},
    'CH2_singlet': {'C': 1, 'H': 2},
    'CH3': {'C': 1, 'H': 3},
    'CH4': {'C': 1, 'H': 4},
    'NH': {'N': 1, 'H': 1},
    'NH2': {'N': 1, 'H': 2},
    'NH3': {'N': 1, 'H': 3},
    'OH': {'O': 1, 'H': 1},
    'H2O': {'O': 1, 'H': 2},
    'HF': {'H': 1, 'F': 1},
    'SiH2': {'Si': 1, 'H': 2},
    'SiH2_triplet': {'Si': 1, 'H': 2},
    'SiH3': {'Si': 1, 'H': 3},
    'SiH4': {'Si': 1, 'H': 4},
    'PH2': {'P': 1, 'H': 2},
    'PH3': {'P': 1, 'H': 3},
    'H2S': {'S': 1, 'H': 2},
    'HCl': {'H': 1, 'Cl': 1},
    'Li2': {'Li': 2},
    'LiF': {'Li': 1, 'F': 1},
    'C2H2': {'C': 2, 'H': 2},
    'C2H4': {'C': 2, 'H': 4},
    'C2H6': {'C': 2, 'H': 6},
    'CN': {'C': 1, 'N': 1},
    'HCN': {'H': 1, 'C': 1, 'N': 1},
    'CO': {'C': 1, 'O': 1},
    'HCO': {'H': 1, 'C': 1, 'O': 1},
    'H2CO': {'H': 2, 'C': 1, 'O': 1},
    'CH3OH': {'C': 1, 'H': 4, 'O': 1},
    'N2': {'N': 2},
    'N2H4': {'N': 2, 'H': 4},
    'NO': {'N': 1, 'O': 1},
    'O2': {'O': 2},
    'H2O2': {'H': 2, 'O': 2},
    'F2': {'F': 2},
    'CO2': {'C': 1, 'O': 2},
    'Na2': {'Na': 2},
    'Si2': {'Si': 2},
    'P2': {'P': 2},
    'S2': {'S': 2},
    'Cl2': {'Cl': 2},
    'NaCl': {'Na': 1, 'Cl': 1},
    'SiO': {'Si': 1, 'O': 1},
    'CS': {'C': 1, 'S': 1},
    'SO': {'S': 1, 'O': 1},
    'ClO': {'Cl': 1, 'O': 1},
    'ClF': {'Cl': 1, 'F': 1},
    'Si2H6': {'Si': 2, 'H': 6},
    'CH3Cl': {'C': 1, 'H': 3, 'Cl': 1},
    'CH3SH': {'C': 1, 'H': 4, 'S': 1},
    'HOCl': {'H': 1, 'O': 1, 'Cl': 1},
    'SO2': {'S': 1, 'O': 2},
}

# Mapping from g2_atomization_reference molecule names to formula names
REFERENCE_TO_FORMULA = {
    'LiH': 'LiH', 'BeH': 'BeH', 'CH': 'CH', 'CH2': 'CH2',
    'CH2_singlet': 'CH2_singlet', 'CH3': 'CH3', 'CH4': 'CH4',
    'NH': 'NH', 'NH2': 'NH2', 'NH3': 'NH3', 'HO': 'OH',
    'H2O': 'H2O', 'HF': 'HF', 'SiH2': 'SiH2',
    'SiH2_triplet': 'SiH2_triplet', 'SiH3': 'SiH3', 'SiH4': 'SiH4',
    'PH2': 'PH2', 'PH3': 'PH3', 'H2S': 'H2S', 'HCl': 'HCl',
    'Li2': 'Li2', 'LiF': 'LiF', 'C2H2': 'C2H2', 'C2H4': 'C2H4',
    'C2H6': 'C2H6', 'CN': 'CN', 'HCN': 'HCN', 'CO': 'CO',
    'HCO': 'HCO', 'H2CO': 'H2CO', 'CH3OH': 'CH3OH', 'N2': 'N2',
    'N2H4': 'N2H4', 'NO': 'NO', 'O2': 'O2', 'H2O2': 'H2O2',
    'F2': 'F2', 'CO2': 'CO2', 'Na2': 'Na2', 'Si2': 'Si2',
    'P2': 'P2', 'S2': 'S2', 'Cl2': 'Cl2', 'NaCl': 'NaCl',
    'SiO': 'SiO', 'CS': 'CS', 'SO': 'SO', 'ClO': 'ClO',
    'ClF': 'ClF', 'Si2H6': 'Si2H6', 'CH3Cl': 'CH3Cl',
    'CH3SH': 'CH3SH', 'HOCl': 'HOCl', 'SO2': 'SO2',
}

# Mapping from g2_reference.py formula names to our 55-molecule set formula names
G2_97_TO_FORMULA = {
    'LiH': 'LiH', 'BeH': 'BeH', 'CH': 'CH', 'CH2_singlet': 'CH2_singlet',
    'CH2': 'CH2', 'CH3': 'CH3', 'CH4': 'CH4', 'NH': 'NH', 'NH2': 'NH2',
    'NH3': 'NH3', 'HO': 'OH', 'H2O': 'H2O', 'HF': 'HF', 'SiH2': 'SiH2',
    'SiH2_triplet': 'SiH2_triplet', 'SiH3': 'SiH3', 'SiH4': 'SiH4',
    'PH2': 'PH2', 'PH3': 'PH3', 'H2S': 'H2S', 'HCl': 'HCl', 'Li2': 'Li2',
    'LiF': 'LiF', 'C2H2': 'C2H2', 'C2H4': 'C2H4', 'C2H6': 'C2H6', 'CN': 'CN',
    'CNH': 'HCN', 'CO': 'CO', 'CHO': 'HCO', 'CH2O': 'H2CO', 'CH4O': 'CH3OH',
    'N2': 'N2', 'N2H4': 'N2H4', 'NO': 'NO', 'O2': 'O2', 'H2O2': 'H2O2',
    'F2': 'F2', 'CO2': 'CO2', 'Na2': 'Na2', 'Si2': 'Si2', 'P2': 'P2',
    'S2': 'S2', 'Cl2': 'Cl2', 'NaCl': 'NaCl', 'SiO': 'SiO', 'CS': 'CS',
    'OS': 'SO', 'OCl': 'ClO', 'FCl': 'ClF', 'Si2H6': 'Si2H6', 'CH3Cl': 'CH3Cl',
    'CH4S': 'CH3SH', 'HOCl': 'HOCl', 'O2S': 'SO2',
}

# Mapping from ccsd_energies.csv molecule names to our formula names (xTC-PP-CCSD(T) data)
CCSD_TO_FORMULA = {
    'LiH': 'LiH', 'BeH': 'BeH', 'CH': 'CH', 'CH2_singlet': 'CH2_singlet',
    'CH2': 'CH2', 'CH3': 'CH3', 'CH4': 'CH4', 'NH': 'NH', 'NH2': 'NH2',
    'NH3': 'NH3', 'HO': 'OH', 'H2O': 'H2O', 'HF': 'HF', 'SiH2': 'SiH2',
    'SiH2_triplet': 'SiH2_triplet', 'SiH3': 'SiH3', 'SiH4': 'SiH4',
    'PH2': 'PH2', 'PH3': 'PH3', 'H2S': 'H2S', 'HCl': 'HCl', 'Li2': 'Li2',
    'LiF': 'LiF', 'C2H2': 'C2H2', 'C2H4': 'C2H4', 'C2H6': 'C2H6', 'CN': 'CN',
    'CNH': 'HCN', 'CO': 'CO', 'CHO': 'HCO', 'CH2O': 'H2CO', 'CH4O': 'CH3OH',
    'N2': 'N2', 'N2H4': 'N2H4', 'NO': 'NO', 'O2': 'O2', 'H2O2': 'H2O2',
    'F2': 'F2', 'CO2': 'CO2', 'Na2': 'Na2', 'Si2': 'Si2', 'P2': 'P2',
    'S2': 'S2', 'Cl2': 'Cl2', 'NaCl': 'NaCl', 'SiO': 'SiO', 'CS': 'CS',
    'OS': 'SO', 'OCl': 'ClO', 'FCl': 'ClF', 'Si2H6': 'Si2H6', 'CH3Cl': 'CH3Cl',
    'CH4S': 'CH3SH', 'HOCl': 'HOCl', 'O2S': 'SO2',
}

# Mapping from molpro_energies.csv system names to our formula names
# Molpro uses lowercase with underscores, e.g., "c_h4" -> "CH4"
MOLPRO_TO_FORMULA = {
    'li_h': 'LiH', 'be_h': 'BeH', 'c_h': 'CH', 'c_h2_singlet': 'CH2_singlet',
    'c_h2_triplet': 'CH2', 'c_h3': 'CH3', 'c_h4': 'CH4', 'n_h': 'NH', 'n_h2': 'NH2',
    'n_h3': 'NH3', 'o_h': 'OH', 'h2_o': 'H2O', 'h_f': 'HF', 'si_h2_singlet': 'SiH2',
    'si_h2_triplet': 'SiH2_triplet', 'si_h3': 'SiH3', 'si_h4': 'SiH4',
    'p_h2': 'PH2', 'p_h3': 'PH3', 'h2_s': 'H2S', 'h_cl': 'HCl', 'li2': 'Li2',
    'li_f': 'LiF', 'c2_h2': 'C2H2', 'c2_h4': 'C2H4', 'c2_h6': 'C2H6', 'c_n': 'CN',
    'h_c_n': 'HCN', 'c_o': 'CO', 'h_c_o': 'HCO', 'h2_c_o': 'H2CO', 'h3_c_o_h': 'CH3OH',
    'n2': 'N2', 'n2_h4': 'N2H4', 'n_o': 'NO', 'o2': 'O2', 'h2_o2': 'H2O2',
    'f2': 'F2', 'c_o2': 'CO2', 'na2': 'Na2', 'si2': 'Si2', 'p2': 'P2',
    's2': 'S2', 'cl2': 'Cl2', 'na_cl': 'NaCl', 'si_o': 'SiO', 'c_s': 'CS',
    's_o': 'SO', 'cl_o': 'ClO', 'cl_f': 'ClF', 'si2_h6': 'Si2H6', 'c_h3_cl': 'CH3Cl',
    'h3_c_s_h': 'CH3SH', 'h_o_cl': 'HOCl', 's_o2': 'SO2',
    # Atoms
    'h': 'H', 'li': 'Li', 'be': 'Be', 'c': 'C', 'n': 'N', 'o': 'O', 'f': 'F',
    'na': 'Na', 'si': 'Si', 'p': 'P', 's': 'S', 'cl': 'Cl',
}


# Mapping from ccsd_energies.csv method column names to short slugs for output columns
METHOD_COLUMNS = {
    'xTC-CCSD(T)': 'xtc_ccsdt',
    'md_xTC-CCSD(T)': 'md_xtc_ccsdt',
    'xTC-DC-CCSDT': 'xtc_dc_ccsdt',
    'md_xTC-DC-CCSDT': 'md_xtc_dc_ccsdt',
    'xTC-FCIQMC': 'xtc_fciqmc',
    'md_xTC-FCIQMC': 'md_xtc_fciqmc',
}

# Non-md methods eligible for QZ HF correction
NON_MD_METHOD_COLUMNS = {
    'xTC-CCSD(T)': 'xtc_ccsdt',
    'xTC-DC-CCSDT': 'xtc_dc_ccsdt',
    'xTC-FCIQMC': 'xtc_fciqmc',
}


def load_f12_data():
    """Load F12 reference data and convert to DataFrame."""
    f12_df = pd.DataFrame(G2_97_ATOMIZATION)
    f12_df['our_formula'] = f12_df['formula'].map(G2_97_TO_FORMULA)
    f12_df = f12_df[f12_df['our_formula'].notna()].copy()
    f12_df['E_f12_kcal'] = f12_df['E_f12'] * KJ_TO_KCAL
    f12_df['E_exp_f12_kcal'] = f12_df['E_exp'] * KJ_TO_KCAL
    return f12_df[['our_formula', 'E_f12_kcal', 'E_exp_f12_kcal']].rename(
        columns={'our_formula': 'formula', 'E_f12_kcal': 'De_f12', 'E_exp_f12_kcal': 'Exp_f12'}
    )


def load_ccsd_data(ccsd_csv_path: Path, basis_short: str = 'avdz',
                   method_col: str = 'xTC-CCSD(T)', col_slug: str = 'xtc_ccsdt'):
    """
    Load energies from ccsd_energies.csv for a given method column and compute
    atomization energies.

    Args:
        ccsd_csv_path: Path to ccsd_energies.csv
        basis_short: Basis set tag (e.g. 'avdz', 'avtz')
        method_col: Column name in the CSV to read energies from
        col_slug: Short slug for naming the output De column
    """
    df = pd.read_csv(ccsd_csv_path)
    df = df[df['basis_short'] == basis_short].copy()

    if method_col not in df.columns:
        print(f"Warning: Column '{method_col}' not found in {ccsd_csv_path}")
        return pd.DataFrame()

    if len(df) == 0:
        print(f"Warning: No data found for basis {basis_short}")
        return pd.DataFrame()

    # Drop rows where the method column is empty / NaN
    df[method_col] = pd.to_numeric(df[method_col], errors='coerce')
    df = df[df[method_col].notna()].copy()

    if len(df) == 0:
        return pd.DataFrame()

    # Build atomic energies dictionary
    atoms = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl']
    atomic_energies = {}
    for _, row in df.iterrows():
        if row['molecule'] in atoms:
            atomic_energies[row['molecule']] = row[method_col]

    # H atom: exact value for single electron (no correlation)
    atomic_energies['H'] = -0.5

    # Compute atomization energies
    results = []
    for _, row in df.iterrows():
        ccsd_mol = row['molecule']
        our_formula = CCSD_TO_FORMULA.get(ccsd_mol)

        if our_formula is None or our_formula not in FORMULA_TO_ATOMS:
            continue

        atoms_dict = FORMULA_TO_ATOMS[our_formula]
        mol_energy = row[method_col]

        # Skip unphysical energies
        if mol_energy > 0 or mol_energy < -1e8:
            continue

        total_atomic = 0.0
        missing_atom = False
        for atom, count in atoms_dict.items():
            if atom not in atomic_energies:
                missing_atom = True
                break
            total_atomic += count * atomic_energies[atom]

        if missing_atom:
            continue

        De = total_atomic - mol_energy
        De_kcal = De * HA_TO_KCAL

        results.append({
            'formula': our_formula,
            f'De_{col_slug}_{basis_short}': De_kcal,
        })

    return pd.DataFrame(results)


def load_ccsd_data_qz_corrected(ccsd_csv_path: Path, basis_short: str,
                                method_col: str, col_slug: str):
    """
    Load xTC energies and apply QZ HF correction to compute atomization energies.

    For each system (atom or molecule):
        E_corr = E_method(basis) - E_xTC_HF(basis) + E_xTC_HF(avqz)

    Then De is computed from corrected energies as usual.

    Args:
        ccsd_csv_path: Path to ccsd_energies.csv
        basis_short: Basis set tag to correct (e.g. 'avdz', 'avtz')
        method_col: Method column name (e.g. 'xTC-CCSD(T)')
        col_slug: Short slug for naming output column

    Returns:
        DataFrame with formula and De_{col_slug}_{basis_short}_QZ_corr column (kcal/mol)
    """
    full_df = pd.read_csv(ccsd_csv_path)

    df_basis = full_df[full_df['basis_short'] == basis_short].copy()
    df_qz = full_df[full_df['basis_short'] == 'avqz'].copy()

    if len(df_basis) == 0 or len(df_qz) == 0:
        return pd.DataFrame()

    for col in [method_col, 'xTC-HF']:
        if col in df_basis.columns:
            df_basis[col] = pd.to_numeric(df_basis[col], errors='coerce')
        if col in df_qz.columns:
            df_qz[col] = pd.to_numeric(df_qz[col], errors='coerce')

    if method_col not in df_basis.columns or 'xTC-HF' not in df_basis.columns:
        return pd.DataFrame()

    # Build lookup: molecule -> xTC-HF energy at avqz
    hf_qz = {}
    for _, row in df_qz.iterrows():
        if pd.notna(row.get('xTC-HF')):
            hf_qz[row['molecule']] = row['xTC-HF']

    # Build corrected atomic energies
    atoms = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl']
    atomic_energies_corr = {}
    for _, row in df_basis.iterrows():
        mol = row['molecule']
        if mol not in atoms:
            continue
        e_method = row[method_col]
        e_hf_basis = row['xTC-HF']
        e_hf_qz = hf_qz.get(mol)
        if pd.isna(e_method) or pd.isna(e_hf_basis) or e_hf_qz is None:
            continue
        atomic_energies_corr[mol] = e_method - e_hf_basis + e_hf_qz

    # H atom: single electron, no correlation or HF correction
    atomic_energies_corr['H'] = -0.5

    # Compute corrected atomization energies
    results = []
    for _, row in df_basis.iterrows():
        mol = row['molecule']
        our_formula = CCSD_TO_FORMULA.get(mol)
        if our_formula is None or our_formula not in FORMULA_TO_ATOMS:
            continue

        e_method = row[method_col]
        e_hf_basis = row['xTC-HF']
        e_hf_qz = hf_qz.get(mol)
        if pd.isna(e_method) or pd.isna(e_hf_basis) or e_hf_qz is None:
            continue

        # Skip unphysical energies
        if e_method > 0 or e_method < -1e8:
            continue

        e_corr = e_method - e_hf_basis + e_hf_qz

        atoms_dict = FORMULA_TO_ATOMS[our_formula]
        total_atomic = 0.0
        missing = False
        for atom, count in atoms_dict.items():
            if atom not in atomic_energies_corr:
                missing = True
                break
            total_atomic += count * atomic_energies_corr[atom]

        if missing:
            continue

        De = total_atomic - e_corr
        De_kcal = De * HA_TO_KCAL
        results.append({
            'formula': our_formula,
            f'De_{col_slug}_{basis_short}_QZ_corr': De_kcal,
        })

    return pd.DataFrame(results)


def load_molpro_ccsdt_data(data_dir: Path, basis: str = '3z'):
    """
    Load CCSDT energies from molpro_energies.csv and compute atomization energies.
    
    Args:
        data_dir: Directory containing molpro_energies.csv
        basis: Basis set to use (2z, 3z, 4z, 5z, 2dz, 3dz, 4dz, 5dz)
    
    Returns:
        DataFrame with formula and De_ccsdt_{basis} column (kcal/mol)
    """
    molpro_csv = data_dir / 'umrigar_molpro_reference.csv'
    if not molpro_csv.exists():
        print(f"Warning: {molpro_csv} not found")
        return pd.DataFrame()
    
    df = pd.read_csv(molpro_csv)
    # Strip whitespace from method column (has trailing space in "CCSD ")
    df['method'] = df['method'].str.strip()
    
    # Filter for CCSDT method and specified basis
    df = df[(df['method'] == 'CCSDT') & (df['basis'] == basis)].copy()
    
    if len(df) == 0:
        print(f"Warning: No CCSDT data found for basis {basis}")
        return pd.DataFrame()
    
    # Map system names to our formula names
    df['formula'] = df['system'].map(MOLPRO_TO_FORMULA)
    df = df[df['formula'].notna()].copy()
    
    # Build atomic energies dictionary
    atoms = ['H', 'Li', 'Be', 'C', 'N', 'O', 'F', 'Na', 'Si', 'P', 'S', 'Cl']
    atomic_energies = {}
    for _, row in df.iterrows():
        if row['formula'] in atoms:
            atomic_energies[row['formula']] = row['energy']
    
    # Compute atomization energies for molecules
    results = []
    for _, row in df.iterrows():
        formula = row['formula']
        if formula not in FORMULA_TO_ATOMS:
            continue
        
        atoms_dict = FORMULA_TO_ATOMS[formula]
        mol_energy = row['energy']
        
        total_atomic = 0.0
        missing_atom = False
        for atom, count in atoms_dict.items():
            if atom not in atomic_energies:
                missing_atom = True
                break
            total_atomic += count * atomic_energies[atom]
        
        if missing_atom:
            continue
        
        De = total_atomic - mol_energy
        De_kcal = De * HA_TO_KCAL
        
        results.append({
            'formula': formula,
            f'De_ccsdt_{basis}': De_kcal,
        })
    
    return pd.DataFrame(results)


def load_data(data_dir: Path):
    """Load atomization reference data."""
    atomization_ref_df = pd.read_csv(data_dir / 'umrigar_ae_reference.csv')
    return atomization_ref_df


def merge_with_corrections(reference_df: pd.DataFrame,
                          ccsd_csv_path: Path = None,
                          data_dir: Path = None):
    """Build merged atomization DataFrame from reference data and other sources."""
    reference_df = reference_df.copy()
    reference_df['formula'] = reference_df['molecule'].map(REFERENCE_TO_FORMULA)
    
    merged = reference_df[['formula', 'SHCI_De', 'ZPE', 'SR_SO', 'CV', 'experiment', 'D0_exp2', 'D0_calc']].copy()
    merged = merged[merged['formula'].notna()].reset_index(drop=True)
    
    merged = merged.rename(columns={
        'SHCI_De': 'De_shci_ref',
        'D0_calc': 'D0_shci_ref',
        'D0_exp2': 'D0_shci_pbe',
    })
    
    # Derive raw De for SHCI+PBE reference (invert corrections)
    if 'D0_shci_pbe' in merged.columns:
        merged['De_shci_pbe'] = merged['D0_shci_pbe'] - merged['ZPE'] - merged['SR_SO'] - merged['CV']
    
    # Load and merge F12 reference data (CC F12 CBS from Knizia paper)
    f12_df = load_f12_data()
    merged = merged.merge(f12_df, on='formula', how='left')
    
    if 'De_f12' in merged.columns:
        merged['D0_f12'] = merged['De_f12'] + merged['ZPE'] + merged['SR_SO'] + merged['CV']
    if 'Exp_f12' in merged.columns:
        merged['D0_exp_f12'] = merged['Exp_f12'] + merged['ZPE'] + merged['SR_SO'] + merged['CV']
    
    # Load and merge xTC atomization energies from ccsd_energies.csv for all methods
    if ccsd_csv_path is None:
        ccsd_csv_path = Path(__file__).parent.parent / 'data' / 'ccsd_energies.csv'

    if ccsd_csv_path.exists():
        # Load raw xTC atomization energies for all methods and bases
        for basis in ['avdz', 'avtz', 'avqz', 'pvdz', 'pvtz', 'pvqz']:
            for method_col, col_slug in METHOD_COLUMNS.items():
                ccsd_df = load_ccsd_data(ccsd_csv_path, basis,
                                         method_col=method_col, col_slug=col_slug)
                if len(ccsd_df) > 0:
                    merged = merged.merge(ccsd_df, on='formula', how='left')
                    de_col = f'De_{col_slug}_{basis}'
                    d0_col = f'D0_{col_slug}_{basis}'
                    if de_col in merged.columns:
                        merged[d0_col] = merged[de_col] + merged['ZPE'] + merged['SR_SO'] + merged['CV']

        # Load QZ-HF-corrected xTC atomization energies (non-md methods only)
        for basis in ['avdz', 'avtz', 'pvdz', 'pvtz']:
            for method_col, col_slug in NON_MD_METHOD_COLUMNS.items():
                ccsd_df = load_ccsd_data_qz_corrected(
                    ccsd_csv_path, basis,
                    method_col=method_col, col_slug=col_slug)
                if len(ccsd_df) > 0:
                    merged = merged.merge(ccsd_df, on='formula', how='left')
                    de_col = f'De_{col_slug}_{basis}_QZ_corr'
                    d0_col = f'D0_{col_slug}_{basis}_QZ_corr'
                    if de_col in merged.columns:
                        merged[d0_col] = merged[de_col] + merged['ZPE'] + merged['SR_SO'] + merged['CV']
    
    # Load and merge CCSDT atomization energies from molpro_energies.csv for different basis sets
    if data_dir is None:
        data_dir = Path(__file__).parent.parent / 'reference_data'
    
    for basis in ['2z', '3z', '4z', '5z']:
        ccsdt_df = load_molpro_ccsdt_data(data_dir, basis)
        if len(ccsdt_df) > 0:
            merged = merged.merge(ccsdt_df, on='formula', how='left')
            de_col = f'De_ccsdt_{basis}'
            d0_col = f'D0_ccsdt_{basis}'
            if de_col in merged.columns:
                merged[d0_col] = merged[de_col] + merged['ZPE'] + merged['SR_SO'] + merged['CV']
    
    return merged


def print_summary(merged_df: pd.DataFrame):
    """Print a summary of the merged data."""
    print(f"\nMerged {len(merged_df)} molecules with all data sources")
    print("\nColumns in output:")
    for col in merged_df.columns:
        non_null = merged_df[col].notna().sum()
        print(f"  {col}: {non_null}/{len(merged_df)} non-null")


def main():
    parser = argparse.ArgumentParser(description='Build merged atomization energy data')
    parser.add_argument('--data-dir', type=Path,
                        default=Path(__file__).parent.parent / 'reference_data',
                        help='Directory containing reference CSV files')
    parser.add_argument('--ccsd-csv', type=Path, default=None,
                        help='Path to ccsd_energies.csv (default: data/ccsd_energies.csv)')
    parser.add_argument('--output', type=Path, default=None,
                        help='Output CSV file (default: atomization_comparison.csv in data-dir)')
    args = parser.parse_args()
    
    data_dir = args.data_dir
    output_csv = args.output or (data_dir / 'atomization_comparison.csv')
    ccsd_csv_path = args.ccsd_csv or (Path(__file__).parent.parent / 'data' / 'ccsd_energies.csv')
    
    # Load data
    print("Loading data...")
    reference_df = load_data(data_dir)
    print(f"  Loaded {len(reference_df)} molecules from reference data")
    
    # Merge with correction data and other sources
    print("Merging with corrections and other data sources...")
    merged_df = merge_with_corrections(reference_df, ccsd_csv_path, data_dir)
    
    print_summary(merged_df)
    
    # Save merged data
    merged_df.to_csv(output_csv, index=False)
    print(f"\nSaved merged data to {output_csv}")


if __name__ == '__main__':
    main()

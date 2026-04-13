#!/usr/bin/env python3
"""
Compute molecular atomization energies from ΛCCSD(T) energy data.

Atomization energy = Σ(atomic energies) - molecular energy

Usage:
    python calculate_atomization_energies.py
"""

import argparse
import csv
from pathlib import Path

# Hartree-to-kcal/mol conversion factor
HARTREE_TO_KCAL = 627.5094740631

# Atom name mapping (name in file -> element symbol)
ATOM_NAMES = {
    'h': 'H',
    'c': 'C',
    'n': 'N',
    'o': 'O',
    'f': 'F',
    'li': 'Li',
    'be': 'Be',
    'na': 'Na',
    'si': 'Si',
    'p': 'P',
    's': 'S',
    'cl': 'Cl',
}

# Molecule composition map: molecule name -> {element: atom_count}
# Defined according to species names used in ccsd_t_energies.csv
MOLECULE_COMPOSITION = {
    # Diatomic molecules
    'li_h': {'Li': 1, 'H': 1},
    'be_h': {'Be': 1, 'H': 1},
    'c_h': {'C': 1, 'H': 1},
    'n_h': {'N': 1, 'H': 1},
    'o_h': {'O': 1, 'H': 1},
    'h_f': {'H': 1, 'F': 1},
    'h_cl': {'H': 1, 'Cl': 1},
    'li_f': {'Li': 1, 'F': 1},
    'li2': {'Li': 2},
    'na2': {'Na': 2},
    'na_cl': {'Na': 1, 'Cl': 1},
    'c_n': {'C': 1, 'N': 1},
    'c_o': {'C': 1, 'O': 1},
    'c_s': {'C': 1, 'S': 1},
    'n2': {'N': 2},
    'n_o': {'N': 1, 'O': 1},
    'o2': {'O': 2},
    'f2': {'F': 2},
    'si2': {'Si': 2},
    'p2': {'P': 2},
    's2': {'S': 2},
    'cl2': {'Cl': 2},
    'si_o': {'Si': 1, 'O': 1},
    's_o': {'S': 1, 'O': 1},
    'cl_o': {'Cl': 1, 'O': 1},
    'cl_f': {'Cl': 1, 'F': 1},
    
    # Triatomic and larger molecules
    'c_h2_singlet': {'C': 1, 'H': 2},
    'c_h2_triplet': {'C': 1, 'H': 2},
    'n_h2': {'N': 1, 'H': 2},
    'h2_o': {'H': 2, 'O': 1},
    'h2_s': {'H': 2, 'S': 1},
    'si_h2_singlet': {'Si': 1, 'H': 2},
    'si_h2_triplet': {'Si': 1, 'H': 2},
    'p_h2': {'P': 1, 'H': 2},
    'c_o2': {'C': 1, 'O': 2},
    's_o2': {'S': 1, 'O': 2},
    'h2_o2': {'H': 2, 'O': 2},
    
    'c_h3': {'C': 1, 'H': 3},
    'n_h3': {'N': 1, 'H': 3},
    'si_h3': {'Si': 1, 'H': 3},
    'p_h3': {'P': 1, 'H': 3},
    
    'c_h4': {'C': 1, 'H': 4},
    'si_h4': {'Si': 1, 'H': 4},
    
    'c2_h2': {'C': 2, 'H': 2},
    'c2_h4': {'C': 2, 'H': 4},
    'c2_h6': {'C': 2, 'H': 6},
    'n2_h4': {'N': 2, 'H': 4},
    'si2_h6': {'Si': 2, 'H': 6},
    
    'h_c_n': {'H': 1, 'C': 1, 'N': 1},
    'h_c_o': {'H': 1, 'C': 1, 'O': 1},
    'h2_c_o': {'H': 2, 'C': 1, 'O': 1},
    'h3_c_o_h': {'H': 4, 'C': 1, 'O': 1},  # CH3OH
    'h3_c_s_h': {'H': 4, 'C': 1, 'S': 1},  # CH3SH
    'c_h3_cl': {'C': 1, 'H': 3, 'Cl': 1},  # CH3Cl
    'h_o_cl': {'H': 1, 'O': 1, 'Cl': 1},   # HOCl
}


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Compute molecular atomization energies from ΛCCSD(T) energy data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-i', '--input_dirs',
        nargs='+',
        default=[str(Path(__file__).parent.resolve())],
        help='List of folders containing ccsd_t_energies.csv (default: script directory)'
    )
    parser.add_argument(
        '-f', '--filename',
        type=str,
        default='ccsd_t_energies.csv',
        help='Energy CSV filename (default: ccsd_t_energies.csv)'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=str((Path(__file__).parent / '../analysis/atomization_energies.csv').resolve()),
        help='Output CSV path'
    )
    parser.add_argument(
        '--H_05',
        action='store_true',
        help='If set, force H atomic energy to -0.5 Ha when computing atomization energies'
    )
    return parser.parse_args()


def read_ccsd_t_energies(filepath):
    """Read CCSD(T) energy data from CSV."""
    energies = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            species = row['species']
            energy = float(row['energy_hartree'])
            energies[species] = energy
    return energies


def calculate_atomization_energy(molecule, mol_energy, atom_energies):
    """
    Compute atomization energy.
    
    Args:
        molecule: Molecule name
        mol_energy: Molecular energy (Hartree)
        atom_energies: Atomic energy dict {element_symbol: energy}
    
    Returns:
        Atomization energy (Hartree), atomization energy (kcal/mol), formula string
    """
    if molecule not in MOLECULE_COMPOSITION:
        return None, None, None
    
    composition = MOLECULE_COMPOSITION[molecule]
    
    # Sum atomic energies
    sum_atom_energies = 0.0
    composition_str_parts = []
    
    for element, count in sorted(composition.items()):
        atom_key = [k for k, v in ATOM_NAMES.items() if v == element]
        if not atom_key:
            print(f"  Warning: missing atomic energy mapping for element {element}")
            return None, None, None
        atom_key = atom_key[0]
        
        if atom_key not in atom_energies:
            print(f"  Warning: missing energy data for atom {element} ({atom_key})")
            return None, None, None
        
        sum_atom_energies += atom_energies[atom_key] * count
        composition_str_parts.append(f"{element}{count}" if count > 1 else element)
    
    composition_str = ''.join(composition_str_parts)
    
    # Atomization energy = sum(atomic energies) - molecular energy
    ae_hartree = sum_atom_energies - mol_energy
    ae_kcal = ae_hartree * HARTREE_TO_KCAL
    
    return ae_hartree, ae_kcal, composition_str


def main():
    args = parse_args()
    output_csv = Path(args.output).resolve()
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    all_results = []
    print("-" * 80)
    
    for input_dir in args.input_dirs:
        input_dir_path = Path(input_dir).resolve()
        input_csv = input_dir_path / args.filename

        if not input_csv.exists():
            print(f"Skipping {input_dir_path}: {args.filename} not found")
            continue

        print(f"Reading data: {input_csv}")
        energies = read_ccsd_t_energies(input_csv)
        print(f"Read energies for {len(energies)} species")

        # Extract atomic energies
        atom_energies = {}
        for atom_key, element in ATOM_NAMES.items():
            if atom_key in energies:
                atom_energies[atom_key] = energies[atom_key]

        # Optional: force H energy to -0.5 Ha for sensitivity/manual override
        if getattr(args, 'H_05', False):
            atom_energies['h'] = -0.5
            print("  Note: --H_05 enabled, H energy forced to -0.5 Ha")

        print(f"Computing atomization energies... ({input_dir_path.name})")
        print("-" * 80)

        for species, mol_energy in energies.items():
            # Skip atoms
            if species in ATOM_NAMES:
                continue

            ae_hartree, ae_kcal, composition = calculate_atomization_energy(
                species, mol_energy, atom_energies
            )

            if ae_hartree is not None:
                all_results.append({
                    'basis_folder': input_dir_path.name,
                    'species': species,
                    'formula': composition,
                    'mol_energy_hartree': mol_energy,
                    'atomization_energy_hartree': ae_hartree,
                    'atomization_energy_kcal_mol': ae_kcal,
                })
                print(f"{species:20s} ({composition:8s}): {ae_kcal:10.2f} kcal/mol")
            else:
                print(f"{species:20s}: cannot compute (missing composition or atomic energies)")

    # Sort and write CSV
    all_results.sort(key=lambda x: (x['basis_folder'], x['species']))

    with open(output_csv, 'w', newline='') as f:
        fieldnames = [
            'basis_folder',
            'species',
            'formula',
            'mol_energy_hartree',
            'atomization_energy_hartree',
            'atomization_energy_kcal_mol',
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    print("-" * 80)
    print(f"\nResults saved to: {output_csv}")
    print(f"Computed atomization energies for {len(all_results)} molecules")


if __name__ == "__main__":
    main()

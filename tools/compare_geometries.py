#!/usr/bin/env python3
"""
Compare _G2_GEOM_RAW in g2_structures.py vs g2_structures_from_UMRIGAR.py.

Checks for structural differences:
  - element composition (symbols multiset)
  - multiplicity
  - charge
  - geometry (interatomic distance matrix, tolerant to rotation/translation)

Usage:
    python3 tools/compare_geometries.py
"""

import sys
import os
import re
from collections import Counter

import numpy as np

# Make src/ importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from g2_structures import _G2_GEOM_RAW as RAW_NEW
from g2_structures_from_UMRIGAR import _G2_GEOM_RAW as RAW_REF

COORD_TOL = 1e-3   # Bohr - tolerance for distance-matrix comparison


def normalize_key(key: str) -> str:
    """Lowercase + strip underscores for fuzzy matching between naming styles."""
    return re.sub(r"[_\-]", "", key).lower()


def distance_matrix(coords: np.ndarray) -> np.ndarray:
    """Return sorted upper-triangle of pairwise distances (rotation-invariant)."""
    diff = coords[:, None, :] - coords[None, :, :]
    dmat = np.sqrt((diff ** 2).sum(axis=-1))
    n = len(coords)
    idx = np.triu_indices(n, k=1)
    return np.sort(dmat[idx])


def compare_entry(key_new, d_new, key_ref, d_ref):
    """
    Compare two geometry dict entries. Returns list of difference strings.
    """
    diffs = []

    sym_new = Counter(d_new["symbols"])
    sym_ref = Counter(d_ref["symbols"])
    if sym_new != sym_ref:
        diffs.append(
            f"  symbols: {dict(sym_new)} (new)  vs  {dict(sym_ref)} (ref)"
        )

    mult_new = d_new.get("multiplicity", None)
    mult_ref = d_ref.get("multiplicity", None)
    if mult_new != mult_ref:
        diffs.append(f"  multiplicity: {mult_new} (new)  vs  {mult_ref} (ref)")

    chg_new = d_new.get("charge", 0)
    chg_ref = d_ref.get("charge", 0)
    if chg_new != chg_ref:
        diffs.append(f"  charge: {chg_new} (new)  vs  {chg_ref} (ref)")

    # Geometry comparison only makes sense when same elements
    if sym_new == sym_ref and len(d_new["symbols"]) > 1:
        dm_new = distance_matrix(np.array(d_new["coords_bohr"]))
        dm_ref = distance_matrix(np.array(d_ref["coords_bohr"]))
        if dm_new.shape == dm_ref.shape:
            max_diff = np.max(np.abs(dm_new - dm_ref))
            if max_diff > COORD_TOL:
                diffs.append(
                    f"  geometry: max |Δd| = {max_diff:.4f} Bohr"
                )
        else:
            diffs.append(
                f"  geometry: atom count mismatch "
                f"({len(d_new['symbols'])} vs {len(d_ref['symbols'])})"
            )

    return diffs


def main():
    # Build normalized-key → (original_key, data) maps
    norm_new = {normalize_key(k): (k, v) for k, v in RAW_NEW.items()}
    norm_ref = {normalize_key(k): (k, v) for k, v in RAW_REF.items()}

    all_norm = sorted(set(norm_new) | set(norm_ref))

    only_new = []
    only_ref = []
    different = []
    matching = []

    for nk in all_norm:
        in_new = nk in norm_new
        in_ref = nk in norm_ref

        if in_new and not in_ref:
            only_new.append(norm_new[nk][0])
        elif in_ref and not in_new:
            only_ref.append(norm_ref[nk][0])
        else:
            key_new, d_new = norm_new[nk]
            key_ref, d_ref = norm_ref[nk]
            diffs = compare_entry(key_new, d_new, key_ref, d_ref)
            if diffs:
                different.append((key_new, key_ref, diffs))
            else:
                matching.append(key_new)

    # ── Report ────────────────────────────────────────────────────────────────
    print(f"{'='*60}")
    print(f"Matched entries : {len(matching) + len(different)}")
    print(f"Identical       : {len(matching)}")
    print(f"Different       : {len(different)}")
    print(f"Only in new     : {len(only_new)}")
    print(f"Only in ref     : {len(only_ref)}")
    print(f"{'='*60}\n")

    if different:
        print("── Structural differences ──────────────────────────────────")
        for key_new, key_ref, diffs in different:
            print(f"\n[{key_new}]  (ref key: {key_ref})")
            for d in diffs:
                print(d)

    if only_new:
        print("\n── Only in g2_structures.py ────────────────────────────────")
        for k in sorted(only_new):
            print(f"  {k}")

    if only_ref:
        print("\n── Only in g2_structures_from_UMRIGAR.py ───────────────────")
        for k in sorted(only_ref):
            print(f"  {k}")

    if not different and not only_new and not only_ref:
        print("No differences found — the two files are structurally identical.")


if __name__ == "__main__":
    main()

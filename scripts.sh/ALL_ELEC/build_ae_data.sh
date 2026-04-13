#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if [[ ! -d "$REPO_ROOT/.venv" ]]; then
    echo "Error: .venv not found. Run: python3 -m venv .venv && .venv/bin/pip install -r requirements.txt"
    exit 1
fi

mapfile -t OUTPUT_DIRS < <(
    find "$REPO_ROOT/data/ALL_ELEC/RAW" -type f -name '*.csv' -printf '%h\n' | sort -u
)

if [[ ${#OUTPUT_DIRS[@]} -gt 0 ]]; then
    "$REPO_ROOT/.venv/bin/python" "$REPO_ROOT/tools/ALL_ELEC/calculate_atomization_energies.py" \
        -i "${OUTPUT_DIRS[@]}" \
        -o "$REPO_ROOT/data/ALL_ELEC/atomization_energies_H05.csv" \
        --H_05
fi

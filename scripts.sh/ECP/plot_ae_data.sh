#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if [[ ! -d "$REPO_ROOT/.venv" ]]; then
    echo "Error: .venv not found. Run: python3 -m venv .venv && .venv/bin/pip install -r requirements.txt"
    exit 1
fi

for basis in pvdz pvtz pvqz avdz avtz avqz; do
    echo "--- basis: $basis ---"
    "$REPO_ROOT/.venv/bin/python" "$REPO_ROOT/tools/ECP/atomization_plots.py" \
        "$REPO_ROOT/data/ECP/atomization_comparison.csv" \
        --output-dir "$REPO_ROOT/figures/ECP" \
        --html-output-dir "$REPO_ROOT/docs/figures/ECP" \
        --no-show \
        --basis "$basis" \
        "$@"
done

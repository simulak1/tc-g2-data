#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ ! -d "$REPO_ROOT/.venv" ]]; then
    echo "Error: .venv not found. Run: python3 -m venv .venv && .venv/bin/pip install -r requirements.txt"
    exit 1
fi

"$REPO_ROOT/.venv/bin/python" "$REPO_ROOT/tools/make_latex_table.py" \
    --data-dir  "$REPO_ROOT/data" \
    --output-dir "$REPO_ROOT/tables" \
    "$@"

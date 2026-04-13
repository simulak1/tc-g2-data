#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PLOT_SCRIPT="$REPO_ROOT/tools/ALL_ELEC/plot_atomization_energy_discrepancy.py"
PNG_ROOT="$REPO_ROOT/figures/ALL_ELEC"
HTML_ROOT="$REPO_ROOT/docs/figures/ALL_ELEC"
CSV_ROOT="$REPO_ROOT/docs/figures/ALL_ELEC"

### MIX BASIS
basis_folder=("NacV-aug-ano-pV(D+d)Z" "NacV-aug-ano-pV(T+d)Z" "NacV-aug-ano-pV(Q+d)Z" "NacV-aug-ano-pV(5+d)Z" "NacV-aug-cc-pV(D+d)Z" "NacV-aug-cc-pV(T+d)Z" "NacV-aug-cc-pV(Q+d)Z" "NacV-aug-cc-pV(5+d)Z")
# basis_folder=("NacV-aug-ano-pV(Q+d)Z" "NacV-aug-ano-pV(5+d)Z" "NacV-aug-cc-pV(Q+d)Z" "NacV-aug-cc-pV(5+d)Z")
# python3 "$PLOT_SCRIPT" "$REPO_ROOT/data/ALL_ELEC/atomization_energies_H05.csv" \
#   -o atomization_H05/EXP/Mix/Prune \
#   --png-output-root "$PNG_ROOT" \
#   --html-output-root "$HTML_ROOT" \
#   --csv-output-root "$CSV_ROOT" \
#   --basis-folder "${basis_folder[@]}" \
#   --basis-order "${basis_folder[@]}" \
#   --stats

python3 "$PLOT_SCRIPT" "$REPO_ROOT/data/ALL_ELEC/atomization_energies_H05.csv" \
  -o atomization_H05/PBE/Mix/Prune \
  --png-output-root "$PNG_ROOT" \
  --html-output-root "$HTML_ROOT" \
  --csv-output-root "$CSV_ROOT" \
  --basis-folder "${basis_folder[@]}" \
  --basis-order "${basis_folder[@]}" \
  --stats --theo-ref "PBE"

#### Ref corrected
basis_folder=("aug-ano-pV(D+d)Z_refcorr" "aug-ano-pV(T+d)Z_refcorr" "aug-ano-pV(Q+d)Z_refcorr" "aug-cc-pV(D+d)Z_refcorr" "aug-cc-pV(T+d)Z_refcorr" "aug-cc-pV(Q+d)Z_refcorr")
# python3 "$PLOT_SCRIPT" "$REPO_ROOT/data/ALL_ELEC/atomization_energies_H05.csv" \
#   -o atomization_H05/EXP/Ref-corr \
#   --png-output-root "$PNG_ROOT" \
#   --html-output-root "$HTML_ROOT" \
#   --csv-output-root "$CSV_ROOT" \
#   --basis-folder "${basis_folder[@]}" \
#   --basis-order "${basis_folder[@]}" \
#   --stats
python3 "$PLOT_SCRIPT" "$REPO_ROOT/data/ALL_ELEC/atomization_energies_H05.csv" \
  -o atomization_H05/PBE/Ref-corr \
  --png-output-root "$PNG_ROOT" \
  --html-output-root "$HTML_ROOT" \
  --csv-output-root "$CSV_ROOT" \
  --basis-folder "${basis_folder[@]}" \
  --basis-order "${basis_folder[@]}" \
  --stats --theo-ref "PBE"

# ### SVD
basis_folder=("aug-cc-pV(D+d)Z_svd_proj" "aug-cc-pV(T+d)Z_svd_proj" "aug-cc-pV(Q+d)Z_svd_proj" "aug-ano-pV(D+d)Z_svd_proj" "aug-ano-pV(T+d)Z_svd_proj" "aug-ano-pV(Q+d)Z_svd_proj")
# python3 "$PLOT_SCRIPT" "$REPO_ROOT/data/ALL_ELEC/atomization_energies_H05.csv" \
#  -o atomization_H05/EXP/SVD\
#  --png-output-root "$PNG_ROOT" \
#  --html-output-root "$HTML_ROOT" \
#  --csv-output-root "$CSV_ROOT" \
#  --basis-folder "${basis_folder[@]}"\
#  --basis-order "${basis_folder[@]}"\
#  --stats
python3 "$PLOT_SCRIPT" "$REPO_ROOT/data/ALL_ELEC/atomization_energies_H05.csv" \
 -o atomization_H05/PBE/SVD\
 --png-output-root "$PNG_ROOT" \
 --html-output-root "$HTML_ROOT" \
 --csv-output-root "$CSV_ROOT" \
 --basis-folder "${basis_folder[@]}"\
 --basis-order "${basis_folder[@]}"\
 --theo-ref "PBE"\
 --stats

#### F12
basis_folder=("F12_aug-cc-pV(D+d)Z" "F12_aug-cc-pV(T+d)Z" "F12_aug-cc-pV(Q+d)Z" "F12_aug-cc-pV(5+d)Z")
# python3 "$PLOT_SCRIPT" "$REPO_ROOT/data/ALL_ELEC/F12_atomization_energies_H05.csv" \
#  -o atomization_H05/EXP/F12\
#  --png-output-root "$PNG_ROOT" \
#  --html-output-root "$HTML_ROOT" \
#  --csv-output-root "$CSV_ROOT" \
#  --basis-folder "${basis_folder[@]}"\
#  --basis-order "${basis_folder[@]}"\
#  --stats
python3 "$PLOT_SCRIPT" "$REPO_ROOT/data/ALL_ELEC/F12_atomization_energies_H05.csv" \
 -o atomization_H05/PBE/F12\
 --png-output-root "$PNG_ROOT" \
 --html-output-root "$HTML_ROOT" \
 --csv-output-root "$CSV_ROOT" \
 --basis-folder "${basis_folder[@]}"\
 --basis-order "${basis_folder[@]}"\
 --theo-ref "PBE"\
 --stats

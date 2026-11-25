#!/bin/bash
# Reproduce every synthetic result. Run from anywhere; paths derive from this script.
#
#   bash scripts/run_all_studyd.sh
#
# Override the interpreter with CLDYN_PY. If your matplotlib dependencies live in a
# separate site directory, point CLDYN_SITE at it.
#
# The real-data analyses are NOT here: GSE173958 and GSE161363 must be retrieved first.
# See docs/data_provenance.md and section 6 of the README.
set -euo pipefail
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
ROOT="${CLDYN_ROOT:-$(cd "$_HERE/.." && pwd)}"
PY="${CLDYN_PY:-python}"
cd "$ROOT"
mkdir -p results/tables results/figures

echo "== Milestone 1: bounds, scenarios, design (~80 s)"
$PY -I experiments/run_all.py > results/run_all.log 2>&1

echo "== Milestone 2: lineage simulator, mechanism signatures and recovery"
$PY -I experiments/lineage_identifiability/run_e8.py > results/e8.log 2>&1
$PY -I experiments/lineage_identifiability/run_e9.py > results/e9.log 2>&1

echo "== Milestone 3: occupancy vs switching rate, sampling design"
$PY -I experiments/lineage_identifiability/run_e10_focus.py > results/e10focus.log 2>&1
$PY -I experiments/sampling_design/run_e11.py > results/e11.log 2>&1

echo "== Milestone 4: Fisher audit, full-tree information, calibration, capacity law"
$PY -I experiments/sampling_design/fisher_audit.py > results/e12.log 2>&1
$PY -I experiments/sampling_design/window_recompute.py > results/e12b.log 2>&1
$PY -I experiments/sampling_design/full_tree_information.py --mc > results/e13.log 2>&1
$PY -I experiments/calibration/run_e14.py > results/e14.log 2>&1          # ~1 h on 12 cores
$PY -I experiments/sampling_design/run_e15b.py > results/e15b.log 2>&1
$PY -I experiments/sampling_design/capacity_verification.py results/tables > results/e16.log 2>&1

echo "== Figures"
if [ -n "${CLDYN_SITE:-}" ]; then export PYTHONPATH="$CLDYN_SITE"; fi
$PY experiments/figures.py
$PY experiments/figures_m2.py results/tables results/figures
$PY experiments/figures_m3.py results/tables results/figures
$PY experiments/figures_m4.py results/tables results/figures

echo "== done; tables in results/tables, figures in results/figures"

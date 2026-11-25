#!/bin/bash
# Phase 4: held-out prediction (estimand group D).
set -u
# Paths are derived from this script's location so the repository runs anywhere.
# Override with CLDYN_BENCH (the benchmarks/ directory) and CLDYN_ROOT (the repo root).
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
B="${CLDYN_BENCH:-$(cd "$_HERE/../.." && pwd)}"
S="${CLDYN_ROOT:-$(cd "$B/.." && pwd)}"
PY="${CLDYN_PY:-python}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ls -d "$B"/data/sim_holdout/*/ 2>/dev/null | sed 's:/$::' > "$B/data/manifests/holdout.txt"
wc -l "$B/data/manifests/holdout.txt"
"${CLDYN_RSCRIPT:-Rscript}" "$B/src/adapters/run_competitors.R" \
  "$B" "$B/data/manifests/holdout.txt" true "$B/results/metrics/holdout_comp.csv" 1 \
  > "$B/logs/holdout_comp.log" 2>&1 &
cd "$S" && "$PY" -I "$B/src/adapters/run_studyd.py" \
  "$S" "$B/data/manifests/holdout.txt" true "$B/results/metrics/holdout_sd.csv" \
  > "$B/logs/holdout_sd.log" 2>&1 &
wait
echo PHASE4_DONE

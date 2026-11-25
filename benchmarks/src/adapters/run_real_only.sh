#!/bin/bash
# Real-data benchmark only, on the corrected trees.
set -u
# Paths are derived from this script's location so the repository runs anywhere.
# Override with CLDYN_BENCH (the benchmarks/ directory) and CLDYN_ROOT (the repo root).
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
B="${CLDYN_BENCH:-$(cd "$_HERE/../.." && pwd)}"
S="${CLDYN_ROOT:-$(cd "$B/.." && pwd)}"
PY="${CLDYN_PY:-python}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
real_job() {
  a=$1

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  [ -s "$B/data/manifests/real_${a}.txt" ] || { echo "skip $a"; return; }
  "${CLDYN_RSCRIPT:-Rscript}" "$B/src/adapters/run_competitors.R" \
    "$B" "$B/data/manifests/real_${a}.txt" unit "$B/results/metrics/real_${a}_comp.csv" 1 \
    > "$B/logs/real_${a}_comp.log" 2>&1
  cd "$S" && "$PY" -I "$B/src/adapters/run_studyd.py" \
    "$S" "$B/data/manifests/real_${a}.txt" unit "$B/results/metrics/real_${a}_sd.csv" \
    > "$B/logs/real_${a}_sd.log" 2>&1
  echo "done real $a"
}
export -f real_job
printf '%s\n' quinn_5k simeonov_M2 quinn_100k quinn_10k simeonov_M1 quinn_30k \
  | xargs -P 6 -I{} bash -c 'real_job {}'
echo REAL_ONLY_DONE

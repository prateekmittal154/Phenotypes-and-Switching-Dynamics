#!/bin/bash
# worker_gt.sh <scenario> <variant> <engine>   engine = comp | sd
set -u
# Paths are derived from this script's location so the repository runs anywhere.
# Override with CLDYN_BENCH (the benchmarks/ directory) and CLDYN_ROOT (the repo root).
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
B="${CLDYN_BENCH:-$(cd "$_HERE/../.." && pwd)}"
S="${CLDYN_ROOT:-$(cd "$B/.." && pwd)}"
PY="${CLDYN_PY:-python}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
sc=$1; v=$2; eng=$3
if [ "$eng" = "comp" ]; then
  "${CLDYN_RSCRIPT:-Rscript}" "$B/src/adapters/run_competitors.R" \
    "$B" "$B/data/manifests/list_${sc}.txt" "$v" "$B/results/metrics/gt_${sc}_${v}.csv" 1 \
    > "$B/logs/gt_${sc}_${v}.log" 2>&1
else
  cd "$S" && "$PY" -I "$B/src/adapters/run_studyd.py" \
    "$S" "$B/data/manifests/list_${sc}.txt" "$v" "$B/results/metrics/sd_${sc}_${v}.csv" \
    > "$B/logs/sd_${sc}_${v}.log" 2>&1
fi
echo "done $sc $v $eng"

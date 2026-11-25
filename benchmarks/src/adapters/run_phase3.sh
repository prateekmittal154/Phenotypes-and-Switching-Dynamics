#!/bin/bash
# Phase 3: reconstruction-uncertainty benchmark (IV).
set -u
# Paths are derived from this script's location so the repository runs anywhere.
# Override with CLDYN_BENCH (the benchmarks/ directory) and CLDYN_ROOT (the repo root).
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
B="${CLDYN_BENCH:-$(cd "$_HERE/../.." && pwd)}"
S="${CLDYN_ROOT:-$(cd "$B/.." && pwd)}"
PY="${CLDYN_PY:-python}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ls -d "$B"/data/sim_recon/R_rep*/ 2>/dev/null | sed 's:/$::' > "$B/data/manifests/recon.txt"
wc -l "$B/data/manifests/recon.txt"

rjob() {
  v=$1

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  "${CLDYN_RSCRIPT:-Rscript}" "$B/src/adapters/run_competitors.R" \
    "$B" "$B/data/manifests/recon.txt" "$v" "$B/results/metrics/recon_${v}_comp.csv" 1 \
    > "$B/logs/recon_${v}_comp.log" 2>&1
  cd "$S" && "$PY" -I "$B/src/adapters/run_studyd.py" \
    "$S" "$B/data/manifests/recon.txt" "$v" "$B/results/metrics/recon_${v}_sd.csv" \
    > "$B/logs/recon_${v}_sd.log" 2>&1
  echo "done recon $v"
}
export -f rjob
printf '%s\n' true unit greedy nj upgma | xargs -P 5 -I{} bash -c 'rjob {}'
echo PHASE3_DONE

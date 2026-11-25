#!/bin/bash
# Real-data benchmark: Simeonov (macsGESTALT) and Quinn (Cassiopeia), per clone.
set -u
# Paths are derived from this script's location so the repository runs anywhere.
# Override with CLDYN_BENCH (the benchmarks/ directory) and CLDYN_ROOT (the repo root).
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
B="${CLDYN_BENCH:-$(cd "$_HERE/../.." && pwd)}"
S="${CLDYN_ROOT:-$(cd "$B/.." && pwd)}"
PY="${CLDYN_PY:-python}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# one manifest per animal so the jobs parallelise and animals stay separable
for a in simeonov_M1 simeonov_M2 quinn_5k quinn_10k quinn_30k quinn_100k; do
  ls -d "$B"/data/real_trees/${a}__*/ 2>/dev/null | sed 's:/$::' > "$B/data/manifests/real_${a}.txt"
done
wc -l "$B"/data/manifests/real_*.txt

run_one() {
  a=$1

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  [ -s "$B/data/manifests/real_${a}.txt" ] || { echo "skip $a (empty)"; return; }
  "${CLDYN_RSCRIPT:-Rscript}" "$B/src/adapters/run_competitors.R" \
     "$B" "$B/data/manifests/real_${a}.txt" unit \
     "$B/results/metrics/real_${a}_comp.csv" 1 > "$B/logs/real_${a}_comp.log" 2>&1
  cd "$S" && "$PY" -I "$B/src/adapters/run_studyd.py" \
     "$S" "$B/data/manifests/real_${a}.txt" unit \
     "$B/results/metrics/real_${a}_sd.csv" > "$B/logs/real_${a}_sd.log" 2>&1
  echo "done $a"
}
export -f run_one
printf '%s\n' simeonov_M1 simeonov_M2 quinn_5k quinn_10k quinn_30k quinn_100k \
  | xargs -P 6 -I{} bash -c 'run_one {}'
echo REAL_DONE

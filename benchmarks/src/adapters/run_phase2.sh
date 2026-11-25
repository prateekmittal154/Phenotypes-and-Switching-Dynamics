#!/bin/bash
# Phase 2: tree-geometry benchmark (II) and real-data benchmark (V).
set -u
# Paths are derived from this script's location so the repository runs anywhere.
# Override with CLDYN_BENCH (the benchmarks/ directory) and CLDYN_ROOT (the repo root).
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
B="${CLDYN_BENCH:-$(cd "$_HERE/../.." && pwd)}"
S="${CLDYN_ROOT:-$(cd "$B/.." && pwd)}"
PY="${CLDYN_PY:-python}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

for g in G_shallow G_deep G_clustered G_dispersed G_starlike; do
  ls -d "$B"/data/sim_geom/${g}_rep*/ 2>/dev/null | sed 's:/$::' > "$B/data/manifests/geom_${g}.txt"
done
for a in simeonov_M1 simeonov_M2 quinn_5k quinn_10k quinn_30k quinn_100k; do
  ls -d "$B"/data/real_trees/${a}__*/ 2>/dev/null | sed 's:/$::' > "$B/data/manifests/real_${a}.txt"
done

geom_job() {
  g=$1

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  "${CLDYN_RSCRIPT:-Rscript}" "$B/src/adapters/run_competitors.R" \
    "$B" "$B/data/manifests/geom_${g}.txt" true "$B/results/metrics/geom_${g}_comp.csv" 1 \
    > "$B/logs/geom_${g}_comp.log" 2>&1
  cd "$S" && "$PY" -I "$B/src/adapters/run_studyd.py" \
    "$S" "$B/data/manifests/geom_${g}.txt" true "$B/results/metrics/geom_${g}_sd.csv" \
    > "$B/logs/geom_${g}_sd.log" 2>&1
  echo "done geom $g"
}
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
export -f geom_job real_job
printf '%s\n' G_shallow G_deep G_clustered G_dispersed G_starlike \
  | xargs -P 5 -I{} bash -c 'geom_job {}'
echo GEOM_DONE
printf '%s\n' simeonov_M1 simeonov_M2 quinn_5k quinn_10k quinn_30k quinn_100k \
  | xargs -P 6 -I{} bash -c 'real_job {}'
echo PHASE2_DONE

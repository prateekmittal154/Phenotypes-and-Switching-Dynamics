#!/bin/bash
# Finalise: run every analysis, rebuild figures and tables, write checksums.
set -u
# Paths are derived from this script's location so the repository runs anywhere.
# Override with CLDYN_BENCH (the benchmarks/ directory) and CLDYN_ROOT (the repo root).
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
B="${CLDYN_BENCH:-$(cd "$_HERE/.." && pwd)}"
S="${CLDYN_ROOT:-$(cd "$B/.." && pwd)}"
PY="${CLDYN_PY:-python}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY="$PY"
cd "$S"

echo "== analyses =="
for a in ground_truth geometry reconstruction real_data holdout; do
  if [ -f "$B/src/metrics/analyse_$a.py" ]; then
    echo "-- $a"
    $PY -I "$B/src/metrics/analyse_$a.py" "$B" > "$B/logs/analyse_$a.out" 2>&1 \
      && echo "   ok" || echo "   FAILED (see logs/analyse_$a.out)"
  fi
done

echo "== static tables =="
$PY -I "$B/src/metrics/make_tableB1.py" "$B" > /dev/null && echo "   B1 ok"
$PY -I "$B/src/metrics/make_tableB6.py" "$B" > /dev/null && echo "   B6 ok"
"${CLDYN_RSCRIPT:-Rscript}" "$B/src/metrics/software_table.R" "$B" > /dev/null 2>&1 \
  && echo "   B2 ok"

echo "== figures =="
$PY -I "$B/src/metrics/make_bench_figures.py" "$B" 2>&1 | grep -E "wrote|failed|skipping"

echo "== checksums =="
mkdir -p "$B/checksums"
cd "$B"
find results protocols reports src environments/install_deps.R README.md -type f \
  ! -name '*.pyc' 2>/dev/null | sort | xargs sha256sum > checksums/BENCHMARK_SHA256.txt
echo "   $(wc -l < checksums/BENCHMARK_SHA256.txt) files checksummed"
echo "FINALISE_DONE"

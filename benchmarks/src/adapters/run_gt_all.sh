#!/bin/bash

# Paths are derived from this script's location so the repository runs anywhere.
# Override with CLDYN_BENCH (the benchmarks/ directory) and CLDYN_ROOT (the repo root).
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
B="${CLDYN_BENCH:-$(cd "$_HERE/../.." && pwd)}"
S="${CLDYN_ROOT:-$(cd "$B/.." && pwd)}"
PY="${CLDYN_PY:-python}"

for v in true unit; do for sc in S1 S2 S3 S4 S5 S6 S7 S8; do echo "$sc $v comp"; done; done \
  | xargs -P 10 -n 3 "$B/src/adapters/worker_gt.sh"
echo COMPETITORS_DONE
for v in true unit; do for sc in S1 S2 S3 S4 S5 S6 S7 S8; do echo "$sc $v sd"; done; done \
  | xargs -P 10 -n 3 "$B/src/adapters/worker_gt.sh"
echo ALL_GT_DONE

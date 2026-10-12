#!/usr/bin/env bash
set -euo pipefail
benchmark_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
analysis_root="$(cd -- "$benchmark_dir/../.." && pwd)"
benchmark_python="$analysis_root/.runtime-tools/venv/bin/python"
if [[ ! -x "$benchmark_python" ]]; then
  echo "Run python3 scripts/benchmark/install_tools.py first (see scripts/benchmark/README.md)." >&2
  exit 1
fi
if [[ "${1:-}" == "--pilot" ]]; then
  shift
  exec "$benchmark_python" "$benchmark_dir/benchmark_direct.py" \
    --output "$analysis_root/data/benchmark/runtime_direct_pilot" \
    --pdbs 1A9N,1DUQ,1FFZ,1KQS --runs 3 --workers 4 --resume "$@"
fi
exec "$benchmark_python" "$benchmark_dir/benchmark_direct.py" \
  --output "$analysis_root/data/benchmark/runtime_direct_all" --workers 4 --resume "$@"

#!/usr/bin/env bash
# Regenerate every PARSE-RNA-Analysis figure by executing its notebook in place.
set -euo pipefail
cd "$(dirname "$0")/../notebooks"
PY="${PARSE_PYTHON:-python}"
fail=0
for nb in */*.ipynb; do
    echo ">>> $nb"
    if ! "$PY" -m jupyter nbconvert --to notebook --execute --inplace "$nb" >/dev/null 2>&1; then
        echo "    FAILED: $nb"; fail=1
    fi
done
echo
if [ "$fail" -eq 0 ]; then
    echo "All figures regenerated into ../figures/"
else
    echo "Some notebooks failed (see above)."; exit 1
fi

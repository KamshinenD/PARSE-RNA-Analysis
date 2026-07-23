#!/usr/bin/env bash
# Fetch the full data/ from NRDStor (HCC) into this repo's data/.
# Requires an HCC account with access to the yesselmanlab allocation.
#
# Usage:
#   bash scripts/fetch_data.sh <your_hcc_username>
#
set -euo pipefail

USER_ID="${1:-${HCC_USER:-}}"
if [ -z "$USER_ID" ]; then
    echo "Usage: bash scripts/fetch_data.sh <your_hcc_username>" >&2
    exit 1
fi

SRC="${USER_ID}@swan.unl.edu:/mnt/nrdstor/yesselmanlab/dewan/PARSE-data/data/"
DEST="$(cd "$(dirname "$0")/.." && pwd)/data/"
mkdir -p "$DEST"

echo "Syncing $SRC -> $DEST"
rsync -avz --partial --progress "$SRC" "$DEST"
echo "Done. All data is in $DEST"

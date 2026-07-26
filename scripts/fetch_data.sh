#!/usr/bin/env bash
# Fetch the archived data/ for PARSE-RNA-Analysis.
#
# Downloads data.zip from Google Drive and unzips it into this repo's data/.
# The archive contains a top-level data/ folder, so it lands in the right place.
#
# Usage:
#   bash scripts/fetch_data.sh
#
set -euo pipefail

# Google Drive file id for data.zip (share link: .../file/d/<FILE_ID>/view).
FILE_ID="${PARSE_DATA_FILE_ID:-1U8nDJfSeA5H59TzdkIEYEDe5-kxcwZcw}"

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ZIP="$REPO_ROOT/data.zip"
DEST="$REPO_ROOT/data"

if [ -d "$DEST" ] && [ -n "$(ls -A "$DEST" 2>/dev/null)" ]; then
    echo "data/ already exists and is non-empty — nothing to do."
    echo "Remove $DEST to re-fetch."
    exit 0
fi

echo "Downloading data.zip from Google Drive..."
if command -v gdown >/dev/null 2>&1; then
    gdown --id "$FILE_ID" -O "$ZIP"
else
    # Fallback: curl through Drive's large-file virus-scan gate. The first
    # request returns an interstitial carrying a per-download uuid; the second
    # request replays it with confirm=t to get the real bytes.
    BASE="https://drive.usercontent.google.com/download?id=${FILE_ID}&export=download"
    UUID="$(curl -sL "$BASE" \
        | grep -o 'name="uuid" value="[^"]*"' | sed 's/.*value="//;s/"//')"
    curl -L "${BASE}&confirm=t&uuid=${UUID}" -o "$ZIP"
fi

echo "Unzipping into $DEST ..."
unzip -q -o "$ZIP" -d "$REPO_ROOT"
rm -f "$ZIP"
echo "Done. All data is in $DEST"

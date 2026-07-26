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

echo "Downloading data.zip from Google Drive..."
if command -v gdown >/dev/null 2>&1; then
    gdown "https://drive.google.com/uc?id=${FILE_ID}" -O "$ZIP"
else
    # Fallback: curl through Drive's large-file virus-scan gate. The first
    # request returns an interstitial carrying a per-download uuid; the second
    # request replays it with confirm=t to get the real bytes.
    BASE="https://drive.usercontent.google.com/download?id=${FILE_ID}&export=download"
    UUID="$(curl -sL "$BASE" \
        | grep -o 'name="uuid" value="[^"]*"' | sed 's/.*value="//;s/"//')"
    curl -L "${BASE}&confirm=t&uuid=${UUID}" -o "$ZIP"
fi

if ! unzip -tq "$ZIP" >/dev/null 2>&1; then
    echo "Downloaded file is not a valid zip (got $(wc -c <"$ZIP") bytes)." >&2
    echo "Check that the Drive link is shared 'Anyone with the link', or that" >&2
    echo "the daily download quota was not exceeded, then retry." >&2
    rm -f "$ZIP"
    exit 1
fi

# -n = never overwrite: only files missing from data/ are extracted, so any
# tables/JSONs regenerated locally (Levels 2-3) are left untouched.
echo "Unzipping into $DEST (existing files are kept) ..."
unzip -qn "$ZIP" -d "$REPO_ROOT"
rm -f "$ZIP"
echo "Done. All data is in $DEST"

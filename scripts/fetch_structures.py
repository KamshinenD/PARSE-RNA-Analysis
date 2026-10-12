#!/usr/bin/env python
"""Download the mmCIF files the analysis runs on into data/structures/.

The PARSE engine downloads structures by PDB id on its own, but the other tools
in the tool comparison and the runtime benchmark (DSSR, FR3D, RNAView, barnaba,
RNApolis, BPNET) read local files. This fetches them from RCSB once, so every
script reads the same copies from data/structures/<ID>.cif.

Usage (from the repo root):
    python scripts/fetch_structures.py              # the 2,697 non-redundant RNAs
    python scripts/fetch_structures.py --all        # every RNA entry in the PDB (~10k)
    python scripts/fetch_structures.py --ids 1GID,1EHZ

The non-redundant set is `unique_pdb_id` in
data/reference/reference_sets/uniqueRNAS.csv; --all uses `pdb_id` in
all_pdb_rnas.csv (regenerate that list with scripts/fetch_all_pdb_rnas.py).
Files already present are skipped, so an interrupted run can simply be repeated.
RCSB occasionally re-issues an entry, so a fresh download may differ slightly from
an older copy; the benchmark records each file's checksum.
"""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
import time
import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
SETS = ROOT / "data" / "reference" / "reference_sets"
OUT = ROOT / "data" / "structures"
URL = "https://files.rcsb.org/download/{}.cif"


def read_ids(path, column):
    with open(path, newline="") as handle:
        return sorted({row[column].strip().upper() for row in csv.DictReader(handle)
                       if row[column].strip()})


def fetch(pdb, out, retries=3):
    """Download one entry; returns (id, 'ok' | 'present' | error text)."""
    target = out / f"{pdb}.cif"
    if target.exists() and target.stat().st_size > 0:
        return pdb, "present"
    part = target.with_suffix(".cif.part")
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(URL.format(pdb), timeout=120) as response:
                part.write_bytes(response.read())
            if not part.read_bytes()[:5] == b"data_":
                raise ValueError("not an mmCIF file")
            part.replace(target)           # never leaves a half-written <ID>.cif
            return pdb, "ok"
        except urllib.error.HTTPError as error:
            if error.code == 404:          # obsolete or withdrawn: retrying won't help
                part.unlink(missing_ok=True)
                return pdb, "not found (404)"
            last = f"HTTP {error.code}"
        except (OSError, ValueError) as error:
            last = f"{type(error).__name__}: {error}"
        time.sleep(2 * attempt)
    part.unlink(missing_ok=True)
    return pdb, last


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--all", action="store_true",
                       help="every RNA-containing PDB entry (all_pdb_rnas.csv), not just the unique set")
    group.add_argument("--ids", help="comma-separated PDB ids instead of a list file")
    ap.add_argument("--output", type=Path, default=OUT, help=f"default: {OUT.relative_to(ROOT)}")
    ap.add_argument("--workers", type=int, default=4, help="parallel downloads (default 4)")
    args = ap.parse_args()

    if args.ids:
        ids = sorted({i.strip().upper() for i in args.ids.split(",") if i.strip()})
    elif args.all:
        ids = read_ids(SETS / "all_pdb_rnas.csv", "pdb_id")
    else:
        ids = read_ids(SETS / "uniqueRNAS.csv", "unique_pdb_id")
    args.output.mkdir(parents=True, exist_ok=True)
    print(f"{len(ids):,} structures -> {args.output}", flush=True)

    counts, failed = {"ok": 0, "present": 0}, []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        for n, (pdb, status) in enumerate(pool.map(lambda i: fetch(i, args.output), ids), 1):
            if status in counts:
                counts[status] += 1
            else:
                failed.append((pdb, status))
            if n % 100 == 0 or n == len(ids):
                print(f"  {n:,}/{len(ids):,}  downloaded {counts['ok']:,}, "
                      f"already present {counts['present']:,}, failed {len(failed)}", flush=True)
    for pdb, status in failed:
        print(f"  FAILED {pdb}: {status}", file=sys.stderr)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

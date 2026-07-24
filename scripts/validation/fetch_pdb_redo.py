"""Download PDB-REDO re-refined coordinates for the validation set.

For each candidate PDB id, fetch ``<id>_final.pdb`` from pdb-redo.eu into
``data/validation/redo_files/<ID>/<ID>.pdb``. IDs without a PDB-REDO entry
(HTTP 404) are skipped — the ids that succeed ARE the validation set that
compare_pdb_redo_cpp.py then scores against the original RCSB coordinates.

Candidates default to the non-redundant RNA set
(reference_sets/uniqueRNAS.csv), so the validation set is reproducible from a
shipped input. Note PDB-REDO is a living database: re-fetching may add/refresh a
few entries versus the archived comparison.
"""

import argparse
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
DEFAULT_LIST = REPO / "data" / "reference" / "reference_sets" / "uniqueRNAS.csv"
DEFAULT_OUT = REPO / "data" / "validation" / "redo_files"
URL = "https://pdb-redo.eu/db/{p}/{p}_final.pdb"


def candidate_ids(list_path: Path) -> list[str]:
    """4-char upper PDB ids from a CSV (unique_pdb_id col) or a plain id list."""
    if list_path.suffix.lower() == ".csv":
        df = pd.read_csv(list_path)
        col = "unique_pdb_id" if "unique_pdb_id" in df.columns else df.columns[0]
        raw = df[col].dropna().astype(str)
    else:
        raw = [ln.strip() for ln in list_path.read_text().splitlines()
               if ln.strip() and not ln.startswith("#")]
    return sorted({str(x).upper()[:4] for x in raw})


def fetch_one(pdb: str, out_dir: Path) -> tuple[str, str]:
    dest = out_dir / pdb / f"{pdb}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return pdb, "cached"
    req = Request(URL.format(p=pdb.lower()), headers={"User-Agent": "PARSE-RNA-Analysis"})
    try:
        with urlopen(req, timeout=60) as r:
            data = r.read()
    except HTTPError as e:
        return pdb, "no_redo" if e.code == 404 else f"http_{e.code}"
    except (URLError, TimeoutError, OSError):
        return pdb, "error"
    if not data.startswith(b"HEADER") and b"ATOM" not in data[:2000]:
        return pdb, "bad"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    return pdb, "ok"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb-list", type=Path, default=DEFAULT_LIST,
                    help="CSV (unique_pdb_id col) or plain id list of candidates")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    ids = candidate_ids(args.pdb_list)
    if args.limit:
        ids = ids[:args.limit]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Candidates: {len(ids)} from {args.pdb_list}\nOutput: {args.out_dir}/")

    counts: dict[str, int] = {}
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for i, (_, status) in enumerate(
                ex.map(lambda p: fetch_one(p, args.out_dir), ids), 1):
            counts[status] = counts.get(status, 0) + 1
            if i % 100 == 0 or i == len(ids):
                print(f"  [{i}/{len(ids)}] {counts}  {i/(time.time()-t0):.1f}/s", flush=True)

    have = counts.get("ok", 0) + counts.get("cached", 0)
    print(f"\nDone: {counts}")
    print(f"  REDO coordinate files available in {args.out_dir}: {have}")


if __name__ == "__main__":
    main()

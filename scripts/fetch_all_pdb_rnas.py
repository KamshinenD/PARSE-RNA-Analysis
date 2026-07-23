#!/usr/bin/env python
"""Fetch a CSV of EVERY RNA-containing structure in the PDB (redundant included).

Uses the public RCSB APIs — no clustering, so this is the full deposited set:
  1. Search API: all entries with >= 1 RNA polymer entity.
  2. Data API (GraphQL): batch-fetch per-entry metadata.

Usage:
    python fetch_all_pdb_rnas.py                 # contains RNA (incl. protein-RNA complexes)
    python fetch_all_pdb_rnas.py --rna-only      # RNA-only (no protein entity)
    python fetch_all_pdb_rnas.py -o all_pdb_rnas.csv

Output columns: pdb_id, title, method, resolution, deposit_date, release_date,
                n_rna_entities, polymer_composition, n_polymer_residues
"""
import argparse
import csv
import json
import sys
import time
import urllib.request

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"


def _post(url, payload):
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def search_ids(rna_only: bool):
    """Return all PDB IDs containing >=1 RNA entity (optionally RNA-only)."""
    rna_term = {
        "type": "terminal", "service": "text",
        "parameters": {
            "attribute": "rcsb_entry_info.polymer_entity_count_RNA",
            "operator": "greater_or_equal", "value": 1},
    }
    query = rna_term
    if rna_only:
        no_protein = {
            "type": "terminal", "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                "operator": "equals", "value": 0},
        }
        query = {"type": "group", "logical_operator": "and",
                 "nodes": [rna_term, no_protein]}
    payload = {
        "query": query,
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    data = _post(SEARCH_URL, payload)
    return [hit["identifier"] for hit in data.get("result_set", [])]


GQL = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info {
      resolution_combined
      polymer_entity_count_RNA
      polymer_composition
      deposited_polymer_monomer_count
    }
    rcsb_accession_info { deposit_date initial_release_date }
  }
}"""


def fetch_meta(ids, batch=250):
    rows = []
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        data = _post(GRAPHQL_URL, {"query": GQL, "variables": {"ids": chunk}})
        for e in data["data"]["entries"] or []:
            info = e.get("rcsb_entry_info") or {}
            acc = e.get("rcsb_accession_info") or {}
            res = info.get("resolution_combined") or []
            rows.append({
                "pdb_id": e["rcsb_id"],
                "title": (e.get("struct") or {}).get("title", ""),
                "method": ";".join(m["method"] for m in (e.get("exptl") or [])),
                "resolution": res[0] if res else "",
                "deposit_date": acc.get("deposit_date", ""),
                "release_date": acc.get("initial_release_date", ""),
                "n_rna_entities": info.get("polymer_entity_count_RNA", ""),
                "polymer_composition": info.get("polymer_composition", ""),
                "n_polymer_residues": info.get("deposited_polymer_monomer_count", ""),
            })
        print(f"  metadata {min(i+batch, len(ids))}/{len(ids)}", file=sys.stderr)
        time.sleep(0.2)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="all_pdb_rnas.csv")
    ap.add_argument("--rna-only", action="store_true",
                    help="only RNA-only structures (exclude protein-RNA complexes)")
    args = ap.parse_args()

    print("Querying RCSB search API for RNA-containing entries ...", file=sys.stderr)
    ids = search_ids(args.rna_only)
    print(f"  {len(ids):,} entries", file=sys.stderr)

    rows = fetch_meta(ids)
    cols = ["pdb_id", "title", "method", "resolution", "deposit_date",
            "release_date", "n_rna_entities", "polymer_composition", "n_polymer_residues"]
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows):,} rows -> {args.output}")


if __name__ == "__main__":
    main()

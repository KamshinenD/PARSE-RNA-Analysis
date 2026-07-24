# PARSE-RNA-Analysis

Notebooks, figures, and the full data-generation pipeline for the PARSE
manuscript. Everything here is reproducible from the public PARSE-RNA engine and
the archived data — from viewing a finished figure down to regenerating every
feature from raw structures.

## Structure

```
notebooks/
  figures/               main-text figure notebooks
  supplemental_figures/  supplemental figure notebooks
figures/                 rendered PNGs
scripts/
  data_gen/              feature + scoring-table regeneration
  validation/            PDB-REDO validation
data/                    fetched from the archive (not in git)
```

## Setup

```bash
pip install -r requirements.txt
```

## Reproducibility levels

Each level rebuilds the inputs of the one above it. Pick the deepest one you need.

### Level 0 — View the figures
The rendered PNGs live in `figures/`, and each executed notebook embeds its
output. Nothing to install or fetch.

### Level 1 — Re-run the figures from the archived data
The `data/` directory is not tracked in git (see [Data](#data)). Fetch it, then
re-execute every notebook:

```bash
bash scripts/fetch_data.sh <your_hcc_username>   # rsync the archive into data/
bash scripts/run_all.sh                          # execute notebooks/*/*.ipynb in place
```

Every figure is rebuilt into `figures/`. Set `PARSE_PYTHON` to choose the
interpreter.

### Level 2 — Rebuild the scoring tables from the feature tables
The scoring tables in `data/reference/scoring_tables/` are derived from the
per-pair / per-residue feature CSVs in `data/reference/high_quality_features/`
(and `low_quality_features/` for the penalty weights). Regenerate them with the
self-contained `data_gen` scripts (no engine needed):

```bash
cd scripts/data_gen
python build_prosco_distributions.py       # -> prosco_distributions.json
python build_z_tables.py                    # -> z_tables.json
python build_backbone_distributions.py      # -> backbone_prosco_distributions.json
python calculate_penalty_weights.py         # -> penalty_weights.json
```

### Level 3 — Regenerate the feature tables from raw structures
The feature CSVs come straight from the PARSE C++ engine. Build the engine,
then re-extract:

```bash
# 1. Build the engine (C++20, CMake >= 3.28)
git clone https://github.com/KamshinenD/PARSE-RNA
cmake -S PARSE-RNA -B PARSE-RNA/build -DCMAKE_BUILD_TYPE=Release
cmake --build PARSE-RNA/build -j
export PARSE_ENGINE="$PWD/PARSE-RNA/build/parse"     # or put `parse` on PATH

# 2. Re-extract every feature table (the engine downloads each mmCIF by id)
cd scripts/data_gen
python extract_pair_features.py        # -> pairs.csv, hbonds.csv, torsions.csv
python extract_backbone_torsions.py    # -> torsions_all.csv
```

The commands above regenerate the GOOD set. The penalty weights also need the
POOR set's pair/hbond features:

```bash
python extract_pair_features.py \
    --reference-set ../../data/reference/reference_sets/reference_set_poor.json \
    --output-dir    ../../data/reference/low_quality_features
```

Then rebuild the tables (Level 2) and the figures (Level 1). Both extractors take
`--limit N` for a quick smoke test and `--engine /path/to/parse` to override
`PARSE_ENGINE`.

#### Validation data (PDB-REDO)
The PDB-REDO comparison behind the validation figure is regenerated the same
way — fetch the re-refined coordinates, then rescore original vs REDO with the
engine:

```bash
cd scripts/validation
python fetch_pdb_redo.py          # uniqueRNAS.csv -> data/validation/redo_files/
python compare_pdb_redo_cpp.py    # score original (RCSB) vs REDO -> the two JSONs
```

`fetch_pdb_redo.py` downloads `<id>_final.pdb` from pdb-redo.eu for every
candidate in `reference_sets/uniqueRNAS.csv` that has a PDB-REDO entry (the rest
are skipped — the ones that succeed *are* the validation set).
`compare_pdb_redo_cpp.py` then writes `pdb_redo_comparison_all.json` and
`pdb_redo_pairs_comparison.json`. PDB-REDO and RCSB are living databases, so a
re-fetch may add or refresh a few entries versus the archived comparison.

## Data

`data/` is not tracked in git (see `.gitignore`), so a fresh clone can **view**
the figures (Level 0) but needs the archive fetched for anything below that
(Levels 1–3). The fetch is a single command — see [Level 1](#level-1--re-run-the-figures-from-the-archived-data):

```bash
bash scripts/fetch_data.sh <your_hcc_username>
```

It rsyncs `/mnt/nrdstor/yesselmanlab/dewan/PARSE-data/data/` (NRDStor, HCC —
requires a yesselmanlab allocation) into `data/`; files already present are
skipped. The archive will later be moved to Zenodo.

### What's in `data/reference/`

| subdir | contents |
|---|---|
| `reference_sets/` | curated GOOD / POOR PDB sets + the RCSB RNA universe |
| `high_quality_features/` | per-pair / per-hbond / per-residue features (GOOD set) |
| `low_quality_features/` | same, for the POOR set (penalty-weight training) |
| `all_uniques_features/` | features over the full non-redundant set |
| `scoring_tables/` | ProSco / Z′ / backbone / penalty-weight tables |
| `scores/`, `metadata/` | per-structure PARSE scores and PDB metadata |

The engine is public at https://github.com/KamshinenD/PARSE-RNA.

# PARSE-RNA-Analysis

Notebooks, figures, and the full data-generation pipeline for the PARSE
manuscript. Everything here is reproducible from the public PARSE-RNA engine and
the archived data — from viewing a finished figure down to regenerating every
feature from raw structures.

## Structure

```
notebooks/
  figures/               figure3.ipynb, figure4.ipynb          (main-text figures)
  supplemental_figures/  figS1..figS5.ipynb                    (supplemental figures)
figures/                 rendered PNGs (also embedded in the executed notebooks)
scripts/
  run_all.sh             re-execute every notebook -> figures/
  fetch_data.sh          rsync the archived data/ into place
  fetch_all_pdb_rnas.py  RCSB query -> reference_sets/all_pdb_rnas.csv
  data_gen/              scoring-table + feature regeneration (see Tier 2 / 3)
    extract_pair_features.py, extract_backbone_torsions.py     (features from the engine)
    build_prosco_distributions.py, build_z_tables.py,
    build_backbone_distributions.py, calculate_penalty_weights.py  (tables from features)
    _engine.py           thin wrapper around the PARSE `parse` binary
    _vendor/             small self-contained helpers (no external deps)
data/                    NOT in git — fetched from the archive (see Data)
requirements.txt
```

## Setup

```bash
pip install -r requirements.txt
```

## Reproducibility tiers

Each tier rebuilds the inputs of the one above it. Pick the deepest one you need.

### Tier 0 — View the figures
The rendered PNGs live in `figures/`, and each executed notebook embeds its
output. Nothing to install or fetch.

### Tier 1 — Re-run the figures from the archived data
Fetch the data (below), then:

```bash
bash scripts/run_all.sh          # executes notebooks/*/*.ipynb in place
```

Every figure is rebuilt into `figures/`. Set `PARSE_PYTHON` to choose the
interpreter.

### Tier 2 — Rebuild the scoring tables from the feature tables
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

### Tier 3 — Regenerate the feature tables from raw structures
The feature CSVs come straight from the PARSE C++ engine — no Python
pair-finder is involved. Build the engine, then re-extract:

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

Then rebuild the tables (Tier 2) and the figures (Tier 1). Both extractors take
`--limit N` for a quick smoke test and `--engine /path/to/parse` to override
`PARSE_ENGINE`.

## Data

Large data files are not tracked in git (see `.gitignore`). A fresh clone can
still **view** the figures (Tier 0); everything below that needs `data/`.

The archived `data/` currently lives on NRDStor (HCC). With a yesselmanlab
allocation, fetch it all in one command:

```bash
bash scripts/fetch_data.sh <your_hcc_username>
```

It rsyncs `/mnt/nrdstor/yesselmanlab/dewan/PARSE-data/data/` into `data/`;
files already present are skipped. (A public Zenodo mirror will replace NRDStor
at publication — only `scripts/fetch_data.sh` changes.)

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

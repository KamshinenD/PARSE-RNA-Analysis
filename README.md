# PARSE-RNA-Analysis

Notebooks, figures, and the full data-generation pipeline for the PARSE
manuscript. Everything here is reproducible from the public PARSE-RNA engine and
the archived data - from viewing a finished figure down to regenerating every
feature from raw structures.

## Structure

```
notebooks/
  figures/               main-text figure notebooks
  supplemental_figures/  supplemental figure notebooks
figures/
  main/                  rendered main-text PNGs
  supplemental/          rendered supplemental PNGs
scripts/
  data_gen/              feature + scoring-table regeneration
  validation/            PDB-REDO validation + tool comparison
  benchmark/             runtime benchmark of PARSE and six other annotators
  fetch_data.sh          download the data archive into data/
  fetch_structures.py    download the mmCIF structures into data/structures/
data/                    fetched from the archive (not in git)
  structures/            mmCIF files, downloaded separately (see Data)
  pymol/                 PyMOL sessions and the images the figure notebooks read
  reference/, validation/, comparison/, benchmark/   analysis inputs and outputs
```

## Setup

```bash
git clone https://github.com/KamshinenD/PARSE-RNA-Analysis.git
cd PARSE-RNA-Analysis
pip install -r requirements.txt
```

## Reproducing the analysis

### Level 0: View the figures
The rendered PNGs live in `figures/`, and each executed notebook embeds its
output. Nothing to install or fetch.

### Level 1: Re-run the figures from the archived data
The `data/` directory is not tracked in git (see [Data](#data)). Fetch it:

```bash
bash scripts/fetch_data.sh
```

Then open the notebooks under `notebooks/figures/` and
`notebooks/supplemental_figures/` and run each one; every figure is written into
`figures/main/` or `figures/supplemental/`.

### Level 2: Rebuild the scoring tables from the feature tables
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

### Level 3: Regenerate the feature tables from raw structures
The feature CSVs come straight from the PARSE C++ engine. Build the engine,
then re-extract:

Run everything below from the repo root (`PARSE-RNA-Analysis/`); the engine is
cloned as a sibling of `scripts/` and `data/`.

```bash
# 1. Build the engine (C++20, CMake >= 3.28): from the repo root
git clone https://github.com/KamshinenD/PARSE-RNA
cmake -S PARSE-RNA -B PARSE-RNA/build -DCMAKE_BUILD_TYPE=Release
cmake --build PARSE-RNA/build -j
export PARSE_ENGINE="$PWD/PARSE-RNA/build/parse"     # or put `parse` on PATH

# 2. Re-extract every feature table (the engine reads data/structures/<ID>.cif;
#    an entry that is not there is downloaded into data/structures/engine_cache/)
cd scripts/data_gen
python extract_pair_features.py        # -> pairs.csv, hbonds.csv, torsions.csv
python extract_backbone_torsions.py    # -> torsions_all.csv
```

The commands above regenerate the curated reference set. The penalty weights also
need the non-curated set's pair/hbond features:

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
way - fetch the re-refined coordinates, then rescore original vs REDO with the
engine:

```bash
cd scripts/validation             # from the repo root
python fetch_pdb_redo.py          # uniqueRNAS.csv -> data/validation/redo_files/
python compare_pdb_redo_cpp.py    # score original (RCSB) vs REDO -> the two JSONs
```

`fetch_pdb_redo.py` downloads `<id>_final.pdb` from pdb-redo.eu for every
candidate in `reference_sets/uniqueRNAS.csv` that has a PDB-REDO entry (the rest
are skipped - the ones that succeed *are* the validation set).
`compare_pdb_redo_cpp.py` then writes `pdb_redo_comparison_all.json` and
`pdb_redo_pairs_comparison.json`. PDB-REDO and RCSB are living databases, so a
re-fetch may add or refresh a few entries versus the archived comparison.

## Data

`data/` is not tracked in git (see `.gitignore`), so a fresh clone can **view**
the figures (Level 0) but needs the archive fetched for anything below that
(Levels 1-3). The fetch is a single command - see [Level 1](#level-1--re-run-the-figures-from-the-archived-data):

```bash
bash scripts/fetch_data.sh
```

It downloads `data.zip` (~400 MB) from Google Drive and unzips it into `data/`.
Extraction is an update rather than a plain overwrite: a file is replaced only
when the archived copy is newer than the one on disk, anything missing is added,
and tables regenerated locally (Levels 2-3) are newer and so survive. Re-running
after the archive is refreshed therefore picks up the new data without
clobbering your own. The archive will later be moved to Zenodo.

### Structures (`data/structures/`)

The mmCIF files are not in the archive (12 GB); they are downloaded straight from
RCSB. Every script reads structures from here; when the engine is given an entry
that is not here it downloads it into `data/structures/engine_cache/`, so nothing
is read from or written to anywhere outside the repository. Fetch them with:

```bash
python scripts/fetch_structures.py            # the 2,697 non-redundant RNAs (~12 GB)
python scripts/fetch_structures.py --all      # every RNA entry in the PDB (~10k)
```

The non-redundant set is the `unique_pdb_id` column of
`data/reference/reference_sets/uniqueRNAS.csv`; `--all` uses
`all_pdb_rnas.csv`, which `scripts/fetch_all_pdb_rnas.py` regenerates from the RCSB
search API. Existing files are skipped, so an interrupted download can be repeated.
RCSB occasionally re-issues an entry, so a new download can differ slightly from
the copy used in the paper; the runtime benchmark records each file's checksum.

### What's in `data/reference/`

| subdir | contents |
|---|---|
| `reference_sets/` | the curated reference set, the non-curated set, the non-redundant list (`uniqueRNAS.csv`) and the RCSB RNA universe |
| `high_quality_features/` | per-pair / per-hbond / per-residue features (curated set) |
| `low_quality_features/` | same, for the non-curated set (penalty-weight derivation) |
| `all_uniques_features/` | features over the full non-redundant set |
| `scoring_tables/` | ProSco / Z′ / backbone / penalty-weight tables |
| `scores/`, `metadata/` | per-structure PARSE scores and PDB metadata; per-pair scores (`pair_scores_uniques_parse.jsonl`, one row per selected pair, each non-redundant entry scored once) are regenerated by `scripts/reference/score_pairs_uniques.py`; nucleotide counts in `metadata_cache.json` are recomputed from the local structure files by `scripts/reference/count_nucleotides.py` (counts every residue of an RNA/DNA polymer entity, so modified, locked and L-nucleotides are included) |
| `resolved_again/` | the re-solved cohort (`resolve_strict.csv`, built by `scripts/validation/build_resolved_cohort.py`) and the per-molecule / per-structure tables behind Figure 5; see its `PROVENANCE.txt` |
| `engine_archive/` | engine output the re-solved cohort builder reads: per-residue and per-pair rows for the non-redundant set and for the extra members of each BGSU equivalence class (read through `scripts/validation/_archive.py`) |

### PyMOL sessions and images (`data/pymol/`)

| subdir | contents |
|---|---|
| `sessions/` | every PyMOL session (`.pse`) behind the structure panels, each saved in its final camera view; `README.txt` maps sessions to panels |
| `renders/` | the ray-traced images the figure notebooks read (Figures 3, 5 and 6, Supplementary Figure S2), plus `7D_labels.json`, the label anchors for Figure 6D |

The notebooks read only `renders/`, never the sessions. The images are rendered
from the sessions by a local script that is not part of the repository.

The engine is public at https://github.com/KamshinenD/PARSE-RNA.

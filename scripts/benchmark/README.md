# Direct command-line runtime comparison

Run from `PARSE-RNA-Analysis`. All seven tools are installed/configured in
`.runtime-tools/`, with the existing licensed DSSR executable. The structures are
read from `data/structures/`; download them first with
`python scripts/fetch_structures.py` (see the main README).

```bash
# Full original Figure 2 timing cohort: 2,096 original mmCIF files.
bash scripts/benchmark/run_all.sh
# On macOS, keep the machine awake:
caffeinate -i bash scripts/benchmark/run_all.sh
```

The runner now uses `benchmark_direct.py`. It times each tool's native command
on the **original, unchanged mmCIF file**. There is no model filtering, chain
renaming, synthetic PDB conversion, grouped annotation, or extra Python worker.
The previous workflow benchmark was stopped and superseded; its pilot and
`runtime_all/` preparation are not publication runtime results. Existing Figure 2
plots and the archived timing CSV remain unchanged.

Results are saved in `data/benchmark/runtime_direct_all/`. Repeating the same
command resumes completed trials. Use an otherwise idle machine. The full run
takes about a day; each trial has a 600-second timeout. `run_all.sh` benchmarks 4
structures at once (`--workers 4`); within a structure, tools run one at a time.
Keep the worker count at or below half the machine's performance cores (the
reference run used 4 on an 8-performance-core Apple M2 Pro). Add a
`.metadata_never_index` file to the output folder on macOS so Spotlight does not
index the logs while the run is timing.

```bash
# Four representative structures: 48, 104, 500, and 2,882 nucleotides.
bash scripts/benchmark/run_all.sh --pilot
# Optional deterministic size-stratified subset:
bash scripts/benchmark/run_all.sh --limit 40 --output data/benchmark/runtime_direct_40
```

## Timing protocol

The stopwatch starts immediately before launching the tool and stops when that
process exits. File parsing and computation performed by the tool itself are
included. Python tools naturally include their own interpreter/import startup.
Our runner's setup, input hashing, symlink creation, output checking, logging,
and summaries are excluded. No structures are downloaded during benchmarking.

Fresh working directories contain symlinks to the original input files so that
tool-generated output stays separate from the source data. These symlinks do
not change any input bytes. There is no up-front input transformation or copy
phase. Original-file checksums are recorded as each structure is benchmarked.

The commands (paths abbreviated) are:

| Mode | Direct command |
| --- | --- |
| PARSE with scoring | `parse input.cif --no-download --cache-dir WORK` |
| PARSE without scoring | same command plus `--no-score` |
| DSSR | `x3dna-dssr -i=input.cif --pair-only` |
| FR3D-Python | `python -m fr3d.classifiers.NA_pairwise_interactions -c basepair -i WORK -o WORK input.cif` |
| RNAView | `rnaview --cif input.cif` |
| barnaba | `barnaba ANNOTATE --pdb input.cif -o annotation` |
| RNApolis | `annotator input.cif --json annotation.json` |
| BPNET | `bpnet input.cif` |

Despite barnaba's option name, its installed MDTraj reader supports mmCIF
natively. RNApolis uses its normal CLI, including secondary-structure output;
its CBC solver dependency is installed. Different programs perform different
additional work and may handle models differently by default. These measurements
represent the stated commands on identical files, not identical internal kernels.
Record the commands in the manuscript.

PARSE and DSSR use the original Figure 2 command options and three measured
repeats. Instead of warm-up runs, each input file is read in full (its checksum),
untimed, just before the structure's first trial, so no tool pays for the disk
read (`--warmups N` restores untimed warm-up runs). Every repeat uses seeded
random tool order. Common BLAS/OpenMP thread settings are capped at one, so each
tool process uses one core; running 4 structures at once changed the median
times of a one-at-a-time pilot by -7% to +6% (one RNAView time on a 48-nt
structure, 0.08 s, by +30%).
The new summary reports medians; the old Figure 2 table used arithmetic means.
Raw repeat times are saved so either statistic can be calculated consistently
for all tools. Rerunning the same commands does not guarantee the old numeric
timings: software builds, machine load, and measurement details can differ.
Do not mix old and new timing protocols in a single speedup calculation.

Native failures, unsupported input, missing/incomplete output, nonzero exits,
and timeouts are recorded without changing the input or retrying through a
conversion pipeline. A tool/structure is eligible for speedups only if every
measured repeat succeeds with consistent annotation output. Valid zero-pair
outputs are allowed. Report coverage alongside speedups.

`speedup = competitor median / PARSE median` for the same structure. Values
above one mean PARSE is faster. Both PARSE modes are supplied as baselines.
Output hashes check repeatability within each tool; pair counts and formats
are not an accuracy comparison across tools.

## Installation and custom inputs

Requirements: macOS or Linux, Python 3.12, Git, CMake, make, GNU C/C++ and
Fortran, the PARSE engine source (`PARSE-RNA/`, cloned as in the main README's Level 3, or `--parse-source`), the mmCIF inputs in
`data/structures/`, and a licensed DSSR executable. On macOS, build prerequisites can be installed with
`brew install python@3.12 cmake gcc`.

```bash
python3 scripts/benchmark/install_tools.py \
  --python /path/to/python3.12 --dssr /path/to/x3dna-dssr
python3 scripts/benchmark/install_tools.py --check
bash scripts/benchmark/run_all.sh --pilot
```

The installer uses pinned Python packages and pinned RNAView/BPNET/FR3D
revisions. NumPy 2.2.6 supplies the API needed by barnaba 0.1.9. PuLP is pinned to
3.3.2, whose bundled CBC is the solver RNApolis uses: PuLP 4 made `LpVariable`
bounds keyword-only, which crashes RNApolis 0.16.17 whenever it calls its
dot-bracket solver (1FFZ and 1KQS in the pilot; small structures need no solve). Short reference-directory aliases accommodate
legacy path-length limits without modifying tool algorithms. Build logs and
configuration stay in `.runtime-tools/`. DSSR is not redistributed.

```bash
.runtime-tools/venv/bin/python scripts/benchmark/benchmark_direct.py \
  --input data/structures \
  --cohort data/benchmark/dssr_PARSE_time_benchmark.csv \
  --output data/benchmark/runtime_direct_custom --runs 3 --warmups 1 --resume
```

The cohort CSV needs `PDB,Nucleotides` columns and each ID's local `PDB.cif`.
The runner verifies files exist immediately, then begins timing without a
structure-transformation phase.

## Outputs and verification

- `trials.jsonl`: status, direct-command wall time, command, original-input hash,
  annotation count/hash, and log path for each trial.
- `per_structure.csv`: median/range for complete, stable repeat sets.
- `speedups.csv`: matched ratios against both PARSE modes.
- `coverage.csv`: eligible sample counts and median speedup versus full PARSE.
- `metadata.json`: protocol, settings, environment, versions, and code/binary hashes.
- `manifest.json`: original input paths, sizes, modification times, and cohort.
- `logs/`: native standard output/error, including failure details.

Changed code, environment, arguments, or inputs require a new output directory.
Trials are flushed individually; a partial final record after interruption is
preserved and retried. A normal interruption also terminates the active tool's
process group.

```bash
.runtime-tools/venv/bin/python -m unittest discover \
  -s scripts/benchmark -p 'test_*.py' -v
```

Tests cover the direct command contract, unchanged input bytes, timing boundaries,
output validation, failure exclusion, matched ratios, interrupted records, and
process cleanup. The direct pilot exercises the actual installed commands.

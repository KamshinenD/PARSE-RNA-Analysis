# PARSE-RNA-Analysis

Notebooks and figures for the PARSE manuscript. Each notebook is self-contained
and produces one figure.

## Structure

```
notebooks/
  figures/               main-figure panels (Figures 3, 4)
  supplemental_figures/  supplementary-figure panels (S1–S5)
figures/    the generated PNG for each notebook
data/       inputs the notebooks read
scripts/    run_all.sh, fetch_all_pdb_rnas.py, data_gen/
```

## Running

Regenerate every figure:

```bash
bash scripts/run_all.sh
```

## Note

Large data files are not tracked in git (see `.gitignore`). A fresh clone can
view the figures — they are embedded in the executed notebooks and saved in
`figures/` — but re-running a notebook needs its data present in `data/`.

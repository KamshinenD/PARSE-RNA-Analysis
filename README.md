# PARSE-RNA-Analysis

Notebooks and figures for the PARSE manuscript. Each notebook is self-contained
and produces one figure.

## Structure

```
notebooks/
  figures/
  supplemental_figures/
figures/
data/
scripts/
```

## Running

Regenerate every figure:

```bash
bash scripts/run_all.sh
```

## Data

Large data files are not tracked in git (see `.gitignore`). A fresh clone can
view the figures — they are embedded in the executed notebooks and saved in
`figures/` — but re-running a notebook needs its data in `data/`.

The full `data/` lives on NRDStor (HCC):

```
/mnt/nrdstor/yesselmanlab/dewan/PARSE-data/data/
```

With an HCC account (yesselmanlab allocation), fetch it all into `data/` with one
command:

```bash
bash scripts/fetch_data.sh <your_hcc_username>
```

It rsyncs the full `data/` from NRDStor; files already present in the clone are
skipped, so only the large files download.

# kranocyte

Characterization of kranocyte and FAP populations in mouse skeletal muscle by scRNA-seq. The
project combines one in-house dataset with four public datasets, all reprocessed uniformly.

## Contents

| Directory | What it holds |
|---|---|
| `workflow/` | FASTQ to matrices: download, STAR index and STARsolo for the five datasets |
| `src/kranocyte/` | the project package: QC, normalization, characterization, DEGs, ontology, interactions, trajectories |
| `notebooks/` | the analysis, numbered: `0*` per-dataset processing, `1AA` characterization, `2*` pseudotemporal dynamics, `3*` interactions, `4` integration, `9` deposit preparation |
| `config/` | configuration and the `samples.tsv` of each dataset |
| `data/` | `raw`, `external`, `interim`, `processed` — outside git |
| `results/`, `figures/` | outputs — outside git |

## Installation

You need Docker or Apptainer configured with NVIDIA support. CellBender uses the GPU (it runs on
CPU, but much more slowly), and STARsolo runs in a container.

```bash
conda env create -f environment.yml      # or: conda-lock install --name kranocyte conda-lock.yml
conda activate kranocyte
pip install -e .
```

> **Warning: the project is pinned to Python 3.11.** This is not a preference. It is the version
> that produced the results, and eight of the dependencies **no longer publish releases
> compatible with 3.11**: `liana` (last 1.10.0, current 2.1.1 — a major version change, with a
> different API), `cellrank` (2.0.7 / 2.4.0), `scanpy` (1.11.5 / 1.13.0), `anndata`
> (0.12.19 / 0.13.4), `scib-metrics`, `numpy`, `scipy` and `networkx`. Installing without the
> upper bounds in `environment.yml` resolves to a different environment from the one that
> produced the figures. **The authoritative record of the environment is the lockfile**, not
> `environment.yml`.

`scranPY` is not on PyPI and is installed from GitHub. **Do not confuse it with `scranpy`**, which
is a different package and is also used. `normalization.py` documents this.

## How to reproduce

There are three entry points, depending on what you want to redo.

### 1. From the FASTQ files — everything

Reproduces the full analysis, including the QC and normalization decisions.

```bash
cd workflow
cp ../config/config.example.env ../config/config.env   # then edit DATA_ROOT
bash 00_setup_common.sh        # GRCm39 reference from NCBI and STAR index
bash 01_oprescu_2020.sh        # ... up to 05_inhouse_ARAUZO_03.sh
```

Then run the notebooks in order, from `0*` to `4`.

You need disk space: the public FASTQ files are **about 790 GB** and must be downloaded again,
because they are not kept. The `--sjdbOverhang 149` index is not on local disk either. The
in-house dataset FASTQ files are on disk, and so is their index (`--sjdbOverhang 89`, because
their R2 reads are 90 nt long).

### 2. From the processed objects — notebook 1 onwards

This is the normal path to reproduce the figures and tables without touching a single FASTQ file.

```bash
# download the Zenodo data record (DOI in PROJECT_CONTEXT.md) to data/processed/
python scripts/make_manifest.py --check     # check against the deposit manifest
jupyter lab notebooks/1AA_cell_type_characterization.ipynb
```

The deposited objects carry `X` (scran + log1p) plus the `counts_raw`, `counts_cellbender`,
`spliced`, `unspliced` and `ambiguous` layers. Only three layers were removed — `norm_scran_log1p`
(identical to `X`), `norm_log1p` and `norm_cpm_log1p` — none of which any notebook reads;
`docs/H5AD_SLIMMING_SPEC.md` gives the per-object table. `prepare_cpdb_input()` needs no argument
change: it detects the missing layer, falls back to `expm1(X)` when `uns['normalization']` records
scran, and warns otherwise.

**What cannot be redone from here** is the comparison between normalizations, since the two
alternative normalized layers are gone. The ambient contamination audit and the velocyto sum
check remain possible, because `counts_raw` and the velocyto layers are kept.

### 3. Browse the results only

The Zenodo code record carries the notebooks **with their outputs**, the final tables and the
paper figures, with their manifest. No installation is needed.

## Where outputs go

- `results/tables/` — tables, with the date in the file name.
- `figures/paper/` — the paper figures (versioned).
- `figures/exploratory/<fecha>/` — everything else (outside git).
- `deposits/` — the bundles produced by notebook 9.

`data/`, `results/` and `figures/` are outside git on purpose. Git keeps only the minimal
reproducible set. The artifacts live in the Zenodo deposits, with checksum manifests that let you
check that a new run gives the same result.

## Data

| Dataset | Source | Platform |
|---|---|---|
| in-house ARAUZO_03 | ENA (accession in preparation) | Illumina NovaSeq 6000 |
| Oprescu 2020 | GSE138826 | Illumina NovaSeq 6000 |
| Song 2023 | GSE215922 | Illumina NovaSeq 6000 |
| Southerland 2023 | GSE227075 | DNBSEQ-G400 |
| Nicoletti 2023 | GSE221736 | NextSeq 500 |

All five were quantified with STARsolo 2.7.11b against the same reference, NCBI RefSeq GRCm39
`GCF_000001635.27`. The matrices published by the original authors are **not** identical to
these: the uniform reprocessing is deliberate.

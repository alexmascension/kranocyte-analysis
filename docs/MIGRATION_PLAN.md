# Restructuring the `kranocito` repository

Detailed outline of the migration to a canonical structure. It includes the classification of
the 4,498 files in the repository and the script that carries it out.

Files that accompany this document:

- [`docs/migration_map.tsv`](migration_map.tsv) — one row per file: `current_path`, `size_bytes`, `modified`, `sha256`, `action`, `reason` and `new_path`.
- [`migrate.sh`](../migrate.sh) — runs the migration by reading that map. Dry run by default.
- [`docs/inventory.tsv`](inventory.tsv) — the raw inventory from which everything above is derived.

---

## 1. Starting point

4,498 files, 363.0 GB. The problem is not the volume. The problem is that the current
structure does not distinguish the irreplaceable from the rebuildable: the in-house FASTQ files,
the intermediate STARsolo matrices and the final objects sit together in the same `data/`, and
the per-dataset folders mix data from two generations of the analysis (v1, from 2021–2023, and
v2, from 2026). A new user cannot tell what can be deleted without losing anything.

Resulting classification:

| `action` | Files | GB | Meaning |
|---|---|---|---|
| `KEEP_raw` | 18 | 28.41 | irreplaceable raw data |
| `KEEP_derived` | 39 | 103.14 | h5ad consumed by the notebooks |
| `KEEP_result` | 2,998 | 3.93 | tables and figures |
| `KEEP_code` / `KEEP_meta` | 71 | 1.72 | code, scripts and metadata |
| `REGEN` | 1,144 | 193.40 | rebuildable from FASTQ or download |
| `ARCHIVE` | 22 | 0.80 | v1 and abandoned notebooks |
| `DELETE` | 206 | 31.59 | orphaned, duplicated or v1-only |

`REGEN` **does not mean delete**. It means that, if space is ever needed, this tier is the first
to be sacrificed, because `workflow/` rebuilds it. The new structure puts each tier in its own
directory, so the path of a file declares its retention policy.

---

## 2. Proposed structure

```text
kranocito/
├── README.md                    # what it is, how to reproduce it, where the data lives
├── CHANGELOG.md
├── TODOS.md
├── pyproject.toml               # makes the `kranocito` package installable
├── environment.yml              # human-readable environment
├── conda-lock.yml               # exact environment
├── .gitignore
│
├── config/
│   ├── config.example.env       # paths; the real config.env is not versioned
│   └── samples/                 # samples.tsv of each dataset, now IN git
│       ├── inhouse_ARAUZO_03.tsv
│       ├── oprescu_2020.tsv  ...
│
├── workflow/                    # FASTQ -> matrices. Everything upstream of the notebooks
│   ├── README.md
│   ├── lib/common.sh
│   ├── 00_setup_common.sh
│   ├── 01_oprescu_2020.sh  ...  04_nicoletti_2023.sh
│   └── 05_inhouse_ARAUZO_03.sh  # formerly run_preprocessing_starsolo.sh
│
├── src/kranocito/               # formerly src/pyfuncs, now an installable package
│   ├── io.py  qc.py  normalization.py  processing.py
│   ├── characterization.py  degs.py  ontology.py
│   ├── interactions.py  trajectory.py  plot_functions.py
│   ├── cell_types.py  dropletQC.py  general.py
│   └── config.py                # replaces common_vars.py
│
├── notebooks/                   # formerly notebooks_v2, without the suffix
│   ├── 0AA_processing_inhouse_dataset.ipynb
│   ├── 0NI / 0OP / 0SO / 0SU_processing_*.ipynb
│   ├── 1AA_cell_type_characterization.ipynb
│   ├── 2AA / 2OP_pseudotemporal_dynamics.ipynb
│   ├── 3AA / 3NI / 3SO / 3SU_interactions.ipynb
│   ├── 4_integration_and_LR_interactions.ipynb
│   ├── 9_prepare_submissions.ipynb   # slimming, manifests, CELLxGENE
│   └── rendered/                # the .html files with outputs
│
├── data/
│   ├── raw/                     # in-house FASTQ. Never regenerated. -> ENA
│   ├── external/                # reference, GTF, whitelists, orthologs. Downloadable
│   ├── interim/                 # STARsolo, CellBender, BAM, logs. Rebuildable
│   └── processed/               # final h5ad per dataset. Kept
│
├── results/
│   ├── tables/                  # xlsx and csv, with the date in the name
│   ├── legacy/                  # v1 results, pending review
│   └── manifest.tsv             # hash and provenance of everything not versioned
│
├── figures/
│   ├── paper/                   # the paper figures. Versioned in git
│   ├── exploratory/<fecha>/     # the rest. Outside git
│   └── legacy/                  # v1 figures
│
├── docs/                        # MIGRATION_PLAN.md, H5AD_SLIMMING_SPEC.md and their tables
├── submission/ena/              # templates, md5 and checklist
└── scripts/
    ├── make_manifest.py         # writes and checks the deposit MANIFEST.tsv
    ├── capture_env.sh
    └── reload_cloudflare.sh
```

### Resulting distribution

| Destination | Files | GB |
|---|---|---|
| `data/interim` | 1,080 | 156.29 |
| `data/processed` | 39 | 103.14 |
| `data/external` | 53 | 36.86 |
| `data/raw` | 18 | 28.41 |
| `notebooks/rendered` | 13 | 1.72 |
| `figures/exploratory` | 2,700 | 1.60 |
| `results/legacy` | 123 | 1.22 |
| `results/tables` | 59 | 1.01 |
| `figures/legacy` | 136 | 0.35 |
| `src/kranocito` | 15 | 0.00 |
| `config/samples` | 5 | 0.00 |
| individual files (13 notebooks, 8 under `workflow/` including `lib/common.sh`, 2 in `submission/ena`, `scripts/reload_cloudflare.sh`, `README.md`, `CHANGELOG.md`, `TODOS.md`, `.gitignore`, `config/config.example.env`) | 29 | 0.00 |
| **Total moved** | **4,270** | **330.6** |

### Why this structure

**`data/` in four tiers (`raw`, `external`, `interim`, `processed`).** This is the Cookiecutter
Data Science convention, which is immediately recognizable. Here it matches the retention
classification exactly: `raw` is all of `KEEP_raw`, `processed` is all of `KEEP_derived`, and
`external` and `interim` are the downloadable and the recomputable. The path of a file is enough
to know whether it can be deleted. Today that requires reading the notebooks.

**`workflow/` separate from `notebooks/`.** The step from FASTQ to matrices is bash, containers
and a cluster. The downstream analysis is Python in Jupyter. These are two different execution
environments and they should look like it. The in-house dataset script becomes
`05_inhouse_ARAUZO_03.sh` and sits next to the four public ones, with the same conventions. That
makes it evident that the five datasets were processed the same way.

**`src/kranocito/` as an installable package.** It removes the `sys.path.append('..')` from the
header of each notebook, which is what currently ties the notebooks to their position in the tree.

**`config/samples/` versioned.** The `samples.tsv` files are metadata, not data: they define what
each GSM is and which condition it belongs to. Today they live under `data/`, which is excluded
from git entirely, so **they are not versioned anywhere**. This is the most serious loss of
traceability in the current repository, and simply moving them fixes it.

**`figures/` split into `paper` and `exploratory`.** This allows versioning the paper figures
without dragging along the 2,700 exploratory ones.

**`results/tables/` with the date in the name, not in the folder.** Fourteen date folders force
you to open them one by one to find the latest version of a table.
`1AA_df_degs_2026-11-18.xlsx` says it on its own, and sorts by name.

---

## 3. What is removed

| Source | Files | GB | `reason` |
|---|---|---|---|
| `data/giordani` | 12 | 16.78 | dataset referenced only by v1 |
| `data/processed` | 31 | 8.99 | dataset referenced only by v1 |
| `data/lung_atlas.h5ad` | 1 | 1.02 | dataset unrelated to the project, no references |
| `data/oprescu` | 5 | 0.91 | dataset referenced only by v1 |
| `data/fan` | 9 | 0.87 | dataset referenced only by v1 |
| `data/demicheli` | 4 | 0.62 | dataset referenced only by v1 |
| `data/mouse_ref` | 15 | 0.52 | dataset referenced only by v1 |
| `data/oprescu-WOT` | 15 | 0.49 | dataset referenced only by v1 |
| `data/leinroth` | 6 | 0.39 | dataset referenced only by v1 |
| `data/scott` | 42 | 0.36 | dataset referenced only by v1 |
| `data/demicheli_mouse` | 5 | 0.36 | dataset referenced only by v1 |
| `data/proietti_2021` | 3 | 0.22 | dataset referenced only by v1 |
| `data/dellorso` | 6 | 0.02 | dataset referenced only by v1 |
| `results/2026-10-24` | 1 | 0.01 | content identical to that of 2026-11-18 (verified part by part) |
| `data/magic_blast` | 13 | 0.01 | dataset referenced only by v1 |
| `data/ronzoni` | 2 | 0.00 | dataset referenced only by v1 |
| `.ipynb_checkpoints/`, `__pycache__/`, `results/DEGs` | 32 | 0.00 | transient artefact |
| `.cache/fonts` | 4 | 0.00 | regenerable font cache |
| **Total** | **206** | **31.59** | |

The per-directory rows above pool the `wget` duplicates (`.1`/`.2` suffixes) with the rest of
their directory: 30 files and 0.48 GB, of which 27 are in `data/scott`.

The fourteen old dataset directories (`giordani`, `processed`, `oprescu`, `fan`,
`demicheli`, `demicheli_mouse`, `leinroth`, `scott`, `proietti_2021`, `dellorso`, `ronzoni`,
`mouse_ref`, `magic_blast`, `oprescu-WOT`) add up to 30.56 GB and **no v2 notebook opens
them**. This was checked by resolving every read and write path of the 13 current notebooks:
they only touch `data/ARAUZO_03/` and `data/public_scRNAseq/datasets/`.

Two specific cases:

- `data/lung_atlas.h5ad` (1.02 GB) is not cited in any notebook or in `pyfuncs`. It is a test
  dataset unrelated to the project.
- `results/2026-10-24/1AA_df_degs.xlsx` is **identical in content** to the 18 November one: its
  16 internal parts match hash for hash, and only the ZIP timestamps differ. In contrast,
  `1AA_df_GOEA.xlsx` did change between the two dates (10 of 16 parts differ), so both versions
  of the ontology table are kept.

**No phase of the script runs `rm`.** What is marked for deletion is moved to `_TO_DELETE/`,
keeping its structure, so it can be reviewed before deleting it by hand.

---

## 4. Retiring v1

The 12 notebooks in `src/notebooks_v1/` and the 4 in `src/src_forgotten_nbs/` are moved aside to
`_ARCHIVE/`. Before that, pin them in the git history:

```bash
git checkout -b legacy/v1 && git push -u origin legacy/v1
git checkout answer-to-reviewers
```

This way the v1 code remains recoverable from a branch without taking up space in the working tree.

Some analyses exist in v1 and **have no equivalent in v2**. Decide on them one by one before
archiving:

| v1 notebook | Content | Covered by v2 |
|---|---|---|
| `4C_..._waddington-OT` | Waddington-OT on Oprescu | partially, by `2OP` |
| `4D_..._cellrank2` | CellRank2 on Oprescu | partially, by `2OP` |
| `6_Supersetting_krano_dataset` | building the superset | no |
| `7_cell_to_cell_interactions` | CellPhoneDB | replaced by LIANA in `3*`/`4` |
| `9_exporting_cellxgene_krano_dataset` | export to the viewer | no |

The cellxgene export and the *supersetting* have no replacement in v2. If they are still needed,
port them before archiving the folder.

`results/legacy/` and `figures/legacy/` collect the outputs of that period without deleting them,
so the decision can be made without haste.

---

## 5. Code changes the migration forces

Moving files breaks paths. These are the exact points:

1. **`src/pyfuncs/common_vars.py`** has `BASE_DIR = "/mnt/data/Proyectos/kranocito"` hard-coded.
   It moves to `src/kranocito/config.py`, which resolves the root from the `KRANOCITO_ROOT`
   environment variable, falling back to the package location.
2. **`config.env`** sets `DATA_ROOT="/scratch/alexmascension/kranocito/data/public_scRNAseq"`,
   the HPC path where the preprocessing ran. It is renamed to `config.example.env` and the real
   one is no longer versioned, which is what the workflow README itself already asks for.
3. **`src/notebooks_v1/...`** also contains
   `program_dir = /media/seth/SETH_DATA/SETH_Alex/Programs/`, a third machine. It goes with v1.
4. **Header of the 13 notebooks**: `sys.path.append('..')` plus
   `from pyfuncs.X import ...` becomes `from kranocito.X import ...` without touching `sys.path`.
5. **`DATA_DIR` of each processing notebook** changes from
   `{BASE_DIR}/data/public_scRNAseq/datasets/<GSE largo>` to `{BASE_DIR}/data/interim/<dataset>`,
   and reading `samples.tsv` changes to `{BASE_DIR}/config/samples/<dataset>.tsv`.
6. **`FIG_DIR`** changes from `figures/{TODAY}/` to `figures/exploratory/{TODAY}/`.
7. **`savefig()`** in `plot_functions.py` defaults to `fig_dir="../figures/"`, relative to the
   notebook. It will be resolved from the package root instead.

---

## 6. Execution order

```bash
# 0. Safety net
cd /mnt/data/Proyectos/kranocito
git status                      # the tree must be clean
git checkout -b legacy/v1 && git push -u origin legacy/v1
git checkout answer-to-reviewers

# 1. Dry run. Touches nothing; writes the details to migration_move_*.log
bash migrate.sh --phase=move

# 2. Move. Source and destination are on the same file system,
#    so these are instant renames: the 330 GB are not copied.
bash migrate.sh --phase=move --apply

# 3. Set aside v1 and abandoned material
bash migrate.sh --phase=archive --apply

# 4. Set aside what is marked for deletion, to _TO_DELETE/
bash migrate.sh --phase=delete --apply

# 5. Review _TO_DELETE/ and only then delete by hand
du -sh _TO_DELETE/ && ls _TO_DELETE/data/

# 6. Apply the code changes from section 5 and check
#    that notebook 0AA starts
```

Each run leaves a `migration_<fase>_<fecha>.log` with the details and an
`undo_<fase>_<fecha>.sh` that reverts the moves **in the reverse order** in which they were made
(tested on a test repository: a round trip restores the original structure). The script refuses
to act if the git tree has uncommitted changes, skips files whose destination already exists, and
can be rerun without additional effects.

---

## 7. Space

| | GB |
|---|---|
| Now | 363.0 |
| After setting aside `DELETE` (31.59) and `ARCHIVE` (0.80) | 330.6 |
| If `data/interim/` (156.29) is also purged | 174.3 |
| If all of `REGEN` is purged — `data/interim/` **and** `data/external/` (193.40 together) | 137.2 |

The last two rows are only advisable once the in-house FASTQ files are deposited in ENA and
confirmed. From that point on, `workflow/` rebuilds `data/interim/` and re-downloads
`data/external/`, and `data/raw/` is backed up off the machine. Rebuilding `data/interim/` also
means re-downloading the ≈790 GB of public FASTQ, so purging it is a one-way door in practice.

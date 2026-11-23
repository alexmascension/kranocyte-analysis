# PROJECT_CONTEXT

Project-specific decisions and constraints. The sections follow the order required by
`guidelines/CODING_GUIDELINES.md` §2. A section that does not apply is kept with `None.`

## 1. Purpose

Characterize the kranocyte and FAP populations of mouse skeletal muscle by scRNA-seq. The project
integrates one in-house control dataset with four public datasets of injury, ischemia and
denervation, all reprocessed uniformly from FASTQ. Details in `README.md`.

## 2. Stack

- **Python 3.11** — pinned, not chosen. Eight dependencies no longer publish compatible releases
  (see the warning in `README.md`). The authoritative record of the environment is the
  **lockfile**.
- Conda environment `kranocyte`, described in `environment.yml` and pinned in the lock.
  Creation: `conda env create -f environment.yml && pip install -e .`
- Main stack: scanpy, anndata, triku, harmonypy, scvelo, cellrank, scFates, palantir, phate,
  liana, CellPhoneDB, gseapy, scranPY + scranpy, rpy2 + bioconductor-scran.
- Containerized external tools: **STARsolo 2.7.11b** (`gcfntnu/star:2.7.11b`, via
  Apptainer) and **CellBender 0.3.2** (`us.gcr.io/broad-dsde-methods/cellbender:0.3.2`, via
  `docker run --gpus all`). Both are recorded in `workflow/containers.tsv` (`HPC-ENV-002`). The
  sha256 of the STARsolo `.sif` is recorded. The CellBender image digest **is pending**, because
  there is no local copy to compute it from.
- The public datasets were preprocessed on the DIPC Hyperion cluster
  (`guidelines/hpc_sites/DIPC_HYPERION.md`).

## 3. Layout

> **Status:** the table describes the target structure of `docs/MIGRATION_PLAN.md`, which has
> **not been applied yet**. Today the package is in `src/pyfuncs/`, the notebooks in
> `src/notebooks_v2/` and the data in `data/ARAUZO_03/` and `data/public_scRNAseq/`.
> `pyfuncs.config.processed_sources()` works with both layouts.

| Path | Contents | Do the guidelines apply? |
|---|---|---|
| `src/kranocyte/` | project package | yes |
| `workflow/` | bash: FASTQ -> matrices | yes (`shellcheck`) |
| `notebooks/` | analysis, numbered | yes, with the exceptions in §22 |
| `config/` | configuration and `samples.tsv` | yes |
| `data/raw` | in-house FASTQ; never regenerated | read-only (`CODE-PROV-001`) |
| `data/external` | reference, GTF, whitelists, orthologs; downloadable | — |
| `data/interim` | STARsolo, CellBender, BAM, logs; rebuildable | — |
| `data/processed` | final h5ad per dataset; kept | — |
| `results/`, `figures/`, `deposits/` | outputs | — |
| `_ARCHIVE/` | retired v1 notebooks | exempt |

## 4. Domain vocabulary

Population names and their palette are declared in `src/kranocyte/cell_types.py`, the single
source: `DICT_MARKERS_*`, `DICT_RENAMING`, `PALETTE_CELL_TYPE` and their integration variants.

Underlying decision: what the project called *kranocytes* was renamed as FAP subpopulations
(`Krano_A` -> `FAP.5`, `Krano_B`/`Krano_C` -> `FAP.6`). The krano/FAP split blurs in the data and
in the lab experiments, so the integration is done a posteriori.

## 5. Data model and invariants

- Atomic unit: one cell from one library (`gsm`), with `condition` and `replicate` taken from
  `config/samples/<dataset>.tsv`.
- `X` of the processed objects is the **scran-normalized, log1p** matrix. The
  `norm_scran_log1p` layer is identical to `X` and is removed on deposit.
- `layers['counts_cellbender']` holds the decontaminated counts. It is the default
  `counts_layer` of `degs`, `characterization` and `interactions`.
- The `spliced`/`unspliced`/`ambiguous` layers come from Velocyto **without decontamination**, on
  purpose: the soup does not have the same spliced/unspliced composition (`io.py`).
- The `*_FAPs*.h5ad` objects are the reprocessed FAP subset, not a view.

## 6. Conventions

- Genome and annotation: **NCBI RefSeq GRCm39, `GCF_000001635.27`**, the unfiltered GTF. The same
  for all five datasets.
- `--sjdbOverhang`: **89** for the in-house dataset (90 nt R2) and **149** for the public ones.
  These are two distinct STAR indices, not an error.
- Quantification: `--soloFeatures GeneFull Velocyto`, `--soloMultiMappers EM`,
  `EmptyDrops_CR`. The notebooks load `GeneFull/raw/UniqueAndMult-EM.mtx`.
- Normalization: scran size factors by deconvolution (Lun et al.), then log1p.
- Seed: `SEED` in `src/kranocyte/config.py`.

## 7. Interfaces

- `data/processed/<dataset>/*_processed*.h5ad` is the **re-entry point boundary**: notebooks
  `1AA` onwards read nothing upstream of it. The Zenodo data deposit publishes exactly those
  objects, slimmed.
- `config/samples/<dataset>.tsv` is consumed by both `workflow/` and the `0*` notebooks.
- `workflow/containers.tsv` is consumed by the workflow README and by notebook 9.

## 8. Overrides

Adopted guidelines and their version: *to be pinned when the guidelines are copied to `guidelines/`.*

| Rule | What is done instead | Reason |
|---|---|---|
| `CODE-NB-003` | Figures go to `figures/paper/` and `figures/exploratory/<fecha>/`, not to `results/figures/<prefijo>_<nombre>.svg` | The project separates paper figures from exploratory ones so that only the former are versioned |
| `CODE-PROV-002` | `config/samples/*.tsv`, `workflow/containers.tsv` and the bundle `MANIFEST.tsv` files are versioned | Exception foreseen by the rule itself: small derived tables, generated by script |

## 9. Decisions log

- **2026-09-13** — Split between the v1 (inherited) and v2 (current) annotations. Reanalysis from
  FASTQ of the four public datasets with STARsolo and a common reference.
- **2026-11-20** — On deposit, the h5ad files are slimmed from 48.65 GB to 32.88 GB (upper bound).
  Only three layers are removed: `norm_scran_log1p` (byte-for-byte identical to `X`), `norm_log1p`
  and `norm_cpm_log1p` (produced by no notebook and consumed by none). `counts_raw` and the
  velocyto layers (`spliced`, `unspliced`, `ambiguous`) are **kept in all eleven objects**: there
  is room for them, and keeping them lets a third party re-audit the ambient correction and rerun
  velocyto. `kranocyte.submissions.SlimSpec` holds the per-object policy.
- **2026-11-20** — Versioning policy: git keeps only the minimal reproducible set. Results and
  figures go to two Zenodo deposits with a checksum manifest.
- **2026-11-20** — Fixed a regression in `workflow/lib/common.sh`: the completion check looked for
  `Solo.out/Gene/raw/matrix.mtx`, which `--soloFeatures GeneFull Velocyto` does not produce. It
  now checks `Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx`, the matrix the notebooks load. The
  versioned copy had never been run end to end.

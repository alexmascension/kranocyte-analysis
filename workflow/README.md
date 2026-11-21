# FASTQ → matrices preprocessing (STARsolo) — five mouse skeletal muscle datasets

This folder (`workflow/`, formerly `src/dataset_download_and_preprocess_scripts/`)
contains the scripts that take the 10x Genomics Single Cell 3′ v3/v3.1 FASTQ files
to the STARsolo count matrices. The five datasets of the project, the in-house
one and the four public ones, are processed:

- with the **same reference**: NCBI RefSeq GRCm39, assembly `GCF_000001635.27`
  (genomic FASTA + unfiltered GTF);
- with the **same STAR version**: STARsolo **2.7.11b**, image
  `gcfntnu/star:2.7.11b`;
- with the same conventions for configuration, container, output layout and
  restart.

There is a single structural difference, and it is deliberate: the in-house
dataset uses **its own STAR index with `--sjdbOverhang 89`** (90 nt R2); the
public ones share one with `--sjdbOverhang 149`. It is explained in the section
[sjdbOverhang](#sjdboverhang-two-distinct-star-indices).

## Included datasets

| Script | Dataset id (`DATASET_ID`) | Source | Samples | Conditions |
|---|---|---|---|---|
| `01_oprescu_2020.sh` | `oprescu_2020_GSE138826` | GEO GSE138826 · BioProject PRJNA577458 · SRA SRP225546 | 7 (GSM4120117–GSM4120123) | The seven time points: uninjured, 0.5, 2, 3.5, 5, 10 and 21 dpi |
| `02_song_2023.sh` | `song_2023_GSE215922` | GEO GSE215922 · PRJNA891295 · SRP403061 | 2: GSM6647486 (sham), GSM6647487 (HLI, day 14) | Sham and HLI day 14 |
| `03_southerland_2023.sh` | `southerland_2023_GSE227075` | GEO GSE227075 · PRJNA943182 · SRP426606 | 4 (GSM7091134–GSM7091137) | C57BL/6 sham and HLI day 1, two replicates each |
| `04_nicoletti_2023.sh` | `nicoletti_2023_GSE221736` | GEO GSE221736 · PRJNA915640 · SRP414775 | 8 (GSM6893974–GSM6893981) | Non-DEN, DEN2d, DEN5d and DEN15d, two replicates each |
| `05_inhouse_ARAUZO_03.sh` | `inhouse_ARAUZO_03` | In-house, sequencing run ARAUZO_03 (local FASTQ). ENA submission in preparation (`submission/ena/`); **no accession yet** | 2: AZ7845, AZ7846 | CTRL, replicates 1 and 2 |

The samples of each dataset are defined in the `SAMPLES` array of its script,
with the format `sample|experiment|condition|replicate`. Each run writes that
list as `DATA_ROOT/datasets/<id>/samples.tsv`. The copies in
`config/samples/<id>.tsv` reproduce that content. If you edit a sample, the
source the scripts use is the `SAMPLES` array.

### Per-dataset details (reads and size)

The figures in this table come from two sources: the STAR `Log.out` files kept
on disk (number of runs per GSM) and the ENA file API, queried on
2026-11-19 (bases per spot, `fastq_bytes`).

| Dataset | Runs per GSM | Platform (ENA) | Bases per spot (R1+R2) | FASTQ.gz in ENA |
|---|---|---|---|---|
| Oprescu 2020 | 1 | Illumina NovaSeq 6000 | 302 | 309.4 GB |
| Song 2023 | 16 (GSM6647486) and 12 (GSM6647487) | Illumina NovaSeq 6000 | 300 | 164.0 GB |
| Southerland 2023 | 1 | DNBSEQ-G400 | 200 (GSM7091134, GSM7091137) and 128 (GSM7091135, GSM7091136) | 190.2 GB |
| Nicoletti 2023 | 8 | NextSeq 500 | 119 | 126.8 GB |
| **Total public** | | | | **790.4 GB** |
| ARAUZO_03 (in-house) | 4 FASTQ pairs per sample: 2 flowcells × 2 lanes (renamed L001–L004; see below) | Illumina NovaSeq 6000 (instruments A01001 and A00500, according to the FASTQ headers) | R1 = 28 nt, R2 = 90 nt (measured in the FASTQ) | 27 GB (16 files) |

ENA only gives the total bases per spot, not the R1/R2 split. The most likely
reading is 2×151 (Oprescu), 2×150 (Song), 2×100 or 28+100 (Southerland) and 28+91
(Nicoletti), but this is an **inference** that cannot be checked against the FASTQ
because they are no longer on disk. In the in-house dataset, AZ7845 has 235.96 M
read pairs, a figure that matches between the facility sheet
(`ARAUZO_03_Sample_Stats.xls`) and `Log.final.out`.

**Sequencer and lane structure of the in-house dataset.** Both samples were
sequenced on **Illumina NovaSeq 6000** (exact term from the ENA controlled
vocabulary). The FASTQ headers give the instruments `A01001` and `A00500`; the
`A` prefix corresponds to NovaSeq 6000. The `L001`–`L004` files are **not** four
lanes of a single flowcell, but **two flowcells with two lanes each**. Renaming
them to `L001`–`L004` flattened that structure:

| Files | Header (instrument:run:flowcell) | Original lane | `MACHINE` in `ARAUZO_03.xls` |
|---|---|---|---|
| `_L001_` | `A01001:472:HH3GYDMXY` | 1 | 02NB |
| `_L002_` | `A01001:472:HH3GYDMXY` | 2 | 02NB |
| `_L003_` | `A00500:534:HH5LLDMXY` | 1 | 01NB |
| `_L004_` | `A00500:534:HH5LLDMXY` | 2 | 01NB |

`data/ARAUZO_03/ARAUZO_03.xls` confirms this with 8 rows: 2 flowcells × 2 lanes
× 2 samples. The structure makes no difference to STARsolo, because the four
pairs of a sample are passed together. It does matter for describing the runs
when submitting to ENA.

**Chemistry (well-founded inference): 10x Genomics Single Cell 3′ v3.1 Dual
Index.** It rests on three facts:

- `ARAUZO_03.xls` assigns dual-index sample indices: SI-TT-E8 (AZ7845)
  and SI-TT-E9 (AZ7846), from the 10x TT family.
- R1 is exactly 28 nt: 16 nt CB + 12 nt UMI, the 10x 3′ v3/v3.1
  structure.
- The 10x v3 whitelist (`3M-february-2018.txt`) is the one that was used.

The commercial name of the kit does not appear literally in the files.

### How the FASTQ files are obtained

**Public.** The scripts resolve, through the ENA file-report API, all the runs
of each SRA experiment (SRX). ENA is the associated European archive that
mirrors the SRA data. The scripts download the paired FASTQ files already
compressed, verify their MD5 and pass all the runs of a GSM to STARsolo
together, in a single invocation. This matters especially for Song (12–16 runs
per GSM) and Nicoletti (8 runs per GSM): it is equivalent to merging the runs
logically, without concatenating files.

**In-house.** The FASTQ files are local: `<sample>_L00N_R{1,2}_001.fastq.gz`,
with four pairs per sample (2 flowcells × 2 lanes, renamed `L001`–`L004`;
see the table above). On each run they are checked against `md5sums.txt`,
which lives in the same directory and plays the role of the ENA MD5s. The four
pairs of a sample are passed to STARsolo together, like the runs of a public
GSM.

## Requirements

- Linux with Bash 4 or later.
- Docker or Apptainer.
- `curl`, `gzip`, `md5sum`, `awk`, `find` and standard GNU utilities.
- Enough disk space and memory for the FASTQ files, the matrices and the
  GRCm39 STAR indices. Values measured in this project:
  - each STAR index takes ≈25 GiB (`SA` alone is 22.3 GB);
  - the public FASTQ files total 790 GB in ENA;
  - the in-house dataset BAMs take 16.6 GB (AZ7845) and 15.3 GB (AZ7846).
  - The public datasets do **not** produce BAMs.

  Check the free space before downloading the full collection.

## Configuration

`config.env` holds local paths and **is not versioned**. Only the template
`config/config.example.env` is versioned:

```bash
cp config/config.example.env config/config.env
nano config/config.env            # define DATA_ROOT, at least
export CONFIG_FILE="$PWD/config/config.env"
```

Configuration points:

- **`CONFIG_FILE`.** By default, `lib/common.sh` looks for `config.env` in the
  parent folder of `lib/`, i.e. `workflow/config.env`. With the configuration in
  `config/`, export `CONFIG_FILE` or create `workflow/config.env`; both are
  local, unversioned files.
- **`DATA_ROOT`.** The scripts never infer the storage location from the working
  directory. `DATA_ROOT` must be defined in `config.env` or exported in the
  shell, and it must be an absolute path.
- **Precedence.** `config.env` is loaded with plain assignments, so any variable
  it defines **overrides** the one exported in the shell. To override something
  from the shell, do not define it in `config.env`.
- **`SJDB_OVERHANG=149`.** This is the default for the shared index of the
  public datasets, chosen for the longest cDNA reads in the collection
  (≈150 nt). It can be changed before building the index. **It does not affect
  `05_inhouse_ARAUZO_03.sh`**, which sets 89 internally.
- **`THREADS`** (default 16). All existing runs used 16.
- **`CURL_RETRIES`** (default 8).
- **`STAR_IMAGE`** (default `gcfntnu/star:2.7.11b`).
- **Optional variables for the in-house dataset:**
  - `INHOUSE_FASTQ_DIR`: absolute path of the directory with the 16 FASTQ files
    and `md5sums.txt`. Default: `DATA_ROOT/datasets/inhouse_ARAUZO_03/fastq`. In
    the new repository structure the FASTQ files live in
    `data/raw/inhouse_ARAUZO_03/`.
  - `INHOUSE_STAR_INDEX_DIR`: absolute path of an already built sjdbOverhang=89
    index, to reuse it instead of
    `DATA_ROOT/common/STAR/mm39_sjdbOverhang_89`. The script checks
    `sjdbOverhang` in `genomeParameters.txt` and aborts if it is not 89.

### Container

`CONTAINER_ENGINE=auto` looks for Docker first and, if it is absent, uses
Apptainer. To force the cluster runtime:

```bash
CONTAINER_ENGINE="apptainer"
```

Both paths use the same image, `gcfntnu/star:2.7.11b`. With Apptainer, the
first run does an `apptainer pull` of `docker://gcfntnu/star:2.7.11b`, converts
it once into `DATA_ROOT/common/containers/star_2.7.11b.sif` and reuses that SIF
for the indices and all the samples. That first run is `00_setup_common.sh`, or
`05_inhouse_ARAUZO_03.sh` if it is launched on its own. `STAR_SIF` changes the
location to an absolute path.

On clusters whose compute nodes have no internet access, run
`00_setup_common.sh` on a login node with internet before submitting the
dataset jobs. If Apptainer is provided through modules, load its module before
running any script. Docker runs with `-u $(id -u):$(id -g)`, so the outputs
belong to the user.

## STARsolo settings

### Common to the five datasets

| Parameter | Value | Reason |
|---|---|---|
| `--soloType` | `CB_UMI_Simple` | 10x 3′: CB and UMI in R1, cDNA in R2 |
| `--readFilesIn` | `R2 R1` (comma-separated list of all runs/lanes) | STAR aligns the first file (cDNA) and reads CB/UMI from the second |
| `--soloCBstart/--soloCBlen` | 1 / 16 | 10x v3/v3.1 cell barcode, bases 1–16 |
| `--soloUMIstart/--soloUMIlen` | 17 / 12 | 10x v3/v3.1 UMI, bases 17–28 |
| `--soloCBwhitelist` | `3M-february-2018.txt` | 10x v3 whitelist (6,794,880 codes of 16 nt; validated on download) |
| `--soloCBmatchWLtype` | `1MM_multi_Nbase_pseudocounts` | Cell Ranger ≥3-style CB correction: allows 1 mismatch and N bases, and resolves CBs compatible with several whitelist codes using pseudocounts |
| `--clipAdapterType` | `CellRanger4` | Trims the TSO at the 5′ end and the poly-A tail at the 3′ end of the cDNA, like Cell Ranger ≥4 |
| `--soloCellFilter` | `EmptyDrops_CR` | Cell calling equivalent to that of Cell Ranger ≥3 (EmptyDrops). Produces `filtered/` in addition to `raw/`, so the filtered matrix is comparable with Cell Ranger's. The full `raw/` matrix is kept for the downstream analysis |
| `--soloMultiMappers` | `EM` | Distributes multimapping reads among genes by EM. `raw/matrix.mtx` holds only the unique reads and `raw/UniqueAndMult-EM.mtx` includes the multimappers |
| `--readFilesCommand` | `zcat` | Compressed FASTQ |
| `--runThreadN` | `THREADS` (16 in all existing runs) | Does not change the matrices |

### Differences between in-house and public (all inherited from the original runs)

| Parameter | Public (`lib/common.sh`, 01–04) | In-house (`05_inhouse_ARAUZO_03.sh`) |
|---|---|---|
| STAR index | `common/STAR/mm39_sjdbOverhang_149` | `common/STAR/mm39_sjdbOverhang_89` |
| `--soloBarcodeReadLength` | `0` | not set (STAR default) |
| `--soloFeatures` | `Gene GeneFull Velocyto` (according to the 21 `Log.out` files; see [Discrepancies](#known-and-pending-discrepancies)) | `Gene GeneFull SJ Velocyto` |
| `--outSAMtype` | `None` (no BAM) | `BAM SortedByCoordinate` |
| `--outSAMattributes` | — | `NH HI nM AS CR UR CB UB GX GN sS sQ sM` |

**`--soloBarcodeReadLength 0` is deliberate for the public datasets.** With the
default value, STAR requires R1 to be exactly CB+UMI (28 nt). Several archived
runs keep R1 beyond the first 28 CB+UMI bases: according to ENA, Oprescu, Song
and part of Southerland have R1 of ≈100–150 nt (inference, see above). With `0`,
STARsolo still reads CB and UMI at the given positions without rejecting those
longer barcode reads. In the in-house dataset, R1 is exactly 28 nt (checked in
the FASTQ), so the original run did not pass this parameter and the script keeps
it that way. To make the behavior explicit, `05` aborts if R1 is not 28 nt.

**Features.**

- `Gene` counts exonic reads, like Cell Ranger.
- `GeneFull` counts exons plus introns (pre-mRNA).
- `Velocyto` gives the spliced/unspliced/ambiguous matrices for RNA velocity.
- `SJ` (in-house only) gives per-cell splice junction counts.

The in-house BAM carries the corrected and uncorrected cell and UMI tags
(`CB/UB`, `CR/UR`) and the gene tags (`GX/GN`), which allow reanalyzing the
reads without realigning.

### Reference

The common reference is the **unfiltered** NCBI RefSeq GRCm39 GTF, not the mm10
or Cell Ranger reference used in the original publications. This is
intentional, so that the reprocessing is uniform, but the published matrices
will not match the ones in this repository numerically.

| File | Source | Observed MD5 (copies on disk) |
|---|---|---|
| `GCF_000001635.27_GRCm39_genomic.fna` | `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/` | `dfb4be77c2031b87848c0de819763b90` |
| `GCF_000001635.27_GRCm39_genomic.gtf` | same | `3a05e4fd87ce5894c6ba8a2c5f3932e0` |
| `3M-february-2018.txt` | `https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz` | `3d36d0a4021fd292b265e2b5e72aaaf3` |

The in-house dataset copy (`data/ARAUZO_03/database/`, downloaded in February
2024) and the public datasets copy (`data/public_scRNAseq/common/reference/`,
downloaded in August 2026) have **identical** MD5s, and so does the whitelist.
Both indices therefore start from the same sequence and the same annotation.
NCBI may update the GTF within the same assembly. If a new download does not
reproduce these MD5s, the reference is no longer the same.

## sjdbOverhang: two distinct STAR indices

STAR recommends `--sjdbOverhang` = maximum read length − 1. The value is fixed
in the index when it is built, so each value implies a different index.

| Index | `--sjdbOverhang` | Used by | Rationale |
|---|---|---|---|
| `mm39_sjdbOverhang_89` | 89 | `05_inhouse_ARAUZO_03.sh` | In-house R2 = 90 nt (measured; the original run computed it with the `awk` line in its script) |
| `mm39_sjdbOverhang_149` | 149 | `01`–`04` | Longest cDNA reads in the collection ≈150 nt (Oprescu, Song) |

**This is neither an error nor an inconsistency:** they are two distinct
indices. All other `genomeGenerate` options are the same in both: same FASTA,
same GTF and same STAR 2.7.11b. The public datasets with shorter R2 (Southerland
≈100 nt and Nicoletti ≈91 nt, inferred from ENA) also use the 149 index. An
overhang larger than read length − 1 is valid in STAR; what was not done was to
build one index per dataset.

`05_inhouse_ARAUZO_03.sh` ignores `SJDB_OVERHANG` from `config.env` (it records
this in the log), builds or reuses the 89 index and checks `sjdbOverhang 89` in
`genomeParameters.txt`.

## Execution order

```bash
# 0. Configuration (once)
cp config/config.example.env config/config.env && nano config/config.env
export CONFIG_FILE="$PWD/config/config.env"

# 1. Public: common inputs + 149 index (once)
bash workflow/00_setup_common.sh

# 2. Public: any dataset, in any order
bash workflow/01_oprescu_2020.sh
bash workflow/02_song_2023.sh
bash workflow/03_southerland_2023.sh
bash workflow/04_nicoletti_2023.sh

# 3. In-house: does not depend on 00. Downloads reference and whitelist if missing
#    (idempotent, the same files as 00), prepares the container,
#    builds or reuses the 89 index and runs STARsolo
bash workflow/05_inhouse_ARAUZO_03.sh
```

Without arguments, each dataset script processes all its samples. To download
and process only some of them, pass them as positional arguments:

```bash
bash workflow/02_song_2023.sh GSM6647486
bash workflow/04_nicoletti_2023.sh GSM6893980 GSM6893981
bash workflow/05_inhouse_ARAUZO_03.sh AZ7846
```

An identifier that does not belong to the dataset aborts before doing anything.

## Output layout

```text
DATA_ROOT/
├── common/
│   ├── reference/GRCm39_NCBI_GCF_000001635.27/   # .fna(.gz), .gtf(.gz)
│   ├── STAR/
│   │   ├── mm39_sjdbOverhang_149/                # public (00_setup_common.sh)
│   │   └── mm39_sjdbOverhang_89/                 # in-house (05), unless INHOUSE_STAR_INDEX_DIR
│   ├── containers/star_2.7.11b.sif               # Apptainer only
│   └── whitelists/3M-february-2018.txt(.gz)
└── datasets/
    ├── <public dataset_id>/
    │   ├── samples.tsv
    │   ├── downloads/<GSM>/
    │   │   ├── ena_file_report.tsv
    │   │   └── fastq/*_{1,2}.fastq.gz
    │   └── results_STAR/<GSM>/
    │       ├── Log.final.out, Log.out, Log.progress.out, SJ.out.tab
    │       └── Solo.out/{Gene,GeneFull,Velocyto}/{raw,filtered}/
    └── inhouse_ARAUZO_03/
        ├── samples.tsv
        ├── fastq/                                # default; or INHOUSE_FASTQ_DIR
        │   ├── <sample>_L00{1..4}_R{1,2}_001.fastq.gz
        │   └── md5sums.txt
        └── results_STAR/<sample>/
            ├── Aligned.sortedByCoord.out.bam
            ├── Log.final.out, Log.out, Log.progress.out, SJ.out.tab
            └── Solo.out/{Gene,GeneFull,SJ,Velocyto}/{raw,filtered}/   # SJ raw/ only
```

Downstream, the notebooks left other files inside each sample's `Solo.out/`:
`adata_raw.h5ad`, `adata_cellbender_velocyto.h5ad` and `CellBender/`.
**These scripts do not produce them.** They do not interfere with the
completed-sample detection either.

**Mapping to the new repository structure.** This layout is the scripts' own;
`DATA_ROOT` can be on scratch. The migration map places its products as follows:

| Pipeline product | Destination in the repository |
|---|---|
| `common/reference/`, `common/whitelists/`, `common/containers/` | `data/external/` |
| In-house FASTQ and `md5sums.txt` | `data/raw/inhouse_ARAUZO_03/` |
| In-house `results_STAR/<sample>/` | `data/interim/inhouse_ARAUZO_03/<sample>/` |
| `samples.tsv` | `config/samples/<id>.tsv` |

The short repository identifiers (`oprescu_2020`, …) do not carry the GSE suffix
that `DATASET_ID` uses in the scripts (`oprescu_2020_GSE138826`, …). To work
directly on the repository, point `INHOUSE_FASTQ_DIR` and
`INHOUSE_STAR_INDEX_DIR` to those paths.

## Restart behavior

All scripts can be rerun:

- **Downloads.** Partial downloads carry a `.part` suffix and `curl` resumes
  them (`--continue-at -`, `CURL_RETRIES` retries). They are only renamed to the
  final name once complete.
- **FASTQ verification.** Already complete FASTQ files are checked again on
  every run: against the ENA MD5 for the public datasets and against
  `md5sums.txt` for the in-house one. A wrong MD5 aborts; move the file aside
  before retrying. For the in-house dataset no re-download is possible: you must
  restore the original. Verifying the 8 FASTQ files of an in-house sample
  (≈13 GB) takes ≈1.5 min.
- **Completed samples.** A sample with complete output is skipped.
  - Public: `Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx`,
    `Solo.out/Velocyto/raw/spliced.mtx` and `Log.final.out` — the same matrix the
    notebooks load (see discrepancy 1 below; the original check named
    `Solo.out/Gene/raw/matrix.mtx`, which `--soloFeatures GeneFull Velocyto`
    never produces).
  - In-house: in addition, `SJ/raw/matrix.mtx` and
    `Aligned.sortedByCoord.out.bam`.
- **Partial outputs.** A STARsolo results directory that is not empty but is
  incomplete **stops the script**, so that a new run is not mixed with a partial
  output. The same happens if STAR fails midway. Move that directory aside and
  run again.
- **STAR indices.** A complete index (`Genome`, `SA`, `SAindex`,
  `genomeParameters.txt`) is reused. An index directory that is not empty but
  incomplete aborts; move it aside before rebuilding. `05` also checks that the
  index `sjdbOverhang` is 89.
- **Apptainer image.** It is generated as `*.part.sif` and only renamed if it is
  not empty. An existing but empty SIF aborts.
- **`samples.tsv`.** It is rewritten atomically on every run from
  `SAMPLES`.

## What is on disk today and what must be re-downloaded

Verified on 2026-11-19 in `/mnt/data/Proyectos/kranocyte/data/`:

| Item | On disk? | Where | To reproduce |
|---|---|---|---|
| Public FASTQ (21 GSM) | **No** (deleted after processing) | — | Re-download from ENA with 01–04: **≈790 GB** (309.4 + 164.0 + 190.2 + 126.8) |
| STAR 149 index (public) | **No** (`public_scRNAseq/common/STAR/` empty) | — | Rebuild with `00_setup_common.sh` |
| GRCm39 reference (fna + gtf, also .gz) | Yes | `public_scRNAseq/common/reference/…` and `ARAUZO_03/database/` (identical MD5s) | Nothing as long as it is kept. The migration map marks it as regenerable; if it is deleted, 00/05 re-download it from NCBI |
| 10x v3 whitelist | Yes (two identical copies) | `public_scRNAseq/common/whitelists/`, `ARAUZO_03/CB_whitelist/` | Nothing |
| `star_2.7.11b.sif` image | Yes (228 MB) | `public_scRNAseq/common/containers/` | Nothing; with Docker, run `docker pull` |
| Public STARsolo results | Yes, 21/21 GSM (`Log.*`, `SJ.out.tab`, `Solo.out/`) | `public_scRNAseq/datasets/<id>/results_STAR/` | — |
| In-house FASTQ | Yes, 16 files, 27 GB + `md5sums.txt` | `ARAUZO_03/FASTQ_renamed/` → `data/raw/inhouse_ARAUZO_03/` | Nothing. Irreplaceable until they are in ENA |
| STAR 89 index (in-house) | **Yes** (`genomeParameters.txt`: `sjdbOverhang 89`; `Log.out`: STAR 2.7.11b) | `ARAUZO_03/database/mm39_STAR/` → `data/external/mm39_STAR/` according to the map | Reuse with `INHOUSE_STAR_INDEX_DIR`, or rebuild with 05 |
| In-house STARsolo results (including BAM) | Yes | `ARAUZO_03/results_STAR/AZ784{5,6}/` | — |

In short: **reproducing the public datasets requires re-downloading ≈790 GB
from ENA and rebuilding the 149 index**. The in-house dataset can be reproduced
today without downloading anything, as long as its FASTQ files and its 89 index
are kept. If the 89 index is deleted, it has to be rebuilt.

## Known and pending discrepancies

Found during unification. Only item 1 has been fixed in the code (in `lib/common.sh`, shared by
00–04); items 2 onwards are recorded but left as they are:

1. **The versioned `lib/common.sh` is a regression of the one that was run, and is FIXED.**
   The HPC copy (`/data/alexmascension/SCRATCH_COPY/`, 21 Aug, 17,231 B) is
   **internally consistent**: it passes `--soloFeatures Gene GeneFull Velocyto` and
   checks `Solo.out/Gene/raw/matrix.mtx`. It is the version that produced the
   matrices. The versioned copy is later and differs in three points:

   | | HPC (run) | git (before the fix) |
   |---|---|---|
   | `--soloFeatures` | `Gene GeneFull Velocyto` | `GeneFull Velocyto`, without updating the check |
   | Download | resume loop, 30 attempts, verifies the `.part` | a single `curl --retry 8` call |
   | `normalise_ena_url` | respects `ftp://` | rewrites everything to `https://` |

   **The versioned version was never run end to end.** With it,
   `run_starsolo` never found `Solo.out/Gene/raw/matrix.mtx`: it neither skipped
   an already completed sample (it always reprocessed) nor passed the final
   check, aborting with "STARsolo did not produce the expected output" after
   hours of compute.

   **The fix applied is not to restore `Gene`.** The notebooks load
   `GeneFull/raw/UniqueAndMult-EM.mtx`: `pyfuncs.io.load_full_adata` defaults to
   `gene_dir="GeneFull"` and `matrix_file="UniqueAndMult-EM.mtx"`, and **none**
   of the ten calls in the five `0*` notebooks overrides them. So
   `--soloFeatures GeneFull Velocyto` is correct, and what was wrong was the
   check. `completed_matrix` becomes
   `Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx`, verified present and non-empty
   in the 23 samples already processed. The HPC resume loop is restored as well,
   which matters with ~790 GB to download.

   Side effect: the analysis **does not use** the `Solo.out/Gene/` directories on
   disk (21 GB in total), and they can be removed.
2. **`config.env` versioned with an HPC path**
   (`DATA_ROOT="/scratch/alexmascension/kranocito/data/public_scRNAseq"`,
   `CONTAINER_ENGINE="apptainer"`). It becomes `config/config.example.env` with
   an empty `DATA_ROOT`; the real `config.env` leaves version control
   (`git rm --cached` + `.gitignore`). The previous README mentioned a
   `config.env.example` that did not exist in the repository.
3. **The original in-house script would no longer work as is**
   (`data/ARAUZO_03/run_preprocessing_starsolo.sh`):
   - it uses an absolute `BASEDIR`;
   - it mounts `$PWD/CB_whitelist` and `$PWD/results_STAR`, so it had to be
     launched from `data/ARAUZO_03/`;
   - it reads the FASTQ files from `20230623/FASTQ_renamed`, today an empty
     directory (the FASTQ files are in `FASTQ_renamed/`);
   - it also downloads the GFF, which is not used.

   `05_inhouse_ARAUZO_03.sh` replaces it with the same STAR parameters,
   compared argument by argument with the original `Log.out` files of the index
   and of AZ7845; only the mount paths inside the container change.
4. **`data/ARAUZO_03/command_run.sh`** (nf-core/scrnaseq with alevin and GRCm38)
   is a **superseded** pipeline. It produced none of the current matrices and
   is archived.
5. **Name of the 89 index in the migration map.** The map moves it to
   `data/external/mm39_STAR/`, a name that does not state the overhang.
   `data/external/STAR/mm39_sjdbOverhang_89/` is recommended, consistent with the
   pipeline layout and so that it is not confused with the 149 index.
6. **ENA accession of the in-house dataset.** It cannot be determined from the
   files: the submission is in preparation. Once it exists, add it to the
   datasets table.

## What the cluster keeps

The copy at `/data/alexmascension/SCRATCH_COPY/kranocyte` on Hyperion holds
what is no longer available locally and avoids rebuilding it:

| Artifact | Where | Why it matters |
|---|---|---|
| STAR index `mm39_sjdbOverhang_149` | `data/public_scRNAseq/common/STAR/` | No need to rebuild it; it is not available locally |
| `star_2.7.11b.sif` | `data/public_scRNAseq/common/containers/` | It is **the exact image** that was run. sha256 `c3c99cfff54ff21d849385af7651c7ecdade7985bef5cafa84fc9d6069b74708`, identical to the local copy. Stronger than pinning a tag |
| GRCm39 reference (FASTA + GTF) | `data/public_scRNAseq/common/reference/` | Downloadable anyway |
| `results_STAR` of the four public datasets | `data/public_scRNAseq/datasets/*/` | 87 GB; backup of `data/interim` |

**FASTQ: none.** Neither in `/data` nor in `/scratch` (that directory no longer
exists). Reproducing the four public datasets from scratch requires
re-downloading the ~790 GB from ENA.

`workflow/containers.tsv` records the project's two images, in line with
`HPC-ENV-002`. The STARsolo one carries its sha256 (the `.sif` is on disk). The
CellBender one (`us.gcr.io/broad-dsde-methods/cellbender:0.3.2`) runs with
`docker run`, there is no `.sif`, and **its digest is still pending**:

```bash
docker image inspect us.gcr.io/broad-dsde-methods/cellbender:0.3.2 \
  --format '{{index .RepoDigests 0}}'
```

Before redistributing the STARsolo `.sif` in a deposit, its terms must be
resolved: its `org.opencontainers.image.licenses` label is **empty** and the
image is built on top of `anaconda/docker-images`.

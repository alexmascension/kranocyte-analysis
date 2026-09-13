# Public mouse skeletal-muscle scRNA-seq preprocessing

These scripts download and uniformly reprocess selected public 10x Genomics
Single Cell 3′ v3/v3.1 datasets with STARsolo 2.7.11b against the same NCBI
GRCm39 reference.

## Included datasets

| Script | GEO series | Samples |
|---|---|---|
| `01_oprescu_2020.sh` | GSE138826 | All seven time points: non-injured, 0.5, 2, 3.5, 5, 10 and 21 dpi |
| `02_song_2023.sh` | GSE215922 | GSM6647486 (Sham) and GSM6647487 (HLI, day 14) |
| `03_southerland_2023.sh` | GSE227075 | GSM7091134–GSM7091137: C57BL/6 sham and HLI day 1, two replicates each |
| `04_nicoletti_2023.sh` | GSE221736 | GSM6893974–GSM6893981: Non-DEN, DEN2d, DEN5d and DEN15d, two replicates each |

The dataset scripts resolve all runs belonging to each SRA experiment (SRX)
through the ENA file-report API. ENA is the European partner archive mirroring
the SRA data. The scripts download its already-compressed paired FASTQs, verify
their MD5 checksums, and pass all runs for one GSM together to STARsolo. This is
particularly important for Song (12–16 runs per GSM) and Nicoletti (8 runs per
GSM).

## Requirements

- Linux with Bash 4 or newer
- Docker or Apptainer
- `curl`, `gzip`, `md5sum`, `awk`, `find` and standard GNU utilities
- Sufficient disk space and memory for the FASTQs, count matrices and the
  GRCm39 STAR index. No BAM is generated. The full collection is still large;
  check the available space before downloading every dataset.

## Configuration

Copy the example configuration and edit the absolute data path:

```bash
cd dataset_download_and_preprocess_scripts
cp config.env.example config.env
nano config.env
```

The scripts never infer storage from the current working directory. `DATA_ROOT`
must be explicitly defined in `config.env` or exported in the shell. The default
index uses `SJDB_OVERHANG=149`, chosen for the longest approximately 150-nt cDNA
reads in this collection. It can be changed before building the index.

`CONTAINER_ENGINE=auto` detects Docker first and otherwise uses Apptainer. To
force the cluster runtime, set:

```bash
CONTAINER_ENGINE="apptainer"
```

Both paths use the same image, `gcfntnu/star:2.7.11b`. With Apptainer,
`00_setup_common.sh` pulls `docker://gcfntnu/star:2.7.11b`, converts it once to
`DATA_ROOT/common/containers/star_2.7.11b.sif`, and reuses that SIF for the
index and all samples. `STAR_SIF` can override this location with an absolute
path. On clusters whose compute nodes have no internet access, run
`00_setup_common.sh` on an internet-enabled login node before submitting the
dataset jobs. If Apptainer is exposed through the modules system, load its
module before running any script.

## Run order

First download the common inputs and build the STAR index once:

```bash
bash 00_setup_common.sh
```

Then run any dataset script:

```bash
bash 01_oprescu_2020.sh
bash 02_song_2023.sh
bash 03_southerland_2023.sh
bash 04_nicoletti_2023.sh
```

With no arguments, a dataset script processes every listed sample. To download
and process only selected GSMs, provide them as positional arguments:

```bash
bash 02_song_2023.sh GSM6647486
bash 04_nicoletti_2023.sh GSM6893980 GSM6893981
```

## Output layout

```text
DATA_ROOT/
├── common/
│   ├── reference/GRCm39_NCBI_GCF_000001635.27/
│   ├── STAR/mm39_sjdbOverhang_149/
│   ├── containers/star_2.7.11b.sif   # Apptainer only
│   └── whitelists/
└── datasets/
    └── <dataset_id>/
        ├── samples.tsv
        ├── downloads/<GSM>/
        │   ├── ena_file_report.tsv
        │   └── fastq/*.fastq.gz
        └── results_STAR/<GSM>/
            ├── Log.final.out
            └── Solo.out/
```

## STARsolo settings

The counting settings reproduce the original pipeline structure:

- 10x 3′ v3/v3.1: CB bases 1–16 and UMI bases 17–28
- v3 whitelist `3M-february-2018.txt`
- `1MM_multi_Nbase_pseudocounts` barcode matching
- `EmptyDrops_CR` cell calling
- `Gene` expression and `Velocyto` spliced/unspliced/ambiguous matrices
- EM assignment of multimappers
- no SAM/BAM output (`--outSAMtype None`)

`--soloBarcodeReadLength 0` is deliberate: several archived runs retain R1
beyond the first 28 CB+UMI bases. STARsolo still reads the CB and UMI from their
specified positions without rejecting those longer barcode reads.

The common reference is the unfiltered NCBI RefSeq GRCm39 GTF, not the mm10 or
Cell Ranger reference used in the original publications. This is intentional
for uniform reprocessing, but published count matrices will therefore not be
numerically identical.

## Restart behaviour

- Partial downloads use a `.part` suffix and are resumed by `curl`.
- Completed FASTQs are rechecked against ENA MD5 values.
- A sample with complete Gene and Velocyto matrices plus `Log.final.out` is
  skipped.
- A non-empty but incomplete STARsolo result directory stops the script rather
  than mixing a new run with partial output. Move that directory aside and rerun.

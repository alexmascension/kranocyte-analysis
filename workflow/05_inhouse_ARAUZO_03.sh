#!/usr/bin/env bash

# In-house dataset ARAUZO_03 (two CTRL replicates, AZ7845 and AZ7846).
#
# Same conventions as 01..04: reads config.env through lib/common.sh, uses the
# same container detection (Docker/Apptainer), the same output layout
# (DATA_ROOT/datasets/<dataset_id>/...) and the same restart behaviour.
#
# Differences with the public datasets, all deliberate and inherited from the
# original run (data/ARAUZO_03/run_preprocessing_starsolo.sh, STAR 2.7.11b,
# August 2025):
#   * FASTQs are local (no ENA download). They are checked against md5sums.txt.
#   * R2 is 90 nt, so this dataset uses its OWN STAR index built with
#     --sjdbOverhang 89 (the public index uses 149). The value is fixed here and
#     ignores SJDB_OVERHANG from config.env.
#   * R1 is exactly 28 nt (16 CB + 12 UMI), so --soloBarcodeReadLength is left
#     at its STAR default, as in the original run.
#   * --soloFeatures Gene GeneFull SJ Velocyto and a coordinate-sorted BAM with
#     CR/UR/CB/UB/GX/GN/sS/sQ/sM attributes are produced, as in the original run.

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

DATASET_ID="inhouse_ARAUZO_03"
# Same record format as the public scripts: sample|experiment|condition|replicate
SAMPLES=(
  'AZ7845|AZ7845|CTRL|1'
  'AZ7846|AZ7846|CTRL|2'
)

# Fixed properties of this dataset. Do not move them to config.env: changing
# any of them means the run is no longer the published preprocessing.
readonly INHOUSE_SJDB_OVERHANG=89      # R2 = 90 nt -> 90 - 1
readonly INHOUSE_R1_LENGTH=28          # 16-nt CB + 12-nt UMI, nothing else
readonly INHOUSE_R2_LENGTH=90
readonly INHOUSE_EXPECTED_LANES=4      # L001..L004 per sample
readonly INHOUSE_MD5_FILE_NAME="md5sums.txt"

# Same sources as 00_setup_common.sh. Keep both scripts in sync.
GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
GTF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz"
WHITELIST_URL="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz"

# config.env (already sourced by lib/common.sh) sets SJDB_OVERHANG for the
# public index. Remember it only to report it, then force the in-house value so
# that initialise_config() resolves STAR_INDEX_DIR to .../mm39_sjdbOverhang_89.
CONFIGURED_SJDB_OVERHANG="${SJDB_OVERHANG:-149}"
SJDB_OVERHANG="${INHOUSE_SJDB_OVERHANG}"

# Optional settings, read from config.env or the environment:
#   INHOUSE_FASTQ_DIR       absolute path to the directory holding
#                           <sample>_L00N_R{1,2}_001.fastq.gz and md5sums.txt.
#                           Default: DATA_ROOT/datasets/inhouse_ARAUZO_03/fastq
#   INHOUSE_STAR_INDEX_DIR  absolute path to an existing sjdbOverhang=89 index
#                           to reuse instead of DATA_ROOT/common/STAR/mm39_sjdbOverhang_89

ensure_inhouse_reference() {
  mkdir -p -- "${REFERENCE_DIR}" "${WHITELIST_DIR}"

  download_file "${GENOME_URL}" "${GENOME_GZ}"
  download_file "${GTF_URL}" "${GTF_GZ}"
  download_file "${WHITELIST_URL}" "${CB_WHITELIST_GZ}"

  decompress_gzip_once "${GENOME_GZ}" "${GENOME_FASTA}"
  decompress_gzip_once "${GTF_GZ}" "${GTF_FILE}"
  decompress_gzip_once "${CB_WHITELIST_GZ}" "${CB_WHITELIST}"

  grep -q '^>' "${GENOME_FASTA}" || die "The genome FASTA is invalid: ${GENOME_FASTA}"
  awk '
    length($0) != 16 || $0 !~ /^[ACGT]+$/ {bad=1; exit}
    END {if (bad || NR < 1000000) exit 1}
  ' "${CB_WHITELIST}" || die "The 10x v3 barcode whitelist failed validation"
}

index_sjdb_overhang() {
  awk -F'\t' '$1 == "sjdbOverhang" {print $2; exit}' "$1/genomeParameters.txt"
}

star_index_complete() {
  [[ -s "${STAR_INDEX_DIR}/Genome" && -s "${STAR_INDEX_DIR}/SA" && \
     -s "${STAR_INDEX_DIR}/SAindex" && -s "${STAR_INDEX_DIR}/genomeParameters.txt" ]]
}

check_index_overhang() {
  local observed
  observed="$(index_sjdb_overhang "${STAR_INDEX_DIR}")"
  [[ "${observed}" == "${INHOUSE_SJDB_OVERHANG}" ]] || \
    die "STAR index ${STAR_INDEX_DIR} was built with sjdbOverhang=${observed:-unknown}; ${DATASET_ID} requires ${INHOUSE_SJDB_OVERHANG}"
}

ensure_inhouse_star_index() {
  mkdir -p -- "${STAR_INDEX_DIR}"

  if star_index_complete; then
    check_index_overhang
    log "STAR index already complete: ${STAR_INDEX_DIR} (sjdbOverhang=${INHOUSE_SJDB_OVERHANG})"
    return 0
  fi

  if [[ -n "$(find "${STAR_INDEX_DIR}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    die "The STAR index directory is non-empty but incomplete: ${STAR_INDEX_DIR}. Move it aside before rebuilding."
  fi

  log "Building GRCm39 STAR index with sjdbOverhang=${INHOUSE_SJDB_OVERHANG} for ${DATASET_ID}"
  run_star_container /index \
    "${REFERENCE_DIR}:/reference:ro" \
    "${STAR_INDEX_DIR}:/index" \
    -- \
    STAR \
      --runThreadN "${THREADS}" \
      --runMode genomeGenerate \
      --genomeDir /index \
      --genomeFastaFiles "/reference/$(basename -- "${GENOME_FASTA}")" \
      --sjdbGTFfile "/reference/$(basename -- "${GTF_FILE}")" \
      --sjdbOverhang "${INHOUSE_SJDB_OVERHANG}"

  star_index_complete || die "STAR did not produce a complete index at ${STAR_INDEX_DIR}"
  check_index_overhang
}

collect_inhouse_lanes() {
  local sample="$1"
  local fastq_dir="$2"
  local index expected_r2

  mapfile -t R1_FILES < <(find "${fastq_dir}" -maxdepth 1 -type f \
    -name "${sample}_L[0-9][0-9][0-9]_R1_001.fastq.gz" -print | sort)
  mapfile -t R2_FILES < <(find "${fastq_dir}" -maxdepth 1 -type f \
    -name "${sample}_L[0-9][0-9][0-9]_R2_001.fastq.gz" -print | sort)

  (( ${#R1_FILES[@]} > 0 )) || die "No Read 1 FASTQs for ${sample} in ${fastq_dir}"
  [[ "${#R1_FILES[@]}" -eq "${#R2_FILES[@]}" ]] || \
    die "Different number of R1 and R2 FASTQs for ${sample} in ${fastq_dir}"
  [[ "${#R1_FILES[@]}" -eq "${INHOUSE_EXPECTED_LANES}" ]] || \
    die "${sample}: expected ${INHOUSE_EXPECTED_LANES} lanes, found ${#R1_FILES[@]}"

  for index in "${!R1_FILES[@]}"; do
    expected_r2="${R1_FILES[${index}]%_R1_001.fastq.gz}_R2_001.fastq.gz"
    [[ "${R2_FILES[${index}]}" == "${expected_r2}" ]] || \
      die "FASTQ pairing mismatch: ${R1_FILES[${index}]} vs ${R2_FILES[${index}]}"
  done
}

# Equivalent of the ENA MD5 check of the public datasets: every FASTQ used for
# a sample must be listed in md5sums.txt and match it, on every run.
verify_inhouse_fastqs() {
  local sample="$1"
  local fastq_dir="$2"
  local md5_file="${fastq_dir}/${INHOUSE_MD5_FILE_NAME}"
  local file filename expected

  [[ -s "${md5_file}" ]] || die "Missing checksum file: ${md5_file}"
  collect_inhouse_lanes "${sample}" "${fastq_dir}"

  for file in "${R1_FILES[@]}" "${R2_FILES[@]}"; do
    filename="$(basename -- "${file}")"
    expected="$(awk -v f="${filename}" '{name=$2; sub(/^\*/, "", name); sub(/^.*\//, "", name)} name == f {print $1; exit}' "${md5_file}")"
    [[ "${expected}" =~ ^[0-9a-f]{32}$ ]] || \
      die "No MD5 entry for ${filename} in ${md5_file}"
    md5_matches "${file}" "${expected}" || \
      die "Existing FASTQ has the wrong MD5: ${file}. Do not rerun until the original file is restored."
  done
  log "${sample}: ${#R1_FILES[@]} FASTQ pair(s) verified against ${INHOUSE_MD5_FILE_NAME}"
}

inhouse_output_complete() {
  local result_dir="$1"
  [[ -s "${result_dir}/Solo.out/Gene/raw/matrix.mtx" && \
     -s "${result_dir}/Solo.out/GeneFull/raw/matrix.mtx" && \
     -s "${result_dir}/Solo.out/SJ/raw/matrix.mtx" && \
     -s "${result_dir}/Solo.out/Velocyto/raw/spliced.mtx" && \
     -s "${result_dir}/Aligned.sortedByCoord.out.bam" && \
     -s "${result_dir}/Log.final.out" ]]
}

run_inhouse_starsolo() {
  local sample="$1"
  local fastq_dir="$2"
  local result_dir="$3"
  local r1_length r2_length r1_csv r2_csv

  collect_inhouse_lanes "${sample}" "${fastq_dir}"
  r1_length="$(first_read_length "${R1_FILES[0]}")"
  r2_length="$(first_read_length "${R2_FILES[0]}")"
  # Without --soloBarcodeReadLength 0, STAR requires R1 == CB + UMI length.
  [[ "${r1_length}" -eq "${INHOUSE_R1_LENGTH}" ]] || \
    die "${sample}: R1 is ${r1_length} nt; this dataset expects exactly ${INHOUSE_R1_LENGTH} nt (16-nt CB + 12-nt UMI)"
  [[ "${r2_length}" -eq "${INHOUSE_R2_LENGTH}" ]] || \
    die "${sample}: R2 is ${r2_length} nt; the sjdbOverhang=${INHOUSE_SJDB_OVERHANG} index assumes ${INHOUSE_R2_LENGTH} nt"
  log "${sample}: ${#R1_FILES[@]} FASTQ pair(s); first lane R1=${r1_length} nt, R2=${r2_length} nt"

  if inhouse_output_complete "${result_dir}"; then
    log "STARsolo output already complete; skipping ${sample}"
    return 0
  fi

  if [[ -d "${result_dir}" ]] && [[ -n "$(find "${result_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    die "Partial/non-empty STARsolo output exists at ${result_dir}. Move it aside before rerunning ${sample}."
  fi
  mkdir -p -- "${result_dir}"

  r1_csv="$(join_container_fastqs "${R1_FILES[@]}")"
  r2_csv="$(join_container_fastqs "${R2_FILES[@]}")"

  # Parameters identical to the original run (see the header). Only the
  # container mount points and the output prefix differ; they do not change
  # the results.
  log "Running STARsolo for ${sample}"
  run_star_container /results \
    "${STAR_INDEX_DIR}:/index:ro" \
    "${WHITELIST_DIR}:/whitelists:ro" \
    "${fastq_dir}:/fastq:ro" \
    "${result_dir}:/results" \
    -- \
    STAR \
      --runThreadN "${THREADS}" \
      --genomeDir /index \
      --readFilesIn "${r2_csv}" "${r1_csv}" \
      --readFilesCommand zcat \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist /whitelists/3M-february-2018.txt \
      --soloCBstart 1 --soloCBlen 16 \
      --soloUMIstart 17 --soloUMIlen 12 \
      --clipAdapterType CellRanger4 \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloCellFilter EmptyDrops_CR \
      --soloFeatures Gene GeneFull SJ Velocyto \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
      --soloMultiMappers EM \
      --outFileNamePrefix /results/

  inhouse_output_complete "${result_dir}" || \
    die "STARsolo did not produce the expected output for ${sample}"
  log "Completed STARsolo for ${sample}"
}

main() {
  local -a requested=("$@")
  local dataset_dir record sample experiment condition replicate result_dir

  initialise_config

  if [[ -n "${INHOUSE_STAR_INDEX_DIR:-}" ]]; then
    [[ "${INHOUSE_STAR_INDEX_DIR}" = /* ]] || \
      die "INHOUSE_STAR_INDEX_DIR must be an absolute path: ${INHOUSE_STAR_INDEX_DIR}"
    STAR_INDEX_DIR="${INHOUSE_STAR_INDEX_DIR}"
    export STAR_INDEX_DIR
  fi
  INHOUSE_FASTQ_DIR="${INHOUSE_FASTQ_DIR:-${DATA_ROOT}/datasets/${DATASET_ID}/fastq}"
  [[ "${INHOUSE_FASTQ_DIR}" = /* ]] || \
    die "INHOUSE_FASTQ_DIR must be an absolute path: ${INHOUSE_FASTQ_DIR}"

  if [[ "${CONFIGURED_SJDB_OVERHANG}" != "${INHOUSE_SJDB_OVERHANG}" ]]; then
    log "SJDB_OVERHANG=${CONFIGURED_SJDB_OVERHANG} from the configuration applies to the public datasets only; ${DATASET_ID} uses ${INHOUSE_SJDB_OVERHANG}"
  fi

  require_cmd curl
  require_cmd gzip
  require_cmd md5sum
  require_cmd awk
  require_cmd find
  validate_requested_samples "${requested[@]}"
  [[ -d "${INHOUSE_FASTQ_DIR}" ]] || \
    die "FASTQ directory not found: ${INHOUSE_FASTQ_DIR}. Set INHOUSE_FASTQ_DIR in config.env."

  ensure_inhouse_reference
  prepare_star_container
  ensure_inhouse_star_index
  check_common_inputs

  dataset_dir="${DATA_ROOT}/datasets/${DATASET_ID}"
  mkdir -p -- "${dataset_dir}"
  write_sample_manifest "${dataset_dir}"

  for record in "${SAMPLES[@]}"; do
    IFS='|' read -r sample experiment condition replicate <<< "${record}"
    sample_requested "${sample}" "${requested[@]}" || continue

    log "Dataset=${DATASET_ID}; sample=${sample}; condition=${condition}; replicate=${replicate}; experiment=${experiment}"
    result_dir="${dataset_dir}/results_STAR/${sample}"
    verify_inhouse_fastqs "${sample}" "${INHOUSE_FASTQ_DIR}"
    run_inhouse_starsolo "${sample}" "${INHOUSE_FASTQ_DIR}" "${result_dir}"
  done
}

# With no arguments, process both samples. To process a subset:
#   bash 05_inhouse_ARAUZO_03.sh AZ7846
main "$@"

#!/usr/bin/env bash

# Shared functions for downloading public 10x 3' v3/v3.1 data and running
# STARsolo. This file is sourced by the executable scripts in the parent folder.

set -Eeuo pipefail

PIPELINE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${PIPELINE_DIR}/config.env}"

if [[ -f "${CONFIG_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${CONFIG_FILE}"
fi

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

log() {
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

initialise_config() {
  [[ -n "${DATA_ROOT:-}" ]] || \
    die "Set DATA_ROOT in ${CONFIG_FILE} or export it before running the script"
  [[ "${DATA_ROOT}" = /* ]] || die "DATA_ROOT must be an absolute path: ${DATA_ROOT}"

  THREADS="${THREADS:-16}"
  CONTAINER_ENGINE="${CONTAINER_ENGINE:-auto}"
  STAR_IMAGE="${STAR_IMAGE:-gcfntnu/star:2.7.11b}"
  SJDB_OVERHANG="${SJDB_OVERHANG:-149}"
  CURL_RETRIES="${CURL_RETRIES:-8}"

  [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
  case "${CONTAINER_ENGINE}" in
    auto|docker|apptainer) ;;
    *) die "CONTAINER_ENGINE must be auto, docker or apptainer: ${CONTAINER_ENGINE}" ;;
  esac
  [[ "${SJDB_OVERHANG}" =~ ^[1-9][0-9]*$ ]] || die "SJDB_OVERHANG must be a positive integer"
  [[ "${CURL_RETRIES}" =~ ^[0-9]+$ ]] || die "CURL_RETRIES must be a non-negative integer"

  COMMON_DIR="${DATA_ROOT}/common"
  CONTAINER_DIR="${COMMON_DIR}/containers"
  REFERENCE_DIR="${COMMON_DIR}/reference/GRCm39_NCBI_GCF_000001635.27"
  WHITELIST_DIR="${COMMON_DIR}/whitelists"
  STAR_INDEX_DIR="${COMMON_DIR}/STAR/mm39_sjdbOverhang_${SJDB_OVERHANG}"
  STAR_SIF="${STAR_SIF:-${CONTAINER_DIR}/star_2.7.11b.sif}"
  [[ "${STAR_SIF}" = /* ]] || die "STAR_SIF must be an absolute path: ${STAR_SIF}"

  GENOME_GZ="${REFERENCE_DIR}/GCF_000001635.27_GRCm39_genomic.fna.gz"
  GTF_GZ="${REFERENCE_DIR}/GCF_000001635.27_GRCm39_genomic.gtf.gz"
  GENOME_FASTA="${GENOME_GZ%.gz}"
  GTF_FILE="${GTF_GZ%.gz}"
  CB_WHITELIST_GZ="${WHITELIST_DIR}/3M-february-2018.txt.gz"
  CB_WHITELIST="${CB_WHITELIST_GZ%.gz}"

  export THREADS CONTAINER_ENGINE STAR_IMAGE STAR_SIF SJDB_OVERHANG CURL_RETRIES
  export COMMON_DIR CONTAINER_DIR REFERENCE_DIR WHITELIST_DIR STAR_INDEX_DIR
  export GENOME_GZ GTF_GZ GENOME_FASTA GTF_FILE CB_WHITELIST_GZ CB_WHITELIST
}

resolve_container_engine() {
  if [[ -n "${RESOLVED_CONTAINER_ENGINE:-}" ]]; then
    return 0
  fi

  case "${CONTAINER_ENGINE}" in
    auto)
      if command -v docker >/dev/null 2>&1; then
        RESOLVED_CONTAINER_ENGINE="docker"
      elif command -v apptainer >/dev/null 2>&1; then
        RESOLVED_CONTAINER_ENGINE="apptainer"
      else
        die "Neither Docker nor Apptainer is available. Load the cluster Apptainer module or install Docker."
      fi
      ;;
    docker)
      require_cmd docker
      RESOLVED_CONTAINER_ENGINE="docker"
      ;;
    apptainer)
      require_cmd apptainer
      RESOLVED_CONTAINER_ENGINE="apptainer"
      ;;
  esac

  export RESOLVED_CONTAINER_ENGINE
  log "Using container engine: ${RESOLVED_CONTAINER_ENGINE}"
}

star_docker_reference() {
  printf '%s\n' "${STAR_IMAGE#docker://}"
}

star_apptainer_reference() {
  if [[ "${STAR_IMAGE}" == docker://* ]]; then
    printf '%s\n' "${STAR_IMAGE}"
  else
    printf 'docker://%s\n' "${STAR_IMAGE}"
  fi
}

prepare_star_container() {
  local partial_sif

  resolve_container_engine
  case "${RESOLVED_CONTAINER_ENGINE}" in
    docker)
      if docker image inspect "$(star_docker_reference)" >/dev/null 2>&1; then
        log "Docker image already present: $(star_docker_reference)"
      else
        log "Pulling Docker image $(star_docker_reference)"
        docker pull "$(star_docker_reference)"
      fi
      ;;
    apptainer)
      if [[ -s "${STAR_SIF}" ]]; then
        log "Apptainer image already present: ${STAR_SIF}"
        return 0
      fi
      [[ ! -e "${STAR_SIF}" ]] || die "Apptainer image exists but is empty or invalid: ${STAR_SIF}"

      mkdir -p -- "$(dirname -- "${STAR_SIF}")"
      partial_sif="${STAR_SIF%.sif}.part.sif"
      log "Converting $(star_apptainer_reference) to reusable SIF: ${STAR_SIF}"
      apptainer pull --force "${partial_sif}" "$(star_apptainer_reference)"
      [[ -s "${partial_sif}" ]] || die "Apptainer did not create a valid SIF: ${partial_sif}"
      mv -- "${partial_sif}" "${STAR_SIF}"
      ;;
  esac
}

run_star_container() {
  local container_workdir="$1"
  shift
  local -a bind_specs=()
  local -a runtime_args=()
  local bind_spec

  while (( $# > 0 )) && [[ "$1" != "--" ]]; do
    bind_specs+=("$1")
    shift
  done
  (( $# > 0 )) || die "Internal error: run_star_container is missing the -- separator"
  shift
  (( $# > 0 )) || die "Internal error: run_star_container received no command"

  resolve_container_engine
  case "${RESOLVED_CONTAINER_ENGINE}" in
    docker)
      for bind_spec in "${bind_specs[@]}"; do
        runtime_args+=(--volume "${bind_spec}")
      done
      docker run --rm -u "$(id -u):$(id -g)" \
        "${runtime_args[@]}" \
        --workdir "${container_workdir}" \
        "$(star_docker_reference)" \
        "$@"
      ;;
    apptainer)
      for bind_spec in "${bind_specs[@]}"; do
        runtime_args+=(--bind "${bind_spec}")
      done
      apptainer exec --cleanenv \
        "${runtime_args[@]}" \
        --pwd "${container_workdir}" \
        "${STAR_SIF}" \
        "$@"
      ;;
  esac
}

download_file() {
  local url="$1"
  local destination="$2"
  local partial="${destination}.part"

  mkdir -p -- "$(dirname -- "${destination}")"

  if [[ -s "${destination}" ]]; then
    log "Already downloaded: ${destination}"
    return 0
  fi
  log "Downloading ${url}"

  curl --fail --location \
    --retry "${CURL_RETRIES}" \
    --retry-delay 5  \
    --continue-at - --output "${partial}" "${url}"
  [[ -s "${partial}" ]] || die "Downloaded file is empty: ${partial}"
  mv -- "${partial}" "${destination}"
}

decompress_gzip_once() {
  local compressed="$1"
  local uncompressed="$2"
  local partial="${uncompressed}.part"

  if [[ -s "${uncompressed}" ]]; then
    log "Already decompressed: ${uncompressed}"
    return 0
  fi

  gzip -t -- "${compressed}" || die "Invalid gzip file: ${compressed}"
  log "Decompressing ${compressed}"
  gzip -dc -- "${compressed}" > "${partial}"
  [[ -s "${partial}" ]] || die "Decompressed file is empty: ${partial}"
  mv -- "${partial}" "${uncompressed}"
}

check_common_inputs() {
  local required
  for required in "${GENOME_FASTA}" "${GTF_FILE}" "${CB_WHITELIST}"; do
    [[ -s "${required}" ]] || die "Missing common input: ${required}. Run 00_setup_common.sh first."
  done

  for required in Genome SA SAindex genomeParameters.txt; do
    [[ -s "${STAR_INDEX_DIR}/${required}" ]] || \
      die "Incomplete STAR index at ${STAR_INDEX_DIR}. Run 00_setup_common.sh first."
  done
}

normalise_ena_url() {
  local remote="$1"
  case "${remote}" in
    https://*) printf '%s\n' "${remote}" ;;
    http://*)  printf 'https://%s\n' "${remote#http://}" ;;
    ftp://*)   printf 'https://%s\n' "${remote#ftp://}" ;;
    *)         printf 'https://%s\n' "${remote}" ;;
  esac
}

md5_matches() {
  local file="$1"
  local expected="$2"
  local observed
  observed="$(md5sum -- "${file}" | awk '{print $1}')"
  [[ "${observed}" == "${expected}" ]]
}

download_ena_fastq() {
  local remote="$1"
  local expected_md5="$2"
  local fastq_dir="$3"
  local url filename destination partial

  url="$(normalise_ena_url "${remote}")"
  filename="$(basename -- "${remote}")"
  destination="${fastq_dir}/${filename}"
  partial="${destination}.part"

  if [[ -s "${destination}" ]]; then
    if md5_matches "${destination}" "${expected_md5}"; then
      log "FASTQ already present and verified: ${filename}"
      return 0
    fi
    die "Existing FASTQ has the wrong MD5: ${destination}. Move it aside and rerun."
  fi

  log "Downloading ${filename}"
  curl --fail --location --show-error \
    --retry "${CURL_RETRIES}" --retry-delay 5 \
    --continue-at - --output "${partial}" "${url}"

  [[ -s "${partial}" ]] || die "Downloaded FASTQ is empty: ${partial}"
  md5_matches "${partial}" "${expected_md5}" || \
    die "MD5 verification failed for ${partial}. Keep the .part file for inspection or move it aside before retrying."
  mv -- "${partial}" "${destination}"
}

download_experiment_fastqs() {
  local experiment="$1"
  local sample_dir="$2"
  local fastq_dir="${sample_dir}/fastq"
  local report="${sample_dir}/ena_file_report.tsv"
  local report_tmp="${report}.part"
  local api="https://www.ebi.ac.uk/ena/portal/api/filereport"
  local run remote_paths remote_md5 remote_bytes layout
  local -a paths md5s
  local index row_count=0

  mkdir -p -- "${fastq_dir}"

  log "Resolving all SRA/ENA runs for ${experiment}"
  curl --fail --location --show-error \
    --retry "${CURL_RETRIES}" --retry-delay 5  \
    --get "${api}" \
    --data-urlencode "accession=${experiment}" \
    --data-urlencode 'result=read_run' \
    --data-urlencode 'fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout' \
    --data-urlencode 'format=tsv' \
    --output "${report_tmp}"

  [[ -s "${report_tmp}" ]] || die "ENA returned an empty report for ${experiment}"
  grep -q $'^run_accession\t' "${report_tmp}" || \
    die "Unexpected ENA report format for ${experiment}: ${report_tmp}"
  mv -- "${report_tmp}" "${report}"

  while IFS=$'\t' read -r run remote_paths remote_md5 remote_bytes layout; do
    [[ -n "${run}" ]] || continue
    layout="${layout%$'\r'}"
    row_count=$((row_count + 1))
    [[ "${layout}" == "PAIRED" ]] || die "${run} is not marked as paired-end in ENA (${layout})"
    [[ -n "${remote_paths}" && -n "${remote_md5}" ]] || \
      die "ENA does not provide FASTQ paths and checksums for ${run}"

    IFS=';' read -r -a paths <<< "${remote_paths}"
    IFS=';' read -r -a md5s <<< "${remote_md5}"
    [[ "${#paths[@]}" -eq 2 ]] || \
      die "Expected exactly two paired FASTQs for ${run}; ENA reported ${#paths[@]}"
    [[ "${#paths[@]}" -eq "${#md5s[@]}" ]] || \
      die "FASTQ/MD5 field count mismatch for ${run}"

    for index in "${!paths[@]}"; do
      download_ena_fastq "${paths[${index}]}" "${md5s[${index}]}" "${fastq_dir}"
    done
  done < <(tail -n +2 -- "${report}")

  (( row_count > 0 )) || die "No sequencing runs were returned for ${experiment}"
  log "Resolved and verified ${row_count} run(s) for ${experiment}"
}

first_read_length() {
  local fastq="$1"
  local length
  length="$({ gzip -cd -- "${fastq}" 2>/dev/null || true; } | awk 'NR == 2 {print length($0); exit}')"
  [[ "${length}" =~ ^[1-9][0-9]*$ ]] || die "Could not determine read length from ${fastq}"
  printf '%s\n' "${length}"
}

collect_fastq_pairs() {
  local fastq_dir="$1"
  mapfile -t R1_FILES < <(find "${fastq_dir}" -maxdepth 1 -type f -name '*_1.fastq.gz' -print | sort)
  mapfile -t R2_FILES < <(find "${fastq_dir}" -maxdepth 1 -type f -name '*_2.fastq.gz' -print | sort)

  (( ${#R1_FILES[@]} > 0 )) || die "No Read 1 FASTQs found in ${fastq_dir}"
  [[ "${#R1_FILES[@]}" -eq "${#R2_FILES[@]}" ]] || \
    die "Different number of R1 and R2 FASTQs in ${fastq_dir}"

  local index r1_stem r2_stem
  for index in "${!R1_FILES[@]}"; do
    r1_stem="$(basename -- "${R1_FILES[${index}]}" _1.fastq.gz)"
    r2_stem="$(basename -- "${R2_FILES[${index}]}" _2.fastq.gz)"
    [[ "${r1_stem}" == "${r2_stem}" ]] || \
      die "FASTQ pairing mismatch: ${R1_FILES[${index}]} vs ${R2_FILES[${index}]}"
  done
}

join_container_fastqs() {
  local file joined="" separator=""
  for file in "$@"; do
    joined+="${separator}/fastq/$(basename -- "${file}")"
    separator=','
  done
  printf '%s\n' "${joined}"
}

run_starsolo() {
  local gsm="$1"
  local fastq_dir="$2"
  local result_dir="$3"
  local completed_matrix="${result_dir}/Solo.out/Gene/raw/matrix.mtx"
  local completed_velocyto="${result_dir}/Solo.out/Velocyto/raw/spliced.mtx"
  local r1_length r2_length r1_csv r2_csv

  collect_fastq_pairs "${fastq_dir}"
  r1_length="$(first_read_length "${R1_FILES[0]}")"
  r2_length="$(first_read_length "${R2_FILES[0]}")"
  (( r1_length >= 28 )) || \
    die "${gsm}: R1 is ${r1_length} nt, shorter than the 16-nt CB + 12-nt UMI required by 10x v3"
  log "${gsm}: ${#R1_FILES[@]} FASTQ pair(s); first run R1=${r1_length} nt, R2=${r2_length} nt"

  if [[ -s "${completed_matrix}" && -s "${completed_velocyto}" && \
        -s "${result_dir}/Log.final.out" ]]; then
    log "STARsolo output already complete; skipping ${gsm}"
    return 0
  fi

  if [[ -d "${result_dir}" ]] && [[ -n "$(find "${result_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    die "Partial/non-empty STARsolo output exists at ${result_dir}. Move it aside before rerunning ${gsm}."
  fi
  mkdir -p -- "${result_dir}"

  r1_csv="$(join_container_fastqs "${R1_FILES[@]}")"
  r2_csv="$(join_container_fastqs "${R2_FILES[@]}")"

  log "Running STARsolo for ${gsm}"
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
      --soloBarcodeReadLength 0 \
      --clipAdapterType CellRanger4 \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloCellFilter EmptyDrops_CR \
      --soloFeatures GeneFull Velocyto \
      --outSAMtype None \
      --soloMultiMappers EM \
      --outFileNamePrefix /results/

  [[ -s "${completed_matrix}" && -s "${completed_velocyto}" && \
     -s "${result_dir}/Log.final.out" ]] || \
    die "STARsolo did not produce the expected output for ${gsm}"
  log "Completed STARsolo for ${gsm}"
}

sample_requested() {
  local gsm="$1"
  shift
  local requested
  (( $# == 0 )) && return 0
  for requested in "$@"; do
    [[ "${requested}" == "${gsm}" ]] && return 0
  done
  return 1
}

validate_requested_samples() {
  local -a requested=("$@")
  local requested_gsm record gsm found
  (( ${#requested[@]} == 0 )) && return 0

  for requested_gsm in "${requested[@]}"; do
    found=0
    for record in "${SAMPLES[@]}"; do
      IFS='|' read -r gsm _ <<< "${record}"
      if [[ "${requested_gsm}" == "${gsm}" ]]; then
        found=1
        break
      fi
    done
    (( found == 1 )) || die "${requested_gsm} is not part of ${DATASET_ID}"
  done
}

write_sample_manifest() {
  local dataset_dir="$1"
  local manifest="${dataset_dir}/samples.tsv"
  local temporary="${manifest}.part"
  local record gsm experiment condition replicate

  {
    printf 'gsm\texperiment\tcondition\treplicate\n'
    for record in "${SAMPLES[@]}"; do
      IFS='|' read -r gsm experiment condition replicate <<< "${record}"
      printf '%s\t%s\t%s\t%s\n' "${gsm}" "${experiment}" "${condition}" "${replicate}"
    done
  } > "${temporary}"
  mv -- "${temporary}" "${manifest}"
}

run_dataset() {
  local -a requested=("$@")

  initialise_config
  local dataset_dir="${DATA_ROOT}/datasets/${DATASET_ID}"
  local record gsm experiment condition replicate sample_dir result_dir
  require_cmd curl
  require_cmd gzip
  require_cmd md5sum
  require_cmd awk
  require_cmd find
  check_common_inputs
  prepare_star_container
  validate_requested_samples "${requested[@]}"

  mkdir -p -- "${dataset_dir}"
  write_sample_manifest "${dataset_dir}"

  for record in "${SAMPLES[@]}"; do
    IFS='|' read -r gsm experiment condition replicate <<< "${record}"
    sample_requested "${gsm}" "${requested[@]}" || continue

    log "Dataset=${DATASET_ID}; sample=${gsm}; condition=${condition}; replicate=${replicate}; experiment=${experiment}"
    sample_dir="${dataset_dir}/downloads/${gsm}"
    result_dir="${dataset_dir}/results_STAR/${gsm}"
    download_experiment_fastqs "${experiment}" "${sample_dir}"
    run_starsolo "${gsm}" "${sample_dir}/fastq" "${result_dir}"
  done
}

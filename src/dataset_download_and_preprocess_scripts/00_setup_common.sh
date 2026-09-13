#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

initialise_config
require_cmd curl
require_cmd gzip
require_cmd awk

GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
GTF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz"
WHITELIST_URL="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz"

mkdir -p -- "${REFERENCE_DIR}" "${WHITELIST_DIR}" "${STAR_INDEX_DIR}"

download_file "${GENOME_URL}" "${GENOME_GZ}"
download_file "${GTF_URL}" "${GTF_GZ}"
download_file "${WHITELIST_URL}" "${CB_WHITELIST_GZ}"

decompress_gzip_once "${GENOME_GZ}" "${GENOME_FASTA}"
decompress_gzip_once "${GTF_GZ}" "${GTF_FILE}"
decompress_gzip_once "${CB_WHITELIST_GZ}" "${CB_WHITELIST}"

grep -q '^>' "${GENOME_FASTA}" || die "The downloaded genome FASTA is invalid"
awk '
  length($0) != 16 || $0 !~ /^[ACGT]+$/ {bad=1; exit}
  END {if (bad || NR < 1000000) exit 1}
' "${CB_WHITELIST}" || die "The 10x v3 barcode whitelist failed validation"

if [[ -s "${STAR_INDEX_DIR}/Genome" && -s "${STAR_INDEX_DIR}/SA" && \
      -s "${STAR_INDEX_DIR}/SAindex" && -s "${STAR_INDEX_DIR}/genomeParameters.txt" ]]; then
  log "STAR index already complete: ${STAR_INDEX_DIR}"
  exit 0
fi

if [[ -n "$(find "${STAR_INDEX_DIR}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  die "The STAR index directory is non-empty but incomplete: ${STAR_INDEX_DIR}. Move it aside before rebuilding."
fi

prepare_star_container

log "Building GRCm39 STAR index with sjdbOverhang=${SJDB_OVERHANG}"
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
    --sjdbOverhang "${SJDB_OVERHANG}"

check_common_inputs
log "Common reference setup completed"

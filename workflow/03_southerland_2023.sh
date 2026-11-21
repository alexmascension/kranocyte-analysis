#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

DATASET_ID="southerland_2023_GSE227075"
SAMPLES=(
  'GSM7091134|SRX19633161|C57BL6_sham|1'
  'GSM7091135|SRX19633162|C57BL6_sham|2'
  'GSM7091136|SRX19633163|C57BL6_HLI_day1|1'
  'GSM7091137|SRX19633164|C57BL6_HLI_day1|2'
)

run_dataset "$@"


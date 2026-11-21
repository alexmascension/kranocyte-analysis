#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

DATASET_ID="oprescu_2020_GSE138826"
SAMPLES=(
  'GSM4120117|SRX6988624|non_injured|1'
  'GSM4120118|SRX6988625|0.5_dpi|1'
  'GSM4120119|SRX6988626|2_dpi|1'
  'GSM4120120|SRX6988627|3.5_dpi|1'
  'GSM4120121|SRX6988628|5_dpi|1'
  'GSM4120122|SRX6988629|10_dpi|1'
  'GSM4120123|SRX6988630|21_dpi|1'
)

# With no arguments, process every sample. To process a subset, pass one or
# more GSM accessions, e.g.:
#   bash 01_oprescu_2020.sh GSM4120117 GSM4120118
run_dataset "$@"


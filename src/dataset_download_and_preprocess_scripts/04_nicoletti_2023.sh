#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

DATASET_ID="nicoletti_2023_GSE221736"
SAMPLES=(
  'GSM6893980|SRX18848366|non_DEN|1'
  'GSM6893981|SRX18848367|non_DEN|2'
  'GSM6893976|SRX18848362|DEN_2d|1'
  'GSM6893977|SRX18848363|DEN_2d|2'
  'GSM6893978|SRX18848364|DEN_5d|1'
  'GSM6893979|SRX18848365|DEN_5d|2'
  'GSM6893974|SRX18848360|DEN_15d|1'
  'GSM6893975|SRX18848361|DEN_15d|2'
)

# Every scRNA-seq replicate contains eight SRA runs. They are merged logically
# by passing all paired FASTQs to one STARsolo invocation per GSM.
run_dataset "$@"


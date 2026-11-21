#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

DATASET_ID="song_2023_GSE215922"
SAMPLES=(
  'GSM6647486|SRX17918111|sham|1'
  'GSM6647487|SRX17918113|HLI_14d|1'
)

# Each GSM contains many sequencing runs. All runs returned for its SRX are
# downloaded, paired and supplied together to a single STARsolo invocation.
run_dataset "$@"


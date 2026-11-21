#!/usr/bin/env bash
# capture_env.sh - freeze the environment that produced the current results.
#
# Run this BEFORE touching anything else, with the analysis environment active:
#
#   conda activate <environment-name>
#   bash capture_env.sh
#
# It leaves four files in env_snapshot_<date>/ :
#   environment.full.yml     full conda dump, with exact builds
#   environment.nobuild.yml  the same without builds, more portable across OSes
#   pip-freeze.txt           everything installed by pip, with exact versions
#   versions.txt             version of python, of R and of the key packages
#
# That dump, and not environment.yml, is the authoritative record of the environment.

set -euo pipefail

[[ -n "${CONDA_DEFAULT_ENV:-}" ]] || { echo "Activate the conda environment first." >&2; exit 1; }

OUT="env_snapshot_$(date +%Y%m%d)"
mkdir -p "$OUT"
echo "environment: $CONDA_DEFAULT_ENV  ->  $OUT/"

conda env export                 > "$OUT/environment.full.yml"
conda env export --no-builds     > "$OUT/environment.nobuild.yml"
pip freeze                       > "$OUT/pip-freeze.txt"

{
  echo "conda environment : $CONDA_DEFAULT_ENV"
  echo "date              : $(date -Iseconds)"
  echo "python            : $(python -c 'import sys; print(sys.version.split()[0])')"
  echo
  echo "--- key packages ---"
  python - <<'PY'
import importlib.metadata as md
for p in ["scanpy","anndata","numpy","scipy","pandas","scvelo","cellrank","liana",
          "triku","scib-metrics","scFates","palantir","phate","harmonypy",
          "cellbender","cellphonedb","doubletdetection","scranpy","rpy2","matplotlib"]:
    try: print(f"{p:18s} {md.version(p)}")
    except md.PackageNotFoundError: print(f"{p:18s} NOT INSTALLED")
PY
  echo
  echo "--- scranPY (git, no version on PyPI) ---"
  pip show scranPY 2>/dev/null | sed -n '1,4p' || echo "scranPY not installed by pip"
  echo
  echo "--- R ---"
  R --version 2>/dev/null | head -1 || echo "R not available on PATH"
  echo
  echo "--- container images: pin by digest, not by tag ---"
  for img in gcfntnu/star:2.7.11b us.gcr.io/broad-dsde-methods/cellbender:0.3.2; do
    docker image inspect "$img" --format '{{.RepoDigests}}' 2>/dev/null \
      | sed "s|^|$img -> |" || echo "$img -> not present locally"
  done
} > "$OUT/versions.txt"

echo "done. Review $OUT/versions.txt and commit it."
#!/usr/bin/env bash
# capture_env.sh - congela el entorno que produjo los resultados actuales.
#
# Ejecutar ANTES de tocar nada mas, con el entorno de analisis activado:
#
#   conda activate <nombre-del-entorno>
#   bash capture_env.sh
#
# Deja cuatro ficheros en env_snapshot_<fecha>/ :
#   environment.full.yml   volcado completo de conda, con builds exactos
#   environment.nobuild.yml  igual pero sin builds, mas portable entre SOs
#   pip-freeze.txt         todo lo instalado por pip, con versiones exactas
#   versions.txt           version de python, de R y de los paquetes clave
#
# Ese volcado, y no environment.yml, es el registro autoritativo del entorno.

set -euo pipefail

[[ -n "${CONDA_DEFAULT_ENV:-}" ]] || { echo "Activa primero el entorno de conda." >&2; exit 1; }

OUT="env_snapshot_$(date +%Y%m%d)"
mkdir -p "$OUT"
echo "entorno: $CONDA_DEFAULT_ENV  ->  $OUT/"

conda env export                 > "$OUT/environment.full.yml"
conda env export --no-builds     > "$OUT/environment.nobuild.yml"
pip freeze                       > "$OUT/pip-freeze.txt"

{
  echo "entorno conda : $CONDA_DEFAULT_ENV"
  echo "fecha         : $(date -Iseconds)"
  echo "python        : $(python -c 'import sys; print(sys.version.split()[0])')"
  echo
  echo "--- paquetes clave ---"
  python - <<'PY'
import importlib.metadata as md
for p in ["scanpy","anndata","numpy","scipy","pandas","scvelo","cellrank","liana",
          "triku","scib-metrics","scFates","palantir","phate","harmonypy",
          "cellbender","cellphonedb","doubletdetection","scranpy","rpy2","matplotlib"]:
    try: print(f"{p:18s} {md.version(p)}")
    except md.PackageNotFoundError: print(f"{p:18s} NO INSTALADO")
PY
  echo
  echo "--- scranPY (git, sin version en PyPI) ---"
  pip show scranPY 2>/dev/null | sed -n '1,4p' || echo "scranPY no instalado por pip"
  echo
  echo "--- R ---"
  R --version 2>/dev/null | head -1 || echo "R no disponible en el PATH"
  echo
  echo "--- imagenes de contenedor: fijar por digest, no por tag ---"
  for img in gcfntnu/star:2.7.11b us.gcr.io/broad-dsde-methods/cellbender:0.3.2; do
    docker image inspect "$img" --format '{{.RepoDigests}}' 2>/dev/null \
      | sed "s|^|$img -> |" || echo "$img -> no presente en local"
  done
} > "$OUT/versions.txt"

echo "listo. Revisa $OUT/versions.txt y confirmalo en git."

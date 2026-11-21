#!/usr/bin/env bash
# migrate.sh - reestructuracion del repositorio kranocito
#
# Por defecto NO modifica nada: imprime lo que haria. Para ejecutar, anadir --apply.
# Las tres fases son independientes y se lanzan por separado.
#
#   bash migrate.sh --phase=move                 # ensayo de los movimientos
#   bash migrate.sh --phase=move    --apply      # mover a la estructura nueva
#   bash migrate.sh --phase=archive --apply      # apartar v1 y notebooks abandonados
#   bash migrate.sh --phase=delete  --apply      # apartar lo marcado para borrar
#
# Ninguna fase ejecuta `rm`. La fase `delete` traslada a _TO_DELETE/ conservando la
# estructura, para poder revisar antes de borrar a mano; `archive` traslada a _ARCHIVE/.
# Cada ejecucion escribe un .log con el detalle y un undo_*.sh que revierte los
# movimientos en orden inverso. Los registros se escriben en el directorio desde el
# que se invoca el script, nunca dentro del repositorio.
#
# Variables: REPO (raiz del repositorio), MAP (ruta de migration_map.tsv).

set -euo pipefail

REPO="${REPO:-/mnt/data/Proyectos/kranocito}"
OUTDIR="$PWD"
MAP="${MAP:-$OUTDIR/migration_map.tsv}"
PHASE="move"; APPLY=0

for arg in "$@"; do
  case "$arg" in
    --phase=*) PHASE="${arg#*=}" ;;
    --apply)   APPLY=1 ;;
    -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
    *) echo "argumento no reconocido: $arg" >&2; exit 2 ;;
  esac
done

MAP="$(readlink -f "$MAP")"
[[ -f "$MAP" ]] || { echo "no encuentro el mapa: $MAP" >&2; exit 1; }
[[ -d "$REPO" ]] || { echo "no encuentro el repositorio: $REPO" >&2; exit 1; }

case "$PHASE" in
  move)    SELECT='($5 ~ /^KEEP/ || $5 == "REGEN") && $7 != "-"'; PREFIX="" ;;
  archive) SELECT='$5 == "ARCHIVE"';                             PREFIX="_ARCHIVE/" ;;
  delete)  SELECT='$5 == "DELETE"';                              PREFIX="_TO_DELETE/" ;;
  *) echo "fase no valida: $PHASE (move|archive|delete)" >&2; exit 2 ;;
esac

TS="$(date +%Y%m%d_%H%M%S)"
LOG="$OUTDIR/migration_${PHASE}_${TS}.log"
UNDO="$OUTDIR/undo_${PHASE}_${TS}.sh"

cd "$REPO"

# Un arbol de git sucio y una migracion masiva no se mezclan.
if [[ $APPLY -eq 1 ]] && [[ -n "$(git status --porcelain 2>/dev/null)" ]]; then
  echo "ERROR: el arbol de git tiene cambios sin confirmar. Confirmalos o guardalos antes." >&2
  exit 1
fi

: > "$LOG"
UNDO_BODY="$(mktemp)"
trap 'rm -f "$UNDO_BODY"' EXIT

n_ok=0; n_skip=0; n_miss=0; n_same=0; bytes=0
while IFS=$'\t' read -r src size mtime sha accion motivo dst; do
  # (el mapa no contiene celdas vacias: migration_map.tsv usa "-" como marcador)
  [[ "$src" == "ruta_actual" ]] && continue
  [[ -n "$PREFIX" ]] && dst="${PREFIX}${src}"
  [[ "$dst" == "-" || -z "$dst" ]] && continue

  if [[ ! -e "$src" ]]; then
    if [[ -e "$dst" ]]; then n_skip=$((n_skip+1))
    else n_miss=$((n_miss+1)); echo "AUSENTE            $src" >> "$LOG"; fi
    continue
  fi
  if [[ "$src" == "$dst" ]]; then
    echo "SIN_CAMBIO         $src" >> "$LOG"; n_same=$((n_same+1)); continue
  fi
  if [[ -e "$dst" ]]; then
    echo "DESTINO_YA_EXISTE  $src -> $dst" >> "$LOG"; n_skip=$((n_skip+1)); continue
  fi

  echo "MOVER              $src -> $dst" >> "$LOG"
  n_ok=$((n_ok+1)); bytes=$((bytes+size))
  if [[ $APPLY -eq 1 ]]; then
    mkdir -p "$(dirname "$dst")"
    mv -n -- "$src" "$dst"
    printf 'mkdir -p %q && mv -n -- %q %q\n' "$(dirname "$src")" "$dst" "$src" >> "$UNDO_BODY"
  fi
done < <(awk -F'\t' "NR>1 && ( $SELECT )" "$MAP")

# El undo deshace en ORDEN INVERSO al de los movimientos: si un destino quedo dentro
# de un directorio creado por un movimiento anterior, revertirlo antes evita dejar
# directorios intermedios ocupados.
{ echo "#!/usr/bin/env bash"; echo "set -euo pipefail"; echo "cd \"$REPO\""; } > "$UNDO"
tac "$UNDO_BODY" >> "$UNDO"
chmod +x "$UNDO"
awk -v f="$PHASE" -v a="$APPLY" -v n="$n_ok" -v b="$bytes" -v s="$n_skip" -v m="$n_miss" -v q="$n_same" \
    -v l="$LOG" -v u="$UNDO" 'BEGIN{
  printf "fase=%s  aplicado=%s\n", f, (a==1 ? "si" : "NO (ensayo)");
  printf "  a mover:   %d ficheros (%.2f GB)\n", n, b/1e9;
  printf "  sin cambio: %d (ya estaban en su sitio)\n", q;
  printf "  omitidos:  %d (ya movidos o destino ocupado)\n", s;
  printf "  ausentes:  %d\n", m;
  printf "  detalle:   %s\n  deshacer:  %s\n", l, u; }'

if [[ $APPLY -eq 1 ]]; then
  find . -type d -empty -not -path './.git/*' -not -path './_TO_DELETE/*' -not -path './_ARCHIVE/*' -delete 2>/dev/null || true
  echo "  directorios vacios eliminados."
fi

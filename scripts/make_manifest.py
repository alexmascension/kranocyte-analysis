#!/usr/bin/env python3
"""
Genera results/manifest.tsv: el registro de todo lo que NO esta en git.

Para cada fichero de data/, results/ y figures/ anota ruta, tamano, fecha, sha256
y, cuando se puede deducir, el notebook que lo produjo. Es lo que mantiene la
trazabilidad entre un commit y las salidas que genero, ahora que los datos y las
figuras exploratorias viven fuera del control de versiones.

    python scripts/make_manifest.py              # escribe results/manifest.tsv
    python scripts/make_manifest.py --check      # falla si algo cambio sin anotarse

El hash de los ficheros de mas de HASH_LIMIT bytes se omite (columna '-'): hashear
cientos de GB en cada ejecucion no compensa. Para esos se registra tamano y fecha.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

HASH_LIMIT = 100 * 1024**2
TRACKED = ("data", "results", "figures")
SKIP = {"__pycache__", ".ipynb_checkpoints", "_TO_DELETE", "_ARCHIVE"}
SKIP_FILES = {"results/manifest.tsv"}   # the manifest does not record itself

#: filename prefix -> notebook that produces it
PRODUCERS = {
    "0AA": "notebooks/0AA_processing_inhouse_dataset.ipynb",
    "0NI": "notebooks/0NI_processing_nicoletti_2023.ipynb",
    "0OP": "notebooks/0OP_processing_oprescu_2020.ipynb",
    "0SO": "notebooks/0SO_processing_song_2023.ipynb",
    "0SU": "notebooks/0SU_processing_southerland_2023.ipynb",
    "1AA": "notebooks/1AA_cell_type_characterization.ipynb",
    "2AA": "notebooks/2AA_pseudotemporal_dynamics.ipynb",
    "2OP": "notebooks/2OP_pseudotemporal_dynamics.ipynb",
    "3AA": "notebooks/3AA_interactions.ipynb",
    "3NI": "notebooks/3NI_interactions.ipynb",
    "3SO": "notebooks/3SO_interactions.ipynb",
    "3SU": "notebooks/3SU_interactions.ipynb",
    "4_":  "notebooks/4_integration_and_LR_interactions.ipynb",
}


def repo_root() -> Path:
    env = os.environ.get("KRANOCITO_ROOT")
    if env:
        return Path(env).expanduser()
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "pyproject.toml").is_file():
            return parent
    return Path.cwd()


def producer(name: str) -> str:
    for prefix, nb in PRODUCERS.items():
        if name.startswith(prefix):
            return nb
    return "-"


def sha256(path: Path, size: int) -> str:
    if size > HASH_LIMIT:
        return "-"
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
    except OSError:
        return "ERROR"
    return digest.hexdigest()


def git_commit(root: Path) -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=root,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "-"


def collect(root: Path) -> list[tuple[str, int, str]]:
    found = []
    for top in TRACKED:
        base = root / top
        if not base.is_dir():
            continue
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = [d for d in dirnames if d not in SKIP]
            for name in filenames:
                full = Path(dirpath) / name
                if not full.is_file() or full.is_symlink():
                    continue
                try:
                    stat = full.stat()
                except OSError:
                    continue
                rel = str(full.relative_to(root))
                if rel in SKIP_FILES:
                    continue
                stamp = datetime.fromtimestamp(stat.st_mtime, timezone.utc).strftime("%Y-%m-%d")
                found.append((rel, stat.st_size, stamp))
    return sorted(found)


def build(root: Path) -> list[str]:
    entries = collect(root)
    with ThreadPoolExecutor(16) as pool:
        hashes = list(pool.map(lambda e: sha256(root / e[0], e[1]), entries))
    commit = git_commit(root)
    lines = ["path\tbytes\tdate\tsha256\tnotebook\tcommit"]
    for (rel, size, stamp), digest in zip(entries, hashes):
        lines.append(f"{rel}\t{size}\t{stamp}\t{digest}\t{producer(Path(rel).name)}\t{commit}")
    return lines


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="compare against the existing manifest and exit 1 if it differs")
    args = parser.parse_args()

    root = repo_root()
    target = root / "results" / "manifest.tsv"
    lines = build(root)

    if args.check:
        if not target.is_file():
            print(f"{target} does not exist", file=sys.stderr)
            return 1
        old = target.read_text().splitlines()
        # the commit changes on every run; everything else is compared
        strip = lambda rows: [r.rsplit("\t", 1)[0] for r in rows]
        if strip(old) == strip(lines):
            print(f"manifest up to date: {len(lines) - 1} files")
            return 0
        a, b = set(strip(old)), set(strip(lines))
        print(f"manifest differs: {len(b - a)} new or changed, {len(a - b)} missing",
              file=sys.stderr)
        for row in sorted(b - a)[:10]:
            print("  +", row.split("\t")[0], file=sys.stderr)
        for row in sorted(a - b)[:10]:
            print("  -", row.split("\t")[0], file=sys.stderr)
        return 1

    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines) + "\n")
    total = sum(int(r.split("\t")[1]) for r in lines[1:])
    print(f"{target}: {len(lines) - 1} files, {total / 1e9:.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

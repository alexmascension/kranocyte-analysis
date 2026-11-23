"""
Project paths and constants.

Replaces ``pyfuncs.common_vars``, which hard-coded ``BASE_DIR`` to one machine's
absolute path. The root is now resolved in this order:

1. the ``KRANOCYTE_ROOT`` environment variable, if set;
2. ``config/config.env`` in the repository, if present (key ``KRANOCYTE_ROOT``);
3. the working directory and its parents (for notebooks and installed packages);
4. the position of this package, walking up until the repository root is found.

The same code therefore works on the workstation and on the HPC without edits, which
is where the ``/scratch/alexmascension/...`` path in the old ``config.env`` came from.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = [
    "ROOT", "DATA_RAW", "DATA_EXTERNAL", "DATA_INTERIM", "DATA_PROCESSED",
    "CONFIG", "SAMPLES", "RESULTS", "TABLES", "FIGURES", "FIGURES_PAPER",
    "SEED", "CELLBENDER_FIXED_ARGS", "DATASETS",
    "layout", "processed_sources", "dataset_dir", "samples_file", "figure_dir",
]


def _abs(path: Path) -> Path:
    """Absolute, normalised path, independent of how ``resolve()`` behaves."""
    return Path(os.path.abspath(os.path.expanduser(str(path))))


def _from_env_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("KRANOCYTE_ROOT"):
            _, _, value = line.partition("=")
            return value.strip().strip('"').strip("'") or None
    return None


def _discover_root() -> Path:
    env = os.environ.get("KRANOCYTE_ROOT")
    if env:
        return _abs(Path(env).expanduser())

    cwd = Path.cwd().resolve()
    here = Path(__file__).resolve()
    candidates = dict.fromkeys((cwd, *cwd.parents, *here.parents))
    for parent in candidates:
        value = _from_env_file(parent / "config" / "config.env")
        if value:
            return _abs(Path(value).expanduser())
        # Repository-root markers. BOTH layouts are accepted: the current one
        # (src/<package>, src/notebooks_v2) and the migration target
        # (pyproject.toml + notebooks/), so the same code works before and after.
        if (parent / "pyproject.toml").is_file() and (parent / "notebooks").is_dir():
            return _abs(parent)
        if (parent / ".git").exists() and (parent / "src").is_dir():
            return _abs(parent)

    raise RuntimeError(
        "Cannot find the project root. Set KRANOCYTE_ROOT or create "
        "config/config.env from config/config.example.env, then start Python "
        "from the project directory or one of its subdirectories."
    )


ROOT: Path = _discover_root()

CONFIG         = ROOT / "config"
SAMPLES        = CONFIG / "samples"
DATA_RAW       = ROOT / "data" / "raw"
DATA_EXTERNAL  = ROOT / "data" / "external"
DATA_INTERIM   = ROOT / "data" / "interim"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS        = ROOT / "results"
TABLES         = RESULTS / "tables"
FIGURES        = ROOT / "figures"
FIGURES_PAPER  = FIGURES / "paper"

SEED = 0

# CellBender runs in a container; --cuda requires a GPU visible to that container.
CELLBENDER_FIXED_ARGS = [
    "--epochs", "150",
    "--fpr", "0.05",
    "--learning-rate", "1e-4",
    "--low-count-threshold", "5",
    "--z-dim", "64",
    "--z-layers", "512",
    "--empty-drop-training-fraction", "0.2",
    "--projected-ambient-count-threshold", "0.1",
    "--checkpoint-mins", "15",
    "--cuda",
]

#: Short dataset identifiers, as they end up under data/.
DATASETS = (
    "inhouse_ARAUZO_03",
    "oprescu_2020",
    "song_2023",
    "southerland_2023",
    "nicoletti_2023",
)


def layout() -> str:
    """``"migrated"`` if data/processed/<dataset>/ exists, ``"current"`` otherwise."""
    if (DATA_PROCESSED / "inhouse_ARAUZO_03").is_dir():
        return "migrated"
    return "current"


def processed_sources() -> list[Path]:
    """
    The ``*_processed*.h5ad`` objects of the five datasets, in either repository layout.

    Before the migration they live in ``data/ARAUZO_03/processed_adatas/`` and
    ``data/public_scRNAseq/datasets/<GSE>/processed_adatas/``; afterwards, in
    ``data/processed/<dataset>/``. This function exists so notebook 9 does not depend
    on which of the two is on disk.

    Returns
    -------
    list of pathlib.Path
        Sorted paths to the processed objects.

    Raises
    ------
    FileNotFoundError
        If no processed object is found under either layout.
    """
    if layout() == "migrated":
        found = sorted(DATA_PROCESSED.glob("*/*_processed*.h5ad"))
    else:
        found = sorted(
            list((ROOT / "data" / "ARAUZO_03" / "processed_adatas").glob("*_processed*.h5ad"))
            + list((ROOT / "data" / "public_scRNAseq" / "datasets").glob("*/processed_adatas/*_processed*.h5ad"))
        )
    if not found:
        raise FileNotFoundError(
            f"No processed objects found under {ROOT} (detected layout: {layout()})"
        )
    return found


def dataset_dir(dataset: str, tier: str = "interim") -> Path:
    """Directory of one dataset at the requested tier (raw, interim or processed)."""
    tiers = {"raw": DATA_RAW, "interim": DATA_INTERIM, "processed": DATA_PROCESSED}
    if tier not in tiers:
        raise ValueError(f"tier must be one of {sorted(tiers)}, not {tier!r}")
    if dataset not in DATASETS:
        raise ValueError(f"Unknown dataset: {dataset!r}. Known: {DATASETS}")
    return tiers[tier] / dataset


def samples_file(dataset: str) -> Path:
    """Sample table of the dataset, now under version control."""
    return SAMPLES / f"{dataset}.tsv"


def figure_dir(day: str | None = None, paper: bool = False) -> Path:
    """Figure destination: figures/paper/ or figures/exploratory/<date>/."""
    if paper:
        return FIGURES_PAPER
    from datetime import date
    return FIGURES / "exploratory" / (day or str(date.today()))

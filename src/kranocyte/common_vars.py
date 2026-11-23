"""
Backwards-compatible shim.

Everything here now lives in :mod:`kranocyte.config`, which resolves the project root
instead of hard-coding one machine's absolute path. Kept so that existing notebooks
can keep importing ``SEED`` and ``CELLBENDER_FIXED_ARGS`` from here.
"""

from .config import ROOT, SEED, CELLBENDER_FIXED_ARGS

#: Deprecated: use :data:`kranocyte.config.ROOT`.
BASE_DIR = str(ROOT)

__all__ = ["BASE_DIR", "ROOT", "SEED", "CELLBENDER_FIXED_ARGS"]

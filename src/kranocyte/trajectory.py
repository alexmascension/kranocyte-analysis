"""
Trajectory inference: RNA velocity (scVelo), directed graph (PAGA) and
PHATE embedding.

Replaces velocity_check.py, which can be deleted.

WHAT IT ANSWERS AND WHAT IT DOES NOT
------------------------------------
With scVelo + PAGA you can ask "does A go towards B?". You cannot ask "what
fraction of A ends up in B?": those are fate probabilities and require
Palantir or CellRank. Keep this in mind when writing conclusions.

PHATE is included here as an EMBEDDING, for visual inspection, not as an
inference method. Graphs are for inference; embeddings are for plotting.
Reducing to two dimensions and then measuring distances in them invalidates
whatever is measured.

THE UNDERLYING PROBLEM
----------------------
These methods ALWAYS return a trajectory. On a homogeneous cluster without
any real transition, PAGA draws a graph with edges and scVelo draws coherent
arrows. The figure is not evidence; the figure is the default output of the
method. This module is therefore built around three checks that do not
depend on external labels:

  1. STABILITY. Does the edge survive resampling of cells?
     -> paga_stability()
  2. INTERNAL CONSISTENCY. Does the direction given by PAGA agree with the
     order given by the pseudotime? They are two readings of the same velocity
     field and need not agree.
     -> computed automatically by paga_report()
  3. ROOT SENSITIVITY. If moving the root changes the order of populations,
     you are reading your hypothesis, not the data.
     -> run_velocity(root_key=...) with different roots

If an external ordinal variable is also available (days post-injury, stage,
dose, a proliferation score), pass it as `external_key` and it is used as an
anchor. It is the strongest evidence available, but it is OPTIONAL: the
module does not need it and is not designed around it.

ORIENTATION OF THE TRANSITION MATRIX
------------------------------------
The convention of whether transitions_confidence[i, j] means "from i to j" or
"from j to i" changes between scanpy/scVelo versions and is a classic source
of inverted conclusions. It is not assumed here: paga_report() checks the sign
of each edge against the pseudotime order and warns if the matrix appears to
be transposed.

WARNING
-------
The scVelo calls in this module were written without being able to install
scvelo in the development environment. The orchestration is tested with mocks;
the library calls are not. Verify every output before trusting it.
"""

from __future__ import annotations

import warnings
from typing import Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


# --------------------------------------------------------------------------
# normalisation of the Velocyto layers
# --------------------------------------------------------------------------

def normalize_velocyto_layers(
    adata: ad.AnnData,
    size_factor_key: str = "size_factors",
    layers: Sequence[str] = ("spliced", "unspliced"),
    keep_raw_suffix: str = "_raw",
    verbose: bool = True,
) -> ad.AnnData:
    """
    Scales spliced/unspliced by the SAME size factors that were used for X.

    The Velocyto layers are raw counts on purpose. The documented scVelo
    workflow is filter_and_normalize -> moments, and calling moments directly
    skips normalisation. Using the size factors already computed keeps S, U
    and X on the same scale, which is what is wanted.

    NOTE: this is NOT batch correction. It is depth normalisation. These are
    two different things that were mixed inside the same function in 0A.

    Stores the originals as '<layer>_raw' so that the step can be undone.
    """
    if size_factor_key not in adata.obs:
        raise ValueError(
            f"obs['{size_factor_key}'] not found; the size factors come from "
            f"scran_size_factors() on the full object"
        )
    sf = adata.obs[size_factor_key].to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(sf)) or (sf <= 0).any():
        raise ValueError("Some size factors are non-positive or non-finite")

    inv = sp.diags(1.0 / sf)
    for lay in layers:
        if lay not in adata.layers:
            warnings.warn(f"Layer '{lay}' not found")
            continue
        raw_key = f"{lay}{keep_raw_suffix}"
        if raw_key not in adata.layers:
            adata.layers[raw_key] = adata.layers[lay].copy()
        # sparse @ diags keeps CSR. (Note: in modern scipy `csr / sf` does not
        # densify either -- it returns COO -- so the direct division in 0A was
        # not the problem I thought it was; but @ diags avoids the format
        # change and the cost of rebuilding the COO.)
        adata.layers[lay] = sp.csr_matrix(
            inv @ sp.csr_matrix(adata.layers[raw_key]), dtype=np.float32
        )

    adata.uns.setdefault("velocity", {})["normalization"] = {
        "size_factor_key": size_factor_key, "layers": list(layers),
    }
    if verbose:
        print(f"[velo] {list(layers)} scaled by '{size_factor_key}'; "
              f"originals in '*{keep_raw_suffix}'")
    return adata


def _harmony_su(
    adata: ad.AnnData,
    batch_key: str,
    n_comps: int = 100,
    seed: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Second batch correction: on the S and U matrices themselves (variant
    '2bc' of 0A).

    Decomposes M = S + U, corrects M with harmony in PCA space, and splits the
    corrected M between S and U preserving the spliced fraction R = S/M of
    each cell and gene.

    What was wrong in the original version: Z_corr.T assumed harmonypy 1.x
    (PCs x cells). In 2.0 it comes as (cells x PCs) and that line produces
    garbage or crashes. The orientation is checked here.

    It is kept here so that it can be COMPARED, not because it is recommended:
    applying harmony to S and U separately corrects the same effect twice.
    """
    import harmonypy
    from sklearn.decomposition import PCA

    S = sp.csr_matrix(adata.layers["spliced"], dtype=np.float64)
    U = sp.csr_matrix(adata.layers["unspliced"], dtype=np.float64)
    M = (S + U).tocsr()

    Md = np.asarray(M.todense())
    Sd = np.asarray(S.todense())
    soporte = Md > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        R = np.where(soporte, Sd / Md, 0.0)         # spliced fraction
    del Sd

    M_log = np.log1p(Md)
    del Md

    pca = PCA(n_components=min(n_comps, min(M_log.shape) - 1),
              random_state=seed, svd_solver="arpack")
    Z = pca.fit_transform(M_log)

    ho = harmonypy.run_harmony(Z, adata.obs, [batch_key], max_iter_harmony=30)
    Zc = np.asarray(ho.Z_corr)
    if Zc.shape == (Z.shape[1], adata.n_obs):        # harmonypy 1.x
        Zc = Zc.T
    elif Zc.shape != (adata.n_obs, Z.shape[1]):
        raise ValueError(f"Z_corr has unexpected shape: {Zc.shape}")

    M_log_h = Zc @ pca.components_ + pca.mean_
    M_h = np.expm1(M_log_h)
    np.clip(M_h, 0, None, out=M_h)
    # BUG FIXED. The PCA reconstruction is DENSE: it returns a positive value
    # in entries that were zero in the original. Since R = 0 there, all of
    # that value went entirely to 'unspliced' (M_h * (1 - 0)), turning every
    # zero of the matrix into unprocessed signal. Result: U/S systematically
    # inflated across the whole object and velocity biased towards induction
    # everywhere, in addition to losing sparsity. The reconstruction is
    # restricted to the original support.
    M_h[~soporte] = 0.0
    del M_log, M_log_h, Zc, Z, soporte

    adata.layers["spliced"] = sp.csr_matrix((M_h * R).astype(np.float32))
    adata.layers["unspliced"] = sp.csr_matrix((M_h * (1.0 - R)).astype(np.float32))
    del M_h, R

    if verbose:
        print(f"[velo] 2bc: S and U reconstructed after harmony on M = S + U")
    return adata


# --------------------------------------------------------------------------
# a single velocity run, on the object passed in
# --------------------------------------------------------------------------

def run_velocity(
    adata: ad.AnnData,
    use_rep: str = "X_harmony",
    batch_key: Optional[str] = None,
    correct_su: bool = False,
    normalize_su: bool = True,
    ambient_max: Optional[float] = 0.2,
    ambient_key: str = "ambient_fraction",
    n_neighbors: Optional[int] = None,
    mode: str = "dynamical",
    root_key: Optional[str] = None,
    n_jobs: int = 4,
    seed: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    scVelo on THE OBJECT PASSED IN.

    In 0A this function wrote to the global `adata` in three of its four
    calls, so velocities were computed on all cells and the subset was left
    empty. There is no reference to global variables here.

    correct_su
        False -> '1bc': batch is corrected only in the representation that
                 feeds the neighbour graph (use_rep='X_harmony').
        True  -> '2bc': S and U are corrected as well. For comparison, not default.

    ambient_max
        Excludes from the velocity gene set the genes whose var[ambient_key]
        is above this value. The soup is 100% spliced, so in contaminated
        genes the observed U/S is biased downwards and they produce spurious
        negative velocity. None to disable.

    root_key
        Boolean obs column with the root cells. WITHOUT it, scVelo chooses
        the root itself and the sign of the pseudotime is arbitrary in each
        run: two executions on the same object can give pseudotimes that
        correlate at -1 with nothing having changed. To COMPARE runs
        (variants, subsets, resamplings) the root must be anchored and the
        same one used in all of them. For a single run it does not matter.
    """
    import scvelo as scv

    if use_rep not in adata.obsm:
        raise ValueError(f"obsm['{use_rep}'] not found")
    if n_neighbors is None:
        n_neighbors = max(10, int(0.5 * adata.n_obs ** 0.5))

    if normalize_su:
        normalize_velocyto_layers(adata, verbose=verbose)
    if correct_su:
        if batch_key is None:
            raise ValueError("correct_su=True requires batch_key")
        if adata.obs[batch_key].nunique() < 2:
            raise ValueError(
                f"correct_su=True with a single category of '{batch_key}': "
                f"there is no batch to correct"
            )
        _harmony_su(adata, batch_key, seed=seed, verbose=verbose)

    n_pcs = adata.obsm[use_rep].shape[1]
    scv.pp.moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
    scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
    scv.tl.velocity(adata, mode=mode)

    # ---- filter contaminated genes BEFORE building the graph
    if ambient_max is not None:
        if ambient_key not in adata.var:
            warnings.warn(
                f"var['{ambient_key}'] not found; genes are not filtered by ambient. "
                f"Compute it with qc.ambient_fraction_per_gene() on the concatenated "
                f"object and carry it over to the subset."
            )
        else:
            vg = adata.var["velocity_genes"].to_numpy().astype(bool)
            amb = adata.var[ambient_key].to_numpy(dtype=float)
            drop = vg & (amb > ambient_max)
            adata.var["velocity_genes"] = vg & ~(amb > ambient_max)
            if verbose:
                quitados = adata.var_names[drop].tolist()
                print(f"[velo] {int(drop.sum())} velocity genes excluded by "
                      f"{ambient_key} > {ambient_max}: {quitados[:10]}"
                      f"{' ...' if len(quitados) > 10 else ''}")
                print(f"[velo] {int(adata.var['velocity_genes'].sum())} "
                      f"velocity genes remain")

    scv.tl.velocity_graph(adata, n_jobs=n_jobs)
    if root_key is not None:
        if root_key not in adata.obs:
            raise ValueError(f"obs['{root_key}'] not found")
        scv.tl.velocity_pseudotime(adata, root_key=root_key)
    else:
        scv.tl.velocity_pseudotime(adata)
        if verbose:
            print("[velo] root chosen by scVelo: the SIGN of the pseudotime "
                  "is arbitrary. Do not compare it across runs without root_key.")
    scv.tl.velocity_confidence(adata)

    adata.uns.setdefault("velocity", {})["run"] = {
        "use_rep": use_rep, "correct_su": bool(correct_su),
        "normalize_su": bool(normalize_su), "ambient_max": ambient_max,
        "mode": mode, "n_neighbors": n_neighbors, "root_key": root_key,
        "n_velocity_genes": int(adata.var["velocity_genes"].sum()),
    }
    if verbose:
        print(f"[velo] mean confidence: "
              f"{float(adata.obs['velocity_confidence'].mean()):.3f}")
    return adata


# --------------------------------------------------------------------------
# PAGA: directed graph between populations
# --------------------------------------------------------------------------

def _parche_scvelo_igraph() -> bool:
    """
    scv.tl.paga crashes with modern scipy. It is not the fault of the user or the pipeline.

    scvelo/tools/paga.py does:

        csr_matrix((weights, zip(*edges)), shape=shape)

    A `zip` is a generator. Old scipy materialised it along the way; new
    scipy requires real sequences and finds zero index arrays, hence the
    error:

        ValueError: mismatching number of index arrays for shape;
                    got 0, expected 2

    The patch does NOT touch scVelo's edge logic: it only replaces the
    csr_matrix seen by that module with one that materialises the generator
    before passing it to scipy. Idempotent; returns True if patching was needed.
    """
    import importlib
    import sys

    # NOTE: `import scvelo.tools.paga as _pg` does NOT work. The __init__ of
    # scvelo.tools does `from .paga import paga`, so the package attribute
    # `paga` points to the FUNCTION and shadows the submodule; the
    # `import ... as` form resolves by attribute and returns the function. The
    # actual module is in sys.modules under its full name.
    importlib.import_module("scvelo.tools.paga")
    _pg = sys.modules["scvelo.tools.paga"]
    if not hasattr(_pg, "csr_matrix"):
        raise RuntimeError(
            "scvelo.tools.paga does not expose csr_matrix: your scvelo version is "
            "not the one this patch expects. Check whether scv.tl.paga already works without it."
        )

    if getattr(_pg.csr_matrix, "_patched", False):
        return False

    _csr = sp.csr_matrix

    def csr_matrix(arg1, shape=None, dtype=None, **kw):
        if isinstance(arg1, tuple) and len(arg1) == 2 and not np.isscalar(arg1[0]):
            data, ij = arg1
            if not isinstance(ij, (np.ndarray, list, tuple)):
                arg1 = (np.asarray(data), tuple(np.asarray(x) for x in ij))
        return _csr(arg1, shape=shape, dtype=dtype, **kw)

    csr_matrix._parcheado = True
    _pg.csr_matrix = csr_matrix
    return True


def _cats(adata: ad.AnnData, group_key: str) -> list:
    s = adata.obs[group_key]
    if not isinstance(s.dtype, pd.CategoricalDtype):
        adata.obs[group_key] = s.astype("category")
    elif len(s.cat.categories) != s.nunique():
        # Categories without cells silently break PAGA indices: the matrix
        # comes out with empty rows and the labels become misaligned.
        adata.obs[group_key] = s.cat.remove_unused_categories()
    return list(adata.obs[group_key].cat.categories)


def _dense(M) -> np.ndarray:
    if M is None:
        return None
    return np.asarray(M.todense()) if sp.issparse(M) else np.asarray(M)


def paga_report(
    adata: ad.AnnData,
    group_key: str,
    use_time_prior: str = "velocity_pseudotime",
    external_key: Optional[str] = None,
    min_cells: int = 20,
    verbose: bool = True,
) -> dict:
    """
    Velocity-directed PAGA, returned as NUMBERS.

    Uses scv.tl.paga, not sc.tl.paga(use_rna_velocity=True), which is the
    broken/deprecated path.

    Returns a dict:

      'connectivity'  Symmetric DataFrame (connectivity). Derived from the kNN
                      graph; it states whether two clusters are ADJACENT. It
                      is a statement about the embedding, not about biology
                      or direction. It is what PAGA draws by default and what
                      is misread as lineage.
      'transitions'  Asymmetric DataFrame (transitions). Derived from the
                      velocity. This is the one that carries the arrow.
      'net_flow'    T - T.T, antisymmetric (net flow). Positive at [a, b] =
                      flow goes from a to b. This is the robust reading: it
                      does not depend on the root or on the global pseudotime sign.
      'edges'       table (edges) sorted by |flujo_neto| with both metrics
                      alongside, plus the pseudotime step of each edge.
      'ordering'         mean pseudotime per population (order).
      'coherence'    fraction of edges in which the sign of the net flow
                      agrees with the pseudotime order (consistency).

    Does NOT set a threshold. The appearance of the graph depends heavily on
    it, and it is a plotting decision, not a result: choose it when plotting
    and report it.

    external_key
        Optional ORDINAL column (days, stage, dose, a score). If given, its
        Spearman correlation with the pseudotime is added. Not required.
    """
    import scvelo as scv

    if _parche_scvelo_igraph() and verbose:
        print("[paga] scvelo/scipy compatibility patch applied "
              "(see _parche_scvelo_igraph)")

    cats = _cats(adata, group_key)
    tam = adata.obs[group_key].value_counts()
    peq = tam[tam < min_cells]
    if len(peq) and verbose:
        print(f"[paga] WARNING: {list(peq.index)} have fewer than {min_cells} "
              f"cells ({dict(peq)}). Their edges are noise; do not interpret them.")

    if use_time_prior and use_time_prior not in adata.obs:
        raise ValueError(
            f"obs['{use_time_prior}'] not found; run run_velocity() first"
        )
    scv.tl.paga(adata, groups=group_key, use_time_prior=use_time_prior)

    P = adata.uns["paga"]
    conect = pd.DataFrame(_dense(P.get("connectivities")), index=cats, columns=cats)
    T = _dense(P.get("transitions_confidence"))
    if T is None:
        raise RuntimeError(
            "scv.tl.paga did not return 'transitions_confidence'. Without it only "
            "the undirected graph is available, which says nothing about direction."
        )
    trans = pd.DataFrame(T, index=cats, columns=cats)
    neto = trans - trans.T

    orden = adata.obs.groupby(group_key, observed=True)[use_time_prior].mean()

    filas = []
    for i, a in enumerate(cats):
        for b in cats[i + 1:]:
            filas.append({
                "a": a, "b": b,
                "net_flow": float(neto.loc[a, b]),
                "trans_a_b": float(trans.loc[a, b]),
                "trans_b_a": float(trans.loc[b, a]),
                "connectivity": float(conect.loc[a, b]),
                "delta_pseudotime": float(orden[b] - orden[a]),
                "n_a": int(tam[a]), "n_b": int(tam[b]),
            })
    aristas = pd.DataFrame(filas)
    aristas["direction"] = np.where(aristas["net_flow"] > 0,
                                    aristas["a"] + " -> " + aristas["b"],
                                    aristas["b"] + " -> " + aristas["a"])
    aristas = aristas.reindex(
        aristas["net_flow"].abs().sort_values(ascending=False).index
    ).reset_index(drop=True)

    # ---- internal consistency: does the flow go towards increasing pseudotime?
    vivas = aristas[aristas["net_flow"].abs() > 1e-8]
    if len(vivas):
        coherencia = float(
            (np.sign(vivas["net_flow"]) == np.sign(vivas["delta_pseudotime"])).mean()
        )
    else:
        coherencia = np.nan

    out = {"connectivity": conect, "transitions": trans, "net_flow": neto,
           "edges": aristas, "ordering": orden, "coherence": coherencia}

    if external_key is not None:
        from scipy.stats import spearmanr
        v = adata.obs[external_key]
        v = v.cat.codes.to_numpy(float) if isinstance(v.dtype, pd.CategoricalDtype) \
            else v.to_numpy(float)
        pt = adata.obs[use_time_prior].to_numpy(float)
        m = np.isfinite(v) & np.isfinite(pt)
        out["external_corr"] = float(spearmanr(pt[m], v[m]).statistic)
        out["external_ordering"] = adata.obs.groupby(external_key, observed=True)[
            use_time_prior].mean()

    adata.uns.setdefault("trajectory", {})["paga"] = {
        "group_key": group_key, "use_time_prior": use_time_prior,
        "coherence": coherencia, "n_groups": len(cats),
    }

    if verbose:
        print(f"\n[paga] mean pseudotime per population:")
        print(orden.sort_values().round(3).to_string())
        print(f"\n[paga] edges sorted by |net flow|:")
        print(aristas.head(15).round(4).to_string(index=False))
        print(f"\n[paga] internal consistency: {coherencia:.2f} of the edges "
              f"have flow towards increasing pseudotime.")
        if np.isfinite(coherencia) and coherencia < 0.5:
            print(
                "[paga] SERIOUS WARNING: below 0.5. Either the transition\n"
                "       matrix uses the transposed convention in your\n"
                "       scanpy/scVelo version (then ALL your arrows are\n"
                "       reversed and transposing 'transitions' suffices),\n"
                "       or the velocity field and the pseudotime contradict\n"
                "       each other and there is no trajectory to read. Check one\n"
                "       edge by hand against the velocity_embedding_grid figure\n"
                "       before continuing."
            )
        elif np.isfinite(coherencia) and coherencia < 0.8:
            print("[paga] mediocre consistency: there are edges where PAGA and the "
                  "pseudotime disagree. Inspect them one by one.")
        if "external_corr" in out:
            print(f"\n[paga] Spearman(pseudotime, {external_key}) = "
                  f"{out['external_corr']:+.3f}")
            print(out["external_ordering"].round(3).to_string())
        print("\n[paga] 'connectivity' is UNDIRECTED and comes from the kNN graph: "
              "it indicates adjacency, not lineage. Direction is in 'net_flow'.")
    return out


def paga_net_transition(adata: ad.AnnData, group_key: str,
                        pair: Sequence[str]) -> float:
    """
    Net flow between two populations: T[a,b] - T[b,a].

    Special case of paga_report(); kept because compare_batch_correction uses
    it as a reference comparable across runs (it does not depend on the root).
    """
    import scvelo as scv

    _parche_scvelo_igraph()
    a, b = pair
    if "paga" not in adata.uns or "transitions_confidence" not in adata.uns.get("paga", {}):
        scv.tl.paga(adata, groups=group_key)
    T = _dense(adata.uns["paga"]["transitions_confidence"])
    cats = _cats(adata, group_key)
    if a not in cats or b not in cats:
        raise ValueError(f"{pair} are not both in {group_key}: {cats}")
    i, j = cats.index(a), cats.index(b)
    return float(T[i, j] - T[j, i])


def paga_stability(
    adata: ad.AnnData,
    group_key: str,
    n_boot: int = 10,
    frac: float = 0.8,
    use_time_prior: str = "velocity_pseudotime",
    min_cells: int = 20,
    n_jobs: int = 4,
    seed: int = 0,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Which PAGA edges survive resampling of cells?

    This is the check that makes PAGA reportable. A graph is always produced;
    the question is which of its edges are a property of the data and which
    depend on the particular cells that were sampled.

    Returns, per pair of populations: the net flow in the full object, the
    median and interquartile range across repetitions, and `estabilidad` =
    fraction of repetitions in which the sign matches that of the full
    object. An edge with stability < 0.8 does not support a figure.

    COST AND SCOPE. The recover_dynamics fits of the full object are reused
    and only velocity_graph and paga are recomputed on each subsample. This
    is what makes the runtime affordable, and it means that what is measured
    is the stability OF THE GRAPH AND THE PARTITION, not of the kinetic fit.
    The kinetic fit is another source of variability not captured here.
    """
    import scvelo as scv

    _parche_scvelo_igraph()
    for req in ("velocity", "Ms", "Mu"):
        if req not in adata.layers:
            raise ValueError(
                f"layers['{req}'] is missing; run run_velocity() on this "
                f"object before measuring stability"
            )

    cats = _cats(adata, group_key)
    base = paga_report(adata, group_key, use_time_prior=use_time_prior,
                       min_cells=min_cells, verbose=False)
    neto_full = base["net_flow"]

    rng = np.random.default_rng(seed)
    n_sub = int(round(frac * adata.n_obs))
    muestras = []

    for k in range(n_boot):
        idx = rng.choice(adata.n_obs, size=n_sub, replace=False)
        sub = adata[np.sort(idx)].copy()
        sub.obs[group_key] = sub.obs[group_key].cat.remove_unused_categories()
        if len(sub.obs[group_key].cat.categories) < len(cats):
            faltan = set(cats) - set(sub.obs[group_key].cat.categories)
            warnings.warn(f"repetition {k}: no cells from {faltan}")
        try:
            scv.tl.velocity_graph(sub, n_jobs=n_jobs)
            scv.tl.paga(sub, groups=group_key, use_time_prior=use_time_prior)
            Tk = _dense(sub.uns["paga"]["transitions_confidence"])
            ck = list(sub.obs[group_key].cat.categories)
            Tk = pd.DataFrame(Tk, index=ck, columns=ck)
            muestras.append(Tk - Tk.T)
        except Exception as e:
            warnings.warn(f"repetition {k} failed: {e}")
        del sub
        if verbose:
            print(f"[stab] {k + 1}/{n_boot}", end="\r")

    if not muestras:
        raise RuntimeError("No repetition completed")

    filas = []
    for i, a in enumerate(cats):
        for b in cats[i + 1:]:
            vals = np.array([
                float(m.loc[a, b]) for m in muestras
                if a in m.index and b in m.index
            ])
            vals = vals[np.isfinite(vals)]
            ref = float(neto_full.loc[a, b])
            if vals.size == 0:
                est = np.nan
            elif abs(ref) < 1e-8:
                est = np.nan
            else:
                est = float((np.sign(vals) == np.sign(ref)).mean())
            filas.append({
                "a": a, "b": b,
                "full_flow": ref,
                "boot_median": float(np.median(vals)) if vals.size else np.nan,
                "boot_iqr": float(np.subtract(*np.percentile(vals, [75, 25])))
                            if vals.size > 3 else np.nan,
                "stability": est,
                "n_rep": int(vals.size),
            })

    out = pd.DataFrame(filas)
    out = out.reindex(
        out["full_flow"].abs().sort_values(ascending=False).index
    ).reset_index(drop=True)

    if verbose:
        print(" " * 30, end="\r")
        print(out.round(4).to_string(index=False))
        firmes = out[(out["stability"] >= 0.8) & out["stability"].notna()]
        print(f"\n[stab] {len(firmes)} of {len(out)} edges keep their sign "
              f"in >=80% of {n_boot} resamplings at {int(100*frac)}%.")
        if len(firmes):
            print("[stab] stable: " + ", ".join(
                f"{r.a}-{r.b}" for r in firmes.itertuples()))
        flojas = out[(out["stability"] < 0.8) & out["stability"].notna()]
        if len(flojas):
            print("[stab] NOT stable: " + ", ".join(
                f"{r.a}-{r.b} ({r.estabilidad:.2f})" for r in flojas.itertuples()))
        print("[stab] Note: this measures stability of the graph and the "
              "partition with the kinetic fits held fixed, not of the fit.")
    return out


# --------------------------------------------------------------------------
# comparison of batch-correction variants
# --------------------------------------------------------------------------

def compare_batch_correction(
    adata: ad.AnnData,
    sample_key: str = "gsm",
    group_key: str = "minor_population",
    pair: Optional[Sequence[str]] = None,
    use_rep_joint: str = "X_harmony",
    ambient_max: Optional[float] = 0.2,
    reprocess_fn=None,
    seed: int = 0,
    verbose: bool = True,
):
    """
    Runs the variants and compares them numerically on the SAME cells.

        per_<sample>    each sample alone, not integrated. There is no batch
                        to correct, so this is the reference.
        1bc             joint, harmony only in the graph representation.
        2bc             joint, harmony also on S and U.

    Criterion: if the separate samples agree with each other, that is the
    correct direction. The joint variant that reproduces it is the one to use.
    If the separate samples do NOT agree, the comparison decides nothing
    and the problem lies elsewhere.

    reprocess_fn
        Callable(adata_sub, integrate: bool) -> adata_sub, to recompute
        HVGs/PCA/harmony/neighbours in each variant. Typically a wrapped
        characterization.subset_and_reprocess, or processing.build_embeddings.
        If None, the object is assumed to already contain X_pca/X_harmony.

    Returns (pseudotiempos, resumen).
    """
    from scipy.stats import spearmanr

    muestras = list(adata.obs[sample_key].astype(str).unique())
    if len(muestras) < 2:
        raise ValueError(
            f"Only one category in '{sample_key}': there is nothing to compare"
        )
    if pair is None:
        raise ValueError(
            "pair is required. The criterion of this function is 'does the joint "
            "variant reproduce what the separate samples say?', and the "
            "per-sample runs share NO cells: their pseudotime correlation is "
            "NaN by construction. The only common reference is the net flow "
            "between two populations, which in addition does not depend on the "
            "root. Pass pair=('PopA','PopB')."
        )

    variantes = {}
    for m in muestras:
        variantes[f"per_{m}"] = dict(
            mask=(adata.obs[sample_key].astype(str) == m).to_numpy(),
            use_rep="X_pca", correct_su=False, integrate=False)
    variantes["1bc"] = dict(mask=np.ones(adata.n_obs, bool),
                            use_rep=use_rep_joint, correct_su=False, integrate=True)
    variantes["2bc"] = dict(mask=np.ones(adata.n_obs, bool),
                            use_rep=use_rep_joint, correct_su=True, integrate=True)

    pt, conf, trans = {}, {}, {}
    for nombre, cfg in variantes.items():
        if verbose:
            print(f"\n{'='*60}\n{nombre}  ({int(cfg['mask'].sum())} cells)\n{'='*60}")
        sub = adata[cfg["mask"]].copy()
        if reprocess_fn is not None:
            sub = reprocess_fn(sub, cfg["integrate"])
        run_velocity(sub, use_rep=cfg["use_rep"], batch_key=sample_key,
                     correct_su=cfg["correct_su"], ambient_max=ambient_max,
                     seed=seed, verbose=verbose)
        pt[nombre] = pd.Series(sub.obs["velocity_pseudotime"].to_numpy(),
                               index=sub.obs_names)
        conf[nombre] = float(sub.obs["velocity_confidence"].mean())
        if pair is not None:
            try:
                trans[nombre] = paga_net_transition(sub, group_key, pair)
            except Exception as e:      # a population may be missing from a subset
                trans[nombre] = np.nan
                warnings.warn(f"{nombre}: PAGA failed ({e})")
        del sub

    pseudotiempos = pd.DataFrame(pt)

    nombres = list(pseudotiempos.columns)
    corr = pd.DataFrame(np.nan, index=nombres, columns=nombres, dtype=float)
    for i in nombres:
        for j in nombres:
            if i == j:
                corr.loc[i, j] = 1.0
                continue
            # concat with new keys: pseudotiempos[[i, j]] with i == j would give
            # a DataFrame with duplicated columns and spearmanr would return a matrix
            comun = pd.concat([pseudotiempos[i], pseudotiempos[j]],
                              axis=1, keys=["a", "b"]).dropna()
            if len(comun) > 20:
                corr.loc[i, j] = float(spearmanr(comun["a"], comun["b"]).statistic)

    resumen = pd.DataFrame({
        "mean_confidence": pd.Series(conf),
        "net_flow": pd.Series(trans) if trans else np.nan,
    })
    per = [c for c in nombres if c.startswith("per_")]
    # NOTE: the sign of velocity_pseudotime is ARBITRARY in each run -- it
    # depends on which cell ends up chosen as root. Averaging signed
    # correlations makes +0.9 and -0.8 cancel out to ~0, which is exactly the
    # opposite of what happens (both are strong). |r| is averaged.
    resumen["corr_with_per_sample"] = [
        corr.loc[n, per].drop(index=n, errors="ignore").abs().mean()
        for n in nombres
    ]
    resumen["min_abs_corr"] = [
        corr.loc[n, per].drop(index=n, errors="ignore").abs().min()
        for n in nombres
    ]

    if verbose:
        print("\n" + "="*60)
        print("Spearman correlation of velocity_pseudotime")
        print("="*60)
        print(corr.round(3).to_string())
        print("\n" + resumen.round(3).to_string())
        if len(per) >= 2:
            # NOTE: the per-sample runs share NO cells, so the pseudotime
            # correlation between them is NaN by construction. To know whether
            # they agree, something defined at the POPULATION level, which is
            # shared, must be examined: the net flow of the pair.
            flujos = pd.Series({k: v for k, v in trans.items() if k in per}) \
                if trans else pd.Series(dtype=float)
            if flujos.notna().sum() >= 2:
                de_acuerdo = bool(np.all(np.sign(flujos.dropna()) == np.sign(flujos.dropna().iloc[0]))
                                  and (flujos.dropna().abs() > 1e-6).all())
                print(f"\n[velo] net flow per sample: "
                      f"{flujos.round(4).to_dict()}")
                print(f"[velo] the separate samples "
                      f"{'AGREE' if de_acuerdo else 'do NOT agree'} on the direction.")
                if not de_acuerdo:
                    print("[velo] Without a consistent reference this comparison "
                          "decides nothing: the problem is not the batch "
                          "correction. Check confidence, number of velocity genes "
                          "and whether the subset has sufficient structure.")
                    return pseudotiempos, resumen
            elif pair is None:
                print("\n[velo] Without 'pair' the samples cannot be compared with "
                      "each other: they share no cells, so the pseudotime "
                      "correlation between them is NaN by construction. Pass "
                      "pair=('PopA','PopB') to obtain the reference.")
            if True:
                cand = resumen.drop(index=per)["corr_with_per_sample"]
                mejor = cand.idxmax()
                print(f"[velo] variant most similar to the separate samples: "
                      f"{mejor} ({cand[mejor]:+.3f})")
                print("[velo] corr_con_per_muestra is |r|: it measures whether the "
                      "joint variant reproduces the ORDER of the separate "
                      "samples, not their direction. Direction is given by flujo_neto.")
        if pair is not None and trans:
            print(f"\n[velo] net flow {pair[0]} -> {pair[1]} "
                  f"(positive = from {pair[0]} to {pair[1]}):")
            print(pd.Series(trans).round(4).to_string())
    return pseudotiempos, resumen


# --------------------------------------------------------------------------
# PHATE: embedding for VISUAL INSPECTION
# --------------------------------------------------------------------------

def run_phate(
    adata: ad.AnnData,
    use_rep: str = "X_harmony",
    key_added: str = "X_phate",
    n_components: int = 2,
    knn: Optional[int] = None,
    decay: int = 40,
    t: str | int = "auto",
    n_jobs: int = 1,
    seed: int = 0,
    verbose: bool = True,
    **kwargs,
) -> ad.AnnData:
    """
    PHATE on the SAME representation that feeds the rest of the pipeline.

    PHATE is an EMBEDDING, not an inference method. It is for visual
    inspection: unlike UMAP it better preserves transition distances and
    elongated structures, which is exactly what matters when a continuum is
    suspected. But the coordinates it produces must NOT feed any pseudotime
    or any distance: that would mean reducing to two dimensions and then
    measuring in them, which is the `ndims_rep=2` error in Palantir.
    Inference is performed on the full graph or diffusion map.

    use_rep
        WHY THIS MATTERS MORE THAN ANYTHING: `phate.PHATE().fit_transform(adata.X)`
        -- which is what appears in every tutorial -- runs on the expression
        matrix and therefore IGNORES harmony. With several samples, the batch
        reappears and two parallel clouds that look like biology are
        obtained. Here it is run on `use_rep` (X_harmony by default), which
        is already integrated.

        Since that representation is ALREADY a reduced space, `n_pca=None` is
        passed so that PHATE does not run another PCA on top of it. If a raw
        expression matrix is passed with `use_rep="X"`, then it is advisable
        to let it run its PCA: pass n_pca=100 via kwargs.

    knn
        If not given, the same k as the already-built neighbour graph
        (uns['neighbors']) is used, so that PHATE and the graph view the
        tissue at the same scale. Two methods using neighbourhoods of
        different size is a classic source of "PHATE says one thing and PAGA
        another".

    t
        The diffusion scale, and the parameter that changes the figure the
        most. 'auto' selects it by the knee of the von Neumann entropy. It is
        recorded in uns['trajectory']['phate'] because it is a decision, not a
        constant.

    n_jobs
        1 by default, for reproducibility. As already observed with harmony,
        algorithms that reduce in parallel are not bit-for-bit identical
        across runs. Increase it if speed is preferred and small shifts of
        the embedding are acceptable.
    """
    try:
        import phate
    except ImportError as e:
        raise ImportError("pip install phate") from e

    if use_rep == "X":
        X = adata.X
        X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
        n_pca = kwargs.pop("n_pca", 100)
        if verbose:
            print("[phate] running on X (expression). NOTE: this is NOT "
                  "integrated; with more than one sample, the batch will "
                  "appear in the embedding.")
    else:
        if use_rep not in adata.obsm:
            raise ValueError(
                f"obsm['{use_rep}'] not found. Run build_embeddings() first, or "
                f"pass use_rep='X' if you really want to run on expression."
            )
        X = np.asarray(adata.obsm[use_rep])
        n_pca = kwargs.pop("n_pca", None)      # already a reduced space

    if knn is None:
        knn = adata.uns.get("neighbors", {}).get("params", {}).get("n_neighbors")
        if knn is None:
            knn = max(5, int(0.5 * adata.n_obs ** 0.5))
        if verbose:
            print(f"[phate] knn={knn} (inherited from the neighbour graph)")

    op = phate.PHATE(
        n_components=n_components, knn=knn, decay=decay, t=t,
        n_pca=n_pca, n_jobs=n_jobs, random_state=seed,
        verbose=1 if verbose else 0, **kwargs,
    )

    if n_jobs == 1:
        try:
            from threadpoolctl import threadpool_limits
        except ImportError:
            Y = op.fit_transform(X)
        else:
            with threadpool_limits(limits=1):
                Y = op.fit_transform(X)
    else:
        Y = op.fit_transform(X)

    adata.obsm[key_added] = np.ascontiguousarray(Y, dtype=np.float32)

    t_usado = getattr(op, "t", t)
    adata.uns.setdefault("trajectory", {})["phate"] = {
        "use_rep": use_rep, "knn": int(knn), "decay": decay,
        "t": int(t_usado) if isinstance(t_usado, (int, np.integer)) else str(t_usado),
        "n_pca": n_pca, "n_components": n_components, "seed": seed,
        "n_jobs": n_jobs,
    }

    if verbose:
        print(f"[phate] obsm['{key_added}'] ready; t={t_usado}, knn={knn}, "
              f"on '{use_rep}'")
        pt = "velocity_pseudotime"
        if pt in adata.obs and n_components >= 1:
            from scipy.stats import spearmanr
            v = adata.obs[pt].to_numpy(float)
            m = np.isfinite(v)
            r = float(spearmanr(Y[m, 0], v[m]).statistic)
            print(f"[phate] Spearman(PHATE1, {pt}) = {r:+.3f}. This is DESCRIPTIVE: "
                  f"it tells whether axis 1 of the plot matches the pseudotime, "
                  f"it validates neither of them.")
        print("[phate] Plot with sc.pl.embedding(adata, basis='phate', "
              "color=...). Do not use these coordinates to compute distances "
              "or pseudotimes.")
    return adata


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------

def plot_velocity(adata: ad.AnnData, basis: str = "umap",
                  color: Optional[str] = None, min_mass: float = 2.0,
                  autoscale: bool = False, **kwargs):
    """
    velocity_embedding_grid with the two defaults that matter.

    min_mass
        velocity_embedding_grid draws an arrow in every grid cell, including
        where there are hardly any cells. Those arrows are extrapolation, not
        data, and they are the longest and most conspicuous in the figure.
    autoscale
        True normalises the length per figure, so two figures are no longer
        comparable in magnitude. For comparing variants, False.
    """
    import scvelo as scv
    return scv.pl.velocity_embedding_grid(
        adata, basis=basis, color=color, min_mass=min_mass,
        autoscale=autoscale, **kwargs)
"""
Gene selection, PCA and integration, on the already normalized object.

Order of the calls (two rounds, as triku requires):

    preliminary HVGs (scanpy)   ->  PCA  ->  harmony  ->  neighbors  ->  triku
                                                                          |
                                          final HVGs <-------------------+
                                                |
                                PCA + harmony on the final set
                                                |
                                   X_pca                X_harmony
                              (temporal biology)        (annotation)

The first round exists only to obtain a reasonable neighbor graph: triku
scores genes by how their expression is distributed over that graph, so it
needs one built from something better than all genes. The preliminary HVGs
are a scaffold, not a result.

THE TWO ROUTES
--------------
Both are computed on THE SAME set of HVGs, so that they are comparable. The
only difference is whether harmony is applied:

    X_pca       uncorrected. For composition and dynamics across conditions.
    X_harmony   corrected. For clustering and cell-type annotation.

When sample and condition coincide (one replicate per timepoint), harmony
aligns the timepoints with each other and can merge states that exist in only
one of them. harmony_merge_report() looks for exactly that.
"""

from __future__ import annotations

import warnings
from typing import Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


# --------------------------------------------------------------------------
# 1. Preliminary HVGs (scaffold)
# --------------------------------------------------------------------------

def select_hvgs_preliminary(
    adata: ad.AnnData,
    n_top_genes: int = 3000,
    flavor: str = "seurat",
    batch_key: Optional[str] = None,
    key_added: str = "hvg_prelim",
    verbose: bool = True,
) -> ad.AnnData:
    """
    scanpy HVGs on X, already normalized and log-transformed.

    They serve only to build the graph that triku uses afterwards. Do not use
    them as the final set or report them as a result.

    batch_key
        With batch_key, scanpy selects per sample and combines, which prevents
        a gene variable in a single sample from dominating. Useful, but beware
        if sample == condition: it then penalizes genes that vary BETWEEN
        conditions, which is precisely the signal of interest. With few samples
        it is best left as None, relying on the second round.
    """
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor=flavor,
        batch_key=batch_key, inplace=True,
    )
    adata.var[key_added] = adata.var["highly_variable"].to_numpy()
    adata.uns.setdefault("processing", {})["hvg_prelim"] = {
        "n_top_genes": n_top_genes, "flavor": flavor, "batch_key": batch_key,
        "n_selected": int(adata.var[key_added].sum()),
    }
    if verbose:
        print(f"[hvg] preliminary: {int(adata.var[key_added].sum())} genes "
              f"(flavor={flavor}, batch_key={batch_key})")
    return adata


# --------------------------------------------------------------------------
# 2. PCA + harmony + neighbors
# --------------------------------------------------------------------------

def build_embeddings(
    adata: ad.AnnData,
    hvg_key: str,
    batch_key: str = "gsm",
    n_comps: int = 50,
    n_neighbors: Optional[int] = None,
    metric: str = "correlation",
    pca_key: str = "X_pca",
    harmony_key: str = "X_harmony",
    neighbors_on: str = "harmony",
    skip_harmony: bool = False,
    max_iter_harmony: int = 30,
    deterministic_harmony: bool = True,
    seed: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Computes PCA and its harmony-corrected version on the same set of genes,
    and the neighbor graph on whichever of the two is specified.

    neighbors_on
        'harmony' or 'pca'. Determines which representation feeds the graph
        (and therefore triku and the clustering). BOTH embeddings are stored in
        any case; this only selects which one is used now.

    Stores the graphs under their own neighbors_key ('neighbors_pca' /
    'neighbors_harmony') so that clustering can be run in both spaces without
    overwriting each other, which is what harmony_merge_report() needs.

    deterministic_harmony
        Restricts harmony to one thread so that X_harmony is reproducible. See
        the long comment in the body: with a fixed random_state but several
        threads the result still differs on every run. Set it to False only if
        speed matters more than exact reproducibility.
    """
    if hvg_key not in adata.var:
        raise ValueError(f"var['{hvg_key}'] not found")
    n_hvg = int(adata.var[hvg_key].sum())
    if n_hvg < 50:
        raise ValueError(f"Only {n_hvg} genes in '{hvg_key}'")

    adata.var["highly_variable"] = adata.var[hvg_key].to_numpy()  # scanpy reads this one
    sc.pp.pca(adata, n_comps=min(n_comps, n_hvg - 1, adata.n_obs - 1),
              use_highly_variable=True, random_state=seed)
    if pca_key != "X_pca":
        adata.obsm[pca_key] = adata.obsm["X_pca"].copy()

    if skip_harmony:
        # no batch to correct (e.g. a single sample in the subset):
        # X_harmony = X_pca so that the rest of the code does not need to branch
        adata.obsm[harmony_key] = adata.obsm["X_pca"].copy()
        neighbors_on = "pca"
    else:
      # harmonypy directly, not sce.pp.harmony_integrate: the scanpy wrapper
      # transposes Z_corr, and harmonypy >=2.0 already returns it as (cells x PCs).
      # With that combination the wrapper fails with a confusing shape error.
      import harmonypy

      # DETERMINISM. harmony initializes with sklearn.KMeans, and sklearn.KMeans
      # is NOT bit-for-bit reproducible when running multithreaded: OpenMP
      # floating-point reductions accumulate in a different order depending on
      # how the threads are scheduled. This is documented sklearn behavior, not
      # a harmonypy bug. With a fixed random_state the initialization is the
      # same, but those last-bit differences are amplified over the harmony
      # iterations and produce a different X_harmony on every run —
      # and with it different neighbors and UMAP.
      #
      # It is restricted to one thread ONLY during this call, with threadpoolctl
      # (shipped with sklearn). This is preferable to a global
      # OMP_NUM_THREADS=1: it does not slow down the rest of the pipeline and
      # does not require modifying the environment before the imports.
      def _run():
          try:
              return harmonypy.run_harmony(
                  adata.obsm["X_pca"], adata.obs, [batch_key],
                  max_iter_harmony=max_iter_harmony, random_state=seed,
              )
          except TypeError:
              warnings.warn(
                  "The installed harmonypy does not accept random_state; seeding "
                  "the global numpy RNG as a workaround. Check with two runs that "
                  f"obsm['{harmony_key}'] is identical."
              )
              np.random.seed(seed)
              return harmonypy.run_harmony(
                  adata.obsm["X_pca"], adata.obs, [batch_key],
                  max_iter_harmony=max_iter_harmony,
              )

      np.random.seed(seed)
      if deterministic_harmony:
          try:
              from threadpoolctl import threadpool_limits
          except ImportError:
              warnings.warn(
                  "Without threadpoolctl a single thread cannot be forced in harmony; "
                  "X_harmony may vary between runs. pip install "
                  "threadpoolctl, or export OMP_NUM_THREADS=1 before importing."
              )
              ho = _run()
          else:
              with threadpool_limits(limits=1):
                  ho = _run()
      else:
          ho = _run()
      Z = np.asarray(ho.Z_corr)
      n_pc = adata.obsm["X_pca"].shape[1]
      if Z.shape == (adata.n_obs, n_pc):
          pass
      elif Z.shape == (n_pc, adata.n_obs):
          Z = Z.T
      else:
          raise ValueError(
              f"harmonypy returned Z_corr with shape {Z.shape}; expected "
              f"({adata.n_obs}, {n_pc}) or its transpose"
          )
      adata.obsm[harmony_key] = np.ascontiguousarray(Z, dtype=np.float32)

    if n_neighbors is None:
        n_neighbors = max(10, int(0.5 * adata.n_obs ** 0.5))

    for name, rep in (("pca", pca_key), ("harmony", harmony_key)):
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=rep,
                        metric=metric, random_state=seed,
                        key_added=f"neighbors_{name}")

    # the "active" graph (default .obsp) is the one triku and leiden will use
    src = f"neighbors_{neighbors_on}"
    adata.uns["neighbors"] = dict(adata.uns[src])
    adata.obsp["distances"] = adata.obsp[adata.uns[src]["distances_key"]]
    adata.obsp["connectivities"] = adata.obsp[adata.uns[src]["connectivities_key"]]
    adata.uns["neighbors"]["distances_key"] = "distances"
    adata.uns["neighbors"]["connectivities_key"] = "connectivities"

    adata.uns.setdefault("processing", {})["embeddings"] = {
        "hvg_key": hvg_key, "n_hvg": n_hvg, "batch_key": batch_key,
        "n_comps": int(adata.obsm["X_pca"].shape[1]), "n_neighbors": n_neighbors,
        "metric": metric, "neighbors_on": neighbors_on,
        "seed": seed, "deterministic_harmony": bool(deterministic_harmony),
    }
    if verbose:
        print(f"[emb] PCA on {n_hvg} HVGs -> {adata.obsm['X_pca'].shape[1]} PCs")
        print(f"[emb] harmony by '{batch_key}' -> obsm['{harmony_key}']")
        print(f"[emb] neighbors (k={n_neighbors}, {metric}) in 'neighbors_pca' and "
              f"'neighbors_harmony'; active = {neighbors_on}")
    return adata


# --------------------------------------------------------------------------
# 3. Final HVGs with triku
# --------------------------------------------------------------------------

def select_hvgs_triku(
    adata: ad.AnnData,
    n_features: Optional[int] = None,
    key_added: str = "hvg_triku",
    use_raw: bool = False,
    random_state: int = 0,
    n_procs: Optional[int] = 1,
    verbose: bool = True,
    **kwargs,
) -> ad.AnnData:
    """
    triku on the active neighbor graph. Second round.

    triku scores each gene by how concentrated its expression is on the
    neighbor graph, not by the mean-variance relationship. That is why it
    needs the graph, and why there are two rounds.

    random_state
        NOT OPTIONAL, EVEN THOUGH IT LOOKS IT. triku compares the observed
        distribution against a null built by RANDOMIZING the expression
        (apply_background_correction). Without a seed, every run selects a
        different gene set, and from there PCA, harmony, neighbors, clusters
        and UMAP all change. This is the classic cause of "the UMAP comes out
        different every time the notebook is run".
    n_procs
        1 by default, also for reproducibility: with several processes the
        reduction order can vary. Increase it if speed matters more than
        exact determinism.
    """
    try:
        import triku as tk
    except ImportError as e:
        raise ImportError("pip install triku") from e

    if "neighbors" not in adata.uns:
        raise ValueError("No neighbor graph; run build_embeddings() first")

    tk_kwargs = dict(n_features=n_features, use_raw=use_raw, **kwargs)
    if n_procs is not None:
        tk_kwargs["n_procs"] = n_procs
    try:
        tk.tl.triku(adata, random_state=random_state, **tk_kwargs)
    except TypeError:
        warnings.warn(
            "The installed triku version does not accept random_state and/or "
            "n_procs. Seeding the global numpy RNG as a workaround, but this does "
            "NOT guarantee reproducibility: check with two runs that var['%s'] "
            "is identical." % key_added
        )
        np.random.seed(random_state)
        tk.tl.triku(adata, n_features=n_features, use_raw=use_raw, **kwargs)

    adata.var[key_added] = adata.var["highly_variable"].to_numpy()
    adata.uns.setdefault("processing", {})["hvg_triku"] = {
        "n_features": n_features, "n_selected": int(adata.var[key_added].sum()),
        "random_state": random_state, "n_procs": n_procs,
    }
    if verbose:
        print(f"[hvg] triku: {int(adata.var[key_added].sum())} genes")
    return adata


# --------------------------------------------------------------------------
# 4. Do soup genes leak into the HVGs?
# --------------------------------------------------------------------------

def soup_vs_hvg_report(
    adata: ad.AnnData,
    hvg_keys: Sequence[str] = ("hvg_prelim", "hvg_triku"),
    top_soup: int = 50,
    ambient_key: str = "ambient_expression",
    verbose: bool = True,
):
    """
    Checks whether ambient genes end up selected as variable.

    The question is legitimate because both scanpy and triku correct for
    expression level: if they are well calibrated, a highly expressed gene
    should NOT be selected merely for being highly expressed. So this is a
    test, not an excuse to exclude them beforehand.

    Returns (detalle, resumen):
      detalle : per soup gene, its ambient rank and whether it is in each set
      resumen : per HVG set, how many soup genes it contains, how many would be
                expected by chance, and the hypergeometric p-value

    How to read it: if 'observados' ~ 'esperados' and p > 0.05, the selectors
    are doing their job and nothing needs to be excluded. If there is clear
    enrichment, then yes: remove them from the set before PCA, not after.
    """
    from scipy.stats import hypergeom

    # the soup profile can come from var (one sample) or from varm (concatenated)
    if ambient_key in adata.varm:
        amb = pd.Series(np.asarray(adata.varm[ambient_key]).mean(axis=1),
                        index=adata.var_names, name=ambient_key)
        src = f"varm['{ambient_key}'] (mean of {adata.varm[ambient_key].shape[1]} samples)"
    elif ambient_key in adata.var:
        amb = adata.var[ambient_key].astype(float)
        src = f"var['{ambient_key}']"
    else:
        raise ValueError(f"'{ambient_key}' not found in var or in varm")

    amb = amb.sort_values(ascending=False)
    soup_genes = amb.head(top_soup).index

    detalle = pd.DataFrame({ambient_key: amb.head(top_soup)})
    detalle.insert(0, "rango_soup", np.arange(1, len(soup_genes) + 1))
    for k in hvg_keys:
        if k in adata.var:
            detalle[k] = adata.var.loc[soup_genes, k].to_numpy().astype(bool)

    rows = []
    M = adata.n_vars
    for k in hvg_keys:
        if k not in adata.var:
            warnings.warn(f"var['{k}'] not found")
            continue
        sel = adata.var[k].to_numpy().astype(bool)
        N = int(sel.sum())
        obs = int(adata.var.loc[soup_genes, k].to_numpy().astype(bool).sum())
        esp = top_soup * N / M
        # P(X >= obs) under random selection
        p = float(hypergeom.sf(obs - 1, M, top_soup, N)) if obs > 0 else 1.0
        rows.append({
            "hvg_set": k, "n_hvg": N, "n_genes": M,
            "soup_observados": obs, "soup_esperados": round(esp, 2),
            "enriquecimiento": round(obs / esp, 2) if esp > 0 else np.nan,
            "p_hipergeom": p,
        })
    resumen = pd.DataFrame(rows).set_index("hvg_set")

    if verbose:
        print(f"soup profile from {src}; top {top_soup} genes\n")
        print(resumen.to_string())
        print()
        for k in resumen.index:
            o, e = resumen.loc[k, "soup_observados"], resumen.loc[k, "soup_esperados"]
            p = resumen.loc[k, "p_hipergeom"]
            if p < 0.05 and o > e:
                print(f"[soup] {k}: ENRICHED ({o} vs {e:.1f} expected, "
                      f"p={p:.2g}). Consider excluding them before PCA.")
            else:
                print(f"[soup] {k}: no enrichment ({o} vs {e:.1f}, "
                      f"p={p:.2g}). The selector corrects well for expression; "
                      f"nothing needs to be excluded.")
        sel_cols = [k for k in hvg_keys if k in detalle.columns]
        if sel_cols:
            hit = detalle[detalle[sel_cols].any(axis=1)]
            if len(hit):
                print(f"\nsoup genes selected in at least one set:")
                print(hit.head(25).to_string())
    return detalle, resumen


def exclude_genes_from_hvg(
    adata: ad.AnnData,
    genes: Sequence[str],
    hvg_key: str,
    key_added: Optional[str] = None,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Removes specific genes from an HVG set. Use it ONLY if soup_vs_hvg_report()
    has shown enrichment: excluding blindly discards biology just in
    case.
    """
    key_added = key_added or hvg_key
    present = [g for g in genes if g in adata.var_names]
    missing = set(genes) - set(present)
    if missing:
        warnings.warn(f"{len(missing)} genes are not in var_names: {sorted(missing)[:5]}")

    sel = adata.var[hvg_key].to_numpy().astype(bool).copy()
    idx = adata.var_names.get_indexer(present)
    removed = int(sel[idx].sum())
    sel[idx] = False
    adata.var[key_added] = sel

    if verbose:
        print(f"[hvg] {removed} of {len(present)} genes excluded from '{hvg_key}' "
              f"-> '{key_added}' ({int(sel.sum())} genes)")
    return adata


# --------------------------------------------------------------------------
# 5. What does harmony merge?
# --------------------------------------------------------------------------

def harmony_merge_report(
    adata: ad.AnnData,
    condition_key: str,
    resolution: float = 1.0,
    dominance: float = 0.8,
    seed: int = 0,
    verbose: bool = True,
):
    """
    Clusters in both spaces and looks for clusters that exist without
    correction but dissolve upon integration.

    The failure of harmony when sample == condition does not appear as "the
    signal disappears from the UMAP". It appears in the GRANULARITY: a state
    present in only one timepoint forms its own cluster in X_pca, and harmony
    pushes it towards the equivalent state of the others until leiden merges
    them. The state is not attenuated: it ceases to exist as a category, and
    then no downstream analysis can find it.

    dominance
        A cluster is considered "condition-specific" if this fraction of its
        cells comes from a single condition.

    REQUIREMENT: condition_key must have AT LEAST TWO values. With only one, all
    clusters are trivially 100% pure and the table is noise. That case in fact
    means that batch and condition are not confounded and X_harmony can be
    used without reservation: the two-route dilemma does not arise there.

    Returns (sospechosos, contingencia). 'sospechosos' lists the
    condition-specific clusters in X_pca together with how many X_harmony
    clusters they are split across: if split across several mixed ones,
    harmony has dissolved it.
    """
    for nk in ("neighbors_pca", "neighbors_harmony"):
        if nk not in adata.uns:
            raise ValueError(f"uns['{nk}'] is missing; run build_embeddings() first")
    if condition_key not in adata.obs:
        raise ValueError(f"obs['{condition_key}'] not found")

    n_cond = adata.obs[condition_key].nunique()
    if n_cond < 2:
        raise ValueError(
            f"obs['{condition_key}'] has a single value "
            f"({adata.obs[condition_key].iloc[0]!r}). With a single condition ALL "
            f"clusters are 100% pure by definition and the report says "
            f"nothing.\n"
            f"  This is good news: if batch and condition are not confounded, "
            f"harmony performs legitimate batch correction and there is no "
            f"two-route dilemma. Use X_harmony as is.\n"
            f"  To look for populations exclusive to ONE SAMPLE, "
            f"pass condition_key='<your sample column>'."
        )
    if verbose and n_cond == adata.obs.groupby(condition_key, observed=True).ngroups:
        pass

    for name in ("pca", "harmony"):
        sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_{name}",
                     neighbors_key=f"neighbors_{name}", flavor="igraph",
                     n_iterations=2, directed=False, random_state=seed)

    comp = pd.crosstab(adata.obs["leiden_pca"], adata.obs[condition_key],
                       normalize="index")
    especificos = comp.index[comp.max(axis=1) >= dominance].tolist()

    cont = pd.crosstab(adata.obs["leiden_pca"], adata.obs["leiden_harmony"])

    rows = []
    for c in especificos:
        fila = cont.loc[c]
        fila = fila[fila > 0].sort_values(ascending=False)
        n = int(fila.sum())
        principal = fila.index[0]
        # is the harmony cluster it ends up in mixed?
        mezcla = comp.loc[c].idxmax()
        harm_comp = pd.crosstab(adata.obs["leiden_harmony"],
                                adata.obs[condition_key], normalize="index")
        rows.append({
            "leiden_pca": c, "n_cells_es": n,
            "condition": mezcla,
            "pureza_pca": round(float(comp.loc[c].max()), 3),
            "n_clusters_harmony": int(len(fila)),
            "harmony_principal": principal,
            "pureza_harmony_destino": round(float(harm_comp.loc[principal].max()), 3),
        })
    sospechosos = pd.DataFrame(rows)

    if verbose:
        n_pca = adata.obs["leiden_pca"].nunique()
        n_h = adata.obs["leiden_harmony"].nunique()
        print(f"[harmony] {n_pca} clusters in X_pca, {n_h} in X_harmony "
              f"(resolution={resolution})")
        if sospechosos.empty:
            print(f"[harmony] no X_pca cluster exceeds {dominance:.0%} from "
                  f"a single condition: no candidate condition-specific states.")
        else:
            print(f"[harmony] {len(sospechosos)} condition-specific clusters "
                  f"in X_pca:\n")
            print(sospechosos.to_string(index=False))
            disueltos = sospechosos[sospechosos["pureza_harmony_destino"] < dominance]
            if len(disueltos):
                print(f"\n[harmony] WARNING: {len(disueltos)} of them end up in a "
                      f"MIXED harmony cluster. Inspect their DEGs before accepting "
                      f"the integrated version: it may be a batch artifact that "
                      f"harmony was right to absorb, or the transient state "
                      f"being sought.")
    return sospechosos, cont
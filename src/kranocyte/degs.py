"""
Differential expression between CONDITIONS, within each population.

    res   = deg_by_group(adata_FAP, "cell_type", "condition",
                         reference="non_DEN", test="DEN_2d")
    sets  = deg_sets(res, n_top=100, direction="up")
    plot_upset(sets)
    lfc_concordance(res)
    core_vs_specific(sets)

THE WARNING THAT MUST NOT BE OMITTED
------------------------------------
Comparing conditions while treating CELLS as replicates is the case where
pseudoreplication does the most damage. With one mouse per condition, the test
answers "do these cells differ?" and not "do these conditions differ?", and no
correction is possible: the true n is 1 vs 1. Squair et al. (Nat Commun 2021)
showed that in this scenario p-values are inflated by orders of
magnitude.

With several mice per condition, the correct approach is PSEUDOBULK:
summing counts per (sample x population) and testing between samples. Without
them, the analysis can proceed — but it then RANKS candidate genes and does not
demonstrate differences, and must be reported as such.

AND THE CONFOUNDER SPECIFIC TO THIS DATASET
-------------------------------------------
Depth per cell differs up to 5x between timepoints (D2 ~70,000 median counts
versus ~15,000 in NI). Depth and condition are the same variable in this
design, and the detection rate of a gene depends strongly on depth.
Consequence: lowly expressed genes will appear "up" in the more deeply
sequenced condition without anything having changed.

This is why deg_by_group() requires pts=True and returns the detection
difference, and why depth_report() exists: inspect it BEFORE interpreting the list.
"""

from __future__ import annotations

import warnings
from typing import Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


# --------------------------------------------------------------------------
# 0. the confounder, first of all
# --------------------------------------------------------------------------

def depth_report(adata: ad.AnnData, condition_key: str,
                 group_key: Optional[str] = None,
                 counts_key: str = "total_counts",
                 verbose: bool = True) -> pd.DataFrame:
    """
    Depth and number of cells per condition (and population).

    If the two conditions to be compared differ greatly in counts per cell,
    any DEG list will be contaminated by it: lowly expressed genes are
    detected more often where depth is higher. Check the ratio before
    looking at the genes.
    """
    cols = [condition_key] + ([group_key] if group_key else [])
    g = adata.obs.groupby(cols, observed=True)
    out = pd.DataFrame({
        "n_cells_es": g.size(),
        "median_counts": g[counts_key].median() if counts_key in adata.obs else np.nan,
        "median_genes": (g["n_genes_by_counts"].median()
                          if "n_genes_by_counts" in adata.obs else np.nan),
    })
    if verbose:
        print(out.round(1).to_string())
        if group_key is None and len(out) == 2:
            r = out["median_counts"].max() / max(out["median_counts"].min(), 1)
            print(f"\n[depth] depth ratio between conditions: {r:.1f}x")
            if r > 1.5:
                print("[depth] WARNING: above 1.5x, the detection rate "
                      "already biases the DEG list. Compare with "
                      "deg_by_group(..., downsample=True).")
    return out


# --------------------------------------------------------------------------
# 1. the test, per population
# --------------------------------------------------------------------------

def deg_by_group(
    adata: ad.AnnData,
    group_key: str,
    condition_key: str,
    reference: str,
    test: str,
    groups: Optional[Sequence[str]] = None,
    method: str = "wilcoxon",
    tie_correct: bool = True,
    min_cells: int = 20,
    downsample: bool = False,
    counts_layer: str = "counts_cellbender",
    seed: int = 0,
    verbose: bool = True,
) -> dict:
    """
    Runs the test-vs-reference test WITHIN each population. {population: df}.

    The sign is fixed by construction: positive `logfoldchanges` = UP
    in `test` relative to `reference`. In the earlier loop this was decided by
    the volcano's `cluster=` argument, and choosing the wrong group inverts the
    whole list without anything failing.

    downsample
        Equalizes depth between the two conditions by downsampling counts
        before the test, on `counts_layer`. It has a cost —information is
        discarded— but it is the only way to know whether the list survives the
        depth confounder. Run both versions and compare.

    min_cells
        Populations with fewer cells than this in either condition are
        skipped: with 8 cells no test is meaningful and the list would be noise.
    """
    obs = adata.obs
    if condition_key not in obs or group_key not in obs:
        raise ValueError(f"Missing obs['{condition_key}'] or obs['{group_key}']")

    if groups is None:
        groups = [g for g in obs[group_key].astype(str).unique()]
    salida, resumen = {}, []

    for g in groups:
        m = (obs[group_key].astype(str) == g) & \
            (obs[condition_key].astype(str).isin([reference, test]))
        sub = adata[m.to_numpy()].copy()
        cond = sub.obs[condition_key].astype(str)
        n_ref, n_test = int((cond == reference).sum()), int((cond == test).sum())
        if min(n_ref, n_test) < min_cells:
            if verbose:
                print(f"[deg] {g}: {n_ref} vs {n_test} cells -> skipping "
                      f"(min_cells={min_cells})")
            continue

        if downsample:
            if counts_layer in sub.layers:
                sub.X = sub.layers[counts_layer].copy()
            objetivo = int(np.percentile(
                np.asarray(sub.X.sum(axis=1)).ravel(), 5))
            sc.pp.downsample_counts(sub, counts_per_cell=objetivo,
                                    random_state=seed)
            sc.pp.normalize_total(sub, target_sum=1e4)
            sc.pp.log1p(sub)
            if verbose:
                print(f"[deg] {g}: downsampled to {objetivo} counts/cell")

        sub.obs["_cond"] = pd.Categorical(cond, categories=[reference, test])
        kw = dict(groupby="_cond", groups=[test], reference=reference,
                  method=method, pts=True, key_added="deg")
        if method == "wilcoxon":
            kw["tie_correct"] = tie_correct
        sc.tl.rank_genes_groups(sub, **kw)

        df = sc.get.rank_genes_groups_df(sub, group=test, key="deg")
        # DETECTION difference: the signal of the depth confounder
        pts = sub.uns["deg"].get("pts")
        if pts is not None:
            df = df.merge(
                pts.rename(columns={reference: "pct_ref", test: "pct_test"})
                   .reset_index().rename(columns={"index": "names"}),
                on="names", how="left")
            df["delta_pct"] = df["pct_test"] - df["pct_ref"]
        df["population"] = g
        salida[g] = df
        resumen.append({"population": g, f"n_{reference}": n_ref,
                        f"n_{test}": n_test, "n_genes": len(df)})
        del sub

    if verbose and resumen:
        print(pd.DataFrame(resumen).set_index("population").to_string())
        print(f"\n[deg] POSITIVE logfoldchanges = up in '{test}' "
              f"relative to '{reference}'.")
    return salida


def depth_bias_report(res: Mapping[str, pd.DataFrame], adata: ad.AnnData,
                      layer: str = "counts_cellbender",
                      verbose: bool = True) -> pd.DataFrame:
    """
    Does the DEG list carry the signature of depth bias?

    The signature is unmistakable and can be simulated: with TWO groups of
    IDENTICAL expression but sequenced at 15,000 and 70,000 counts, after
    normalization and log1p ~7% of genes come out "significant", 98% of them
    UP in the deeper group, and the effect concentrates in poorly detected genes
    (mean logFC 0.29 and +42 detection points in low-expression genes,
    versus 0.08 and +2 points in high-expression genes).

    This function reproduces that table on the actual results. If the same
    pattern appears —directional asymmetry and bias concentrated in lowly
    expressed genes—, a large part of the list is depth. If logFC is flat
    with respect to expression level, it is not.
    """
    M = adata.layers[layer] if layer in adata.layers else adata.X
    media = pd.Series(np.asarray(M.mean(axis=0)).ravel(), index=adata.var_names)

    filas = []
    for pop, df in res.items():
        d = df.copy()
        d["media"] = media.reindex(d["names"]).to_numpy()
        d["nivel"] = pd.cut(d["media"], [-.01, .05, .3, 1, 5, 1e9],
                            labels=["~0", "muy bajo", "bajo", "medio", "alto"])
        for niv, sub in d.groupby("nivel", observed=True):
            sig = sub[sub["pvals_adj"] < 0.05]
            filas.append({
                "population": pop, "nivel": str(niv), "n": len(sub),
                "logFC_medio": round(float(sub["logfoldchanges"].mean()), 3),
                "delta_detection": (round(float(sub["delta_pct"].mean()), 3)
                                    if "delta_pct" in sub else np.nan),
                "pct_sig": round(100 * len(sig) / max(len(sub), 1), 1),
                "frac_sig_up": (round(float((sig["logfoldchanges"] > 0).mean()), 2)
                                     if len(sig) else np.nan),
            })
    out = pd.DataFrame(filas)
    if verbose:
        print(out.to_string(index=False))
        alza = out.groupby("population")["frac_sig_up"].mean()
        desbal = alza[(alza > 0.85) | (alza < 0.15)]
        if len(desbal):
            print(f"\n[bias] {list(desbal.index)} have >85% of DEGs in "
                  f"ONE direction only. Real biology is expected to be more "
                  f"balanced; that asymmetry is the signature of technical bias.")
        print("[bias] Also compare 'logFC_medio' across levels: if it grows as "
              "expression decreases, it is depth, not biology.")
    return out


# --------------------------------------------------------------------------
# 2. sets and UpSet
# --------------------------------------------------------------------------

def deg_sets(
    res: Mapping[str, pd.DataFrame],
    n_top: Optional[int] = 100,
    direction: str = "up",
    max_padj: float = 0.05,
    min_lfc: float = 0.5,
    min_delta_pct: Optional[float] = None,
    min_pct_both: Optional[float] = 0.1,
    verbose: bool = True,
) -> dict:
    """
    Gene sets per population, filtered and ranked. {population: set}.

    direction
        'up' or 'down'. ALWAYS SEPARATE. A gene that goes up in FAP.1 and down
        in FAP.3 is not a "shared" gene: it is the opposite. Putting both in the
        same UpSet turns a divergence into an overlap.

    min_pct_both
        THIS is the filter against the depth confounder: it requires the gene
        to be detected in at least this fraction of cells IN BOTH
        conditions. Depth bias concentrates almost entirely in poorly
        detected genes —there the depth difference translates into a
        difference in zeros, and the rank-sum reads mostly zeros—, whereas
        genes well expressed in both groups are barely affected.

        (Correction of an earlier statement: min_delta_pct does NOT protect
        against this. With 5x depth, the artifact PRODUCES large delta_pct;
        requiring it would select artifacts instead of discarding them.)

    min_delta_pct
        Minimum difference in the fraction of cells expressing the gene. Useful
        to retain presence/absence changes rather than magnitude-only
        changes, but it is NOT a depth control. Use it together with
        min_pct_both, never instead of it.
    """
    sets = {}
    for pop, df in res.items():
        d = df[df["pvals_adj"] <= max_padj]
        d = d[d["logfoldchanges"] >= min_lfc] if direction == "up" \
            else d[d["logfoldchanges"] <= -min_lfc]
        if min_pct_both is not None and {"pct_ref", "pct_test"} <= set(d.columns):
            antes = len(d)
            d = d[(d["pct_ref"] >= min_pct_both) & (d["pct_test"] >= min_pct_both)]
            if verbose and antes:
                print(f"[sets] {pop}: {antes} -> {len(d)} after requiring "
                      f"detection >= {min_pct_both} in BOTH conditions")
        if min_delta_pct is not None and "delta_pct" in d:
            antes = len(d)
            d = d[d["delta_pct"].abs() >= min_delta_pct]
            if verbose and antes:
                print(f"[sets] {pop}: {antes} -> {len(d)} after delta_pct "
                      f">= {min_delta_pct}")
        d = d.reindex(d["logfoldchanges"].abs().sort_values(ascending=False).index)
        genes = d["names"].tolist()
        sets[pop] = set(genes[:n_top] if n_top else genes)

    if verbose:
        print(f"\n[sets] direction '{direction}': "
              + ", ".join(f"{k}={len(v)}" for k, v in sets.items()))
        tam = [len(v) for v in sets.values()]
        if tam and max(tam) > 3 * max(min(tam), 1):
            print("[sets] The sets differ widely in size. The UpSet "
                  "favors the large ones: look at the fractions, not only at "
                  "the bars.")
    return sets


def plot_upset(sets: Mapping[str, set], min_subset_size: int = 1,
               sort_by: str = "cardinality", show_counts: bool = True,
               figsize=(11, 6), **kwargs):
    """
    UpSet of the sets. Returns (fig, upsetplot dict).

    How NOT to read it: the "only FAP.3" bar does not mean "gene exclusive to
    FAP.3". It means "gene that entered the top-N of FAP.3 and not that of the
    others", which with truncated lists is very different — it may be at
    rank N+1 in the other four. This is why lfc_concordance(), which does not
    depend on any cutoff, should also be examined.
    """
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_contents

    datos = from_contents({k: sorted(v) for k, v in sets.items()})

    def _dibuja(counts):
        fig = plt.figure(figsize=figsize)
        up = UpSet(datos, min_subset_size=min_subset_size, sort_by=sort_by,
                   show_counts=counts, **kwargs)
        ejes = up.plot(fig=fig)
        fig.canvas.draw()      # the label failure is raised WHEN DRAWING,
        return fig, ejes       # not when building: it must be forced here


    try:
        return _dibuja(show_counts)
    except TypeError as e:
        # upsetplot 0.9 + recent matplotlib versions crash when drawing the
        # count labels. The problem is the library, not the figure.
        if not show_counts:
            raise
        plt.close("all")
        warnings.warn(f"upsetplot failed with show_counts=True ({e}); retrying without "
                      f"the labels. The sizes are in core_vs_specific().")
        return _dibuja(False)


# --------------------------------------------------------------------------
# 3. what does not depend on a cutoff
# --------------------------------------------------------------------------

def lfc_concordance(
    res: Mapping[str, pd.DataFrame],
    restrict: Optional[Sequence[str]] = None,
    max_padj: Optional[float] = None,
    min_abs_lfc: float = 0.0,
    method: str = "pearson",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Correlation of logFC between populations, over ALL shared genes.

    It is the version of the UpSet that does not depend on n_top. The UpSet
    answers "do the lists match?", which is a question about the threshold;
    this answers "is it the same response?", which is the biological question.
    If two populations correlate at 0.8 in logFC but share few genes in the
    top-100, the response is shared and the small overlap was an
    artifact of the cutoff.

    restrict
        PASS THE UNION OF THE DEG SETS. Over the 25,000 genes, the
        correlation is dominated by the thousands that change in no population
        and whose logFC is noise: it comes out low even if the response is
        identical. The useful question is "among the genes that change in SOME
        population, do they change in the same way?".

            union = set().union(*sets.values())
            lfc_concordance(res, restrict=union)

    method
        'pearson' by default: here the MAGNITUDE of the change is informative
        and Spearman discards it. With Spearman, a set of genes that all go up
        by similar amounts is dominated by rank noise and the correlation
        comes out artificially low. Check both if in doubt.
    """
    lfc = {}
    for pop, df in res.items():
        d = df
        if max_padj is not None:
            d = d[d["pvals_adj"] <= max_padj]
        if min_abs_lfc:
            d = d[d["logfoldchanges"].abs() >= min_abs_lfc]
        lfc[pop] = d.set_index("names")["logfoldchanges"]

    M = pd.DataFrame(lfc)
    if restrict is not None:
        M = M.loc[M.index.intersection(list(restrict))]
    elif verbose:
        print("[lfc] no 'restrict': ALL genes are used, and those that "
              "change in no population dilute the correlation. Pass the "
              "union of the DEG sets.")
    corr = M.corr(method=method, min_periods=50)
    if verbose:
        print(f"[lfc] {len(M)} genes; {method} correlation of logFC:")
        print(corr.round(3).to_string())
        vals = corr.where(~np.eye(len(corr), dtype=bool)).stack()
        if len(vals):
            print(f"\n[lfc] mean across pairs: {vals.mean():.3f} "
                  f"(min {vals.min():.3f} between "
                  f"{vals.idxmin()[0]} and {vals.idxmin()[1]})")
            print("[lfc] High = common response to injury. Low = each "
                  "population responds in its own way. This is what the UpSet "
                  "suggests but does not measure.")
    return corr


def query_sets(
    sets: Mapping[str, set],
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
    exact: bool = True,
    n_populations: Optional[int] = None,
    res: Optional[Mapping[str, pd.DataFrame]] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Extracts the genes of a specific cell of the UpSet.

        query_sets(up, ["FAP.4"])                      # ONLY in FAP.4
        query_sets(up, list(up))                       # shared by all
        query_sets(up, ["FAP.1", "FAP.3"])             # only in those two
        query_sets(up, ["FAP.1", "FAP.3"], exact=False)  # in both, possibly more
        query_sets(up, n_populations=3)                # in exactly 3, any

    THE DISTINCTION THAT MATTERS: `exact`
        True  -> the gene is in the populations of `include` AND IN NO OTHER.
                 This is what each UpSet bar shows.
        False -> the gene is at least in those of `include`, and may be in
                 others.

        "Genes shared by FAP.1 and FAP.3" means different things depending on
        which is chosen, and the difference is usually large. exact=True asks
        for "shared by those two and by no other", which is a claim of
        specificity; exact=False, only "present in both". The function
        prints both numbers so they are not confused.

    res
        If the dict from deg_by_group() is passed, adds the logFC of EVERY
        population for the genes found, including those whose set the gene
        did not enter. That is the informative part: a "FAP.4-specific" gene
        whose logFC in FAP.1 is 0.9 —just below the threshold— is not
        specific, it is an artifact of the cutoff. If it is at 0.05 in the
        others, then it is.
    """
    todas = list(sets)
    if n_populations is not None:
        elegidos = [g for g in set().union(*sets.values())
                    if sum(g in v for v in sets.values()) == n_populations]
        etiqueta = f"in exactly {n_populations} populations"
    else:
        if not include:
            raise ValueError("Provide 'include' or 'n_populations'")
        falta = [p for p in include if p not in sets]
        if falta:
            raise ValueError(f"{falta} are not in sets. Available: {todas}")
        dentro = set.intersection(*[sets[p] for p in include])
        fuera = set(exclude) if exclude else (set(todas) - set(include)
                                             if exact else set())
        prohibidos = set().union(*[sets[p] for p in fuera]) if fuera else set()
        elegidos = sorted(dentro - prohibidos)
        etiqueta = (f"only in {list(include)}" if exact and not exclude
                    else f"in {list(include)}"
                         + (f", excluding {list(exclude)}" if exclude else ""))
        if verbose and not exclude:
            otras = set(todas) - set(include)
            en_otras = set().union(*[sets[p] for p in otras]) if otras else set()
            print(f"[query] in {list(include)} (possibly also in others): "
                  f"{len(dentro)} genes")
            print(f"[query] ONLY in {list(include)}: {len(dentro - en_otras)} genes")

    elegidos = sorted(elegidos)
    out = pd.DataFrame({"gene": elegidos})
    out["populations"] = [";".join(sorted(p for p in todas if g in sets[p]))
                          for g in elegidos]
    out["n_populations"] = out["populations"].str.count(";") + 1

    if res is not None:
        for pop in todas:
            if pop not in res:
                continue
            d = res[pop].set_index("names")
            out[f"lfc_{pop}"] = d["logfoldchanges"].reindex(elegidos).round(2).values
        cols_lfc = [c for c in out.columns if c.startswith("lfc_")]
        if cols_lfc:
            out = out.sort_values(cols_lfc[0], ascending=False)

    if verbose:
        print(f"\n[query] {len(out)} genes {etiqueta}")
        if len(out):
            print(", ".join(out["gene"].tolist()[:60])
                  + (" ..." if len(out) > 60 else ""))
    return out.reset_index(drop=True)


def core_vs_specific(sets: Mapping[str, set], verbose: bool = True) -> pd.DataFrame:
    """
    Splits the genes into shared core / intermediate / specific.

    The core (in all populations) is the generic response to injury; the
    specific genes are what distinguishes each population. These are usually
    two different stories and should be passed to GOEA separately — a joint
    analysis is dominated by the core, which is the same for all.
    """
    todos = sorted(set().union(*sets.values())) if sets else []
    n = len(sets)
    filas = []
    for g in todos:
        donde = [k for k, v in sets.items() if g in v]
        filas.append({"gene": g, "n_populations": len(donde),
                      "populations": ";".join(sorted(donde))})
    out = pd.DataFrame(filas).sort_values(["n_populations", "gene"],
                                          ascending=[False, True])
    if verbose and len(out):
        print(out["n_populations"].value_counts().sort_index(ascending=False)
              .rename("n_genes").to_string())
        nucleo = out[out["n_populations"] == n]["gene"].tolist()
        print(f"\n[core] core ({n}/{n} populations), {len(nucleo)} genes:")
        print("  " + ", ".join(nucleo[:40]) + (" ..." if len(nucleo) > 40 else ""))
        unicos = out[out["n_populations"] == 1]
        print(f"\n[core] {len(unicos)} genes in a single population")
        print("[core] Pass core and specific genes to run_goea() SEPARATELY: "
              "together, the core dominates the result and masks what is "
              "specific to each population.")
    return out
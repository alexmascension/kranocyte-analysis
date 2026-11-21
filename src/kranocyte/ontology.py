"""
GO term over-representation analysis (ORA) per population, with gseapy.

    background = build_background(adata, min_cells=10)
    degs       = top_degs(adata, "minor_population", n_top=150)
    res        = run_goea(degs, background, organism="Mouse")
    plot_goea(res)

THE BACKGROUND IS HALF OF THE RESULT
------------------------------------
An over-representation test asks: "of the N genes that COULD have come out,
k from this term did; is that more than expected?". That "could have come
out" is the background, and choosing it badly changes the result more than
any other decision.

The default in almost every tool is the whole genome (~22,000 genes in
mouse). In single-cell data this is wrong, and always in the same direction:
the assay can only detect genes expressed in THE tissue under study, so
comparing against the full genome systematically inflates every
tissue-specific term. In a muscle dataset "muscle contraction" will come out
enriched in every population, including the non-muscle ones.

THE RULE: the background is the set of genes on which the differential
expression TEST WAS RUN. If rank_genes_groups ran on the 25,000 genes that
remained after filter_genes, the background is those 25,000. No more (genes
that could never have come out) and no fewer (genes that did compete).

HVGs AS BACKGROUND? This was asked about, and the answer is: not by default.
HVGs are not "the genes that competed": rank_genes_groups is not restricted to
them. Using them as background answers a different question ("among the
variable genes, is this term over-represented in this population's markers?"),
which is legitimate and more conservative, but is not the one the analyst
believes they are asking. If that question is wanted, build_background(
use_hvg=True) builds it and records it in the name. What must NOT happen is
mixing: an HVG background with DEGs computed on all genes means part of the
foreground is not even in the universe, and the test is ill-defined.

DOUBLE DIPPING
--------------
These DEGs come from comparing clusters that were defined on the same data.
The differential expression p-values are no longer valid as evidence, and the
GO p-values inherit the whole problem. They are useful to RANK hypotheses
about what each population does, not to assert that it does it.

gseapy: MIND WHICH FUNCTION
---------------------------
gp.enrichr() uses the Enrichr web API and IGNORES the custom background. The one
that accepts a custom background is gp.enrich(), which runs Fisher's test
locally. This module always uses gp.enrich().
"""

from __future__ import annotations

import re
import warnings
from typing import Iterable, Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


# The version is deliberately NOT pinned here: Enrichr publishes GO every year,
# and hardcoding a year means silently keeping an outdated annotation.
# latest_go_libraries() asks Enrichr which ones exist NOW and takes the newest.
PREFIJOS_GO = ("GO_Biological_Process", "GO_Molecular_Function",
               "GO_Cellular_Component")

_RE_ANYO = re.compile(r"_(\d{4})$")


def list_libraries(organism: str = "Mouse", pattern: Optional[str] = None,
                   verbose: bool = True) -> list:
    """Libraries available in Enrichr, optionally filtered by substring."""
    import gseapy as gp

    libs = gp.get_library_name(organism=organism)
    if pattern:
        libs = [l for l in libs if pattern.lower() in l.lower()]
    if verbose:
        for l in sorted(libs):
            print(" ", l)
    return sorted(libs)


def latest_go_libraries(organism: str = "Mouse", prefixes: Sequence[str] = PREFIJOS_GO,
                        verbose: bool = True) -> list:
    """
    The MOST RECENT version of each GO branch that Enrichr serves today.

    Querying instead of hardcoding a year has two advantages: the annotation
    does not go stale, and nothing breaks when Enrichr retires a version.

    But record which one was returned. GO changes between versions (terms are
    merged, genes are re-annotated) and two analyses with different versions
    are not comparable. run_goea() stores the resolved names in the result
    for exactly that reason.
    """
    import gseapy as gp

    disponibles = gp.get_library_name(organism=organism)
    elegidas = []
    for pref in prefixes:
        cand = [l for l in disponibles if l.startswith(pref)]
        if not cand:
            warnings.warn(f"Enrichr ({organism}) does not offer any '{pref}*'")
            continue
        def _anyo(l):
            m = _RE_ANYO.search(l)
            return int(m.group(1)) if m else -1
        elegidas.append(max(cand, key=_anyo))
    if verbose:
        print(f"[libs] {organism}: {elegidas}")
    return elegidas


# --------------------------------------------------------------------------
# 1. the background
# --------------------------------------------------------------------------

def build_background(
    adata: ad.AnnData,
    min_cells: int = 10,
    group_key: Optional[str] = None,
    layer: Optional[str] = "counts_cellbender",
    use_hvg: bool = False,
    hvg_key: str = "hvg_triku",
    exclude: Optional[Iterable[str]] = None,
    verbose: bool = True,
) -> list:
    """
    Gene universe for the test. Returns a list of symbols.

    min_cells
        A gene detected in 3 cells out of 9,000 had no real chance of coming
        out as a marker, so including it in the universe only dilutes.

        WHY THIS CRITERION AND NOT TOTAL COUNTS. Trimming by sum(counts)
        looks equivalent but is not: a gene with 100 total counts can be 10
        counts in each of 10 cells of a population of 12, i.e. the perfect
        marker of a rare population. Total counts confound "lowly expressed
        everywhere" with "highly expressed in a few cells", and only the
        second is of interest. Counting CELLS separates those two cases;
        summing counts does not.

    group_key
        If given (e.g. 'minor_population'), a gene enters the background when
        it is detected in >= min_cells cells OF AT LEAST ONE population, rather
        than in the whole object. This is the correct criterion for an ORA of
        cluster markers: the question is "could this gene have come out as a
        marker of ANY population?", and with small populations the global
        threshold removes precisely the markers of rare populations. With
        kranocytes in the object, this is not a subtlety.

    layer
        Counting is done on COUNTS, not on the normalised X, because "expressed
        in N cells" is only meaningful on counts. If the layer is absent, it
        falls back to X with a warning.

    use_hvg
        Changes the question (see the module header). If enabled, make sure
        the DEGs were also computed restricted to HVGs.

    exclude
        Genes to remove from the universe AND, for consistency, from the
        foreground. Globins and anything known to be soup belong here: if a gene
        is in the dataset because of contamination, it should not count as
        "could have come out".
    """
    if use_hvg:
        if hvg_key not in adata.var:
            raise ValueError(f"var['{hvg_key}'] not found")
        genes = adata.var_names[adata.var[hvg_key].to_numpy().astype(bool)]
        origen = f"HVGs ({hvg_key})"
    else:
        if layer and layer in adata.layers:
            M = adata.layers[layer]
            fuente = layer
        else:
            M = adata.X
            fuente = "X"
            if layer:
                warnings.warn(
                    f"Layer '{layer}' not found; counting cells with a non-zero "
                    f"value in X. On log-normalised data the criterion remains "
                    f"valid (zero is zero), but it should be documented."
                )
        if group_key is None:
            n_cel = np.asarray((M > 0).sum(axis=0)).ravel()
            keep = n_cel >= min_cells
            origen = f"detected in >={min_cells} cells on '{fuente}'"
        else:
            if group_key not in adata.obs:
                raise ValueError(f"obs['{group_key}'] not found")
            keep = np.zeros(adata.n_vars, dtype=bool)
            tam = {}
            for grp in adata.obs[group_key].astype(str).unique():
                m = (adata.obs[group_key].astype(str) == grp).to_numpy()
                n_cel = np.asarray((M[m] > 0).sum(axis=0)).ravel()
                keep |= n_cel >= min_cells
                tam[grp] = int(m.sum())
            origen = (f">={min_cells} cells in at least one of the "
                      f"{len(tam)} categories of '{group_key}' on '{fuente}'")
            if verbose:
                peq = {k: v for k, v in tam.items() if v < 2 * min_cells}
                if peq:
                    print(f"[background] WARNING: {peq} have fewer than {2*min_cells} "
                          f"cells, so min_cells={min_cells} is almost their whole "
                          f"size. Lower min_cells or those populations will not "
                          f"contribute any gene to the universe.")
        genes = adata.var_names[keep]

    genes = pd.Index(genes)
    n0 = len(genes)
    if exclude:
        exclude = set(exclude)
        genes = genes[~genes.isin(exclude)]

    if verbose:
        print(f"[background] {len(genes)} genes out of {adata.n_vars} ({origen})")
        if exclude:
            print(f"[background] {n0 - len(genes)} explicitly excluded")
        print("[background] Check that this number matches the genes on "
              "which rank_genes_groups was run. If it does not, the test is "
              "misspecified in one direction or the other.")
    return list(genes)


# --------------------------------------------------------------------------
# 2. the foreground
# --------------------------------------------------------------------------

def top_degs(
    adata: ad.AnnData,
    groupby: str,
    n_top: int = 150,
    key: str = "rank_genes_groups",
    method: str = "wilcoxon",
    tie_correct: bool = True,
    max_padj: float = 0.05,
    min_logfc: float = 0.25,
    min_pct: float = 0.10,
    rank_by: str = "pi",
    exclude: Optional[Iterable[str]] = None,
    recompute: bool = False,
    verbose: bool = True,
) -> dict:
    """
    Top N markers per population, filtered, as a dict {population: [genes]}.

    Taking the first N by score is not enough: among them there may be genes
    with a minute effect or expressed in 2% of cells, and those drag in GO
    terms at random. The filters are applied BEFORE trimming to N.

    tie_correct
        Wilcoxon only, and only when recompute=True. scanpy defaults to
        False; with a matrix full of zeros, ties dominate the ranking and
        not correcting them inflates significance. True by default here.

    min_pct
        Requires rank_genes_groups to have been run with pts=True. If it is
        not available, the filter is skipped with a warning.

    rank_by
        Criterion used for sorting before cutting to n_top:

        'pi'     -log10(padj) * logFC. The pi-value of Xiao et al. (2014):
                 combines evidence and effect size, instead of letting
                 only one of the two decide. This is the default.
        'score'  the z-score returned by scanpy (its native order).
        'logfc'  effect size only. Useful as a control: if the top N changes
                 a lot between 'logfc' and 'pi', the list is being decided by
                 the p-value, i.e. by the cluster size.

        CAVEAT ON 'pi' WITH MANY CELLS. With 9,000 cells the p-values
        underflow to 0 and -log10(p) saturates; here they are clipped to the
        smallest observed positive value so that the product stays finite, but
        this means that dozens of genes share the same -log10(p) and among
        them the LFC alone determines the order. This is not a code defect; the
        p-value simply no longer discriminates at this sample size.

    The dict also reports how many genes each population retained: if one
    ends up with 12 genes and another with 150, their GO results are not
    comparable in power.
    """
    if recompute or key not in adata.uns:
        if verbose:
            print(f"[deg] computing rank_genes_groups ({method}) on '{groupby}'")
        kw = dict(groupby=groupby, method=method, key_added=key, pts=True)
        if method == "wilcoxon":
            # With data this sparse, ties (all the zeros) are
            # everywhere. scanpy defaults to tie_correct=False and that
            # inflates significance; it is enabled here.
            kw["tie_correct"] = tie_correct
        sc.tl.rank_genes_groups(adata, **kw)

    grupos = list(adata.uns[key]["names"].dtype.names)
    exclude = set(exclude or ())
    salida, resumen = {}, []

    for g in grupos:
        df = sc.get.rank_genes_groups_df(adata, group=g, key=key)
        n_ini = len(df)
        df = df[df["names"].notna()]

        if "pvals_adj" in df:
            df = df[df["pvals_adj"] <= max_padj]
        if "logfoldchanges" in df:
            df = df[df["logfoldchanges"] >= min_logfc]

        col_pct = next((c for c in ("pct_nz_group", "pts") if c in df.columns), None)
        if col_pct is not None:
            df = df[df[col_pct] >= min_pct]
        elif verbose and min_pct:
            print(f"[deg] no cell-fraction column; skipping min_pct "
                  f"(re-run with pts=True if it is needed)")

        df = df[~df["names"].isin(exclude)]

        # --- sorting
        if rank_by == "pi":
            col_p = "pvals_adj" if "pvals_adj" in df else "pvals"
            pv = df[col_p].to_numpy(dtype=float)
            positivos = pv[pv > 0]
            suelo = float(positivos.min()) if positivos.size else 1e-300
            pv = np.clip(pv, suelo, 1.0)
            df = df.assign(_pi=-np.log10(pv) * df["logfoldchanges"].to_numpy(float))
            df = df.sort_values("_pi", ascending=False)
        elif rank_by == "logfc":
            df = df.sort_values("logfoldchanges", ascending=False)
        elif rank_by == "score":
            if "scores" in df:
                df = df.sort_values("scores", ascending=False)
        else:
            raise ValueError("rank_by must be 'pi', 'score' or 'logfc'")

        genes = df["names"].tolist()[:n_top]
        salida[g] = genes

        fila = {"population": g, "testados": n_ini,
                "tras_filtros": len(df), "usados": len(genes)}
        if verbose and rank_by != "score" and "scores" in df:
            alt = df.sort_values("scores", ascending=False)["names"].tolist()[:n_top]
            inter = len(set(genes) & set(alt))
            fila["jaccard_vs_score"] = round(
                inter / max(len(set(genes) | set(alt)), 1), 2)
        resumen.append(fila)

    tabla = pd.DataFrame(resumen).set_index("population")
    if verbose:
        print(tabla.to_string())
        pocos = tabla.index[tabla["usados"] < 20].tolist()
        if pocos:
            print(f"[deg] WARNING: {pocos} are left with fewer than 20 genes. Their "
                  f"GOEA will have very little power, and the absence of terms does "
                  f"NOT mean the absence of biology.")
        if "jaccard_vs_score" in tabla:
            j = tabla["jaccard_vs_score"].mean()
            print(f"[deg] Mean Jaccard between the top-{n_top} by '{rank_by}' and "
                  f"by the scanpy score: {j:.2f}. If it is close to 1, the "
                  f"choice of criterion changes nothing and is not worth discussing.")
        desigual = tabla["usados"].max() / max(tabla["usados"].min(), 1)
        if desigual > 3:
            print(f"[deg] The lists differ up to {desigual:.0f}x in size. "
                  f"The number of significant terms scales with the size "
                  f"of the list: do not compare populations by how many terms "
                  f"they yield.")
    return salida


# --------------------------------------------------------------------------
# 3. coverage: do the dataset symbols and the library speak the same language?
# --------------------------------------------------------------------------

def _genes_de_libreria(gs: Mapping[str, Sequence[str]]) -> set:
    out = set()
    for v in gs.values():
        out.update(v)
    return out


def coverage_report(background: Sequence[str], gene_sets: Mapping,
                    nombre: str = "", verbose: bool = True) -> dict:
    """
    Fraction of the background present in the library, as-is and upper-cased.

    This is the check that almost nobody runs and that invalidates entire results:
    if the symbols are mouse RefSeq (Col1a1) and the library is human
    (COL1A1), the overlap is almost nil, gseapy does not fail, and it returns an
    empty table or, worse, one based on the handful of genes that match by chance.
    """
    lib = _genes_de_libreria(gene_sets)
    bg = set(background)
    tal_cual = len(bg & lib) / max(len(bg), 1)
    mayus = len({g.upper() for g in bg} & {g.upper() for g in lib}) / max(len(bg), 1)
    out = {"libreria": nombre, "n_fondo": len(bg), "n_libreria": len(lib),
           "cobertura": round(tal_cual, 3), "cobertura_mayusculas": round(mayus, 3)}
    if verbose:
        print(f"[cov] {nombre}: {tal_cual:.1%} of the background is in the library "
              f"({mayus:.1%} ignoring case)")
        if tal_cual < 0.30:
            print("[cov] WARNING: very low coverage.")
            if mayus > tal_cual + 0.2:
                print("[cov] Ignoring case raises it a lot -> the dataset symbols and "
                      "the library symbols come from different organisms/nomenclatures. "
                      "Check organism= before reading anything.")
    return out


def _resolver_nomenclatura(
    background: Sequence[str],
    gene_dict: Mapping[str, Sequence[str]],
    gene_sets: Mapping[str, Mapping],
    match_case: str = "auto",
    verbose: bool = True,
):
    """
    Puts the dataset symbols and the library symbols in the same alphabet.

    WHY THIS IS NEEDED. Enrichr serves its mouse libraries with symbols in
    UPPER CASE (COL1A1), whereas mouse MGI/RefSeq nomenclature is
    capitalised (Col1a1). Even when organism='Mouse' is requested, the GENES
    come in upper case: the organism determines WHICH annotation is served,
    not how the symbols are written. Without normalisation the overlap is
    0.1% and gseapy does not fail: it returns almost empty tables with the
    few genes that match by chance (those already fully upper case, such as
    mitochondrial genes).

    WHAT THIS DOES NOT FIX. Upper-casing matches the symbols that differ
    only in capitalisation, which are the majority, but not those with
    different names across species: Trp53 does not become TP53. Those genes
    are left out and cannot be recovered without an orthology table.

    match_case: 'auto' | 'upper' | 'exact'.
    """
    bgs = list(dict.fromkeys(background))
    lib = set()
    for gs in gene_sets.values():
        lib |= _genes_de_libreria(gs)

    exacta = len(set(bgs) & lib) / max(len(bgs), 1)
    mayus = len({g.upper() for g in bgs} & {g.upper() for g in lib}) / max(len(bgs), 1)

    if match_case == "auto":
        modo = "upper" if mayus > exacta + 0.10 else "exact"
        if verbose:
            print(f"[case] exact overlap {exacta:.1%}, upper-case "
                  f"{mayus:.1%} -> mode '{modo}'")
            if modo == "upper":
                print("[case] Enrichr serves mouse libraries with symbols "
                      "in UPPER CASE. Both sides are normalised. Note: this does NOT "
                      "match genes with different names across species "
                      "(Trp53 != TP53); those are lost.")
    else:
        modo = match_case

    if modo != "upper":
        return bgs, dict(gene_dict), dict(gene_sets), None, modo

    mapa = {}
    colisiones = []
    for g in bgs:
        u = g.upper()
        if u in mapa and mapa[u] != g:
            colisiones.append((mapa[u], g))
        else:
            mapa[u] = g
    if colisiones and verbose:
        print(f"[case] {len(colisiones)} collisions when upper-casing "
              f"(e.g. {colisiones[:3]}); keeping the first of each pair.")

    bg2 = list(dict.fromkeys(g.upper() for g in bgs))
    gd2 = {k: [g.upper() for g in v] for k, v in gene_dict.items()}
    gs2 = {lib_n: {t: [g.upper() for g in gg] for t, gg in gs.items()}
           for lib_n, gs in gene_sets.items()}
    return bg2, gd2, gs2, mapa, modo


# --------------------------------------------------------------------------
# 4. the test
# --------------------------------------------------------------------------

def run_goea(
    gene_dict: Mapping[str, Sequence[str]],
    background: Sequence[str],
    libraries: Optional[Sequence[str]] = None,
    organism: str = "Mouse",
    gene_sets: Optional[Mapping[str, Mapping]] = None,
    min_set_size: int = 5,
    max_set_size: int = 500,
    match_case: str = "auto",
    global_fdr: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Per-population ORA against a custom background, with gp.enrich() (not gp.enrichr()).

    gene_sets
        Optional, {library_name: {term: [genes]}}. If not given, they are
        downloaded with gp.get_library(), which requires network access. Save
        them and pass them in for repeated runs: that way the analysis is
        reproducible even if Enrichr changes version, which it does.

    min_set_size / max_set_size
        Applied AFTER intersecting each term with the background, which is
        where they make sense: a 300-gene term of which only 4 are in the
        tissue is, for this test, a 4-gene term. Huge terms
        ('biological process') discriminate nothing and tiny ones give
        unstable p-values.

    global_fdr
        gseapy corrects for multiple testing WITHIN each call. Running 8
        populations x 3 libraries gives 24 families of tests, and the correction
        of each one knows nothing of the others. With True, 'padj_global' is added
        with Benjamini-Hochberg over EVERYTHING. Use that column when stating
        "these are the significant terms".
    """
    import gseapy as gp
    from statsmodels.stats.multitest import multipletests

    if type(background) is not list:
        bg = list(dict.fromkeys(background))
    else: 
        bg = background
    bg_set = set(bg)

    if gene_sets is None:
        if libraries is None:
            libraries = latest_go_libraries(organism=organism, verbose=verbose)
        gene_sets = {}
        for lib in libraries:
            if verbose:
                print(f"[goea] downloading {lib} ({organism})...")
            gene_sets[lib] = gp.get_library(name=lib, organism=organism)

    # Enrichr symbols and the dataset symbols may be in different alphabets
    bg, gene_dict, gene_sets, mapa_case, modo_case = _resolver_nomenclatura(
        bg, gene_dict, gene_sets, match_case=match_case, verbose=verbose)
    bg_set = set(bg)

    filas = []
    for lib, gs in gene_sets.items():
        coverage_report(bg, gs, nombre=lib, verbose=verbose)

        # trim the terms to the background: the test universe is the background
        gs_bg = {}
        for term, genes in gs.items():
            inter = [g for g in set(genes) if g in bg_set]
            if min_set_size <= len(inter) <= max_set_size:
                gs_bg[term] = inter
        if verbose:
            print(f"[goea] {lib}: {len(gs_bg)} of {len(gs)} terms with "
                  f"{min_set_size}-{max_set_size} genes within the background")
        if not gs_bg:
            warnings.warn(f"{lib}: no usable term after intersecting "
                          f"with the background. Check the coverage.")
            continue

        for pop, genes in gene_dict.items():
            genes = [g for g in dict.fromkeys(genes) if g in bg_set]
            if len(genes) < 5:
                warnings.warn(f"{pop}/{lib}: only {len(genes)} foreground genes "
                              f"are in the background; skipping it")
                continue
            try:
                enr = gp.enrich(gene_list=genes, gene_sets=gs_bg, background=bg,
                                outdir=None, no_plot=True, verbose=False)
            except Exception as e:
                warnings.warn(f"{pop}/{lib} failed: {e}")
                continue
            if enr is None or enr.results is None or len(enr.results) == 0:
                continue
            df = enr.results.copy()
            df["population"] = pop
            df["libreria"] = lib
            df["n_genes_disponibles"] = len(genes)
            filas.append(df)

    if not filas:
        raise RuntimeError(
            "No test produced results. See the [cov] and [case] lines "
            "above: if the overlap is still low even in upper case, "
            "the problem is not capitalisation but that the symbols are not from "
            "the species/nomenclature of the library."
        )

    res = pd.concat(filas, ignore_index=True)

    ren = {"Term": "termino", "Overlap": "solapamiento", "P-value": "pval",
           "Adjusted P-value": "padj", "Odds Ratio": "odds_ratio",
           "Combined Score": "score", "Genes": "genes"}
    res = res.rename(columns={k: v for k, v in ren.items() if k in res.columns})

    if mapa_case is not None and "genes" in res:
        # return the symbols to the user's nomenclature
        res["genes"] = res["genes"].astype(str).apply(
            lambda s: ";".join(mapa_case.get(x, x) for x in s.split(";")))

    if global_fdr and "pval" in res:
        res["padj_global"] = multipletests(res["pval"].to_numpy(), method="fdr_bh")[1]

    orden = [c for c in ("population", "libreria", "termino", "solapamiento",
                         "odds_ratio", "pval", "padj", "padj_global", "genes",
                         "n_primer_plano") if c in res.columns]
    res = res[orden + [c for c in res.columns if c not in orden]]
    res = res.sort_values(["population", "padj_global" if global_fdr else "padj"])
    # provenance: which GO version was used. GO changes between versions and
    # two analyses with different libraries are not comparable.
    res.attrs["librerias"] = list(gene_sets.keys())
    res.attrs["n_fondo"] = len(bg)
    res.attrs["organism"] = organism
    res.attrs["match_case"] = modo_case

    if verbose:
        col = "padj_global" if global_fdr else "padj"
        sig = res[res[col] < 0.05]
        print(f"\n[goea] {len(sig)} terms with {col} < 0.05, out of {len(res)} tested")
        print(sig.groupby("population").size().to_string()
              if len(sig) else "  (none)")
        print("\n[goea] Remember: the DEGs come from clusters defined on "
              "these same data. This ranks hypotheses; it does not confirm them.")
    return res


# --------------------------------------------------------------------------
# 5. does the conclusion change if the background changes?
# --------------------------------------------------------------------------

def background_sensitivity(
    gene_dict: Mapping[str, Sequence[str]],
    backgrounds: Mapping[str, Sequence[str]],
    alpha: float = 0.05,
    col_p: str = "padj_global",
    verbose: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """
    Runs the same ORA with SEVERAL backgrounds and returns which terms survive
    all of them.

    Debating which background is correct is interesting, but the question that
    really matters is whether the conclusion depends on that choice. If a term
    comes out with all three backgrounds, the debate is theoretical. If it only
    comes out with one, the result is an artefact of a nuisance decision, and that
    is itself the finding.

        backgrounds = {
            "detected_global": build_background(adata, min_cells=10),
            "detected_per_pop": build_background(adata, min_cells=5,
                                                 group_key="minor_population"),
            "hvg": build_background(adata, use_hvg=True),
        }
        rob = background_sensitivity(degs, backgrounds, gene_sets=gs)

    Returns one row per (population, term) with the p-value under each
    background and 'n_fondos_sig' = the number of backgrounds in which it falls below alpha.
    """
    trozos = {}
    for nombre, bg in backgrounds.items():
        if verbose:
            print(f"\n===== background '{nombre}' ({len(set(bg))} genes) =====")
        r = run_goea(gene_dict, bg, verbose=verbose, **kwargs)
        trozos[nombre] = r.set_index(["population", "libreria", "termino"])[col_p]

    tabla = pd.concat(trozos, axis=1)
    tabla.columns = list(trozos.keys())
    tabla["n_fondos_sig"] = (tabla < alpha).sum(axis=1)
    tabla["n_fondos_testado"] = tabla[list(trozos)].notna().sum(axis=1)
    tabla = tabla.sort_values(["n_fondos_sig", list(trozos)[0]],
                              ascending=[False, True])

    if verbose:
        n = len(backgrounds)
        robustos = tabla[tabla["n_fondos_sig"] == n]
        frag = tabla[(tabla["n_fondos_sig"] > 0) & (tabla["n_fondos_sig"] < n)]
        print(f"\n[sens] {len(robustos)} terms significant with ALL {n} "
              f"backgrounds; {len(frag)} only with some.")
        if len(frag):
            print("[sens] Those in the second list depend on a nuisance "
                  "decision. Do not report them without stating which background yields them.")
        print(tabla.head(20).round(4).to_string())
    return tabla


# --------------------------------------------------------------------------
# 6. redundancy: GO terms are not independent
# --------------------------------------------------------------------------

def collapse_redundant(
    res: pd.DataFrame,
    col_p: str = "padj_global",
    alpha: float = 0.05,
    max_jaccard: float = 0.5,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Removes terms that are significant BECAUSE OF THE SAME GENES as a better one.

    Obtaining 1,500 terms does not mean having found 1,500 things. GO is
    a hierarchical graph: 'muscle cell differentiation' is contained in
    'muscle organ development', and if 20 muscle genes trigger both, that
    is ONE result counted twice. The number of significant terms is not a
    measure of anything.

    Strategy (greedy, per population): sort by p-value and keep a term only
    if the gene set that makes it significant (the 'genes' column, which is
    the intersection with the foreground, not the whole term) overlaps less
    than max_jaccard with that of every term already kept.

    This is a poor man's REVIGO: it groups by shared genes, not by semantic
    similarity in the ontology. In exchange it needs nothing external and is
    exactly the criterion that matters here, namely "are they the same genes?".
    """
    if col_p not in res:
        col_p = "padj"
    sig = res[res[col_p] < alpha].copy()
    if not len(sig):
        if verbose:
            print(f"[collapse] no term with {col_p} < {alpha}")
        return sig

    sig["_set"] = sig["genes"].astype(str).apply(lambda x: frozenset(x.split(";")))

    guardados = []
    for pop, sub in sig.groupby("population", sort=False):
        sub = sub.sort_values(col_p)
        elegidos = []
        for i, fila in sub.iterrows():
            s = fila["_set"]
            if not s:
                continue
            redundante = any(
                len(s & t) / max(len(s | t), 1) > max_jaccard for t in elegidos)
            if not redundante:
                elegidos.append(s)
                guardados.append(i)

    out = sig.loc[guardados].drop(columns=["_set"])
    if verbose:
        antes = sig.groupby("population").size()
        despues = out.groupby("population").size()
        tabla = pd.DataFrame({"significativos": antes, "no_redundantes": despues})
        tabla["reduction"] = (1 - tabla["no_redundantes"] / tabla["significativos"]).round(2)
        print(tabla.to_string())
        print(f"\n[collapse] {len(sig)} -> {len(out)} terms "
              f"(gene Jaccard > {max_jaccard} = redundant)")
        print("[collapse] This is NOT a statistical correction: the p-values do not "
              "change. It only stops counting the same finding several times.")
    return out


# --------------------------------------------------------------------------
# 7. plotting
# --------------------------------------------------------------------------

def plot_goea(res: pd.DataFrame, n_top: int = 6, libreria: Optional[str] = None,
              col_p: str = "padj_global", figsize=(11, 8)):
    """
    Dotplot: terms (rows) x populations (columns), size = overlap,
    colour = -log10(p).

    The n_top terms of EACH population are shown, not the global n_top:
    otherwise the population with the longest lists takes over the whole figure.
    """
    import matplotlib.pyplot as plt

    df = res if libreria is None else res[res["libreria"] == libreria]
    if col_p not in df:
        col_p = "padj"

    elegidos = (df.sort_values(col_p).groupby("population").head(n_top)["termino"]
                .drop_duplicates().tolist())
    sub = df[df["termino"].isin(elegidos)].copy()

    sub["k"] = sub["solapamiento"].astype(str).str.split("/").str[0].astype(float)
    sub["mlog"] = -np.log10(sub[col_p].clip(lower=1e-300))

    pops = sorted(sub["population"].unique())
    terms = (sub.groupby("termino")["mlog"].max().sort_values().index.tolist())
    xi = {p: i for i, p in enumerate(pops)}
    yi = {t: i for i, t in enumerate(terms)}

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    sc_ = ax.scatter([xi[p] for p in sub["population"]],
                     [yi[t] for t in sub["termino"]],
                     s=15 * sub["k"], c=sub["mlog"], cmap="magma_r",
                     edgecolors="none")
    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels(pops, rotation=45, ha="right")
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels([t[:70] for t in terms], fontsize=8)
    ax.set(xlabel="", ylabel="")
    ax.grid(alpha=.2)
    fig.colorbar(sc_, ax=ax, label=f"-log10({col_p})", shrink=.5)
    return fig, ax
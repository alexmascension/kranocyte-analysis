"""
Cell-cell communication: CellPhoneDB v5 (statistical method) and LIANA+.

TWO ROUTES, AND THEY ARE NOT INTERCHANGEABLE
--------------------------------------------
CellPhoneDB (sections 1-4): HUMAN database; your symbols must be ortholog-mapped
first, and it expects the normalised matrix WITHOUT log. In exchange it handles
multi-subunit complexes, which for Il6 (IL6 + IL6R + IL6ST) is exactly what matters.

LIANA+ (section 5): NATIVE mouse resource ('mouseconsensus'), no orthology
step, expects the LOG-normalised matrix, and allows running several resources
and seeing what survives. This is the route I recommend starting with.

Note that the scale of the matrix is OPPOSITE between the two. It is the
easiest copy-paste error to make here.

    mapa  = mouse_to_human(adata.var_names)              # once, it is cached
    paths = prepare_cpdb_input(adata, "cell_type", OUT, ortholog_map=mapa)
    res   = run_cpdb_statistical(**paths, cpdb_file_path=DB)
    il    = filter_interactions(res, ["Il6", "Il1b"], ortholog_map=mapa)

CELLPHONEDB IS HUMAN
--------------------
The database consists of HUMAN ligand-receptor pairs, with human symbols. Your
data are from mouse. There is no official murine version, so symbols must be
translated by orthology BEFORE running anything, and that is not a formatting
detail: it is a strong biological assumption. That IL6-IL6R interact in human
does not guarantee that Il6-Il6ra do the same in mouse, and there are human
pairs whose murine ortholog simply does not exist (CXCL8/IL8 is the classic case).

What must NOT be done is to upper-case the symbols and carry on. It works
for most genes (Il6 -> IL6, Il1b -> IL1B) and fails silently exactly where the
name differs between species (Trp53 is not TP53), which are often interesting
genes as well. mouse_to_human() uses Ensembl orthology via pybiomart; the
'upper' mode exists only as a last resort and emits a warning.

DO NOT SUBSET GENES BEFORE RUNNING
---------------------------------
You want to look at Il6 and Il1b, but the statistical method permutes the cell
labels and compares against the null distribution of ALL interactions, and it
also needs the full matrix to evaluate complexes (a multi-subunit receptor
only counts if all subunits are expressed). If you pass it two genes, the
null is computed over two genes and the result means nothing. Run it on
everything and filter the RESULT: that is what filter_interactions() is for.

ON THE P-VALUES
---------------
They come from permuting cell labels, so their minimum resolution is
1/iterations: with 1000 permutations there is no p < 0.001. CellPhoneDB applies
no multiple-testing correction, and there are thousands of pairs x cell-type
combinations. And as always in this pipeline: with 2 mice the cells are
pseudoreplicates, so the test answers "is this mean larger than when shuffling
the labels of THESE cells?", not "does this happen in the mouse?".

WARNING
-------
The cellphonedb calls in this module are NOT tested: the library could not
be installed in the environment where this was written (fbpca does not compile).
The orchestration (orthology, input format, result filtering) has been
tested with mocks. Check the first run by hand.
"""

from __future__ import annotations

import os
import warnings
from typing import Iterable, Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


# --------------------------------------------------------------------------
# 1. mouse -> human orthology
# --------------------------------------------------------------------------

MIRRORS_ENSEMBL = ("http://www.ensembl.org", "http://useast.ensembl.org",
                   "http://asia.ensembl.org")

URL_MGI_HOM = ("https://www.informatics.jax.org/downloads/reports/"
               "HOM_MouseHumanSequence.rpt")


def _orto_biomart(reintentos: int = 4, espera: float = 5.0,
                  mirrors: Sequence[str] = MIRRORS_ENSEMBL,
                  verbose: bool = True) -> pd.DataFrame:
    """
    Orthology from Ensembl with retries and mirror rotation.

    BioMart returns HTTP 429 (too many requests) quite easily, and
    pybiomart makes TWO requests per query (one for the dataset
    configuration and another for the data), so the probability of hitting the
    limit is twice what it seems. It is not an error in your data or in the
    code: the request is retried with increasing waits and, if the mirror keeps
    refusing, the mirror is switched.
    """
    import time

    from pybiomart import Dataset

    ultimo = None
    for host in mirrors:
        for intento in range(reintentos):
            try:
                ds = Dataset(name="mmusculus_gene_ensembl", host=host)
                df = ds.query(attributes=[
                    "external_gene_name",
                    "hsapiens_homolog_associated_gene_name",
                    "hsapiens_homolog_orthology_type",
                ])
                df.columns = ["mouse", "human", "kind"]
                if verbose:
                    print(f"[orto] Ensembl OK ({host})")
                return df
            except Exception as e:
                ultimo = e
                pausa = espera * (2 ** intento)
                if verbose:
                    print(f"[orto] {host} attempt {intento+1}/{reintentos} "
                          f"failed ({type(e).__name__}); waiting {pausa:.0f}s")
                time.sleep(pausa)
    raise RuntimeError(
        f"BioMart does not respond on any of the mirrors {list(mirrors)}. "
        f"Last error: {ultimo}. Try method='mgi', which is a static file "
        f"with no request limit."
    ) from ultimo


def _orto_mgi(url: str = URL_MGI_HOM, verbose: bool = True) -> pd.DataFrame:
    """
    Orthology from the MGI homology report.

    It is a static file, with no API and no request limit, and it is the source
    used by the mouse community itself. Each homology class groups the
    orthologous genes of the two species; we keep the classes that have
    EXACTLY one mouse gene and one human gene, which is the operational
    definition of one-to-one.

    NOT VERIFIED from the environment where this was written: informatics.jax.org
    is blocked there, so I could not check the exact column names of the
    current file. That is why they are detected instead of hard-coded, and
    if detection fails the error tells you which columns it found.
    """
    df = pd.read_csv(url, sep="\t", dtype=str)

    def _buscar(*claves):
        for c in df.columns:
            cl = c.lower()
            if all(k in cl for k in claves):
                return c
        return None

    col_key = _buscar("class", "key") or _buscar("homologene")
    col_org = _buscar("organism") or _buscar("species")
    col_sym = "Symbol" if "Symbol" in df.columns else _buscar("symbol")
    if not all((col_key, col_org, col_sym)):
        raise RuntimeError(
            f"Unrecognised MGI file format. Columns: "
            f"{list(df.columns)}. Detected: key={col_key}, org={col_org}, "
            f"symbol={col_sym}."
        )

    org = df[col_org].str.lower()
    df = df.assign(_esp=np.where(org.str.contains("human"), "human",
                    np.where(org.str.contains("mouse"), "mouse", None)))
    df = df.dropna(subset=["_esp"])

    cuenta = df.groupby([col_key, "_esp"]).size().unstack(fill_value=0)
    claves = cuenta.index[(cuenta.get("mouse", 0) == 1) & (cuenta.get("human", 0) == 1)]
    sub = df[df[col_key].isin(claves)]
    piv = sub.pivot_table(index=col_key, columns="_esp", values=col_sym,
                          aggfunc="first")
    out = piv.reset_index()[["mouse", "human"]].dropna()
    if verbose:
        print(f"[orto] MGI: {len(out)} one-to-one pairs")
    return out


def mouse_to_human(
    genes: Optional[Iterable[str]] = None,
    cache_path: str = "mouse_human_orthologs.tsv",
    method: str = "mgi",
    one_to_one: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Mouse symbol -> human symbol table. DataFrame with ['mouse', 'human'].

    It is cached on disk: the first call downloads, subsequent calls read the
    file. Keep it with the analysis, as orthology annotations change.

    method
        'mgi'      MGI homology report. Static file, with no API and no
                   request limits. It is the default for that reason: BioMart
                   returns 429 easily and leaves the analysis hanging.
        'biomart'  Ensembl via pybiomart, with retries and mirror rotation.
        'upper'    convert to UPPER CASE. Last resort, no network. Matches
                   most symbols and FAILS SILENTLY on those that differ
                   between species (Trp53 != TP53). Warns when used.

    one_to_one
        Discards genes with several orthologs in the other direction. This is
        the conservative choice: with a one-to-many mapping one has to decide which
        human gene receives the murine expression, and that decision would end
        up being made by the row order.
    """
    if os.path.exists(cache_path):
        mapa = pd.read_csv(cache_path, sep="\t", dtype=str)
        if verbose:
            print(f"[orto] read from cache: {cache_path} ({len(mapa)} pairs)")
        return mapa

    if method == "upper":
        if genes is None:
            raise ValueError("method='upper' requires the list of genes")
        mapa = pd.DataFrame({"mouse": list(genes)})
        mapa["human"] = mapa["mouse"].str.upper()
        warnings.warn(
            "Orthology by UPPER CASE. Matches most symbols but "
            "not those that differ between species (Trp53 != TP53), and does not "
            "report which ones were lost. Use it only if you cannot download the real "
            "table, and state it in the methods."
        )
    elif method == "mgi":
        mapa = _orto_mgi(verbose=verbose)
    elif method == "biomart":
        df = _orto_biomart(verbose=verbose)
        df = df.dropna(subset=["mouse", "human"])
        df = df[df["human"] != ""]
        if one_to_one:
            n0 = len(df)
            df = df[df["kind"] == "ortholog_one2one"]
            if verbose:
                print(f"[orto] {len(df)} of {n0} pairs are one-to-one")
        mapa = df[["mouse", "human"]].drop_duplicates()
    else:
        raise ValueError("method must be 'mgi', 'biomart' or 'upper'")

    if one_to_one and method != "upper":
        mapa = mapa[~mapa["mouse"].duplicated(keep=False)]
        mapa = mapa[~mapa["human"].duplicated(keep=False)]

    mapa.to_csv(cache_path, sep="\t", index=False)
    if verbose:
        print(f"[orto] {len(mapa)} pairs saved to {cache_path}")
        if genes is not None:
            g = set(genes)
            print(f"[orto] covers {len(g & set(mapa['mouse']))/max(len(g),1):.1%} "
                  f"of your genes")
            faltan = [x for x in ("Il6", "Il1b") if x in g and x not in set(mapa["mouse"])]
            if faltan:
                print(f"[orto] CAUTION: {faltan} have no one-to-one ortholog in "
                      f"this table. Check it by hand before continuing.")
    return mapa


# --------------------------------------------------------------------------
# 2. prepare the input
# --------------------------------------------------------------------------

def prepare_cpdb_input(
    adata: ad.AnnData,
    group_key: str,
    outdir: str,
    ortholog_map: pd.DataFrame,
    layer: Optional[str] = None,
    log_layer: Optional[str] = "norm_scran_log1p",
    counts_layer: str = "counts_cellbender",
    min_cells: int = 10,
    exclude_groups: Sequence[str] = (),
    verbose: bool = True,
) -> dict:
    """
    Writes meta.tsv and counts.h5ad ready for CellPhoneDB. Returns the paths.

    WHICH MATRIX. CellPhoneDB averages expression per cell type and the
    statistic is the mean of the two means. A mean of log values is NOT
    the log of the mean: log compresses high values, so a ligand highly
    expressed in few cells weighs differently on one scale and on the other. That
    is why the normalised matrix WITHOUT log is supplied here:

      layer      use it if you already have a normalised non-log layer.
      log_layer  otherwise, it is reconstructed with expm1() from the log-normalised layer
                 (exact: it undoes the log1p you applied).
      counts_layer  last resort: library-size normalisation on counts.

    The CellPhoneDB tutorials are not consistent on this and you will see people
    passing log data. If you want to replicate someone, pass their layer via
    `layer=` and record it; but by default it goes without log, which is consistent
    with what the method computes.

    min_cells
        Very small populations give unstable means and permutation does not
        fix that: they are excluded and a message is printed. This is not the same as not having
        the population; it means you cannot say anything about it here.
    """
    os.makedirs(outdir, exist_ok=True)

    # --- cells
    obs = adata.obs[group_key].astype(str)
    tam = obs.value_counts()
    fuera = set(exclude_groups) | set(tam.index[tam < min_cells])
    if fuera and verbose:
        print(f"[cpdb] excluded {sorted(fuera)} "
              f"({dict(tam[list(fuera)])} cells, min_cells={min_cells})")
    keep_cells = ~obs.isin(fuera).to_numpy()

    # --- matrix without log
    if layer is not None:
        M = adata.layers[layer]
        origen = f"layer '{layer}'"
    elif log_layer and log_layer in adata.layers:
        M = adata.layers[log_layer]
        M = M.copy()
        if sp.issparse(M):
            M.data = np.expm1(M.data)
        else:
            M = np.expm1(M)
        origen = f"expm1('{log_layer}')"
    elif log_layer and adata.uns.get("normalization", {}).get("scran") is not None:
        # DEPOSITED objects do not carry the log layer: it was removed because it is
        # identical to X (see kranocyte.submissions.SlimSpec). Falling through to the
        # library-size branch here would be a SILENT change of method -- a different
        # normalisation with no warning. If uns['normalization'] records scran, X is
        # that same log-normalised matrix and expm1(X) is the right thing to use.
        M = adata.X.copy()
        if sp.issparse(M):
            M.data = np.expm1(M.data)
        else:
            M = np.expm1(M)
        origen = f"expm1(X) -- '{log_layer}' absent; X is the scran+log1p matrix"
        if verbose:
            print(f"[cpdb] layer '{log_layer}' is not in the object; using X, which "
                  f"uns['normalization'] records as scran+log1p")
    else:
        if log_layer:
            warnings.warn(
                f"prepare_cpdb_input: log_layer='{log_layer}' was requested but it is not in "
                f"adata.layers and uns['normalization'] does not record scran, so there is no "
                f"way to tell whether X is that matrix. Falling back to library-size "
                f"normalisation over '{counts_layer}', which is NOT the same method. Pass "
                f"layer= or log_layer=None explicitly if that is what you want.",
                stacklevel=2,
            )
        M = adata.layers[counts_layer] if counts_layer in adata.layers else adata.X
        M = sp.csr_matrix(M, dtype=np.float64)
        libs = np.asarray(M.sum(axis=1)).ravel()
        libs[libs == 0] = 1.0
        M = sp.diags(1e4 / libs) @ M
        origen = f"CPM(1e4) on '{counts_layer}'"
    if verbose:
        print(f"[cpdb] expression matrix: {origen} (without log)")

    # --- orthology
    mapa = dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
    en_mapa = np.array([g in mapa for g in adata.var_names])
    if verbose:
        print(f"[cpdb] {en_mapa.sum()} of {adata.n_vars} genes have a human "
              f"ortholog ({en_mapa.mean():.1%})")
        if en_mapa.mean() < 0.4:
            print("[cpdb] WARNING: low coverage. Check the orthology "
                  "table before interpreting anything.")

    sub = ad.AnnData(
        X=sp.csr_matrix(M)[keep_cells][:, en_mapa],
        obs=pd.DataFrame(index=adata.obs_names[keep_cells]),
        var=pd.DataFrame(index=[mapa[g] for g in adata.var_names[en_mapa]]),
    )
    if not sub.var_names.is_unique:
        n_dup = int(sub.var_names.duplicated().sum())
        warnings.warn(f"{n_dup} duplicated human symbols after mapping; "
                      f"keeping the first occurrence of each")
        sub = sub[:, ~sub.var_names.duplicated()].copy()

    # --- check genes of interest
    if verbose:
        for g in ("Il6", "Il1b"):
            if g in mapa:
                print(f"[cpdb] {g} -> {mapa[g]}: "
                      f"{'present' if mapa[g] in sub.var_names else 'ABSENT'}")

    meta = pd.DataFrame({"Cell": sub.obs_names, "cell_type": obs[keep_cells].values})
    meta_path = os.path.join(outdir, "meta.tsv")
    counts_path = os.path.join(outdir, "counts.h5ad")
    meta.to_csv(meta_path, sep="\t", index=False)
    sub.write_h5ad(counts_path)

    if verbose:
        print(f"[cpdb] {sub.n_obs} cells x {sub.n_vars} genes -> {counts_path}")
        print(meta["cell_type"].value_counts().to_string())
    return {"meta_file_path": meta_path, "counts_file_path": counts_path}


# --------------------------------------------------------------------------
# 2b. look at what is INSIDE the database
# --------------------------------------------------------------------------

def inspect_cpdb_database(cpdb_file_path: str, verbose: bool = True) -> dict:
    """
    Opens the CellPhoneDB .zip and returns its tables as {name: DataFrame}.

    It does not use the library API: the .zip contains CSVs and they are read directly. This way
    it works even if the API changes between versions, and you can inspect the database
    BEFORE running anything.

    Four tables matter and they are related to one another:
      gene         gene_name / hgnc_symbol / uniprot / ensembl
      protein      one row per uniprot, with whether it is a receptor, secreted, etc.
      complex      complexes, defined as uniprot_1..uniprot_4
      interaction  the pairs, where each partner is a uniprot OR the name of
                   a complex
    That is why searching for "IL6" in the interaction table is not enough: the partners are not
    symbols.
    """
    import zipfile

    tablas = {}
    with zipfile.ZipFile(cpdb_file_path) as z:
        nombres = [n for n in z.namelist() if n.lower().endswith(".csv")]
        for n in nombres:
            clave = os.path.basename(n).replace(".csv", "")
            try:
                with z.open(n) as fh:
                    tablas[clave] = pd.read_csv(fh, dtype=str)
            except Exception as e:
                warnings.warn(f"Could not read {n}: {e}")
    if verbose:
        print(f"[db] {os.path.basename(cpdb_file_path)}: {len(tablas)} tables")
        for k, v in sorted(tablas.items()):
            print(f"  {k:35s} {v.shape[0]:>7} x {v.shape[1]}")
    return tablas


def _tabla(db: Mapping[str, pd.DataFrame], *claves,
           requiere: Sequence[str] = ()) -> Optional[pd.DataFrame]:
    """
    Looks up a table by name, choosing the BEST candidate rather than the first.

    This matters because the zip contains several tables whose name contains 'gene'
    (gene_input, gene_synonym_to_gene_name, ...) and taking the first one in the
    zip listing is a lottery that depends on the order of the files. Candidates
    are scored: first those that have the required columns, then those with the
    shortest name (gene_input beats gene_synonym_to_gene_name).
    """
    cands = [(k, v) for k, v in db.items()
             if all(c in k.lower() for c in claves)]
    if not cands:
        return None

    def _puntua(kv):
        k, v = kv
        cols = [c.lower() for c in v.columns]
        tiene = sum(any(r in c for c in cols) for r in requiere)
        return (-tiene, len(k))

    cands.sort(key=_puntua)
    return cands[0][1]


def cpdb_gene_coverage(
    db: Mapping[str, pd.DataFrame],
    genes: Sequence[str],
    ortholog_map: Optional[pd.DataFrame] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Which of YOUR genes exist in CellPhoneDB?

    Run this BEFORE running the analysis. If a gene is not in the database, it cannot appear in
    any result, and its absence from the output does not mean "does not interact":
    it means "CellPhoneDB does not cover it". These are very different things and the
    second cannot be read from the results table.

    genes in mouse nomenclature if you pass ortholog_map; otherwise, human.
    """
    gt = _tabla(db, "gene", requiere=("uniprot", "gene_name", "hgnc"))
    if gt is None:
        raise RuntimeError(f"Gene table not found. Tables: {list(db)}")

    cols = [c for c in gt.columns
            if any(k in c.lower() for k in ("gene_name", "hgnc", "symbol"))]
    if not cols:
        raise RuntimeError(f"The gene table has no symbols. Columns: "
                           f"{list(gt.columns)}")
    en_db = set()
    for c in cols:
        en_db |= set(gt[c].dropna().astype(str))

    if ortholog_map is not None:
        mapa = dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
    else:
        mapa = {}

    filas = []
    for g in genes:
        h = mapa.get(g, g if not mapa else None)
        filas.append({"gene": g, "human": h,
                      "has_ortholog": h is not None,
                      "in_cellphonedb": bool(h) and h in en_db})
    out = pd.DataFrame(filas)

    if verbose:
        n = len(out)
        print(f"[db] {out['has_ortholog'].sum()}/{n} with a human ortholog; "
              f"{out['in_cellphonedb'].sum()}/{n} present in CellPhoneDB")
        fuera = out[~out["in_cellphonedb"]]
        if len(fuera) and len(fuera) <= 30:
            print("[db] NOT in the database: "
                  + ", ".join(f"{r.gen}" + ("" if r.tiene_ortologo else " (no ortholog)")
                              for r in fuera.itertuples()))
    return out


def cpdb_interactions_for(
    db: Mapping[str, pd.DataFrame],
    genes: Sequence[str],
    ortholog_map: Optional[pd.DataFrame] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    ALL interactions in the database that involve these genes, regardless
    of your data.

    This is the ceiling: the maximum number of pairs you could possibly recover. If
    IL6 takes part in 6 interactions and your result returns 2, you know that the
    filter was your data; if it returns 6, you have exhausted what the database
    covers. Without this number, the results table cannot be interpreted.

    Resolves the chain symbol -> uniprot -> complexes containing it ->
    interactions in which the uniprot or any of those complexes takes part. Searching
    for the symbol directly in the interaction table does not work: there the
    partners are uniprots and complex names.
    """
    gt = _tabla(db, "gene", requiere=("uniprot", "gene_name", "hgnc"))
    ct = _tabla(db, "complex", requiere=("complex_name", "uniprot_1"))
    it_ = _tabla(db, "interaction", requiere=("partner_a", "partner_b"))
    if gt is None or it_ is None:
        raise RuntimeError(f"Required tables are missing. Tables: {list(db)}")

    if ortholog_map is not None:
        mapa = dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
        humanos = {mapa[g] for g in genes if g in mapa}
        perdidos = [g for g in genes if g not in mapa]
        if perdidos:
            warnings.warn(f"No human ortholog: {perdidos}")
    else:
        humanos = set(genes)

    col_sym = [c for c in gt.columns
               if any(k in c.lower() for k in ("gene_name", "hgnc", "symbol"))]
    col_uni = next((c for c in gt.columns if "uniprot" in c.lower()), None)
    if col_uni is None:
        raise RuntimeError(f"The gene table has no uniprot column: {list(gt.columns)}")

    sel = pd.Series(False, index=gt.index)
    for c in col_sym:
        sel |= gt[c].astype(str).isin(humanos)
    uniprots = set(gt.loc[sel, col_uni].dropna().astype(str))

    # complexes that contain any of those uniprots
    complejos = set()
    if ct is not None:
        cols_u = [c for c in ct.columns if c.lower().startswith("uniprot")]
        col_nom = next((c for c in ct.columns if "complex_name" in c.lower()), None)
        if cols_u and col_nom:
            m = pd.Series(False, index=ct.index)
            for c in cols_u:
                m |= ct[c].astype(str).isin(uniprots)
            complejos = set(ct.loc[m, col_nom].dropna().astype(str))

    partners = uniprots | complejos
    col_pa = next((c for c in it_.columns if c.lower() in ("partner_a", "multidata_1_id", "partner a")), None)
    col_pb = next((c for c in it_.columns if c.lower() in ("partner_b", "multidata_2_id", "partner b")), None)
    if col_pa is None or col_pb is None:
        cand = [c for c in it_.columns if "partner" in c.lower()]
        if len(cand) >= 2:
            col_pa, col_pb = cand[0], cand[1]
        else:
            raise RuntimeError(f"Partner columns not found in the "
                               f"interaction table: {list(it_.columns)}")

    hit = it_[it_[col_pa].astype(str).isin(partners)
              | it_[col_pb].astype(str).isin(partners)].copy()

    # translate uniprot -> symbol so that the table is readable
    u2s = {}
    for c in col_sym:
        u2s.update(dict(zip(gt[col_uni].astype(str), gt[c].astype(str))))
    hit["symbol_a"] = hit[col_pa].astype(str).map(u2s).fillna(hit[col_pa])
    hit["symbol_b"] = hit[col_pb].astype(str).map(u2s).fillna(hit[col_pb])

    if verbose:
        print(f"[db] tables used: gene={gt.shape}, "
              f"complex={None if ct is None else ct.shape}, "
              f"interaction={it_.shape}")
        print(f"[db] {sorted(humanos)} -> {len(uniprots)} uniprots, "
              f"{len(complejos)} complexes containing them")
        if complejos:
            print(f"[db] complexes: {sorted(complejos)[:8]}")
        print(f"[db] {len(hit)} interactions in the database involve them")
        cols_ver = [c for c in ("symbol_a", "symbol_b", "annotation_strategy",
                                "classification", "directionality", "is_integrin")
                    if c in hit.columns]
        if len(hit):
            print(hit[cols_ver].head(25).to_string(index=False))
    return hit


# --------------------------------------------------------------------------
# 3. the statistical method
# --------------------------------------------------------------------------

def download_cpdb_database(target_dir: str, version: str = "v5.0.0",
                           verbose: bool = True) -> str:
    """Downloads the CellPhoneDB database. Returns the path to the .zip."""
    from cellphonedb.utils import db_utils

    os.makedirs(target_dir, exist_ok=True)
    db_utils.download_database(target_dir, version)
    zips = [f for f in os.listdir(target_dir) if f.endswith(".zip")]
    if not zips:
        raise RuntimeError(f"No .zip was downloaded into {target_dir}")
    ruta = os.path.join(target_dir, sorted(zips)[-1])
    if verbose:
        print(f"[cpdb] database: {ruta}")
    return ruta


def run_cpdb_statistical(
    cpdb_file_path: str,
    meta_file_path: str,
    counts_file_path: str,
    output_path: str,
    counts_data: str = "hgnc_symbol",
    iterations: int = 1000,
    threshold: float = 0.1,
    pvalue: float = 0.05,
    threads: int = 4,
    microenvs_file_path: Optional[str] = None,
    score_interactions: bool = True,
    seed: int = 42,
    verbose: bool = True,
) -> dict:
    """
    cpdb_statistical_analysis_method.call(), with the parameters that matter
    made explicit.

    threshold
        MINIMUM fraction of cells of a type that must express the gene for
        it to be considered expressed there. 0.1 is the CellPhoneDB default and it is a
        decision, not a constant: raising it removes interactions supported by
        a minority of cells (which in a heterogeneous population may be
        real) and lowering it fills the result with pairs that depend on three
        cells. Report the value you use.

    iterations
        Number of permutations. It sets the floor of the p-value: with 1000, the minimum
        possible is 0.001 and nothing exists below it. If you are going to compare thousands of
        pairs, that floor is what dominates, not your alpha.

    microenvs_file_path
        Optional TSV (cell_type, microenvironment) that restricts which types
        can interact with each other. If you know that two populations do not share a
        niche, stating it here reduces the number of tests and is more honest than
        filtering it afterwards.
    """
    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

    os.makedirs(output_path, exist_ok=True)
    kwargs = dict(
        cpdb_file_path=cpdb_file_path,
        meta_file_path=meta_file_path,
        counts_file_path=counts_file_path,
        counts_data=counts_data,
        output_path=output_path,
        iterations=iterations,
        threshold=threshold,
        pvalue=pvalue,
        threads=threads,
        debug_seed=seed,
        separator="|",
        microenvs_file_path=microenvs_file_path,
        score_interactions=score_interactions,
    )
    res = cpdb_statistical_analysis_method.call(**kwargs)

    if verbose:
        for k, v in res.items():
            if isinstance(v, pd.DataFrame):
                print(f"[cpdb] {k}: {v.shape}")
        print(f"[cpdb] permutation p-values: floor = {1/iterations:g}. "
              f"CellPhoneDB does NOT correct for multiple testing.")
    return res


# --------------------------------------------------------------------------
# 4. filter the result by specific genes
# --------------------------------------------------------------------------

_META_COLS = ("id_cp_interaction", "interacting_pair", "partner_a", "partner_b",
              "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b",
              "annotation_strategy", "is_integrin", "directionality",
              "classification")


def filter_interactions(
    res: Mapping[str, pd.DataFrame],
    genes: Sequence[str],
    ortholog_map: Optional[pd.DataFrame] = None,
    max_pvalue: float = 0.05,
    include_complexes: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Interactions involving specific genes, in long format.

    genes in MOUSE NOMENCLATURE (Il6, Il1b); they are translated with ortholog_map.
    If no map is passed they are assumed to be human symbols already.

    include_complexes
        In addition to the rows where gene_a/gene_b is exactly your gene, includes
        those whose 'interacting_pair' mentions it: in CellPhoneDB many
        receptors are complexes and the gene appears inside the complex name
        and not in gene_a/gene_b. They are flagged with the 'via_complex' column so that
        you can separate them; they are not the same and conflating them inflates the count.

    Returns one row per (interaction, cell-type pair) with its mean and
    its p-value, already filtered by max_pvalue.
    """
    means = res["means"]
    pvals = res["pvalues"]

    if ortholog_map is not None:
        mapa = dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
        humanos, perdidos = [], []
        for g in genes:
            (humanos.append(mapa[g]) if g in mapa else perdidos.append(g))
        if perdidos:
            warnings.warn(f"No human ortholog: {perdidos}")
    else:
        humanos = list(genes)
    humanos = set(humanos)
    if verbose:
        print(f"[filtro] searching for {sorted(humanos)}")

    meta_cols = [c for c in _META_COLS if c in means.columns]
    par_cols = [c for c in means.columns if c not in meta_cols]

    exacto = (means["gene_a"].isin(humanos) | means["gene_b"].isin(humanos))
    if include_complexes:
        patron = "|".join(sorted(humanos))
        por_nombre = means["interacting_pair"].astype(str).str.contains(
            patron, case=True, regex=True, na=False)
        sel = exacto | por_nombre
    else:
        sel = exacto

    if not sel.any():
        if verbose:
            print("[filtro] no interactions with those genes")
        return pd.DataFrame()

    m = means.loc[sel, meta_cols + par_cols]
    p = pvals.loc[sel, par_cols] if set(par_cols).issubset(pvals.columns) else None

    largo = m.melt(id_vars=meta_cols, var_name="cell_pair", value_name="mean")
    if p is not None:
        pl = pvals.loc[sel, meta_cols + par_cols].melt(
            id_vars=meta_cols, var_name="cell_pair", value_name="pvalue")
        largo = largo.merge(pl, on=meta_cols + ["cell_pair"], how="left")

    largo["via_complex"] = ~(
        largo["gene_a"].isin(humanos) | largo["gene_b"].isin(humanos))

    if "interaction_scores" in res and isinstance(res["interaction_scores"], pd.DataFrame):
        sc = res["interaction_scores"]
        if "interacting_pair" in sc.columns:
            sl = sc[sc["interacting_pair"].isin(largo["interacting_pair"])].melt(
                id_vars=["interacting_pair"], var_name="cell_pair",
                value_name="score")
            largo = largo.merge(sl, on=["interacting_pair", "cell_pair"], how="left")

    if "pvalue" in largo:
        largo = largo[largo["pvalue"] <= max_pvalue]
    largo = largo[largo["mean"] > 0].sort_values(
        ["pvalue", "mean"] if "pvalue" in largo else ["mean"],
        ascending=[True, False] if "pvalue" in largo else [False])

    if verbose:
        print(f"[filtro] {len(largo)} rows (interaction x cell pair) with "
              f"p <= {max_pvalue}")
        n_dir = int((~largo["via_complex"]).sum())
        print(f"[filtro] {n_dir} by direct gene, {len(largo)-n_dir} only by "
              f"appearing in the name of a complex")
        if len(largo):
            print(largo.groupby("interacting_pair").size()
                  .sort_values(ascending=False).head(10).to_string())
    return largo.reset_index(drop=True)


def plot_interactions(largo: pd.DataFrame, n_top: int = 25, figsize=(11, 8)):
    """Dotplot: interaction x cell pair, size = mean, colour = -log10(p)."""
    import matplotlib.pyplot as plt

    if not len(largo):
        raise ValueError("nothing to plot")
    top = (largo.groupby("interacting_pair")["mean"].max()
           .sort_values(ascending=False).head(n_top).index.tolist())
    sub = largo[largo["interacting_pair"].isin(top)].copy()

    pares = sorted(sub["cell_pair"].unique())
    filas = [t for t in top if t in set(sub["interacting_pair"])]
    xi = {p: i for i, p in enumerate(pares)}
    yi = {t: i for i, t in enumerate(filas)}
    col = (-np.log10(sub["pvalue"].clip(lower=1e-4))
           if "pvalue" in sub else sub["mean"])

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    s = ax.scatter([xi[p] for p in sub["cell_pair"]],
                   [yi[t] for t in sub["interacting_pair"]],
                   s=20 + 60 * sub["mean"] / max(sub["mean"].max(), 1e-9),
                   c=col, cmap="magma_r", edgecolors="none")
    ax.set_xticks(range(len(pares)))
    ax.set_xticklabels(pares, rotation=90, fontsize=7)
    ax.set_yticks(range(len(filas)))
    ax.set_yticklabels(filas, fontsize=8)
    ax.grid(alpha=.2)
    fig.colorbar(s, ax=ax, label="-log10(p)" if "pvalue" in sub else "mean",
                 shrink=.5)
    return fig, ax


# ==========================================================================
# LIANA+: several methods and several resources, with consensus
# ==========================================================================
#
# Why LIANA in addition to CellPhoneDB:
#
#   1. NATIVE MOUSE resource ('mouseconsensus'), so the orthology step
#      disappears. With a caveat: part of that resource was built by transferring
#      human annotation, so the assumption is not removed, it is shifted to
#      whoever curated it. Better than doing it yourself by hand; it is not magic.
#   2. Allows running SEVERAL resources and seeing which conclusions survive, which is
#      the only honest way to answer "does it matter which database I use?".
#
# MIND THE MATRIX SCALE. LIANA expects LOG-NORMALISED data, the opposite of
# the standalone CellPhoneDB above, which must be given the matrix without
# log. Your adata.X (norm_scran_log1p) is already fine. What must be set is
# use_raw=False, because LIANA looks at adata.raw by default.

def liana_resources(verbose: bool = True) -> list:
    """Resources that LIANA serves today. Query them, do not guess the list."""
    import liana as li

    rs = list(li.rs.show_resources())
    if verbose:
        print(f"[liana] {len(rs)} resources: {sorted(rs)}")
        murinos = [r for r in rs if "mouse" in r.lower() or "murine" in r.lower()]
        if murinos:
            print(f"[liana] mouse-native: {murinos}")
    return rs


def resource_coverage(
    genes: Sequence[str],
    resources: Sequence[str] = ("mouseconsensus",),
    ortholog_map: Optional[pd.DataFrame] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    In which resources does each gene exist, and with which partners? BEFORE running.

    It is the ceiling of the analysis, just like cpdb_interactions_for(): if Il6
    takes part in 6 pairs of the resource and your run returns 2, the filter was
    your data; if it returns 6, you have exhausted what the resource covers. Without
    this denominator the results table cannot be read, and the absence of
    a pair cannot be distinguished from "it is not annotated".

    ortholog_map
        Only if the resource is HUMAN. With 'mouseconsensus' do NOT translate: the
        symbols are already mouse symbols and translating them would break them.
    """
    import liana as li

    mapa = (dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
            if ortholog_map is not None else {})

    filas = []
    for r in resources:
        df = li.rs.select_resource(r)
        col_l = "ligand" if "ligand" in df.columns else df.columns[0]
        col_r = "receptor" if "receptor" in df.columns else df.columns[1]
        # the subunits of a complex are joined by '_'
        lig = df[col_l].astype(str).str.split("_")
        rec = df[col_r].astype(str).str.split("_")
        for g in genes:
            buscado = mapa.get(g, g)
            como_lig = lig.apply(lambda x: buscado in x)
            como_rec = rec.apply(lambda x: buscado in x)
            sel = como_lig | como_rec
            filas.append({
                "resource": r, "gene": g, "searched": buscado,
                "n_pairs": int(sel.sum()),
                "as_ligand": int(como_lig.sum()),
                "as_receptor": int(como_rec.sum()),
                "partners": "; ".join(sorted(set(
                    df.loc[sel, col_r].astype(str)[como_lig[sel]].tolist()
                    + df.loc[sel, col_l].astype(str)[como_rec[sel]].tolist()
                ))[:12]),
            })
    out = pd.DataFrame(filas)
    if verbose:
        print(out[["resource", "gene", "searched", "n_pairs", "as_ligand",
                   "as_receptor"]].to_string(index=False))
        vacios = out[out["n_pairs"] == 0]
        if len(vacios):
            print(f"\n[liana] WARNING: {list(zip(vacios.gen, vacios.recurso))} do not "
                  f"appear in those resources. Their absence from the result does NOT "
                  f"mean that they do not interact: it means that they are not "
                  f"annotated.")
        for r in out["resource"].unique():
            sub = out[(out.recurso == r) & (out.n_pares > 0)]
            for f in sub.itertuples():
                print(f"\n[liana] {f.gen} in {r}: {f.partners}")
    return out


def run_liana(
    adata: ad.AnnData,
    groupby: str,
    resource: str = "mouseconsensus",
    expr_prop: float = 0.025,
    n_perms: int = 1000,
    use_raw: bool = False,
    key_added: str = "liana_res",
    seed: int = 0,
    n_jobs: int = 4,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    LIANA rank_aggregate: several methods, consensus by rank aggregation.

    TWO RANKINGS, NOT ONE P-VALUE. The output contains:

      magnitude_rank    is this interaction strongly expressed in this
                        pair of cell types?
      specificity_rank  is it MORE SO HERE than in the other pairs?

    They measure different things and conflating them is the classic error. That FAP.3 expresses
    Il1b is not news if everyone expresses it: that is what
    specificity is for. For Il6/Il1b you almost certainly want to sort by
    specificity_rank.

    And the consensus is a RANK AGGREGATION, not a meta-analysis: its p-value
    measures agreement between methods, not biological evidence. Do not report it
    as if it were a p-value of the interaction.

    expr_prop
        Minimum fraction of cells of the type that must express the gene. Same
        role as the CellPhoneDB 'threshold', same status as a decision that
        should be reported.

    n_perms
        Start at 100 to time it. With 9,000 cells and several methods,
        1000 permutations are not free.
    """
    import time

    import liana as li

    if use_raw is False and adata.raw is not None and verbose:
        print("[liana] use_raw=False: adata.X is used. Check that it is "
              "LOG-normalised (LIANA expects it that way, the opposite of the "
              "standalone CellPhoneDB).")

    faltan = [c for c in (groupby,) if c not in adata.obs]
    if faltan:
        raise ValueError(f"Missing obs{faltan}")
    n_grupos = adata.obs[groupby].nunique()
    if verbose:
        print(f"[liana] {resource}, {n_grupos} groups, expr_prop={expr_prop}, "
              f"n_perms={n_perms}")

    t0 = time.time()
    li.mt.rank_aggregate(
        adata, groupby=groupby, resource_name=resource, expr_prop=expr_prop,
        use_raw=use_raw, n_perms=n_perms, seed=seed, n_jobs=n_jobs,
        key_added=key_added, verbose=verbose, inplace=True,
    )
    res = adata.uns[key_added].copy()
    res["resource"] = resource

    if verbose:
        print(f"[liana] {len(res)} interactions in {time.time()-t0:.0f}s")
        print("[liana] sort by 'specificity_rank' for what is specific and by "
              "'magnitude_rank' for what is abundant; they are not the same.")
    return res


# ---------------------------------------------------------------------------
# 5b. Expression context for the LIANA table
# ---------------------------------------------------------------------------
# LIANA returns aggregated scores but NOT the means they are derived from.
# Without them a row cannot be judged: you do not know whether a very low
# `specificity_rank` comes from a receptor well expressed in one population or from a gene
# that is detected in only 30 cells. These functions reconstruct the per-group
# means and the background means, and (this is the important part) VERIFY the
# reconstruction against the columns that LIANA does provide. If my complex
# rule or my logFC definition do not match theirs, the report
# says so instead of leaving you plausible but false columns.

_EPS_CHK = 1e-4


def _matriz_expresion(adata, use_raw=False, layer=None):
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None")
        return adata.raw.X, pd.Index(adata.raw.var_names)
    if layer is not None:
        return adata.layers[layer], pd.Index(adata.var_names)
    return adata.X, pd.Index(adata.var_names)


def _sumas_por_grupo(X, etiquetas, grupos):
    """Sums and sizes per group, for X and for expm1(X)."""
    Xe = X.copy()
    if sp.issparse(Xe):
        Xe.data = np.expm1(Xe.data)
    else:
        Xe = np.expm1(Xe)

    S, Se, n = {}, {}, {}
    for g in grupos:
        m = etiquetas == g
        n[g] = int(m.sum())
        S[g] = np.asarray(X[m].sum(axis=0)).ravel()
        Se[g] = np.asarray(Xe[m].sum(axis=0)).ravel()
    return S, Se, n


def _complejos_por_grupo(M, complejos):
    """
    Value of each complex in each group.

    LIANA reassembles complexes by the subunit with the LOWEST MEAN ('min'
    rule), and then keeps ALL the statistics of that subunit, not
    the minimum of each statistic separately. That is why it is also necessary
    to return which subunit won, in order to obtain its logFC.

    CAUTION: M comes as GENES x GROUPS. Subunits are rows, not columns;
    the minimum is taken over axis=0. Indexing it the other way round gives no error, it gives NaN, which is
    worse.
    """
    medias, elegida = {}, {}
    perdidos = []
    for c in complejos:
        subs = str(c).split("_")
        if not all(s in M.index for s in subs):
            perdidos.append(c)
            continue
        sub_M = M.loc[subs]                       # subunits x groups
        medias[c] = sub_M.min(axis=0)
        elegida[c] = sub_M.idxmin(axis=0)
    return (pd.DataFrame(medias).T, pd.DataFrame(elegida).T, perdidos)


def _lookup(df_ancho, filas, columnas):
    """Vectorised df_ancho.at[fila, col], with NaN where missing."""
    if not len(df_ancho):
        return np.full(len(filas), np.nan)
    i = df_ancho.index.get_indexer(filas)
    j = df_ancho.columns.get_indexer(columnas)
    V = df_ancho.to_numpy()
    # the same lookup serves for means (float) and for subunit names
    # (str); if I allocate float and then insert strings, it breaks.
    relleno = np.nan if V.dtype.kind in "fc" else None
    out = np.full(len(filas), relleno, dtype=V.dtype if V.dtype.kind in "fc"
                  else object)
    ok = (i >= 0) & (j >= 0)
    out[ok] = V[i[ok], j[ok]]
    return out


def annotate_liana(
    res: pd.DataFrame,
    adata: ad.AnnData,
    groupby: str,
    *,
    use_raw: bool = False,
    layer: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Adds to the LIANA table the expression context it lacks.

    New columns
    -----------
    L_source, R_target   mean of the ligand in the sender and of the receptor in the
                         receiver. They are the two pieces from which
                         `lr_means` and `expr_prod` are derived.
    L_target, R_source   the crossed ones. They show reciprocity and
                         autocrine signalling: if R_source is also high, the sender
                         is listening to itself.
    L_global, R_global   mean over ALL cells. The mental denominator
                         for deciding whether 0.6 is high or low in your data.
    lfc_L, lfc_R         logFC of the ligand in the sender and of the receptor in the
                         receiver, each against the rest of the populations.
                         Their mean is exactly `lr_logfc`.
    spec_weight_oe       `spec_weight` divided by its value under uniformity,
                         1/n_populations². It is an observed/expected ratio, NOT an
                         odds ratio: there are no odds anywhere, it is a
                         ratio of fractions. 1 = like the uniform distribution,
                         3 = three times more concentrated in this pair.

    Verification
    ------------
    It is checked that (L_source+R_target)/2 == lr_means, that
    L_source*R_target == expr_prod and that (lfc_L+lfc_R)/2 == lr_logfc. For the
    logFC BOTH possible definitions are tried (difference of means on
    the log matrix, and log2 of the de-logged means) and the one that
    reproduces the LIANA column is kept. If neither matches, it warns and does not give you lfc_L
    or lfc_R, because a badly reconstructed column is worse than none.

    CAUTION with L_global: it is the mean over all cells, not the sum of the
    per-group means. NATMI uses the latter for its denominator, and the two only
    coincide if all populations have the same size.
    """
    out = res.copy()
    col_l = "ligand_complex" if "ligand_complex" in out else "ligand"
    col_r = "receptor_complex" if "receptor_complex" in out else "receptor"

    if groupby not in adata.obs:
        raise ValueError(f"Missing obs['{groupby}']")

    X, var_names = _matriz_expresion(adata, use_raw=use_raw, layer=layer)
    etiquetas = adata.obs[groupby].astype(str).to_numpy()
    grupos = sorted(pd.unique(etiquetas))
    n_pop = len(grupos)
    n_tot = adata.n_obs

    S, Se, n = _sumas_por_grupo(X, etiquetas, grupos)
    tot = np.sum([S[g] for g in grupos], axis=0)
    tot_e = np.sum([Se[g] for g in grupos], axis=0)

    # per-group means (matrix as is, assumed to be log-normalised)
    M = pd.DataFrame({g: S[g] / max(n[g], 1) for g in grupos}, index=var_names)
    glob = pd.Series(tot / n_tot, index=var_names)

    # two logFC candidates, one per definition
    lfcA, lfcB = {}, {}
    for g in grupos:
        resto = max(n_tot - n[g], 1)
        lfcA[g] = S[g] / max(n[g], 1) - (tot - S[g]) / resto
        lfcB[g] = (np.log2(Se[g] / max(n[g], 1) + 1)
                   - np.log2((tot_e - Se[g]) / resto + 1))
    LFC = {"A_dif_medias_log": pd.DataFrame(lfcA, index=var_names),
           "B_log2_medias_crudas": pd.DataFrame(lfcB, index=var_names)}

    # --- complexes
    comps = pd.unique(np.concatenate([out[col_l].astype(str).to_numpy(),
                                      out[col_r].astype(str).to_numpy()]))
    CM, CE, perdidos = _complejos_por_grupo(M, comps)
    if perdidos and verbose:
        print(f"[contexto] {len(perdidos)} complexes with subunits missing "
              f"from var_names (e.g. {perdidos[:5]}); those rows become NaN")
    CG, _, _ = _complejos_por_grupo(glob.to_frame("_g"), comps)

    lig = out[col_l].astype(str).to_numpy()
    rec = out[col_r].astype(str).to_numpy()
    src = out["source"].astype(str).to_numpy()
    tgt = out["target"].astype(str).to_numpy()

    out["L_source"] = _lookup(CM, lig, src)
    out["R_target"] = _lookup(CM, rec, tgt)
    out["L_target"] = _lookup(CM, lig, tgt)
    out["R_source"] = _lookup(CM, rec, src)
    out["L_global"] = _lookup(CG, lig, np.repeat("_g", len(lig)))
    out["R_global"] = _lookup(CG, rec, np.repeat("_g", len(rec)))

    # --- verification of the means
    chks = {}
    if "lr_means" in out:
        chks["lr_means"] = np.nanmax(np.abs(
            (out["L_source"] + out["R_target"]) / 2 - out["lr_means"]))
    if "expr_prod" in out:
        chks["expr_prod"] = np.nanmax(np.abs(
            out["L_source"] * out["R_target"] - out["expr_prod"]))

    # --- logFC: choose the definition that reproduces lr_logfc
    elegido, dif_lfc = None, {}
    if "lr_logfc" in out:
        sub_l = _lookup(CE, lig, src)
        sub_r = _lookup(CE, rec, tgt)
        for nombre, F in LFC.items():
            a = _lookup(F, sub_l, src)
            b = _lookup(F, sub_r, tgt)
            dif_lfc[nombre] = np.nanmax(np.abs((a + b) / 2 - out["lr_logfc"]))
        elegido = min(dif_lfc, key=dif_lfc.get)
        if dif_lfc[elegido] <= _EPS_CHK:
            out["lfc_L"] = _lookup(LFC[elegido], sub_l, src)
            out["lfc_R"] = _lookup(LFC[elegido], sub_r, tgt)
        else:
            warnings.warn(
                f"Could not reproduce 'lr_logfc' with either of the two "
                f"definitions (best residual {dif_lfc[elegido]:.4g} with "
                f"'{elegido}'). Not adding lfc_L/lfc_R: better not to give you "
                f"columns whose origin is unknown.")
            elegido = None

    # --- NATMI observed/expected
    if "spec_weight" in out:
        out["spec_weight_oe"] = out["spec_weight"] * (n_pop ** 2)

    if verbose:
        print(f"[contexto] {n_pop} populations, {n_tot} cells; "
              f"expected spec_weight under uniformity = "
              f"{1/n_pop**2:.4f}")
        for k, v in chks.items():
            estado = "OK" if v <= _EPS_CHK else "MISMATCH"
            print(f"[contexto] check {k:10s} max residual {v:.3g}  {estado}")
        if dif_lfc:
            for k, v in dif_lfc.items():
                print(f"[contexto] check lr_logfc  {k:22s} max residual "
                      f"{v:.3g}")
            if elegido:
                print(f"[contexto] logFC reconstructed with '{elegido}'")
        if any(v > _EPS_CHK for v in chks.values()):
            print("[contexto] WARNING: if lr_means/expr_prod do not match, the "
                  "matrix you are passing here is NOT the same one with which "
                  "you ran LIANA (use_raw, layer, or the wrong adata?). "
                  "Do not trust L_source/R_target until this is fixed.")
    return out


def orden_columnas_liana(res: pd.DataFrame) -> pd.DataFrame:
    """
    Reorders for reading: identity, decision, specificity, context. The
    raw magnitude scores and the CellPhoneDB p go to the end
    because in a cytokine question they add nothing (the p saturates at 0 and the
    magnitude is dominated by the abundant partner).
    """
    delante = ["source", "target", "ligand_complex", "receptor_complex",
               "ligand", "receptor", "role",
               "specificity_rank", "magnitude_rank",
               "lr_logfc", "lfc_L", "lfc_R",
               "spec_weight_oe", "spec_weight", "scaled_weight",
               "L_source", "R_target", "L_global", "R_global",
               "L_target", "R_source"]
    detras = ["lr_means", "expr_prod", "lrscore", "cellphone_pvals",
              "resource"]
    d = [c for c in delante if c in res.columns]
    t = [c for c in detras if c in res.columns]
    resto = [c for c in res.columns if c not in d and c not in t]
    return res[d + resto + t]


def filter_liana(
    res: pd.DataFrame,
    genes: Sequence[str],
    ortholog_map: Optional[pd.DataFrame] = None,
    *,
    max_specificity: float = 0.05,
    min_expr_prod: float = 0.05,
    max_magnitude: Optional[float] = None,
    como: Optional[str] = None,
    sort_by: str = "lr_logfc",
    max_rank: Optional[float] = None,      # backwards compatibility
    rank_col: str = "specificity_rank",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Interactions involving specific genes, with criteria on TWO axes.

    In LIANA the 'ligand_complex'/'receptor_complex' columns contain the
    subunits joined by underscores (IL6R_IL6ST). They must be SPLIT and
    compared element by element: a `contains("IL6")` would also capture
    IL6ST, which is something else. 'role' is set so that you know whether your gene is the
    ligand or a receptor subunit.

    WHY MAGNITUDE IS NOT FILTERED BY RANK
    -------------------------------------
    `magnitude_rank` is a rank against the ENTIRE table, where you compete with
    collagens, fibronectin and Apoe. A cytokine always loses on that rank, and
    yet it can be perfectly well expressed in absolute terms. That
    is why the magnitude floor is `min_expr_prod`, which is on the scale of the
    data, and `max_magnitude` is switched off. Switch it on only if your question is
    open-ended ("what talks most here?") rather than directed at a few genes.

    `specificity_rank` is filtered by rank: there the rank is the metric.
    But it has a tie floor (look at `res[rank_col].value_counts()`), and
    below it there is no ordering. Hence the default ordering is `lr_logfc`
    descending, which is the only specificity column that still
    discriminates within the tie.

    como
        'both' | 'ligand' | 'receptor_subunit'. Filters by the role of your
        genes. 'receptor_subunit' does NOT mean discardable: if you search for Il6 and
        Il11 on Il6ra_Il6st comes up, that is the gp130 family using the same
        receptor, which may be the answer to your question.
    """
    if max_rank is not None:
        warnings.warn("max_rank is deprecated; interpreting it as "
                      "max_specificity. Use max_specificity/min_expr_prod.")
        max_specificity = max_rank

    mapa = (dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
            if ortholog_map is not None else {})
    buscados = {mapa.get(g, g) for g in genes}

    col_l = "ligand_complex" if "ligand_complex" in res else "ligand"
    col_r = "receptor_complex" if "receptor_complex" in res else "receptor"

    subs_l = res[col_l].astype(str).str.split("_")
    subs_r = res[col_r].astype(str).str.split("_")
    es_lig = subs_l.apply(lambda x: bool(buscados & set(x)))
    es_rec = subs_r.apply(lambda x: bool(buscados & set(x)))

    sel = es_lig | es_rec
    out = res[sel].copy()
    out["role"] = np.where(es_lig[sel],
                           np.where(es_rec[sel], "both", "ligand"),
                           "receptor_subunit")
    antes = len(out)

    pasos = []
    if como is not None:
        out = out[out["role"] == como]
        pasos.append(f"como == '{como}': {len(out)}")
    if rank_col in out and max_specificity is not None:
        out = out[out[rank_col] <= max_specificity]
        pasos.append(f"{rank_col} <= {max_specificity}: {len(out)}")
    elif rank_col not in out:
        warnings.warn(f"No '{rank_col}' column; specificity is not filtered")
    if min_expr_prod is not None and "expr_prod" in out:
        out = out[out["expr_prod"] >= min_expr_prod]
        pasos.append(f"expr_prod >= {min_expr_prod}: {len(out)}")
    if max_magnitude is not None and "magnitude_rank" in out:
        out = out[out["magnitude_rank"] <= max_magnitude]
        pasos.append(f"magnitude_rank <= {max_magnitude}: {len(out)}")

    if sort_by in out:
        asc = sort_by.endswith("_rank") or sort_by.startswith("cellphone")
        out = out.sort_values(sort_by, ascending=asc)
    elif verbose:
        print(f"[filtro] no '{sort_by}' column; sorting by {rank_col}. "
              f"Did you forget to run annotate_liana()?")
        if rank_col in out:
            out = out.sort_values(rank_col)

    out = orden_columnas_liana(out)

    if verbose:
        print(f"[filtro] {antes} annotated interactions with "
              f"{sorted(buscados)}")
        for p in pasos:
            print(f"[filtro]   -> {p}")
        if len(out):
            print(out["role"].value_counts().to_string())
            cols = [c for c in ("source", "target", col_l, col_r, "role",
                                "specificity_rank", "lr_logfc", "lfc_L",
                                "lfc_R", "spec_weight_oe", "L_source",
                                "R_target", "L_global", "R_global")
                    if c in out.columns]
            print(out[cols].head(20).round(3).to_string(index=False))
        else:
            print("[filtro] nothing survives. Before relaxing the threshold check "
                  "whether it is due to expr_prod (the pairs are off) or to "
                  "specificity_rank (they are on everywhere): these are "
                  "different conclusions.")
    return out


def compare_resources(
    adata: ad.AnnData,
    groupby: str,
    resources: Sequence[str] = ("mouseconsensus", "cellphonedb",
                                "cellchatdb", "celltalkdb"),
    max_rank: float = 0.05,
    rank_col: str = "specificity_rank",
    verbose: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """
    The same question with several resources: which interactions survive?

    DISTINGUISH TWO THINGS THAT ARE EASILY CONFUSED, as I confused them myself this morning with
    the GOEA backgrounds:

      n_recursos_testado  in how many resources that pair is ANNOTATED
      n_recursos_sig      in how many it falls below max_rank

    A pair absent from CellTalkDB is not "non-significant", it is "not annotated", and
    if you do not separate the two the smallest resource makes everything look
    fragile. The fraction that matters is sig/tested, not sig/total.

    The human resources (cellphonedb, cellchatdb, celltalkdb) are run as
    they are on mouse symbols and will NOT match well unless LIANA
    translates internally. If you see absurdly low coverage in a human resource, use
    li.rs.translate_resource() to ortholog-transfer it, or stick with
    'mouseconsensus' and those that LIANA serves for mouse.
    """
    trozos = []
    for r in resources:
        if verbose:
            print(f"\n===== resource '{r}' =====")
        try:
            res = run_liana(adata, groupby, resource=r,
                            key_added=f"liana_{r}", verbose=verbose, **kwargs)
        except Exception as e:
            warnings.warn(f"{r} failed: {e}")
            continue
        trozos.append(res)

    if not trozos:
        raise RuntimeError("No resource completed")

    todo = pd.concat(trozos, ignore_index=True)
    col_l = "ligand_complex" if "ligand_complex" in todo else "ligand"
    col_r = "receptor_complex" if "receptor_complex" in todo else "receptor"
    clave = ["source", "target", col_l, col_r]

    todo["_sig"] = todo[rank_col] <= max_rank
    tabla = todo.groupby(clave).agg(
        n_recursos_testado=("resource", "nunique"),
        n_recursos_sig=("_sig", "sum"),
        mejor_rank=(rank_col, "min"),
    ).reset_index()
    tabla["fraction"] = (tabla["n_resources_sig"]
                         / tabla["n_resources_tested"].clip(lower=1))
    tabla = tabla.sort_values(["n_resources_sig", "best_rank"],
                              ascending=[False, True])

    if verbose:
        n = todo["resource"].nunique()
        print(f"\n[recursos] {n} resources run")
        print(f"[recursos] pairs annotated in all {n}: "
              f"{int((tabla.n_recursos_testado == n).sum())}")
        print(f"[recursos] significant in all those in which they are "
              f"annotated: {int((tabla.fraccion == 1).sum())}")
        print("[recursos] ALWAYS compare sig/tested, not sig/total: a pair "
              "that exists in only one resource is not fragile, it is exclusive.")
        print(tabla.head(20).round(4).to_string(index=False))
    return tabla


def plot_liana(adata: ad.AnnData, key: str = "liana_res",
               colour: str = "specificity_rank", size: str = "magnitude_rank",
               source_labels=None, target_labels=None, top_n: int = 25,
               **kwargs):
    """
    LIANA dotplot, with specificity in the COLOUR by default.

    The library default usually highlights magnitude, and magnitude rewards
    genes that are abundant everywhere. If what you are looking for is which pair of cell
    types communicates in a particular way, the colour has to be specificity.
    """
    import liana as li

    return li.pl.dotplot(
        adata=adata, uns_key=key, colour=colour, size=size,
        source_labels=source_labels, target_labels=target_labels,
        top_n=top_n, orderby=colour, orderby_ascending=True, **kwargs)


# ---------------------------------------------------------------------------
# 5c. Delta between conditions
# ---------------------------------------------------------------------------

def _clave_liana(df: pd.DataFrame) -> pd.MultiIndex:
    col_l = "ligand_complex" if "ligand_complex" in df else "ligand"
    col_r = "receptor_complex" if "receptor_complex" in df else "receptor"
    return pd.MultiIndex.from_arrays(
        [df["source"].astype(str), df["target"].astype(str),
         df[col_l].astype(str), df[col_r].astype(str)],
        names=["source", "target", "ligand_complex", "receptor_complex"])


def _medias_y_props(adata, groupby, complejos, use_raw=False, layer=None):
    """Means and fraction of expressing cells, per complex and per group."""
    X, var_names = _matriz_expresion(adata, use_raw=use_raw, layer=layer)
    etiquetas = adata.obs[groupby].astype(str).to_numpy()
    grupos = sorted(pd.unique(etiquetas))

    medias, props = {}, {}
    for g in grupos:
        m = etiquetas == g
        sub = X[m]
        ng = max(int(m.sum()), 1)
        medias[g] = np.asarray(sub.mean(axis=0)).ravel()
        if sp.issparse(sub):
            props[g] = np.asarray((sub > 0).sum(axis=0)).ravel() / ng
        else:
            props[g] = (sub > 0).sum(axis=0) / ng
    M = pd.DataFrame(medias, index=var_names)
    P = pd.DataFrame(props, index=var_names)

    CM, _, _ = _complejos_por_grupo(M, complejos)
    # LIANA's expr_prop is applied to EACH subunit: the worst one decides
    CP, _, _ = _complejos_por_grupo(P, complejos)
    return CM, CP


def delta_liana(
    res_a: pd.DataFrame,
    res_b: pd.DataFrame,
    *,
    seleccion: Optional[pd.DataFrame] = None,
    nombres: Sequence[str] = ("A", "B"),
    cols: Sequence[str] = ("expr_prod", "L_source", "R_target", "lr_logfc"),
    adatas: Optional[Sequence[ad.AnnData]] = None,
    groupby: Optional[str] = None,
    expr_prop: float = 0.1,
    use_raw: bool = False,
    layer: Optional[str] = None,
    min_abs_delta: Optional[float] = None,
    top: int = 25,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compares two LIANA runs (one per condition) over a fixed selection.

    THE SELECTION IS PASSED IN, NOT RECOMPUTED
    ------------------------------------------
    `seleccion` is the table you already filtered on the pooled data (the output of
    `filter_liana`). There are deliberately no selection thresholds here: if there
    were, the criterion would live in two places and the argument that the list
    was fixed BEFORE looking at the delta would no longer be verifiable. If you pass
    None all pairs in common are compared, which is useful for exploration but not
    for drawing conclusions.

    RANKS ARE NOT COMPARED
    ----------------------
    `magnitude_rank` and `specificity_rank` are ranks WITHIN each table. Two
    runs have different reference distributions (it is enough that in
    one more interactions pass expr_prop for all ranks to
    shift), so subtracting them measures the size of the table, not the biology.
    That is why `cols` only accepts absolute columns. If you include a '_rank' it warns.

    ABSENT IS NOT ZERO
    ------------------
    A pair may be missing from a table because the gene is not expressed or because it
    did not reach `expr_prop`, which is a fraction of cells, not a mean. These are
    opposite conclusions and the table separates them in 'reason':

      no_expresado      the mean is 0 in that condition
      bajo_expr_prop    there is expression but in fewer than expr_prop of cells
      desconocido       you did not pass me the adatas and I cannot tell

    Distinguishing them requires `adatas=(adata_a, adata_b)` and `groupby`. Without
    them the delta of rows present in only one condition stays NaN, which
    is the honest choice: you do not know whether the change is from 0.6 to 0 or from 0.6 to 0.55.

    And with one replicate per condition this is DESCRIPTIVE. There is no p-value because
    nothing is replicated: the delta and the between-animal variability are
    indistinguishable.
    """
    na, nb = nombres[0], nombres[1]
    malas = [c for c in cols if c.endswith("_rank")]
    if malas:
        warnings.warn(
            f"{malas} are internal ranks of each table and are not comparable "
            f"between runs. Removing them. Use expr_prod / L_source / R_target / "
            f"lr_logfc / spec_weight_oe.")
        cols = [c for c in cols if c not in malas]

    A = res_a.copy(); A.index = _clave_liana(A)
    B = res_b.copy(); B.index = _clave_liana(B)
    for nom, T in ((na, A), (nb, B)):
        if T.index.duplicated().any():
            raise ValueError(f"{nom} has duplicated pairs; did you concatenate "
                             f"resources without separating them?")

    usables = [c for c in cols if c in A.columns and c in B.columns]
    faltan = [c for c in cols if c not in usables]
    if faltan:
        warnings.warn(f"{faltan} are not in both tables (is annotate_liana "
                      f"missing on one of them?); omitting them")
    if not usables:
        raise ValueError("No comparable column")

    if seleccion is not None:
        clave = _clave_liana(seleccion).unique()
    else:
        clave = A.index.union(B.index)
        if verbose:
            print("[delta] no selection: comparing everything. Exploratory, do not "
                  "draw conclusions from this.")

    out = pd.DataFrame(index=clave)
    for c in usables:
        out[f"{c}_{na}"] = A[c].reindex(clave)
        out[f"{c}_{nb}"] = B[c].reindex(clave)

    en_a = clave.isin(A.index)
    en_b = clave.isin(B.index)
    out["status"] = np.where(en_a & en_b, "both_dbs",
                    np.where(en_a, f"only_{na}",
                    np.where(en_b, f"only_{nb}", "in_neither")))
    out["reason"] = np.where(out["status"] == "both_dbs", "", "unknown")

    # --- fill in the missing side, if given the means to do so
    if adatas is not None and groupby is not None:
        comps = pd.unique(np.concatenate([
            clave.get_level_values("ligand_complex").to_numpy(),
            clave.get_level_values("receptor_complex").to_numpy()]))
        lig = clave.get_level_values("ligand_complex").to_numpy()
        rec = clave.get_level_values("receptor_complex").to_numpy()
        src = clave.get_level_values("source").to_numpy()
        tgt = clave.get_level_values("target").to_numpy()

        motivo = out["reason"].to_numpy(dtype=object)
        for nom, adata, presente in ((na, adatas[0], en_a),
                                     (nb, adatas[1], en_b)):
            if presente.all():
                continue
            CM, CP = _medias_y_props(adata, groupby, comps,
                                     use_raw=use_raw, layer=layer)
            Ls = _lookup(CM, lig, src); Rt = _lookup(CM, rec, tgt)
            pL = _lookup(CP, lig, src); pR = _lookup(CP, rec, tgt)
            falta = ~presente
            for col, val in (("L_source", Ls), ("R_target", Rt)):
                if f"{col}_{nom}" in out:
                    v = out[f"{col}_{nom}"].to_numpy(dtype=float)
                    v[falta] = val[falta]
                    out[f"{col}_{nom}"] = v
            if f"expr_prod_{nom}" in out:
                v = out[f"expr_prod_{nom}"].to_numpy(dtype=float)
                v[falta] = (Ls * Rt)[falta]
                out[f"expr_prod_{nom}"] = v
            peor = np.fmin(np.nan_to_num(pL, nan=0.0),
                           np.nan_to_num(pR, nan=0.0))
            cero = np.fmin(np.nan_to_num(Ls, nan=0.0),
                           np.nan_to_num(Rt, nan=0.0)) <= 0
            motivo[falta & cero] = "not_expressed"
            motivo[falta & ~cero & (peor < expr_prop)] = "low_expr_prop"
            motivo[falta & ~cero & (peor >= expr_prop)] = "present_unannotated"
        out["reason"] = motivo

    for c in usables:
        out[f"d_{c}"] = out[f"{c}_{nb}"] - out[f"{c}_{na}"]

    ref = f"d_{usables[0]}"
    out = out.sort_values(ref, key=lambda s: s.abs(), ascending=False)
    if min_abs_delta is not None:
        out = out[out[ref].abs() >= min_abs_delta]

    orden = (["status", "reason"]
             + [f"d_{c}" for c in usables]
             + [f"{c}_{n}" for c in usables for n in (na, nb)])
    out = out[orden].reset_index()

    if verbose:
        print(f"[delta] {len(clave)} pairs compared "
              f"({nb} minus {na}), sorted by |{ref}|")
        print(out["status"].value_counts().to_string())
        m = out.loc[out["reason"] != "", "reason"].value_counts()
        if len(m):
            print(m.to_string())
            if (out["reason"] == "unknown").any():
                print("[delta] pass adatas=(a, b) and groupby to find out whether "
                      "the absent pairs are at zero or only below "
                      "expr_prop; it is not the same.")
        print(out.head(top).round(3).to_string(index=False))
        print("[delta] remember: one replicate per condition. This describes, "
              "it does not test.")
    return out

"""
Loading of STARsolo matrices (GeneFull + Velocyto) with optional CellBender.

Typical usage, in two passes per sample:

    # pass 1: raw matrix -> CellBender input
    a = load_full_adata(solo, include_cellbender=False, load_velocyto=False,
                        keep_raw_layer=False)
    a.layers.clear()
    a.write_h5ad(f"{solo}/adata_raw.h5ad")

    # ... CellBender is run here on that h5ad ...

    # pass 2: CellBender + Velocyto + annotation
    adata = load_full_adata(
        solo,
        include_cellbender=True,
        cellbender_h5=f"{solo}/CellBender/adata_raw_cellbender.h5",
        load_velocyto=True,
        min_cell_probability=0.5,
        gtf_path=GTF,
        gene2ensembl_path=f"{REFERENCE_DIR}/gene2ensembl.gz",
    )

Resulting layers:
    counts_raw          GeneFull exactly as output by STARsolo
    counts_cellbender   decontaminated GeneFull (only if include_cellbender)
    spliced / unspliced / ambiguous     RAW Velocyto, NOT decontaminated

adata.X = counts_cellbender if CellBender is present, otherwise counts_raw.

The Velocyto matrices are deliberately left raw. CellBender models a
gene x droplet matrix and does not preserve the U/S ratio on which RNA
velocity depends; moreover, the soup does not have the same spliced/unspliced
composition as the cells, so rescaling both layers by the same factor does not
correct the bias. Decontamination is used for cell calling, clustering,
embedding and DE; kinetics are computed on the original counts of those cells.

Notes on the CellBender 0.3.x output, learned the hard way:
  - The "full" h5 does NOT contain all input barcodes, only the analysed
    droplets (probable cells + additional). Its set takes precedence, not
    that of GeneFull.
  - If the input h5ad had the IDs in the index rather than in a column,
    CellBender fills var['gene_id'] with 'NA' for all genes. Therefore the
    gene index is chosen by checking, not by guessing (_pick_var_index).
  - It provides var['ambient_expression'] (soup profile: what it is made of)
    and obs['cell_probability'] / ['background_fraction'], which are kept.

Indexing: var_names = gene_id from features.tsv. With an NCBI RefSeq GTF
that gene_id IS the symbol; the stable identifiers (entrez_id, mgi_id)
are added from the GTF with annotate_var_from_gtf().
"""

from __future__ import annotations

import gzip
import os
import re
import urllib.request
import warnings
from typing import Iterable, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


# --------------------------------------------------------------------------
# low-level utilities
# --------------------------------------------------------------------------

def _resolve(path: str) -> str:
    """Returns path or path.gz, whichever exists."""
    if os.path.exists(path):
        return path
    if os.path.exists(path + ".gz"):
        return path + ".gz"
    raise FileNotFoundError(f"File not found: {path} (nor .gz)")


def _make_unique(values, join: str = "-") -> np.ndarray:
    """Like var_names_make_unique, but on a standalone array."""
    out = np.asarray(values, dtype=object).copy()
    seen: dict = {}
    for i, v in enumerate(out):
        if v in seen:
            seen[v] += 1
            out[i] = f"{v}{join}{seen[v]}"
        else:
            seen[v] = 0
    return out.astype(str)


def _read_barcodes(path: str) -> pd.Index:
    bc = pd.read_csv(_resolve(path), sep="\t", header=None, usecols=[0]).iloc[:, 0]
    idx = pd.Index(bc.astype(str))
    # NOTE: pd.Index(Series) inherits the .name of the Series, which with
    # header=None is the integer 0. anndata rejects an index.name that is not str or None.
    idx.name = None
    return idx


def _read_features(path: str) -> pd.DataFrame:
    """STARsolo features.tsv: gene_id, gene_symbol, feature_type."""
    df = pd.read_csv(_resolve(path), sep="\t", header=None, dtype=str)
    ncol = df.shape[1]
    names = ["gene_id", "gene_symbol", "feature_type"][:ncol]
    df.columns = names
    if "gene_symbol" not in df.columns:
        df["gene_symbol"] = df["gene_id"]
    return df


def read_star_matrix(
    mtx_path: str,
    barcodes_path: str,
    features_path: str,
    dtype=np.float32,
) -> ad.AnnData:
    """
    Reads a STARsolo mtx/barcodes/features trio and returns a cells x genes
    AnnData indexed by gene_id.
    """
    X = scipy.io.mmread(_resolve(mtx_path))          # genes x cells
    X = sp.csr_matrix(X.T, dtype=dtype)              # cells x genes

    obs_names = _read_barcodes(barcodes_path)
    var = _read_features(features_path)

    if X.shape != (len(obs_names), len(var)):
        raise ValueError(
            f"Inconsistent dimensions in {mtx_path}: matrix {X.shape}, "
            f"{len(obs_names)} barcodes, {len(var)} features"
        )

    if var["gene_id"].duplicated().any():
        raise ValueError(f"duplicated gene_id in {features_path}")

    var = var.set_index("gene_id")
    var.index.name = None

    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=obs_names), var=var)
    # the symbol may be repeated: it is kept as a column and a unique
    # version is stored separately, in case it is wanted for plots
    adata.var["gene_symbol_unique"] = _make_unique(adata.var["gene_symbol"].values)
    return adata


def _align_indexer(source: pd.Index, target: pd.Index, what: str, strict: bool):
    """Positional indices of target within source. -1 if missing."""
    if source.has_duplicates:
        raise ValueError(f"{what}: the source index has duplicates")
    pos = source.get_indexer(target)
    missing = int((pos < 0).sum())
    if missing:
        msg = f"{what}: {missing}/{len(target)} elements missing from the source matrix"
        if strict:
            raise ValueError(msg)
        warnings.warn(msg + " -> filled with zeros")
    return pos


def _subset_positional(X: sp.spmatrix, row_pos: np.ndarray, col_pos: np.ndarray) -> sp.csr_matrix:
    """Positional subsetting tolerant to -1 (fills with zeros)."""
    rows_ok = row_pos >= 0
    cols_ok = col_pos >= 0

    sub = X[np.where(rows_ok, row_pos, 0)][:, np.where(cols_ok, col_pos, 0)]
    sub = sp.csr_matrix(sub)

    if not rows_ok.all():
        sub = sp.diags(rows_ok.astype(sub.dtype)) @ sub
    if not cols_ok.all():
        sub = sub @ sp.diags(cols_ok.astype(sub.dtype))
    return sp.csr_matrix(sub)


# --------------------------------------------------------------------------
# CellBender
# --------------------------------------------------------------------------

def read_cellbender_h5(path: str) -> ad.AnnData:
    """
    Reads the CellBender 0.3.x output. Uses the official reader if
    available; otherwise falls back to scanpy.
    """
    try:
        from cellbender.remove_background.downstream import anndata_from_h5
        return anndata_from_h5(path)
    except Exception:
        import scanpy as sc
        return sc.read_10x_h5(path)


def inspect_cellbender_h5(path: str, n: int = 5) -> None:
    """
    Diagnostic: shows what is actually in the CellBender output.
    Use it when the alignment of var or obs fails.
    """
    a = read_cellbender_h5(path)
    print(f"{path}\n  shape: {a.shape}")
    print(f"  var_names[:{n}]: {list(a.var_names[:n])}  (unique: {not a.var_names.has_duplicates})")
    print(f"  obs_names[:{n}]: {list(a.obs_names[:n])}  (unique: {not a.obs_names.has_duplicates})")
    print("  var columns:")
    for c in a.var.columns:
        vals = a.var[c].astype(str)
        print(f"    {c:<24} unique={vals.nunique()}/{len(vals)}  e.g.={list(vals[:3])}")
    print("  obs columns:", list(a.obs.columns))


def _pick_var_index(
    adata_cb: ad.AnnData,
    target: pd.Index,
    forced_key: Optional[str] = None,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Chooses which column (or the index itself) of the CellBender output to use
    for matching the GeneFull genes.

    It does not guess: it requires the candidate to be unique and keeps the one
    that recovers the most genes of 'target'. If none is suitable, it prints
    the full table so that the reason is visible.
    """
    tgt = pd.Index(target.astype(str))

    def _cands():
        if forced_key is not None:
            if forced_key == "__index__":
                yield "__index__", pd.Index(adata_cb.var_names.astype(str))
            elif forced_key in adata_cb.var.columns:
                yield forced_key, pd.Index(adata_cb.var[forced_key].astype(str))
            else:
                raise ValueError(
                    f"'{forced_key}' is not in var. Columns: "
                    f"{list(adata_cb.var.columns)}"
                )
            return
        yield "__index__", pd.Index(adata_cb.var_names.astype(str))
        for c in adata_cb.var.columns:
            yield c, pd.Index(adata_cb.var[c].astype(str))

    report, best = [], None
    for name, idx in _cands():
        unique = not idx.has_duplicates
        overlap = int(tgt.isin(idx).sum())
        report.append((name, unique, overlap))
        if unique and overlap and (best is None or overlap > best[1]):
            best = (name, overlap, idx)

    if best is None:
        lines = "\n".join(
            f"    {n:<24} unique={u!s:<5} recovers={o}/{len(tgt)}"
            for n, u, o in report
        )
        raise ValueError(
            "No column of the CellBender output can be used to match the "
            f"GeneFull genes:\n{lines}\n"
            "  Pass cellbender_var_key='<column>' if you know which one it is, or "
            "inspect the file with inspect_cellbender_h5()."
        )

    name, overlap, idx = best
    if verbose:
        print(f"[cellbender] genes matched by '{name}': "
              f"{overlap}/{len(tgt)} GeneFull genes")

    adata_cb = adata_cb.copy()
    if name != "__index__":
        adata_cb.var["cellbender_var_names"] = adata_cb.var_names.astype(str)
        adata_cb.var_names = idx
    return adata_cb


# --------------------------------------------------------------------------
# ambient: diagnostics
# --------------------------------------------------------------------------

def ambient_fraction_per_gene(adata: ad.AnnData) -> pd.Series:
    """
    Fraction of counts removed by CellBender per gene, over the cells
    present in the object. Useful for excluding soup-dominated genes from
    the velocity gene set.
    """
    raw = np.asarray(adata.layers["counts_raw"].sum(axis=0)).ravel()
    cb = np.asarray(adata.layers["counts_cellbender"].sum(axis=0)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(raw > 0, 1.0 - cb / raw, 0.0)
    return pd.Series(np.clip(frac, 0, 1), index=adata.var_names, name="ambient_fraction")



# --------------------------------------------------------------------------
# annotation from the GTF (traceability)
# --------------------------------------------------------------------------

_GTF_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')

GENE2ENSEMBL_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"


def _open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def read_gtf_header(gtf_path: str, max_lines: int = 50) -> dict:
    """
    '#!' header of the GTF. It is the provenance of the annotation: keep it.
    """
    out = {}
    with _open_text(gtf_path) as fh:
        for i, line in enumerate(fh):
            if i >= max_lines or not line.startswith("#"):
                break
            line = line.lstrip("#!").strip()
            if not line:
                continue
            key, _, val = line.partition(" ")
            out[key.strip()] = val.strip()
    return out


def parse_gtf_genes(gtf_path: str, feature: str = "gene") -> pd.DataFrame:
    """
    Extracts one row per gene from the GTF, with the db_xref entries split out.

    Works with NCBI RefSeq GTFs (gene_id = symbol, db_xref with GeneID/MGI)
    and also with Ensembl/GENCODE (gene_id = ENSMUSG..., gene_name = symbol).
    """
    rows = []
    with _open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != feature:
                continue

            attrs, xrefs = {}, []
            for k, v in _GTF_ATTR.findall(f[8]):
                if k == "db_xref":
                    xrefs.append(v)
                elif k not in attrs:          # keep the 1st occurrence
                    attrs[k] = v

            gid = attrs.get("gene_id")
            if not gid:
                continue

            rec = {
                "gene_id": gid,
                # RefSeq uses 'gene', Ensembl/GENCODE use 'gene_name'
                "gene_symbol_gtf": attrs.get("gene") or attrs.get("gene_name", ""),
                "gene_biotype": attrs.get("gene_biotype") or attrs.get("gene_type", ""),
                "description": attrs.get("description", ""),
                "seqname": f[0],
                "start": int(f[3]),
                "end": int(f[4]),
                "strand": f[6],
                "entrez_id": "",
                "mgi_id": "",
                "ensembl_id": "",
            }
            for x in xrefs:
                if x.startswith("GeneID:"):
                    rec["entrez_id"] = x.split(":", 1)[1]
                elif x.startswith("MGI:"):
                    rec["mgi_id"] = x.split(":", 1)[1]      # -> 'MGI:5455983'
                elif x.startswith("Ensembl:"):
                    rec["ensembl_id"] = x.split(":", 1)[1]  # some GTFs do include it

            # Ensembl GTFs: the gene_id itself is already the ENSMUSG
            if not rec["ensembl_id"] and gid.startswith(("ENSMUSG", "ENSG")):
                rec["ensembl_id"] = gid.split(".")[0]

            rows.append(rec)

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"No '{feature}' features found in {gtf_path}")

    dup = df["gene_id"].duplicated()
    if dup.any():
        warnings.warn(
            f"{int(dup.sum())} duplicated gene_id in the GTF (genes on several "
            f"scaffolds); keeping the first occurrence"
        )
        df = df[~dup]
    return df


def annotate_var_from_gtf(
    adata: ad.AnnData,
    gtf_path: str,
    columns: Optional[Sequence[str]] = None,
    overwrite: bool = False,
) -> ad.AnnData:
    """
    Adds to adata.var the identifiers from the GTF used in STAR, matching
    on adata.var_names <-> gene_id. Stores the provenance in adata.uns.

    With a RefSeq GTF this yields 'entrez_id' and 'mgi_id', which are the
    stable IDs. 'ensembl_id' will be empty: RefSeq does not include it (see
    add_ensembl_ids()).
    """
    genes = parse_gtf_genes(gtf_path).set_index("gene_id")
    if columns is not None:
        genes = genes[list(columns)]

    hit = adata.var_names.isin(genes.index)
    if not hit.all():
        warnings.warn(
            f"{int((~hit).sum())}/{adata.n_vars} genes of the matrix are not in "
            f"the GTF. Is it the same file that was used for the STAR index?"
        )

    ann = genes.reindex(adata.var_names)
    for col in ann.columns:
        if col in adata.var.columns and not overwrite:
            continue
        vals = ann[col].values
        adata.var[col] = pd.Categorical(vals) if vals.dtype == object else vals

    adata.uns["annotation"] = {
        "gtf_path": str(gtf_path),
        "gtf_file": os.path.basename(str(gtf_path)),
        **read_gtf_header(gtf_path),
    }
    return adata


def download_gene2ensembl(
    dest: str,
    url: str = GENE2ENSEMBL_URL,
    overwrite: bool = False,
    timeout: int = 120,
) -> str:
    """
    Downloads gene2ensembl.gz from NCBI (~200 MB) if it is not already on disk.

    The download is atomic: it writes to '<dest>.part', checks that the gzip
    opens and only then renames. In this way an interrupted download never
    leaves a truncated file that later fails cryptically midway through
    parsing.
    """
    if os.path.exists(dest) and not overwrite:
        return dest

    parent = os.path.dirname(os.path.abspath(dest))
    os.makedirs(parent, exist_ok=True)
    tmp = dest + ".part"

    print(f"[gene2ensembl] downloading {url}")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp, open(tmp, "wb") as fh:
            total = int(resp.headers.get("Content-Length") or 0)
            done = next_mark = 0
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
                done += len(chunk)
                if total and done >= next_mark:
                    print(f"[gene2ensembl] {100 * done / total:5.1f}% "
                          f"({done / 1e6:.0f}/{total / 1e6:.0f} MB)")
                    next_mark += total // 10

        # integrity: if the .gz is truncated, this fails here rather than inside
        # the pandas loop half an hour later
        with gzip.open(tmp, "rt") as fh:
            fh.readline()

        os.replace(tmp, dest)
        print(f"[gene2ensembl] saved to {dest} ({os.path.getsize(dest) / 1e6:.0f} MB)")
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise

    return dest


def _read_gene2ensembl_header(path: str) -> list:
    """
    Reads the header and strips the leading '#'. NCBI writes it as
    '#tax_id', but some mirrors and versions do not, so it is not assumed.
    """
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    return [h.lstrip("#").strip() for h in header]


def add_ensembl_ids(
    adata: ad.AnnData,
    gene2ensembl_path: str,
    tax_id: int = 10090,          # 10090 = Mus musculus, 9606 = Homo sapiens
    entrez_col: str = "entrez_id",
    chunksize: int = 2_000_000,
    download: bool = True,
    url: str = GENE2ENSEMBL_URL,
) -> ad.AnnData:
    """
    Maps Entrez GeneID -> Ensembl gene ID using the official NCBI file
    (gene2ensembl.gz from https://ftp.ncbi.nlm.nih.gov/gene/DATA/).

    This is the reproducible route: a versioned, offline file rather than an API.
    The mapping is NOT complete nor always 1:1 -- both are recorded in var:
      ensembl_id     first Ensembl ID (or empty if none)
      ensembl_all    all of them, separated by ';'
      ensembl_n      how many (0 = no match)
    """
    if entrez_col not in adata.var.columns:
        raise ValueError(
            f"'{entrez_col}' not found in var; run annotate_var_from_gtf() first"
        )

    if not os.path.exists(gene2ensembl_path):
        if not download:
            raise FileNotFoundError(
                f"File not found: {gene2ensembl_path}, and download=False. "
                f"Download it from {url}"
            )
        download_gene2ensembl(gene2ensembl_path, url=url)

    cols = _read_gene2ensembl_header(gene2ensembl_path)
    tax_col = next((c for c in cols if c.endswith("tax_id")), None)
    for needed in ("GeneID", "Ensembl_gene_identifier"):
        if needed not in cols:
            raise ValueError(
                f"{gene2ensembl_path} does not have the column '{needed}'. "
                f"Columns found: {cols}"
            )
    if tax_col is None:
        raise ValueError(f"Cannot find the tax_id column. Columns: {cols}")

    keep = []
    for chunk in pd.read_csv(
        gene2ensembl_path, sep="\t", chunksize=chunksize, dtype=str,
        header=0, names=cols,
        usecols=[tax_col, "GeneID", "Ensembl_gene_identifier"],
    ):
        chunk = chunk[chunk[tax_col] == str(tax_id)]
        if len(chunk):
            keep.append(chunk[["GeneID", "Ensembl_gene_identifier"]])

    if not keep:
        raise ValueError(f"No rows for tax_id={tax_id} in {gene2ensembl_path}")

    m = pd.concat(keep).drop_duplicates()
    m = m[m["Ensembl_gene_identifier"] != "-"]
    grouped = m.groupby("GeneID")["Ensembl_gene_identifier"].agg(list)

    entrez = adata.var[entrez_col].astype(str)
    hits = entrez.map(grouped)

    adata.var["ensembl_all"] = pd.Categorical(
        [";".join(v) if isinstance(v, list) else "" for v in hits]
    )
    adata.var["ensembl_n"] = np.array(
        [len(v) if isinstance(v, list) else 0 for v in hits], dtype=np.int16
    )
    adata.var["ensembl_id"] = pd.Categorical(
        [v[0] if isinstance(v, list) else "" for v in hits]
    )

    n_map = int((adata.var["ensembl_n"] > 0).sum())
    n_multi = int((adata.var["ensembl_n"] > 1).sum())
    print(
        f"[ensembl] {n_map}/{adata.n_vars} genes mapped "
        f"({100 * n_map / adata.n_vars:.1f}%), {n_multi} with more than one Ensembl ID"
    )
    adata.uns.setdefault("annotation", {})["gene2ensembl"] = os.path.basename(
        str(gene2ensembl_path)
    )
    return adata




# --------------------------------------------------------------------------
# main loading
# --------------------------------------------------------------------------

def load_full_adata(
    STARsolo_path: str,
    raw_counts: bool = True,
    include_cellbender: bool = True,
    load_velocyto: bool = True,
    gene_dir: str = "GeneFull",
    matrix_file: str = "UniqueAndMult-EM.mtx",
    cellbender_h5: Optional[str] = None,
    cellbender_var_key: Optional[str] = None,
    keep_raw_layer: bool = True,
    min_cell_probability: Optional[float] = None,
    restrict_to: Optional[Sequence[str]] = None,
    strict: bool = True,
    check_velocyto_sum: bool = False,
    gtf_path: Optional[str] = None,
    gene2ensembl_path: Optional[str] = None,
    cellbender_meta: Optional[dict] = None,
) -> ad.AnnData:
    """
    Parameters
    ----------
    raw_counts
        True -> 'raw' subdirectory; False -> 'filtered'.
    include_cellbender
        Second pass: uses the CellBender output as X.
    keep_raw_layer
        Also stores pre-CellBender GeneFull in layers['counts_raw'].
        Costs memory but makes it possible to audit what was removed.
    min_cell_probability
        If given and CellBender provides 'cell_probability', filters droplets
        before loading Velocyto. Recommended in raw mode: reduces ~700k
        droplets to a few thousand and saves a great deal of memory.
    restrict_to
        Explicit list of barcodes to keep (alternative to the previous filter).
    gtf_path
        GTF used in STAR. If given, adds to var the stable identifiers
        (entrez_id, mgi_id), biotype and coordinates, and stores the
        provenance of the annotation in uns['annotation'].
    gene2ensembl_path
        NCBI gene2ensembl.gz for mapping entrez_id -> Ensembl. Requires
        gtf_path. Downloaded automatically if it does not exist.
    cellbender_meta
        Free-form dict (version, arguments) stored in uns['cellbender'].
    strict
        If True, fails when a barcode/gene is absent from any matrix.
    check_velocyto_sum
        Checks whether GeneFull ~= spliced + unspliced + ambiguous and stores
        the result in adata.uns['velocyto_sum_check'].
    """
    dirsub = "raw" if raw_counts else "filtered"
    gene_path = f"{STARsolo_path}/{gene_dir}/{dirsub}"
    velo_path = f"{STARsolo_path}/Velocyto/{dirsub}"

    # ---- 1. raw GeneFull ---------------------------------------------------
    adata = read_star_matrix(
        f"{gene_path}/{matrix_file}",
        f"{gene_path}/barcodes.tsv",
        f"{gene_path}/features.tsv",
    )
    if keep_raw_layer or not include_cellbender:
        adata.layers["counts_raw"] = adata.X.copy()

    # ---- 2. CellBender -----------------------------------------------------
    if include_cellbender:
        cb_file = cellbender_h5 or f"{STARsolo_path}/CellBender/adata_raw_cellbender.h5"
        adata_cb = read_cellbender_h5(cb_file)
        adata_cb = _pick_var_index(adata_cb, adata.var_names, cellbender_var_key)

        # NOTE: the "full" CellBender h5 does NOT contain all input barcodes,
        # only the droplets it analysed (probable cells + additional, i.e.
        # --total-droplets-included). The rest are not even evaluated.
        # Hence the CellBender droplets take precedence here, not GeneFull's.
        keep = adata.obs_names.isin(adata_cb.obs_names)
        n_keep = int(keep.sum())
        if n_keep == 0:
            raise ValueError(
                "No GeneFull barcode appears in the CellBender output. "
                "This is usually a different suffix ('-1'), not a data problem.\n"
                f"  GeneFull:   {list(adata.obs_names[:3])}\n"
                f"  CellBender: {list(adata_cb.obs_names[:3])}"
            )
        print(f"[cellbender] {n_keep} droplets analysed out of {adata.n_obs} in GeneFull")
        adata = adata[keep].copy()

        # from here on the obs alignment must be exact
        row_pos = _align_indexer(adata_cb.obs_names, adata.obs_names,
                                 "CellBender obs", strict=True)
        # genes, by contrast, may be missing: CellBender excludes from the
        # analysis those it estimates to have ~0 background and may not return them
        col_pos = _align_indexer(adata_cb.var_names, adata.var_names,
                                 "CellBender var", strict=False)

        cb_X = _subset_positional(
            sp.csr_matrix(adata_cb.X), row_pos, col_pos
        ).astype(adata.X.dtype)

        # Genes absent from the output: filled with the RAW counts, not with
        # zeros. CellBender excluded them for having negligible background,
        # so their decontaminated value is the original one. Setting them to
        # zero would silently delete genes.
        missing_genes = col_pos < 0
        if missing_genes.any():
            warnings.warn(
                f"{int(missing_genes.sum())} genes are not in the CellBender "
                f"output; their raw counts are kept without decontamination"
            )
            # adata.X is still the raw matrix at this point
            cb_X = cb_X + adata.X @ sp.diags(missing_genes.astype(adata.X.dtype))
            adata.var["cellbender_returned"] = ~missing_genes

        adata.layers["counts_cellbender"] = sp.csr_matrix(cb_X)
        adata.X = adata.layers["counts_cellbender"].copy()

        # PER-GENE metrics provided by CellBender:
        #   ambient_expression  soup profile (fraction of the ambient counts
        #                       attributable to each gene; sums to ~1 over genes).
        #                       It tells WHAT the soup is made of.
        #   cellbender_analyzed whether the gene entered the model or was
        #                       excluded for negligible background.
        # Do not confuse ambient_expression with ambient_fraction_per_gene(),
        # which measures WHAT FRACTION OF THAT GENE was removed. They differ:
        # a highly expressed gene can dominate the soup and still lose only 2%.
        cb_var = adata_cb.var.reindex(adata.var_names)
        for col in ("ambient_expression", "cellbender_analyzed"):
            if col in cb_var.columns:
                adata.var[col] = cb_var[col].values

        # per-droplet metrics that CellBender does provide and that are worth keeping
        for col in ("cell_probability", "background_fraction", "cell_size",
                    "droplet_efficiency"):
            if col in adata_cb.obs.columns:
                vals = adata_cb.obs[col].to_numpy()
                out = np.full(adata.n_obs, np.nan, dtype=np.float32)
                ok = row_pos >= 0
                out[ok] = vals[row_pos[ok]]
                adata.obs[col] = out

        if cellbender_meta is not None:
            adata.uns["cellbender"] = dict(cellbender_meta)
        adata.uns.setdefault("cellbender", {})["h5"] = os.path.basename(cb_file)

        del adata_cb

    # ---- 3. droplet filtering BEFORE Velocyto -----------------------------
    if restrict_to is not None:
        keep = adata.obs_names.isin(pd.Index(restrict_to))
        adata = adata[keep].copy()
    elif min_cell_probability is not None:
        if "cell_probability" not in adata.obs:
            raise ValueError("'cell_probability' not found; was CellBender run?")
        adata = adata[adata.obs["cell_probability"] >= min_cell_probability].copy()

    # per-gene ambient: ONLY meaningful on the remaining cells.
    # Computed before filtering it would include empty droplets, where the
    # removed fraction is ~100%, and would inflate the value for all genes.
    if "counts_raw" in adata.layers and "counts_cellbender" in adata.layers:
        adata.var["ambient_fraction"] = ambient_fraction_per_gene(adata).values

        adata.obs["counts_cellbender"] = adata.layers["counts_cellbender"].sum(1).A1
        adata.obs["counts_raw"] = adata.layers["counts_raw"].sum(1).A1

        adata.var["counts_cellbender"] = adata.layers["counts_cellbender"].sum(0).A1
        adata.var["counts_raw"] = adata.layers["counts_raw"].sum(0).A1

    # ---- 4. Velocyto -------------------------------------------------------
    if load_velocyto:
        velo_bc = _read_barcodes(f"{velo_path}/barcodes.tsv")
        velo_var = _read_features(f"{velo_path}/features.tsv")
        velo_genes = pd.Index(velo_var["gene_id"].astype(str))

        row_pos = _align_indexer(velo_bc, adata.obs_names, "Velocyto obs", strict)
        col_pos = _align_indexer(velo_genes, adata.var_names, "Velocyto var", strict)

        for layer in ("spliced", "unspliced", "ambiguous"):
            X = scipy.io.mmread(_resolve(f"{velo_path}/{layer}.mtx"))
            X = sp.csr_matrix(X.T, dtype=adata.X.dtype)
            if X.shape != (len(velo_bc), len(velo_genes)):
                raise ValueError(f"Inconsistent dimensions in {layer}.mtx")
            adata.layers[layer] = _subset_positional(X, row_pos, col_pos)
            del X

        if check_velocyto_sum and "counts_raw" in adata.layers:
            tot = (adata.layers["spliced"]
                   + adata.layers["unspliced"]
                   + adata.layers["ambiguous"])
            diff = abs(tot - adata.layers["counts_raw"]).sum()
            denom = adata.layers["counts_raw"].sum()
            frac = float(diff / denom) if denom else np.nan
            adata.uns["velocyto_sum_check"] = frac
            print(f"[check] |GeneFull - (S+U+A)| / GeneFull = {frac:.4f}")

    # ---- 5. optional annotation -----------------------------------------
    if gtf_path is not None:
        annotate_var_from_gtf(adata, gtf_path)
        if gene2ensembl_path is not None:
            add_ensembl_ids(adata, gene2ensembl_path)
    elif gene2ensembl_path is not None:
        raise ValueError(
            "gene2ensembl_path requires gtf_path: Ensembl IDs are mapped "
            "from entrez_id, and entrez_id comes from the GTF."
        )

    adata.obs_names_make_unique()
    return adata





def _safe_sheet_name(name: str, used: set, max_len: int = 31) -> str:
    """Valid Excel sheet names + avoids duplicates."""
    s = re.sub(r'[:\\/?*\[\]]', '_', str(name))[:max_len] or "Sheet"
    base = s
    i = 1
    while s in used:
        suf = f"_{i}"
        s = base[:max_len - len(suf)] + suf
        i += 1
    used.add(s)
    return s

def save_deg_to_excel_simple(
    df: pd.DataFrame,
    out_path: str,
    group_col: str = 'population',
    cols: tuple[str, ...] = ('gene','lfc','pvals_adj','pvals_adjxlfc')
) -> str:
    """
    Saves an Excel file with one sheet per population.
    - Preserves the order in which populations appear in the df.
    - Exports only the columns listed in `cols`, if they exist.
    """
    if group_col not in df.columns:
        raise ValueError(f"'{group_col}' is not in the DataFrame")

    groups_in_order = df[group_col].dropna().drop_duplicates().tolist()
    used = set()

    with pd.ExcelWriter(out_path) as writer:
        for g in groups_in_order:
            gdf = df[df[group_col] == g]
            sheet = _safe_sheet_name(g, used)
            use_cols = [c for c in cols if c in gdf.columns] or list(gdf.columns)
            gdf[use_cols].to_excel(writer, sheet_name=sheet, index=False)


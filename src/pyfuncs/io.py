"""
Carga de matrices STARsolo (GeneFull + Velocyto) con CellBender opcional.

Uso típico, en dos pases por muestra:

    # pase 1: matriz cruda -> input de CellBender
    a = load_full_adata(solo, include_cellbender=False, load_velocyto=False,
                        keep_raw_layer=False)
    a.layers.clear()
    a.write_h5ad(f"{solo}/adata_raw.h5ad")

    # ... aquí corre CellBender sobre ese h5ad ...

    # pase 2: CellBender + Velocyto + anotación
    adata = load_full_adata(
        solo,
        include_cellbender=True,
        cellbender_h5=f"{solo}/CellBender/adata_raw_cellbender.h5",
        load_velocyto=True,
        min_cell_probability=0.5,
        gtf_path=GTF,
        gene2ensembl_path=f"{REFERENCE_DIR}/gene2ensembl.gz",
    )

Layers resultantes:
    counts_raw          GeneFull tal cual sale de STARsolo
    counts_cellbender   GeneFull descontaminado (solo si include_cellbender)
    spliced / unspliced / ambiguous     Velocyto RAW, SIN descontaminar

adata.X = counts_cellbender si hay CellBender, si no counts_raw.

Las matrices de Velocyto se dejan crudas a propósito. CellBender modela una
matriz gene x droplet y no preserva la relación U/S de la que depende RNA
velocity; además el soup no tiene la misma composición spliced/unspliced que
las células, así que reescalar ambas capas por el mismo factor no corrige el
sesgo. La descontaminación se usa para llamar células, clustering, embedding
y DE; la cinética se calcula sobre los counts originales de esas células.

Notas sobre la salida de CellBender 0.3.x, aprendidas a base de tropezar:
  - El h5 "completo" NO trae todos los barcodes de entrada, solo los droplets
    analizados (probable cells + additional). Manda su conjunto, no el de
    GeneFull.
  - Si el h5ad de entrada tenía los IDs en el índice y no en una columna,
    CellBender rellena var['gene_id'] con 'NA' para todos los genes. Por eso
    el índice de genes se elige comprobando, no adivinando (_pick_var_index).
  - Aporta var['ambient_expression'] (perfil del soup: de qué está hecho) y
    obs['cell_probability'] / ['background_fraction'], que se conservan.

Indexación: var_names = gene_id de features.tsv. Con un GTF de NCBI RefSeq
ese gene_id ES el símbolo; los identificadores estables (entrez_id, mgi_id)
se añaden desde el GTF con annotate_var_from_gtf().
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
# utilidades de bajo nivel
# --------------------------------------------------------------------------

def _resolve(path: str) -> str:
    """Devuelve path o path.gz, lo que exista."""
    if os.path.exists(path):
        return path
    if os.path.exists(path + ".gz"):
        return path + ".gz"
    raise FileNotFoundError(f"No existe {path} (ni .gz)")


def _make_unique(values, join: str = "-") -> np.ndarray:
    """Como var_names_make_unique, pero sobre un array suelto."""
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
    # OJO: pd.Index(Series) hereda el .name de la Series, que con header=None
    # es el entero 0. anndata rechaza un index.name que no sea str o None.
    idx.name = None
    return idx


def _read_features(path: str) -> pd.DataFrame:
    """features.tsv de STARsolo: gene_id, gene_symbol, feature_type."""
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
    Lee un trío mtx/barcodes/features de STARsolo y devuelve un AnnData
    células x genes indexado por gene_id.
    """
    X = scipy.io.mmread(_resolve(mtx_path))          # genes x células
    X = sp.csr_matrix(X.T, dtype=dtype)              # células x genes

    obs_names = _read_barcodes(barcodes_path)
    var = _read_features(features_path)

    if X.shape != (len(obs_names), len(var)):
        raise ValueError(
            f"Dimensiones incoherentes en {mtx_path}: matriz {X.shape}, "
            f"{len(obs_names)} barcodes, {len(var)} features"
        )

    if var["gene_id"].duplicated().any():
        raise ValueError(f"gene_id duplicados en {features_path}")

    var = var.set_index("gene_id")
    var.index.name = None

    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=obs_names), var=var)
    # el símbolo puede repetirse: lo dejamos como columna y guardamos aparte
    # una versión unique, por si quieres usarla para plots
    adata.var["gene_symbol_unique"] = _make_unique(adata.var["gene_symbol"].values)
    return adata


def _align_indexer(source: pd.Index, target: pd.Index, what: str, strict: bool):
    """Índices posicionales de target dentro de source. -1 si falta."""
    if source.has_duplicates:
        raise ValueError(f"{what}: el índice de origen tiene duplicados")
    pos = source.get_indexer(target)
    missing = int((pos < 0).sum())
    if missing:
        msg = f"{what}: {missing}/{len(target)} elementos ausentes en la matriz de origen"
        if strict:
            raise ValueError(msg)
        warnings.warn(msg + " -> se rellenan con ceros")
    return pos


def _subset_positional(X: sp.spmatrix, row_pos: np.ndarray, col_pos: np.ndarray) -> sp.csr_matrix:
    """Subsetting posicional tolerante a -1 (rellena con ceros)."""
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
    Lee la salida de CellBender 0.3.x. Usa el lector oficial si está
    disponible; si no, cae a scanpy.
    """
    try:
        from cellbender.remove_background.downstream import anndata_from_h5
        return anndata_from_h5(path)
    except Exception:
        import scanpy as sc
        return sc.read_10x_h5(path)


def inspect_cellbender_h5(path: str, n: int = 5) -> None:
    """
    Diagnóstico: enseña qué hay realmente en la salida de CellBender.
    Úsalo cuando el alineamiento de var u obs falle.
    """
    a = read_cellbender_h5(path)
    print(f"{path}\n  shape: {a.shape}")
    print(f"  var_names[:{n}]: {list(a.var_names[:n])}  (únicos: {not a.var_names.has_duplicates})")
    print(f"  obs_names[:{n}]: {list(a.obs_names[:n])}  (únicos: {not a.obs_names.has_duplicates})")
    print("  var columns:")
    for c in a.var.columns:
        vals = a.var[c].astype(str)
        print(f"    {c:<24} únicos={vals.nunique()}/{len(vals)}  ej={list(vals[:3])}")
    print("  obs columns:", list(a.obs.columns))


def _pick_var_index(
    adata_cb: ad.AnnData,
    target: pd.Index,
    forced_key: Optional[str] = None,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Elige qué columna (o el propio índice) de la salida de CellBender usar para
    casar con los genes de GeneFull.

    No adivina: exige que la candidata sea única y se queda con la que más
    genes de 'target' recupera. Si ninguna sirve, imprime la tabla completa
    para que se vea por qué.
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
                    f"'{forced_key}' no está en var. Columnas: "
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
            f"    {n:<24} único={u!s:<5} recupera={o}/{len(tgt)}"
            for n, u, o in report
        )
        raise ValueError(
            "Ninguna columna de la salida de CellBender sirve para casar los "
            f"genes de GeneFull:\n{lines}\n"
            "  Pásale cellbender_var_key='<columna>' si sabes cuál es, o mira "
            "el fichero con inspect_cellbender_h5()."
        )

    name, overlap, idx = best
    if verbose:
        print(f"[cellbender] genes casados por '{name}': "
              f"{overlap}/{len(tgt)} genes de GeneFull")

    adata_cb = adata_cb.copy()
    if name != "__index__":
        adata_cb.var["cellbender_var_names"] = adata_cb.var_names.astype(str)
        adata_cb.var_names = idx
    return adata_cb


# --------------------------------------------------------------------------
# ambient: diagnóstico
# --------------------------------------------------------------------------

def ambient_fraction_per_gene(adata: ad.AnnData) -> pd.Series:
    """
    Fracción de counts que CellBender ha eliminado por gen, sobre las células
    presentes en el objeto. Útil para excluir genes dominados por soup del
    conjunto de genes de velocity.
    """
    raw = np.asarray(adata.layers["counts_raw"].sum(axis=0)).ravel()
    cb = np.asarray(adata.layers["counts_cellbender"].sum(axis=0)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(raw > 0, 1.0 - cb / raw, 0.0)
    return pd.Series(np.clip(frac, 0, 1), index=adata.var_names, name="ambient_fraction")



# --------------------------------------------------------------------------
# anotación desde el GTF (trazabilidad)
# --------------------------------------------------------------------------

_GTF_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')

GENE2ENSEMBL_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"


def _open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def read_gtf_header(gtf_path: str, max_lines: int = 50) -> dict:
    """
    Cabecera '#!' del GTF. Es la procedencia de tu anotación: guárdala.
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
    Extrae una fila por gen del GTF, con los db_xref desglosados.

    Funciona con GTFs de NCBI RefSeq (gene_id = símbolo, db_xref con GeneID/MGI)
    y también con Ensembl/GENCODE (gene_id = ENSMUSG..., gene_name = símbolo).
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
                elif k not in attrs:          # nos quedamos con la 1a aparición
                    attrs[k] = v

            gid = attrs.get("gene_id")
            if not gid:
                continue

            rec = {
                "gene_id": gid,
                # RefSeq usa 'gene', Ensembl/GENCODE usan 'gene_name'
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
                    rec["ensembl_id"] = x.split(":", 1)[1]  # algunos GTFs sí lo traen

            # GTFs de Ensembl: el propio gene_id ya es el ENSMUSG
            if not rec["ensembl_id"] and gid.startswith(("ENSMUSG", "ENSG")):
                rec["ensembl_id"] = gid.split(".")[0]

            rows.append(rec)

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"No se han encontrado features '{feature}' en {gtf_path}")

    dup = df["gene_id"].duplicated()
    if dup.any():
        warnings.warn(
            f"{int(dup.sum())} gene_id duplicados en el GTF (genes en varios "
            f"scaffolds); me quedo con la primera aparición"
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
    Añade a adata.var los identificadores del GTF que usaste en STAR, casando
    por adata.var_names <-> gene_id. Guarda la procedencia en adata.uns.

    Con un GTF de RefSeq esto te da 'entrez_id' y 'mgi_id', que son los IDs
    estables. 'ensembl_id' quedará vacío: RefSeq no lo incluye (ver
    add_ensembl_ids()).
    """
    genes = parse_gtf_genes(gtf_path).set_index("gene_id")
    if columns is not None:
        genes = genes[list(columns)]

    hit = adata.var_names.isin(genes.index)
    if not hit.all():
        warnings.warn(
            f"{int((~hit).sum())}/{adata.n_vars} genes de la matriz no están en "
            f"el GTF. ¿Es el mismo fichero que usaste para el índice de STAR?"
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
    Descarga gene2ensembl.gz de NCBI (~200 MB) si no está ya en disco.

    La descarga es atómica: escribe en '<dest>.part', comprueba que el gzip
    se abre y solo entonces renombra. Así una descarga interrumpida nunca
    deja un fichero truncado que luego falle de forma críptica a mitad del
    parseo.
    """
    if os.path.exists(dest) and not overwrite:
        return dest

    parent = os.path.dirname(os.path.abspath(dest))
    os.makedirs(parent, exist_ok=True)
    tmp = dest + ".part"

    print(f"[gene2ensembl] descargando {url}")
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

        # integridad: si el .gz está truncado, esto revienta aquí y no dentro
        # del bucle de pandas media hora después
        with gzip.open(tmp, "rt") as fh:
            fh.readline()

        os.replace(tmp, dest)
        print(f"[gene2ensembl] guardado en {dest} ({os.path.getsize(dest) / 1e6:.0f} MB)")
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise

    return dest


def _read_gene2ensembl_header(path: str) -> list:
    """
    Lee la cabecera y le quita el '#' inicial. NCBI la escribe como
    '#tax_id', pero hay mirrors y versiones que no, así que no lo asumimos.
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
    Mapea Entrez GeneID -> Ensembl gene ID usando el fichero oficial de NCBI
    (gene2ensembl.gz de https://ftp.ncbi.nlm.nih.gov/gene/DATA/).

    Es la ruta reproducible: fichero versionado y offline, en vez de una API.
    El mapeo NO es completo ni siempre 1:1 — se anotan ambas cosas en var:
      ensembl_id     primer Ensembl ID (o vacío si no hay)
      ensembl_all    todos, separados por ';'
      ensembl_n      cuántos (0 = sin correspondencia)
    """
    if entrez_col not in adata.var.columns:
        raise ValueError(
            f"No hay '{entrez_col}' en var; ejecuta antes annotate_var_from_gtf()"
        )

    if not os.path.exists(gene2ensembl_path):
        if not download:
            raise FileNotFoundError(
                f"No existe {gene2ensembl_path} y download=False. "
                f"Descárgalo de {url}"
            )
        download_gene2ensembl(gene2ensembl_path, url=url)

    cols = _read_gene2ensembl_header(gene2ensembl_path)
    tax_col = next((c for c in cols if c.endswith("tax_id")), None)
    for needed in ("GeneID", "Ensembl_gene_identifier"):
        if needed not in cols:
            raise ValueError(
                f"{gene2ensembl_path} no tiene la columna '{needed}'. "
                f"Columnas encontradas: {cols}"
            )
    if tax_col is None:
        raise ValueError(f"No encuentro la columna tax_id. Columnas: {cols}")

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
        raise ValueError(f"Ninguna fila para tax_id={tax_id} en {gene2ensembl_path}")

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
        f"[ensembl] {n_map}/{adata.n_vars} genes mapeados "
        f"({100 * n_map / adata.n_vars:.1f}%), {n_multi} con más de un Ensembl ID"
    )
    adata.uns.setdefault("annotation", {})["gene2ensembl"] = os.path.basename(
        str(gene2ensembl_path)
    )
    return adata




# --------------------------------------------------------------------------
# carga principal
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
        True -> subdirectorio 'raw'; False -> 'filtered'.
    include_cellbender
        Segunda pasada: usa la salida de CellBender como X.
    keep_raw_layer
        Guarda también GeneFull pre-CellBender en layers['counts_raw'].
        Cuesta memoria pero es lo que permite auditar qué se ha quitado.
    min_cell_probability
        Si se indica y CellBender aporta 'cell_probability', filtra droplets
        antes de cargar Velocyto. Recomendado en modo raw: reduce de ~700k
        droplets a unos miles y ahorra muchísima memoria.
    restrict_to
        Lista explícita de barcodes a conservar (alternativa al filtro anterior).
    gtf_path
        GTF que usaste en STAR. Si se pasa, añade a var los identificadores
        estables (entrez_id, mgi_id), biotipo y coordenadas, y guarda la
        procedencia de la anotación en uns['annotation'].
    gene2ensembl_path
        gene2ensembl.gz de NCBI para mapear entrez_id -> Ensembl. Requiere
        gtf_path. Se descarga solo si no existe.
    cellbender_meta
        Dict libre (versión, argumentos) que se guarda en uns['cellbender'].
    strict
        Si True, falla cuando un barcode/gen no aparece en alguna matriz.
    check_velocyto_sum
        Comprueba si GeneFull ~= spliced + unspliced + ambiguous y lo guarda
        en adata.uns['velocyto_sum_check'].
    """
    dirsub = "raw" if raw_counts else "filtered"
    gene_path = f"{STARsolo_path}/{gene_dir}/{dirsub}"
    velo_path = f"{STARsolo_path}/Velocyto/{dirsub}"

    # ---- 1. GeneFull crudo -------------------------------------------------
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

        # OJO: el h5 "completo" de CellBender NO trae todos los barcodes de
        # entrada, solo los droplets que analizó (probable cells + additional,
        # es decir --total-droplets-included). El resto ni se evalúan.
        # Así que aquí mandan los droplets de CellBender, no los de GeneFull.
        keep = adata.obs_names.isin(adata_cb.obs_names)
        n_keep = int(keep.sum())
        if n_keep == 0:
            raise ValueError(
                "Ningún barcode de GeneFull aparece en la salida de CellBender. "
                "Suele ser un sufijo distinto ('-1'), no un problema de datos.\n"
                f"  GeneFull:   {list(adata.obs_names[:3])}\n"
                f"  CellBender: {list(adata_cb.obs_names[:3])}"
            )
        print(f"[cellbender] {n_keep} droplets analizados de {adata.n_obs} en GeneFull")
        adata = adata[keep].copy()

        # a partir de aquí el alineamiento de obs sí debe ser exacto
        row_pos = _align_indexer(adata_cb.obs_names, adata.obs_names,
                                 "CellBender obs", strict=True)
        # los genes, en cambio, pueden faltar: CellBender excluye del análisis
        # los que estima con background ~0 y puede no devolverlos
        col_pos = _align_indexer(adata_cb.var_names, adata.var_names,
                                 "CellBender var", strict=False)

        cb_X = _subset_positional(
            sp.csr_matrix(adata_cb.X), row_pos, col_pos
        ).astype(adata.X.dtype)

        # Genes ausentes en la salida: se rellenan con los counts CRUDOS, no
        # con ceros. CellBender los excluyó por tener background despreciable,
        # o sea que su valor descontaminado es el original. Ponerlos a cero
        # sería borrar genes en silencio.
        missing_genes = col_pos < 0
        if missing_genes.any():
            warnings.warn(
                f"{int(missing_genes.sum())} genes no están en la salida de "
                f"CellBender; se conservan sus counts crudos sin descontaminar"
            )
            # adata.X sigue siendo la matriz cruda en este punto
            cb_X = cb_X + adata.X @ sp.diags(missing_genes.astype(adata.X.dtype))
            adata.var["cellbender_returned"] = ~missing_genes

        adata.layers["counts_cellbender"] = sp.csr_matrix(cb_X)
        adata.X = adata.layers["counts_cellbender"].copy()

        # métricas POR GEN que aporta CellBender:
        #   ambient_expression  perfil del soup (fracción de los counts ambient
        #                       atribuible a cada gen; suma ~1 sobre genes).
        #                       Dice DE QUÉ está hecho el soup.
        #   cellbender_analyzed si el gen entró en el modelo o quedó excluido
        #                       por background despreciable.
        # No confundir ambient_expression con ambient_fraction_per_gene(), que
        # mide QUÉ FRACCIÓN DE ESE GEN se ha eliminado. Son cosas distintas:
        # un gen muy expresado puede dominar el soup y aun así perder un 2%.
        cb_var = adata_cb.var.reindex(adata.var_names)
        for col in ("ambient_expression", "cellbender_analyzed"):
            if col in cb_var.columns:
                adata.var[col] = cb_var[col].values

        # métricas por droplet que CellBender sí aporta y merece la pena guardar
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

    # ---- 3. filtrado de droplets ANTES de Velocyto -------------------------
    if restrict_to is not None:
        keep = adata.obs_names.isin(pd.Index(restrict_to))
        adata = adata[keep].copy()
    elif min_cell_probability is not None:
        if "cell_probability" not in adata.obs:
            raise ValueError("No hay 'cell_probability'; ¿corriste CellBender?")
        adata = adata[adata.obs["cell_probability"] >= min_cell_probability].copy()

    # ambient por gen: SOLO tiene sentido sobre las células que quedan.
    # Calculado antes del filtrado incluiría droplets vacíos, donde la
    # fracción eliminada es ~100%, e inflaría el número para todos los genes.
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
                raise ValueError(f"Dimensiones incoherentes en {layer}.mtx")
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

    # ---- 5. anotación opcional -------------------------------------------
    if gtf_path is not None:
        annotate_var_from_gtf(adata, gtf_path)
        if gene2ensembl_path is not None:
            add_ensembl_ids(adata, gene2ensembl_path)
    elif gene2ensembl_path is not None:
        raise ValueError(
            "gene2ensembl_path necesita gtf_path: los Ensembl IDs se mapean "
            "desde entrez_id, y entrez_id sale del GTF."
        )

    adata.obs_names_make_unique()
    return adata





def _safe_sheet_name(name: str, used: set, max_len: int = 31) -> str:
    """Nombres de hoja válidos en Excel + evita duplicados."""
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
    Guarda un Excel con una hoja por población.
    - Respeta el orden en el que aparecen las poblaciones en el df.
    - Exporta solo las columnas indicadas en `cols` si existen.
    """
    if group_col not in df.columns:
        raise ValueError(f"'{group_col}' no está en el DataFrame")

    groups_in_order = df[group_col].dropna().drop_duplicates().tolist()
    used = set()

    with pd.ExcelWriter(out_path) as writer:
        for g in groups_in_order:
            gdf = df[df[group_col] == g]
            sheet = _safe_sheet_name(g, used)
            use_cols = [c for c in cols if c in gdf.columns] or list(gdf.columns)
            gdf[use_cols].to_excel(writer, sheet_name=sheet, index=False)
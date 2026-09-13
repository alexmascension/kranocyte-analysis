"""
Comunicación célula-célula: CellPhoneDB v5 (método estadístico) y LIANA+.

DOS RUTAS, Y NO SON INTERCAMBIABLES
-----------------------------------
CellPhoneDB (secciones 1-4): base HUMANA, hay que orto-transferir tus símbolos
antes, y espera la matriz normalizada SIN log. A cambio maneja complejos
multi-subunidad, que para Il6 (IL6 + IL6R + IL6ST) es justo lo que importa.

LIANA+ (sección 5): recurso NATIVO de ratón ('mouseconsensus'), sin paso de
ortología, espera la matriz LOG-normalizada, y permite correr varios recursos
y ver qué sobrevive. Es la ruta que recomiendo empezar.

Fíjate en que la escala de la matriz es OPUESTA entre las dos. Es el error de
copiar-pegar más fácil de cometer aquí.

    mapa  = mouse_to_human(adata.var_names)              # una vez, se cachea
    paths = prepare_cpdb_input(adata, "cell_type", OUT, ortholog_map=mapa)
    res   = run_cpdb_statistical(**paths, cpdb_file_path=DB)
    il    = filter_interactions(res, ["Il6", "Il1b"], ortholog_map=mapa)

CELLPHONEDB ES HUMANA
---------------------
La base de datos son pares ligando-receptor HUMANOS, con símbolos humanos. Tus
datos son de ratón. No hay ninguna versión murina oficial, así que hay que
traducir por ortología ANTES de correr nada, y eso no es un detalle de formato:
es una suposición biológica fuerte. Que IL6-IL6R interactúe en humano no
garantiza que Il6-Il6ra haga lo mismo en ratón, y hay pares humanos cuyo
ortólogo murino sencillamente no existe (CXCL8/IL8 es el caso clásico).

Lo que NO hay que hacer es poner los símbolos en mayúsculas y seguir. Funciona
para la mayoría —Il6 -> IL6, Il1b -> IL1B— y falla en silencio justo donde el
nombre difiere entre especies (Trp53 no es TP53), que además son a menudo genes
interesantes. mouse_to_human() usa ortología de Ensembl vía pybiomart; el modo
'upper' existe solo como último recurso y avisa.

NO RECORTES GENES ANTES DE CORRER
---------------------------------
Quieres mirar Il6 e Il1b, pero el método estadístico permuta las etiquetas de
las células y compara contra la distribución nula de TODAS las interacciones,
y además necesita la matriz entera para evaluar los complejos (un receptor de
varias subunidades solo cuenta si todas se expresan). Si le pasas dos genes, el
nulo se calcula sobre dos genes y el resultado no significa nada. Se corre
completo y se filtra el RESULTADO: para eso está filter_interactions().

SOBRE LOS P-VALORES
-------------------
Salen de permutar etiquetas de célula, así que su resolución mínima es
1/iterations: con 1000 permutaciones no existe p < 0.001. CellPhoneDB no aplica
ninguna corrección por test múltiple y son miles de pares x combinaciones de
tipos celulares. Y como siempre en este pipeline: con 2 ratones las células son
pseudorréplicas, así que el test contesta "¿es esta media mayor que barajando
las etiquetas de ESTAS células?", no "¿pasa esto en el ratón?".

AVISO
-----
Las llamadas a cellphonedb de este módulo NO están probadas: la librería no se
pudo instalar en el entorno donde se escribió (fbpca no compila). La
orquestación —ortología, formato de entrada, filtrado del resultado— sí se ha
probado con mocks. Verifica la primera corrida a mano.
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
# 1. ortología ratón -> humano
# --------------------------------------------------------------------------

MIRRORS_ENSEMBL = ("http://www.ensembl.org", "http://useast.ensembl.org",
                   "http://asia.ensembl.org")

URL_MGI_HOM = ("https://www.informatics.jax.org/downloads/reports/"
               "HOM_MouseHumanSequence.rpt")


def _orto_biomart(reintentos: int = 4, espera: float = 5.0,
                  mirrors: Sequence[str] = MIRRORS_ENSEMBL,
                  verbose: bool = True) -> pd.DataFrame:
    """
    Ortología desde Ensembl con reintentos y rotación de espejo.

    BioMart devuelve HTTP 429 (too many requests) con bastante facilidad, y
    pybiomart hace DOS peticiones por consulta —una para la configuración del
    dataset y otra para los datos—, así que la probabilidad de topar con el
    límite es el doble de lo que parece. No es un error de tus datos ni del
    código: se reintenta con espera creciente y, si el espejo sigue negándose,
    se cambia de espejo.
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
                df.columns = ["mouse", "human", "tipo"]
                if verbose:
                    print(f"[orto] Ensembl OK ({host})")
                return df
            except Exception as e:
                ultimo = e
                pausa = espera * (2 ** intento)
                if verbose:
                    print(f"[orto] {host} intento {intento+1}/{reintentos} "
                          f"falló ({type(e).__name__}); espero {pausa:.0f}s")
                time.sleep(pausa)
    raise RuntimeError(
        f"BioMart no responde en ninguno de los espejos {list(mirrors)}. "
        f"Último error: {ultimo}. Prueba method='mgi', que es un fichero "
        f"estático sin límite de peticiones."
    ) from ultimo


def _orto_mgi(url: str = URL_MGI_HOM, verbose: bool = True) -> pd.DataFrame:
    """
    Ortología desde el informe de homología de MGI.

    Es un fichero estático, sin API ni límite de peticiones, y es la fuente que
    usa la propia comunidad de ratón. Cada clase de homología agrupa los genes
    ortólogos de las dos especies; nos quedamos con las clases que tienen
    EXACTAMENTE un gen de ratón y uno humano, que es la definición operativa de
    uno-a-uno.

    NO VERIFICADO desde el entorno donde se escribió esto: informatics.jax.org
    está bloqueado ahí, así que no he podido comprobar los nombres exactos de
    las columnas del fichero actual. Por eso se detectan en vez de fijarse, y
    si la detección falla el error te dice qué columnas ha encontrado.
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
            f"No reconozco el formato del fichero de MGI. Columnas: "
            f"{list(df.columns)}. Detectadas: key={col_key}, org={col_org}, "
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
        print(f"[orto] MGI: {len(out)} parejas uno-a-uno")
    return out


def mouse_to_human(
    genes: Optional[Iterable[str]] = None,
    cache_path: str = "mouse_human_orthologs.tsv",
    method: str = "mgi",
    one_to_one: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Tabla símbolo de ratón -> símbolo humano. DataFrame con ['mouse', 'human'].

    Se cachea en disco: la primera llamada descarga, las siguientes leen el
    fichero. Guárdalo con el análisis — las anotaciones de ortología cambian.

    method
        'mgi'      informe de homología de MGI. Fichero estático, sin API ni
                   límites de peticiones. Es el default por eso: BioMart
                   devuelve 429 con facilidad y deja el análisis colgado.
        'biomart'  Ensembl vía pybiomart, con reintentos y rotación de espejo.
        'upper'    poner en MAYÚSCULAS. Último recurso, sin red. Empareja la
                   mayoría de símbolos y FALLA EN SILENCIO en los que difieren
                   entre especies (Trp53 != TP53). Avisa al usarlo.

    one_to_one
        Descarta genes con varios ortólogos en el otro sentido. Es lo
        conservador: con un mapeo uno-a-muchos hay que decidir qué gen humano
        recibe la expresión del murino, y esa decisión la acabaría tomando el
        orden de las filas.
    """
    if os.path.exists(cache_path):
        mapa = pd.read_csv(cache_path, sep="\t", dtype=str)
        if verbose:
            print(f"[orto] leído de caché: {cache_path} ({len(mapa)} pares)")
        return mapa

    if method == "upper":
        if genes is None:
            raise ValueError("method='upper' necesita la lista de genes")
        mapa = pd.DataFrame({"mouse": list(genes)})
        mapa["human"] = mapa["mouse"].str.upper()
        warnings.warn(
            "Ortología por MAYÚSCULAS. Empareja la mayoría de los símbolos pero "
            "no los que difieren entre especies (Trp53 != TP53), y no avisa de "
            "cuáles ha perdido. Úsalo solo si no puedes bajar la tabla real, y "
            "dilo en los métodos."
        )
    elif method == "mgi":
        mapa = _orto_mgi(verbose=verbose)
    elif method == "biomart":
        df = _orto_biomart(verbose=verbose)
        df = df.dropna(subset=["mouse", "human"])
        df = df[df["human"] != ""]
        if one_to_one:
            n0 = len(df)
            df = df[df["tipo"] == "ortholog_one2one"]
            if verbose:
                print(f"[orto] {len(df)} de {n0} pares son uno-a-uno")
        mapa = df[["mouse", "human"]].drop_duplicates()
    else:
        raise ValueError("method debe ser 'mgi', 'biomart' o 'upper'")

    if one_to_one and method != "upper":
        mapa = mapa[~mapa["mouse"].duplicated(keep=False)]
        mapa = mapa[~mapa["human"].duplicated(keep=False)]

    mapa.to_csv(cache_path, sep="\t", index=False)
    if verbose:
        print(f"[orto] {len(mapa)} pares guardados en {cache_path}")
        if genes is not None:
            g = set(genes)
            print(f"[orto] cubre el {len(g & set(mapa['mouse']))/max(len(g),1):.1%} "
                  f"de tus genes")
            faltan = [x for x in ("Il6", "Il1b") if x in g and x not in set(mapa["mouse"])]
            if faltan:
                print(f"[orto] OJO: {faltan} no tienen ortólogo uno-a-uno en "
                      f"esta tabla. Compruébalo a mano antes de seguir.")
    return mapa


# --------------------------------------------------------------------------
# 2. preparar la entrada
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
    Escribe meta.tsv y counts.h5ad listos para CellPhoneDB. Devuelve las rutas.

    QUÉ MATRIZ. CellPhoneDB promedia la expresión por tipo celular y el
    estadístico es la media de las dos medias. Una media de valores log NO es
    el log de la media: log comprime los valores altos, así que un ligando muy
    expresado en pocas células pesa distinto en una escala y en la otra. Por eso
    aquí se entrega la matriz normalizada SIN log:

      layer      úsala si ya tienes una capa normalizada no-log.
      log_layer  si no, se reconstruye con expm1() sobre la capa log-normalizada
                 (exacto: deshace el log1p que aplicaste tú).
      counts_layer  último recurso: normalización por library size sobre counts.

    Los tutoriales de CellPhoneDB no son consistentes en esto y verás gente
    pasando datos log. Si quieres replicar a alguien, pásale su capa por
    `layer=` y deja constancia; pero por defecto va sin log, que es lo coherente
    con lo que el método calcula.

    min_cells
        Las poblaciones muy pequeñas dan medias inestables y la permutación no
        lo arregla: se excluyen y se avisa. No es lo mismo que no tener la
        población; es que no puedes decir nada de ella aquí.
    """
    os.makedirs(outdir, exist_ok=True)

    # --- células
    obs = adata.obs[group_key].astype(str)
    tam = obs.value_counts()
    fuera = set(exclude_groups) | set(tam.index[tam < min_cells])
    if fuera and verbose:
        print(f"[cpdb] excluidas {sorted(fuera)} "
              f"({dict(tam[list(fuera)])} células, min_cells={min_cells})")
    keep_cells = ~obs.isin(fuera).to_numpy()

    # --- matriz sin log
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
    else:
        M = adata.layers[counts_layer] if counts_layer in adata.layers else adata.X
        M = sp.csr_matrix(M, dtype=np.float64)
        libs = np.asarray(M.sum(axis=1)).ravel()
        libs[libs == 0] = 1.0
        M = sp.diags(1e4 / libs) @ M
        origen = f"CPM(1e4) sobre '{counts_layer}'"
    if verbose:
        print(f"[cpdb] matriz de expresión: {origen} (sin log)")

    # --- ortología
    mapa = dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
    en_mapa = np.array([g in mapa for g in adata.var_names])
    if verbose:
        print(f"[cpdb] {en_mapa.sum()} de {adata.n_vars} genes tienen ortólogo "
              f"humano ({en_mapa.mean():.1%})")
        if en_mapa.mean() < 0.4:
            print("[cpdb] AVISO: cobertura baja. Comprueba la tabla de "
                  "ortología antes de interpretar nada.")

    sub = ad.AnnData(
        X=sp.csr_matrix(M)[keep_cells][:, en_mapa],
        obs=pd.DataFrame(index=adata.obs_names[keep_cells]),
        var=pd.DataFrame(index=[mapa[g] for g in adata.var_names[en_mapa]]),
    )
    if not sub.var_names.is_unique:
        n_dup = int(sub.var_names.duplicated().sum())
        warnings.warn(f"{n_dup} símbolos humanos duplicados tras el mapeo; me "
                      f"quedo con la primera aparición de cada uno")
        sub = sub[:, ~sub.var_names.duplicated()].copy()

    # --- comprobación de genes de interés
    if verbose:
        for g in ("Il6", "Il1b"):
            if g in mapa:
                print(f"[cpdb] {g} -> {mapa[g]}: "
                      f"{'presente' if mapa[g] in sub.var_names else 'AUSENTE'}")

    meta = pd.DataFrame({"Cell": sub.obs_names, "cell_type": obs[keep_cells].values})
    meta_path = os.path.join(outdir, "meta.tsv")
    counts_path = os.path.join(outdir, "counts.h5ad")
    meta.to_csv(meta_path, sep="\t", index=False)
    sub.write_h5ad(counts_path)

    if verbose:
        print(f"[cpdb] {sub.n_obs} células x {sub.n_vars} genes -> {counts_path}")
        print(meta["cell_type"].value_counts().to_string())
    return {"meta_file_path": meta_path, "counts_file_path": counts_path}


# --------------------------------------------------------------------------
# 2b. mirar qué hay DENTRO de la base de datos
# --------------------------------------------------------------------------

def inspect_cpdb_database(cpdb_file_path: str, verbose: bool = True) -> dict:
    """
    Abre el .zip de CellPhoneDB y devuelve sus tablas como {nombre: DataFrame}.

    No usa la API de la librería: el .zip son CSVs y se leen directamente. Así
    esto funciona aunque cambie la API entre versiones, y puedes mirar la base
    ANTES de correr nada.

    Las tablas que importan son cuatro y están relacionadas entre sí:
      gene         gene_name / hgnc_symbol / uniprot / ensembl
      protein      una fila por uniprot, con si es receptor, secretado, etc.
      complex      complejos, definidos como uniprot_1..uniprot_4
      interaction  los pares, donde cada partner es un uniprot O el nombre de
                   un complejo
    Por eso buscar "IL6" en la tabla de interacciones no basta: los partners no
    son símbolos.
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
                warnings.warn(f"No he podido leer {n}: {e}")
    if verbose:
        print(f"[db] {os.path.basename(cpdb_file_path)}: {len(tablas)} tablas")
        for k, v in sorted(tablas.items()):
            print(f"  {k:35s} {v.shape[0]:>7} x {v.shape[1]}")
    return tablas


def _tabla(db: Mapping[str, pd.DataFrame], *claves,
           requiere: Sequence[str] = ()) -> Optional[pd.DataFrame]:
    """
    Busca una tabla por nombre, eligiendo el MEJOR candidato y no el primero.

    Importa porque el zip trae varias tablas cuyo nombre contiene 'gene'
    (gene_input, gene_synonym_to_gene_name, ...) y quedarse con la primera del
    listado del zip es una lotería que depende del orden de los ficheros. Se
    puntúa: primero las que tienen las columnas que hacen falta, luego las de
    nombre más corto (gene_input gana a gene_synonym_to_gene_name).
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
    ¿Cuáles de TUS genes existen en CellPhoneDB?

    Se hace ANTES de correr. Si un gen no está en la base, no puede aparecer en
    ningún resultado, y su ausencia en la salida no significa "no interactúa":
    significa "CellPhoneDB no lo contempla". Son cosas muy distintas y la
    segunda no se puede leer de la tabla de resultados.

    genes en nomenclatura de ratón si pasas ortholog_map; si no, humanos.
    """
    gt = _tabla(db, "gene", requiere=("uniprot", "gene_name", "hgnc"))
    if gt is None:
        raise RuntimeError(f"No encuentro la tabla de genes. Tablas: {list(db)}")

    cols = [c for c in gt.columns
            if any(k in c.lower() for k in ("gene_name", "hgnc", "symbol"))]
    if not cols:
        raise RuntimeError(f"La tabla de genes no trae símbolos. Columnas: "
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
        filas.append({"gen": g, "humano": h,
                      "tiene_ortologo": h is not None,
                      "en_cellphonedb": bool(h) and h in en_db})
    out = pd.DataFrame(filas)

    if verbose:
        n = len(out)
        print(f"[db] {out['tiene_ortologo'].sum()}/{n} con ortólogo humano; "
              f"{out['en_cellphonedb'].sum()}/{n} presentes en CellPhoneDB")
        fuera = out[~out["en_cellphonedb"]]
        if len(fuera) and len(fuera) <= 30:
            print("[db] NO están en la base: "
                  + ", ".join(f"{r.gen}" + ("" if r.tiene_ortologo else " (sin ortólogo)")
                              for r in fuera.itertuples()))
    return out


def cpdb_interactions_for(
    db: Mapping[str, pd.DataFrame],
    genes: Sequence[str],
    ortholog_map: Optional[pd.DataFrame] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    TODAS las interacciones de la base que involucran esos genes, con
    independencia de tus datos.

    Esto es el techo: el máximo de pares que podrías llegar a recuperar. Si
    IL6 participa en 6 interacciones y tu resultado devuelve 2, sabes que el
    filtro han sido tus datos; si devuelve 6, has agotado lo que la base
    contempla. Sin este número, la tabla de resultados no se puede interpretar.

    Resuelve la cadena símbolo -> uniprot -> complejos que lo contienen ->
    interacciones donde participa el uniprot o alguno de esos complejos. Buscar
    el símbolo directamente en la tabla de interacciones no vale: ahí los
    partners son uniprots y nombres de complejo.
    """
    gt = _tabla(db, "gene", requiere=("uniprot", "gene_name", "hgnc"))
    ct = _tabla(db, "complex", requiere=("complex_name", "uniprot_1"))
    it_ = _tabla(db, "interaction", requiere=("partner_a", "partner_b"))
    if gt is None or it_ is None:
        raise RuntimeError(f"Faltan tablas necesarias. Tablas: {list(db)}")

    if ortholog_map is not None:
        mapa = dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
        humanos = {mapa[g] for g in genes if g in mapa}
        perdidos = [g for g in genes if g not in mapa]
        if perdidos:
            warnings.warn(f"Sin ortólogo humano: {perdidos}")
    else:
        humanos = set(genes)

    col_sym = [c for c in gt.columns
               if any(k in c.lower() for k in ("gene_name", "hgnc", "symbol"))]
    col_uni = next((c for c in gt.columns if "uniprot" in c.lower()), None)
    if col_uni is None:
        raise RuntimeError(f"La tabla de genes no trae uniprot: {list(gt.columns)}")

    sel = pd.Series(False, index=gt.index)
    for c in col_sym:
        sel |= gt[c].astype(str).isin(humanos)
    uniprots = set(gt.loc[sel, col_uni].dropna().astype(str))

    # complejos que contienen alguno de esos uniprots
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
            raise RuntimeError(f"No encuentro las columnas de partner en la "
                               f"tabla de interacciones: {list(it_.columns)}")

    hit = it_[it_[col_pa].astype(str).isin(partners)
              | it_[col_pb].astype(str).isin(partners)].copy()

    # traducir uniprot -> símbolo para que la tabla se pueda leer
    u2s = {}
    for c in col_sym:
        u2s.update(dict(zip(gt[col_uni].astype(str), gt[c].astype(str))))
    hit["simbolo_a"] = hit[col_pa].astype(str).map(u2s).fillna(hit[col_pa])
    hit["simbolo_b"] = hit[col_pb].astype(str).map(u2s).fillna(hit[col_pb])

    if verbose:
        print(f"[db] tablas usadas: gene={gt.shape}, "
              f"complex={None if ct is None else ct.shape}, "
              f"interaction={it_.shape}")
        print(f"[db] {sorted(humanos)} -> {len(uniprots)} uniprots, "
              f"{len(complejos)} complejos que los contienen")
        if complejos:
            print(f"[db] complejos: {sorted(complejos)[:8]}")
        print(f"[db] {len(hit)} interacciones en la base los involucran")
        cols_ver = [c for c in ("simbolo_a", "simbolo_b", "annotation_strategy",
                                "classification", "directionality", "is_integrin")
                    if c in hit.columns]
        if len(hit):
            print(hit[cols_ver].head(25).to_string(index=False))
    return hit


# --------------------------------------------------------------------------
# 3. el método estadístico
# --------------------------------------------------------------------------

def download_cpdb_database(target_dir: str, version: str = "v5.0.0",
                           verbose: bool = True) -> str:
    """Baja la base de datos de CellPhoneDB. Devuelve la ruta al .zip."""
    from cellphonedb.utils import db_utils

    os.makedirs(target_dir, exist_ok=True)
    db_utils.download_database(target_dir, version)
    zips = [f for f in os.listdir(target_dir) if f.endswith(".zip")]
    if not zips:
        raise RuntimeError(f"No se ha descargado ningún .zip en {target_dir}")
    ruta = os.path.join(target_dir, sorted(zips)[-1])
    if verbose:
        print(f"[cpdb] base de datos: {ruta}")
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
    cpdb_statistical_analysis_method.call(), con los parámetros que importan
    explícitos.

    threshold
        Fracción MÍNIMA de células de un tipo que deben expresar el gen para
        considerarlo expresado ahí. 0.1 es el default de CellPhoneDB y es una
        decisión, no una constante: subirlo elimina interacciones sostenidas por
        una minoría de células (que en una población heterogénea pueden ser
        reales) y bajarlo llena el resultado de pares dependientes de tres
        células. Declara el valor que uses.

    iterations
        Número de permutaciones. Fija el suelo del p-valor: con 1000, el mínimo
        posible es 0.001 y no existe nada por debajo. Si vas a comparar miles de
        pares, ese suelo es el que manda, no tu alpha.

    microenvs_file_path
        TSV opcional (cell_type, microenvironment) que restringe qué tipos
        pueden interactuar entre sí. Si sabes que dos poblaciones no comparten
        nicho, decirlo aquí reduce el número de tests y es más honesto que
        filtrarlo después.
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
        print(f"[cpdb] p-valores por permutación: suelo = {1/iterations:g}. "
              f"CellPhoneDB NO corrige por test múltiple.")
    return res


# --------------------------------------------------------------------------
# 4. filtrar el resultado por genes concretos
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
    Interacciones que involucran genes concretos, en formato largo.

    genes en NOMENCLATURA DE RATÓN (Il6, Il1b); se traducen con ortholog_map.
    Si no pasas mapa se asume que ya son símbolos humanos.

    include_complexes
        Además de las filas donde gene_a/gene_b es exactamente tu gen, incluye
        aquellas cuyo 'interacting_pair' lo menciona: en CellPhoneDB muchos
        receptores son complejos y el gen aparece dentro del nombre del complejo
        y no en gene_a/gene_b. Se marca con la columna 'via_complejo' para que
        puedas separarlas — no son lo mismo y confundirlas infla el recuento.

    Devuelve una fila por (interacción, par de tipos celulares) con su media y
    su p-valor, ya filtrada por max_pvalue.
    """
    means = res["means"]
    pvals = res["pvalues"]

    if ortholog_map is not None:
        mapa = dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
        humanos, perdidos = [], []
        for g in genes:
            (humanos.append(mapa[g]) if g in mapa else perdidos.append(g))
        if perdidos:
            warnings.warn(f"Sin ortólogo humano: {perdidos}")
    else:
        humanos = list(genes)
    humanos = set(humanos)
    if verbose:
        print(f"[filtro] buscando {sorted(humanos)}")

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
            print("[filtro] ninguna interacción con esos genes")
        return pd.DataFrame()

    m = means.loc[sel, meta_cols + par_cols]
    p = pvals.loc[sel, par_cols] if set(par_cols).issubset(pvals.columns) else None

    largo = m.melt(id_vars=meta_cols, var_name="par_celular", value_name="media")
    if p is not None:
        pl = pvals.loc[sel, meta_cols + par_cols].melt(
            id_vars=meta_cols, var_name="par_celular", value_name="pvalue")
        largo = largo.merge(pl, on=meta_cols + ["par_celular"], how="left")

    largo["via_complejo"] = ~(
        largo["gene_a"].isin(humanos) | largo["gene_b"].isin(humanos))

    if "interaction_scores" in res and isinstance(res["interaction_scores"], pd.DataFrame):
        sc = res["interaction_scores"]
        if "interacting_pair" in sc.columns:
            sl = sc[sc["interacting_pair"].isin(largo["interacting_pair"])].melt(
                id_vars=["interacting_pair"], var_name="par_celular",
                value_name="score")
            largo = largo.merge(sl, on=["interacting_pair", "par_celular"], how="left")

    if "pvalue" in largo:
        largo = largo[largo["pvalue"] <= max_pvalue]
    largo = largo[largo["media"] > 0].sort_values(
        ["pvalue", "media"] if "pvalue" in largo else ["media"],
        ascending=[True, False] if "pvalue" in largo else [False])

    if verbose:
        print(f"[filtro] {len(largo)} filas (interacción x par celular) con "
              f"p <= {max_pvalue}")
        n_dir = int((~largo["via_complejo"]).sum())
        print(f"[filtro] {n_dir} por gen directo, {len(largo)-n_dir} solo por "
              f"aparecer en el nombre de un complejo")
        if len(largo):
            print(largo.groupby("interacting_pair").size()
                  .sort_values(ascending=False).head(10).to_string())
    return largo.reset_index(drop=True)


def plot_interactions(largo: pd.DataFrame, n_top: int = 25, figsize=(11, 8)):
    """Dotplot: interacción x par celular, tamaño = media, color = -log10(p)."""
    import matplotlib.pyplot as plt

    if not len(largo):
        raise ValueError("nada que dibujar")
    top = (largo.groupby("interacting_pair")["media"].max()
           .sort_values(ascending=False).head(n_top).index.tolist())
    sub = largo[largo["interacting_pair"].isin(top)].copy()

    pares = sorted(sub["par_celular"].unique())
    filas = [t for t in top if t in set(sub["interacting_pair"])]
    xi = {p: i for i, p in enumerate(pares)}
    yi = {t: i for i, t in enumerate(filas)}
    col = (-np.log10(sub["pvalue"].clip(lower=1e-4))
           if "pvalue" in sub else sub["media"])

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    s = ax.scatter([xi[p] for p in sub["par_celular"]],
                   [yi[t] for t in sub["interacting_pair"]],
                   s=20 + 60 * sub["media"] / max(sub["media"].max(), 1e-9),
                   c=col, cmap="magma_r", edgecolors="none")
    ax.set_xticks(range(len(pares)))
    ax.set_xticklabels(pares, rotation=90, fontsize=7)
    ax.set_yticks(range(len(filas)))
    ax.set_yticklabels(filas, fontsize=8)
    ax.grid(alpha=.2)
    fig.colorbar(s, ax=ax, label="-log10(p)" if "pvalue" in sub else "media",
                 shrink=.5)
    return fig, ax


# ==========================================================================
# LIANA+: varios métodos y varios recursos, con consenso
# ==========================================================================
#
# Por qué LIANA además de CellPhoneDB:
#
#   1. Recurso NATIVO DE RATÓN ('mouseconsensus'), así que desaparece el paso
#      de ortología. Con matiz: parte de ese recurso se construyó transfiriendo
#      anotación humana, así que la suposición no se elimina, se traslada a
#      quien lo curó. Mejor que hacerlo tú a mano; no es magia.
#   2. Permite correr VARIOS recursos y ver qué conclusiones sobreviven, que es
#      la única forma honesta de responder "¿importa qué base uso?".
#
# OJO CON LA ESCALA DE LA MATRIZ. LIANA espera datos LOG-NORMALIZADOS, al revés
# que el CellPhoneDB standalone de arriba, al que hay que darle la matriz sin
# log. Tu adata.X (norm_scran_log1p) ya vale. Lo que sí hay que poner es
# use_raw=False, porque LIANA mira adata.raw por defecto.

def liana_resources(verbose: bool = True) -> list:
    """Recursos que LIANA sirve hoy. Pregunta, no adivines la lista."""
    import liana as li

    rs = list(li.rs.show_resources())
    if verbose:
        print(f"[liana] {len(rs)} recursos: {sorted(rs)}")
        murinos = [r for r in rs if "mouse" in r.lower() or "murine" in r.lower()]
        if murinos:
            print(f"[liana] nativos de ratón: {murinos}")
    return rs


def resource_coverage(
    genes: Sequence[str],
    resources: Sequence[str] = ("mouseconsensus",),
    ortholog_map: Optional[pd.DataFrame] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    ¿En qué recursos existe cada gen, y con qué partners? ANTES de correr.

    Es el techo del análisis, igual que cpdb_interactions_for(): si Il6
    participa en 6 pares del recurso y tu corrida devuelve 2, el filtro han
    sido tus datos; si devuelve 6, has agotado lo que el recurso contempla. Sin
    este denominador la tabla de resultados no se puede leer, y la ausencia de
    un par no se puede distinguir de "no está anotado".

    ortholog_map
        Solo si el recurso es HUMANO. Con 'mouseconsensus' NO se traduce: los
        símbolos ya son de ratón y traducirlos los rompería.
    """
    import liana as li

    mapa = (dict(zip(ortholog_map["mouse"], ortholog_map["human"]))
            if ortholog_map is not None else {})

    filas = []
    for r in resources:
        df = li.rs.select_resource(r)
        col_l = "ligand" if "ligand" in df.columns else df.columns[0]
        col_r = "receptor" if "receptor" in df.columns else df.columns[1]
        # las subunidades de un complejo vienen unidas por '_'
        lig = df[col_l].astype(str).str.split("_")
        rec = df[col_r].astype(str).str.split("_")
        for g in genes:
            buscado = mapa.get(g, g)
            como_lig = lig.apply(lambda x: buscado in x)
            como_rec = rec.apply(lambda x: buscado in x)
            sel = como_lig | como_rec
            filas.append({
                "recurso": r, "gen": g, "buscado": buscado,
                "n_pares": int(sel.sum()),
                "como_ligando": int(como_lig.sum()),
                "como_receptor": int(como_rec.sum()),
                "partners": "; ".join(sorted(set(
                    df.loc[sel, col_r].astype(str)[como_lig[sel]].tolist()
                    + df.loc[sel, col_l].astype(str)[como_rec[sel]].tolist()
                ))[:12]),
            })
    out = pd.DataFrame(filas)
    if verbose:
        print(out[["recurso", "gen", "buscado", "n_pares", "como_ligando",
                   "como_receptor"]].to_string(index=False))
        vacios = out[out["n_pares"] == 0]
        if len(vacios):
            print(f"\n[liana] AVISO: {list(zip(vacios.gen, vacios.recurso))} no "
                  f"aparecen en esos recursos. Su ausencia del resultado NO "
                  f"significa que no interactúen: significa que no están "
                  f"anotados.")
        for r in out["recurso"].unique():
            sub = out[(out.recurso == r) & (out.n_pares > 0)]
            for f in sub.itertuples():
                print(f"\n[liana] {f.gen} en {r}: {f.partners}")
    return out


def run_liana(
    adata: ad.AnnData,
    groupby: str,
    resource: str = "mouseconsensus",
    expr_prop: float = 0.1,
    n_perms: int = 1000,
    use_raw: bool = False,
    key_added: str = "liana_res",
    seed: int = 0,
    n_jobs: int = 4,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    rank_aggregate de LIANA: varios métodos, consenso por agregación de rangos.

    DOS RANKINGS, NO UN P-VALOR. La salida trae:

      magnitude_rank    ¿está esta interacción fuertemente expresada en este
                        par de tipos celulares?
      specificity_rank  ¿lo está MÁS AQUÍ que en los demás pares?

    Miden cosas distintas y confundirlos es el error clásico. Que FAP.3 exprese
    Il1b no es noticia si lo expresa todo el mundo: para eso está la
    especificidad. Para Il6/Il1b casi seguro quieres ordenar por
    specificity_rank.

    Y el consenso es una AGREGACIÓN DE RANGOS, no un metaanálisis: su p-valor
    mide concordancia entre métodos, no evidencia biológica. No lo reportes
    como si fuera un p-valor de la interacción.

    expr_prop
        Fracción mínima de células del tipo que deben expresar el gen. Mismo
        papel que el 'threshold' de CellPhoneDB, misma condición de decisión
        declarable.

    n_perms
        Empieza en 100 para cronometrar. Con 9.000 células y varios métodos,
        1000 permutaciones no son gratis.
    """
    import time

    import liana as li

    if use_raw is False and adata.raw is not None and verbose:
        print("[liana] use_raw=False: se usa adata.X. Comprueba que está "
              "LOG-normalizada (LIANA la espera así, al revés que el "
              "CellPhoneDB standalone).")

    faltan = [c for c in (groupby,) if c not in adata.obs]
    if faltan:
        raise ValueError(f"No hay obs{faltan}")
    n_grupos = adata.obs[groupby].nunique()
    if verbose:
        print(f"[liana] {resource}, {n_grupos} grupos, expr_prop={expr_prop}, "
              f"n_perms={n_perms}")

    t0 = time.time()
    li.mt.rank_aggregate(
        adata, groupby=groupby, resource_name=resource, expr_prop=expr_prop,
        use_raw=use_raw, n_perms=n_perms, seed=seed, n_jobs=n_jobs,
        key_added=key_added, verbose=verbose, inplace=True,
    )
    res = adata.uns[key_added].copy()
    res["recurso"] = resource

    if verbose:
        print(f"[liana] {len(res)} interacciones en {time.time()-t0:.0f}s")
        print("[liana] ordena por 'specificity_rank' para lo específico y por "
              "'magnitude_rank' para lo abundante; no son lo mismo.")
    return res


# ---------------------------------------------------------------------------
# 5b. Contexto de expresión para la tabla de LIANA
# ---------------------------------------------------------------------------
# LIANA devuelve puntuaciones agregadas pero NO las medias de las que salen.
# Sin ellas no se puede juzgar una fila: no sabes si un `specificity_rank`
# bajísimo viene de un receptor bien expresado en una población o de un gen
# que solo se detecta en 30 células. Estas funciones reconstruyen las medias
# por grupo y las de fondo, y —esto es lo importante— VERIFICAN la
# reconstrucción contra las columnas que sí trae LIANA. Si mi regla de
# complejos o mi definición de logFC no coinciden con las suyas, el informe
# lo dice en vez de dejarte columnas plausibles y falsas.

_EPS_CHK = 1e-4


def _matriz_expresion(adata, use_raw=False, layer=None):
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True pero adata.raw es None")
        return adata.raw.X, pd.Index(adata.raw.var_names)
    if layer is not None:
        return adata.layers[layer], pd.Index(adata.var_names)
    return adata.X, pd.Index(adata.var_names)


def _sumas_por_grupo(X, etiquetas, grupos):
    """Sumas y tamaños por grupo, para X y para expm1(X)."""
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
    Valor de cada complejo en cada grupo.

    LIANA reensambla los complejos por la subunidad de MENOR MEDIA (regla
    'min'), y luego se queda con TODAS las estadísticas de esa subunidad —no
    con el mínimo de cada estadística por separado. Por eso hace falta también
    devolver qué subunidad ganó, para poder sacar su logFC.

    OJO: M viene como GENES x GRUPOS. Las subunidades son filas, no columnas;
    el mínimo va sobre axis=0. Indexarlo al revés no da error, da NaN, que es
    peor.
    """
    medias, elegida = {}, {}
    perdidos = []
    for c in complejos:
        subs = str(c).split("_")
        if not all(s in M.index for s in subs):
            perdidos.append(c)
            continue
        sub_M = M.loc[subs]                       # subunidades x grupos
        medias[c] = sub_M.min(axis=0)
        elegida[c] = sub_M.idxmin(axis=0)
    return (pd.DataFrame(medias).T, pd.DataFrame(elegida).T, perdidos)


def _lookup(df_ancho, filas, columnas):
    """df_ancho.at[fila, col] vectorizado, con NaN donde falte."""
    if not len(df_ancho):
        return np.full(len(filas), np.nan)
    i = df_ancho.index.get_indexer(filas)
    j = df_ancho.columns.get_indexer(columnas)
    V = df_ancho.to_numpy()
    # el mismo lookup sirve para medias (float) y para nombres de subunidad
    # (str); si reservo float y luego meto strings, revienta.
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
    Añade a la tabla de LIANA el contexto de expresión que le falta.

    Columnas nuevas
    ---------------
    L_source, R_target   media del ligando en el emisor y del receptor en el
                         receptor. Son las dos piezas de las que salen
                         `lr_means` y `expr_prod`.
    L_target, R_source   las cruzadas. Sirven para ver reciprocidad y
                         autocrinía: si R_source también es alto, el emisor
                         se está escuchando a sí mismo.
    L_global, R_global   media sobre TODAS las células. El denominador mental
                         para decidir si 0.6 es mucho o poco en tus datos.
    lfc_L, lfc_R         logFC del ligando en el emisor y del receptor en el
                         receptor, cada uno contra el resto de poblaciones.
                         Su media es exactamente `lr_logfc`.
    spec_weight_oe       `spec_weight` dividido por su valor bajo uniformidad,
                         1/n_poblaciones². Es un observado/esperado, NO un
                         odds ratio: no hay odds en ninguna parte, es un
                         cociente de fracciones. 1 = como el reparto uniforme,
                         3 = tres veces más concentrado en este par.

    Verificación
    ------------
    Se comprueba que (L_source+R_target)/2 == lr_means, que
    L_source*R_target == expr_prod y que (lfc_L+lfc_R)/2 == lr_logfc. Para el
    logFC se prueban las DOS definiciones posibles (diferencia de medias sobre
    la matriz log, y log2 de las medias des-logueadas) y se queda la que
    reproduce la columna de LIANA. Si ninguna cuadra, avisa y no te da lfc_L
    ni lfc_R, porque una columna mal reconstruida es peor que ninguna.

    OJO con L_global: es la media sobre todas las células, no la suma de las
    medias por grupo. NATMI usa lo segundo para su denominador, y las dos solo
    coinciden si todas las poblaciones tienen el mismo tamaño.
    """
    out = res.copy()
    col_l = "ligand_complex" if "ligand_complex" in out else "ligand"
    col_r = "receptor_complex" if "receptor_complex" in out else "receptor"

    if groupby not in adata.obs:
        raise ValueError(f"No hay obs['{groupby}']")

    X, var_names = _matriz_expresion(adata, use_raw=use_raw, layer=layer)
    etiquetas = adata.obs[groupby].astype(str).to_numpy()
    grupos = sorted(pd.unique(etiquetas))
    n_pop = len(grupos)
    n_tot = adata.n_obs

    S, Se, n = _sumas_por_grupo(X, etiquetas, grupos)
    tot = np.sum([S[g] for g in grupos], axis=0)
    tot_e = np.sum([Se[g] for g in grupos], axis=0)

    # medias por grupo (matriz tal cual, que se asume log-normalizada)
    M = pd.DataFrame({g: S[g] / max(n[g], 1) for g in grupos}, index=var_names)
    glob = pd.Series(tot / n_tot, index=var_names)

    # dos candidatos de logFC, uno por definición
    lfcA, lfcB = {}, {}
    for g in grupos:
        resto = max(n_tot - n[g], 1)
        lfcA[g] = S[g] / max(n[g], 1) - (tot - S[g]) / resto
        lfcB[g] = (np.log2(Se[g] / max(n[g], 1) + 1)
                   - np.log2((tot_e - Se[g]) / resto + 1))
    LFC = {"A_dif_medias_log": pd.DataFrame(lfcA, index=var_names),
           "B_log2_medias_crudas": pd.DataFrame(lfcB, index=var_names)}

    # --- complejos
    comps = pd.unique(np.concatenate([out[col_l].astype(str).to_numpy(),
                                      out[col_r].astype(str).to_numpy()]))
    CM, CE, perdidos = _complejos_por_grupo(M, comps)
    if perdidos and verbose:
        print(f"[contexto] {len(perdidos)} complejos con subunidades ausentes "
              f"de var_names (ej. {perdidos[:5]}); esas filas van a NaN")
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

    # --- verificación de las medias
    chks = {}
    if "lr_means" in out:
        chks["lr_means"] = np.nanmax(np.abs(
            (out["L_source"] + out["R_target"]) / 2 - out["lr_means"]))
    if "expr_prod" in out:
        chks["expr_prod"] = np.nanmax(np.abs(
            out["L_source"] * out["R_target"] - out["expr_prod"]))

    # --- logFC: elegir la definición que reproduce lr_logfc
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
                f"No consigo reproducir 'lr_logfc' con ninguna de las dos "
                f"definiciones (mejor residuo {dif_lfc[elegido]:.4g} con "
                f"'{elegido}'). No añado lfc_L/lfc_R: prefiero no darte "
                f"columnas que no sé de dónde salen.")
            elegido = None

    # --- observado/esperado de NATMI
    if "spec_weight" in out:
        out["spec_weight_oe"] = out["spec_weight"] * (n_pop ** 2)

    if verbose:
        print(f"[contexto] {n_pop} poblaciones, {n_tot} células; "
              f"esperado de spec_weight bajo uniformidad = "
              f"{1/n_pop**2:.4f}")
        for k, v in chks.items():
            estado = "OK" if v <= _EPS_CHK else "NO CUADRA"
            print(f"[contexto] check {k:10s} residuo máx {v:.3g}  {estado}")
        if dif_lfc:
            for k, v in dif_lfc.items():
                print(f"[contexto] check lr_logfc  {k:22s} residuo máx "
                      f"{v:.3g}")
            if elegido:
                print(f"[contexto] logFC reconstruido con '{elegido}'")
        if any(v > _EPS_CHK for v in chks.values()):
            print("[contexto] AVISO: si lr_means/expr_prod no cuadran, la "
                  "matriz que le estás pasando aquí NO es la misma con la que "
                  "corriste LIANA (¿use_raw, layer, o el adata equivocado?). "
                  "No te fíes de L_source/R_target hasta arreglarlo.")
    return out


def orden_columnas_liana(res: pd.DataFrame) -> pd.DataFrame:
    """
    Reordena para leer: identidad, decisión, especificidad, contexto. Las
    puntuaciones crudas de magnitud y el p de CellPhoneDB se van al final
    porque en una pregunta de citoquinas no aportan (el p satura en 0 y la
    magnitud la domina el partner abundante).
    """
    delante = ["source", "target", "ligand_complex", "receptor_complex",
               "ligand", "receptor", "como",
               "specificity_rank", "magnitude_rank",
               "lr_logfc", "lfc_L", "lfc_R",
               "spec_weight_oe", "spec_weight", "scaled_weight",
               "L_source", "R_target", "L_global", "R_global",
               "L_target", "R_source"]
    detras = ["lr_means", "expr_prod", "lrscore", "cellphone_pvals",
              "recurso"]
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
    max_rank: Optional[float] = None,      # compatibilidad hacia atrás
    rank_col: str = "specificity_rank",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Interacciones que involucran genes concretos, con criterio en DOS ejes.

    En LIANA las columnas 'ligand_complex'/'receptor_complex' traen las
    subunidades unidas por guion bajo (IL6R_IL6ST). Hay que PARTIRLAS y
    comparar elemento a elemento: un `contains("IL6")` capturaría también
    IL6ST, que es otra cosa. Se marca 'como' para que sepas si tu gen es el
    ligando o una subunidad del receptor.

    POR QUÉ LA MAGNITUD NO SE FILTRA POR RANGO
    ------------------------------------------
    `magnitude_rank` es un rango contra la tabla ENTERA, donde compites con
    colágenos, fibronectina y Apoe. Una citoquina pierde ese rango siempre, y
    sin embargo puede estar perfectamente expresada en términos absolutos. Por
    eso el suelo de magnitud es `min_expr_prod`, que está en la escala de los
    datos, y `max_magnitude` viene apagado. Enciéndelo solo si tu pregunta es
    abierta ("¿qué habla más aquí?") en vez de dirigida a unos genes.

    `specificity_rank` sí se filtra por rango: ahí el rango es la métrica.
    Pero tiene un suelo de empate (mira `res[rank_col].value_counts()`), y por
    debajo de él no ordena. De ahí que el orden por defecto sea `lr_logfc`
    descendente, que es la única columna de especificidad que sigue
    discriminando dentro del empate.

    como
        'ambos' | 'ligando' | 'subunidad_receptor'. Filtra por el papel de tus
        genes. 'subunidad_receptor' NO significa descartable: si buscas Il6 y
        sale Il11 sobre Il6ra_Il6st, eso es la familia gp130 usando el mismo
        receptor, que puede ser la respuesta a tu pregunta.
    """
    if max_rank is not None:
        warnings.warn("max_rank está deprecado; lo interpreto como "
                      "max_specificity. Usa max_specificity/min_expr_prod.")
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
    out["como"] = np.where(es_lig[sel],
                           np.where(es_rec[sel], "ambos", "ligando"),
                           "subunidad_receptor")
    antes = len(out)

    pasos = []
    if como is not None:
        out = out[out["como"] == como]
        pasos.append(f"como == '{como}': {len(out)}")
    if rank_col in out and max_specificity is not None:
        out = out[out[rank_col] <= max_specificity]
        pasos.append(f"{rank_col} <= {max_specificity}: {len(out)}")
    elif rank_col not in out:
        warnings.warn(f"No hay columna '{rank_col}'; no filtro especificidad")
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
        print(f"[filtro] no hay '{sort_by}'; ordeno por {rank_col}. "
              f"¿Se te olvidó pasar por annotate_liana()?")
        if rank_col in out:
            out = out.sort_values(rank_col)

    out = orden_columnas_liana(out)

    if verbose:
        print(f"[filtro] {antes} interacciones anotadas con "
              f"{sorted(buscados)}")
        for p in pasos:
            print(f"[filtro]   -> {p}")
        if len(out):
            print(out["como"].value_counts().to_string())
            cols = [c for c in ("source", "target", col_l, col_r, "como",
                                "specificity_rank", "lr_logfc", "lfc_L",
                                "lfc_R", "spec_weight_oe", "L_source",
                                "R_target", "L_global", "R_global")
                    if c in out.columns]
            print(out[cols].head(20).round(3).to_string(index=False))
        else:
            print("[filtro] nada sobrevive. Antes de relajar el umbral mira "
                  "si es por expr_prod (los pares están apagados) o por "
                  "specificity_rank (están encendidos en todas partes): son "
                  "conclusiones distintas.")
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
    La misma pregunta con varios recursos: ¿qué interacciones sobreviven?

    DISTINGUE DOS COSAS QUE SE CONFUNDEN, y me las confundí yo esta mañana con
    los fondos de la GOEA:

      n_recursos_testado  en cuántos recursos está ANOTADO ese par
      n_recursos_sig      en cuántos sale por debajo de max_rank

    Un par ausente de CellTalkDB no es "no significativo", es "no anotado", y
    si no separas las dos cosas el recurso más pequeño hace que todo parezca
    frágil. La fracción que importa es sig/testado, no sig/total.

    Los recursos humanos (cellphonedb, cellchatdb, celltalkdb) se corren tal
    cual sobre símbolos de ratón y NO van a emparejar bien salvo que LIANA
    traduzca por dentro. Si ves coberturas ridículas en un recurso humano, usa
    li.rs.translate_resource() para orto-transferirlo, o quédate con
    'mouseconsensus' y los que LIANA sirva en murino.
    """
    trozos = []
    for r in resources:
        if verbose:
            print(f"\n===== recurso '{r}' =====")
        try:
            res = run_liana(adata, groupby, resource=r,
                            key_added=f"liana_{r}", verbose=verbose, **kwargs)
        except Exception as e:
            warnings.warn(f"{r} falló: {e}")
            continue
        trozos.append(res)

    if not trozos:
        raise RuntimeError("Ningún recurso completó")

    todo = pd.concat(trozos, ignore_index=True)
    col_l = "ligand_complex" if "ligand_complex" in todo else "ligand"
    col_r = "receptor_complex" if "receptor_complex" in todo else "receptor"
    clave = ["source", "target", col_l, col_r]

    todo["_sig"] = todo[rank_col] <= max_rank
    tabla = todo.groupby(clave).agg(
        n_recursos_testado=("recurso", "nunique"),
        n_recursos_sig=("_sig", "sum"),
        mejor_rank=(rank_col, "min"),
    ).reset_index()
    tabla["fraccion"] = (tabla["n_recursos_sig"]
                         / tabla["n_recursos_testado"].clip(lower=1))
    tabla = tabla.sort_values(["n_recursos_sig", "mejor_rank"],
                              ascending=[False, True])

    if verbose:
        n = todo["recurso"].nunique()
        print(f"\n[recursos] {n} recursos corridos")
        print(f"[recursos] pares anotados en los {n}: "
              f"{int((tabla.n_recursos_testado == n).sum())}")
        print(f"[recursos] significativos en todos aquellos en los que están "
              f"anotados: {int((tabla.fraccion == 1).sum())}")
        print("[recursos] Compara SIEMPRE sig/testado, no sig/total: un par "
              "que solo existe en un recurso no es frágil, es exclusivo.")
        print(tabla.head(20).round(4).to_string(index=False))
    return tabla


def plot_liana(adata: ad.AnnData, key: str = "liana_res",
               colour: str = "specificity_rank", size: str = "magnitude_rank",
               source_labels=None, target_labels=None, top_n: int = 25,
               **kwargs):
    """
    dotplot de LIANA, con especificidad en el COLOR por defecto.

    El default de la librería suele destacar magnitud, y la magnitud premia a
    los genes abundantes en todas partes. Si lo que buscas es qué par de tipos
    celulares habla de forma particular, el color tiene que ser especificidad.
    """
    import liana as li

    return li.pl.dotplot(
        adata=adata, uns_key=key, colour=colour, size=size,
        source_labels=source_labels, target_labels=target_labels,
        top_n=top_n, orderby=colour, orderby_ascending=True, **kwargs)


# ---------------------------------------------------------------------------
# 5c. Delta entre condiciones
# ---------------------------------------------------------------------------

def _clave_liana(df: pd.DataFrame) -> pd.MultiIndex:
    col_l = "ligand_complex" if "ligand_complex" in df else "ligand"
    col_r = "receptor_complex" if "receptor_complex" in df else "receptor"
    return pd.MultiIndex.from_arrays(
        [df["source"].astype(str), df["target"].astype(str),
         df[col_l].astype(str), df[col_r].astype(str)],
        names=["source", "target", "ligand_complex", "receptor_complex"])


def _medias_y_props(adata, groupby, complejos, use_raw=False, layer=None):
    """Medias y fracción de células que expresan, por complejo y por grupo."""
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
    # expr_prop de LIANA se aplica a CADA subunidad: la que manda es la peor
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
    Compara dos corridas de LIANA (una por condición) sobre una selección fija.

    LA SELECCIÓN SE PASA, NO SE RECALCULA
    -------------------------------------
    `seleccion` es la tabla que ya filtraste sobre el pooled (la salida de
    `filter_liana`). Aquí no hay umbrales de selección a propósito: si los
    hubiera, tendrías el criterio en dos sitios y el argumento de que la lista
    estaba fija ANTES de mirar el delta dejaría de ser comprobable. Si pasas
    None se comparan todos los pares en común, que sirve para explorar pero no
    para concluir.

    NO SE COMPARAN RANGOS
    ---------------------
    `magnitude_rank` y `specificity_rank` son rangos DENTRO de cada tabla. Dos
    corridas tienen distribuciones de referencia distintas —basta con que en
    una pasen el expr_prop más interacciones para que todos los rangos se
    desplacen— así que restarlos mide el tamaño de la tabla, no la biología.
    Por eso `cols` solo admite columnas absolutas. Si metes un '_rank' avisa.

    AUSENTE NO ES CERO
    ------------------
    Un par puede faltar en una tabla porque el gen no se expresa o porque no
    llegó al `expr_prop`, que es una fracción de células, no una media. Son
    conclusiones opuestas y la tabla las separa en 'motivo':

      no_expresado      la media es 0 en esa condición
      bajo_expr_prop    hay expresión pero en menos del expr_prop de células
      desconocido       no me pasaste los adatas y no puedo saberlo

    Para distinguirlas hace falta `adatas=(adata_a, adata_b)` y `groupby`. Sin
    eso el delta de las filas que solo están en una condición queda en NaN, que
    es lo honesto: no sabes si el cambio es de 0.6 a 0 o de 0.6 a 0.55.

    Y con una réplica por condición esto es DESCRIPTIVO. No hay p-valor porque
    no hay nada que replique: el delta y la variabilidad entre animales son
    indistinguibles.
    """
    na, nb = nombres[0], nombres[1]
    malas = [c for c in cols if c.endswith("_rank")]
    if malas:
        warnings.warn(
            f"{malas} son rangos internos de cada tabla y no son comparables "
            f"entre corridas. Los quito. Usa expr_prod / L_source / R_target / "
            f"lr_logfc / spec_weight_oe.")
        cols = [c for c in cols if c not in malas]

    A = res_a.copy(); A.index = _clave_liana(A)
    B = res_b.copy(); B.index = _clave_liana(B)
    for nom, T in ((na, A), (nb, B)):
        if T.index.duplicated().any():
            raise ValueError(f"{nom} tiene pares duplicados; ¿concatenaste "
                             f"recursos sin separar?")

    usables = [c for c in cols if c in A.columns and c in B.columns]
    faltan = [c for c in cols if c not in usables]
    if faltan:
        warnings.warn(f"{faltan} no están en ambas tablas (¿te falta "
                      f"annotate_liana en alguna?); las omito")
    if not usables:
        raise ValueError("Ninguna columna comparable")

    if seleccion is not None:
        clave = _clave_liana(seleccion).unique()
    else:
        clave = A.index.union(B.index)
        if verbose:
            print("[delta] sin selección: comparo todo. Exploratorio, no "
                  "concluyas de aquí.")

    out = pd.DataFrame(index=clave)
    for c in usables:
        out[f"{c}_{na}"] = A[c].reindex(clave)
        out[f"{c}_{nb}"] = B[c].reindex(clave)

    en_a = clave.isin(A.index)
    en_b = clave.isin(B.index)
    out["estado"] = np.where(en_a & en_b, "ambas",
                    np.where(en_a, f"solo_{na}",
                    np.where(en_b, f"solo_{nb}", "en_ninguna")))
    out["motivo"] = np.where(out["estado"] == "ambas", "", "desconocido")

    # --- rellenar el lado que falta, si me das con qué
    if adatas is not None and groupby is not None:
        comps = pd.unique(np.concatenate([
            clave.get_level_values("ligand_complex").to_numpy(),
            clave.get_level_values("receptor_complex").to_numpy()]))
        lig = clave.get_level_values("ligand_complex").to_numpy()
        rec = clave.get_level_values("receptor_complex").to_numpy()
        src = clave.get_level_values("source").to_numpy()
        tgt = clave.get_level_values("target").to_numpy()

        motivo = out["motivo"].to_numpy(dtype=object)
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
            motivo[falta & cero] = "no_expresado"
            motivo[falta & ~cero & (peor < expr_prop)] = "bajo_expr_prop"
            motivo[falta & ~cero & (peor >= expr_prop)] = "presente_sin_anotar"
        out["motivo"] = motivo

    for c in usables:
        out[f"d_{c}"] = out[f"{c}_{nb}"] - out[f"{c}_{na}"]

    ref = f"d_{usables[0]}"
    out = out.sort_values(ref, key=lambda s: s.abs(), ascending=False)
    if min_abs_delta is not None:
        out = out[out[ref].abs() >= min_abs_delta]

    orden = (["estado", "motivo"]
             + [f"d_{c}" for c in usables]
             + [f"{c}_{n}" for c in usables for n in (na, nb)])
    out = out[orden].reset_index()

    if verbose:
        print(f"[delta] {len(clave)} pares comparados "
              f"({nb} menos {na}), ordenados por |{ref}|")
        print(out["estado"].value_counts().to_string())
        m = out.loc[out["motivo"] != "", "motivo"].value_counts()
        if len(m):
            print(m.to_string())
            if (out["motivo"] == "desconocido").any():
                print("[delta] pásame adatas=(a, b) y groupby para saber si "
                      "los ausentes están a cero o solo por debajo del "
                      "expr_prop; no es lo mismo.")
        print(out.head(top).round(3).to_string(index=False))
        print("[delta] recuerda: una réplica por condición. Esto describe, "
              "no contrasta.")
    return out
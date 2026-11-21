"""
Per-sample QC, prior to concatenation.

Approach: COMPUTE -> INSPECT -> DECIDE -> APPLY, in that order and as
separate steps. No function in this module filters cells on its own or
chooses thresholds for you. Flags are added as boolean columns and the
actual filtering is done at the end, once all samples have been inspected.

    for gsm in samples:
        a = ad.read_h5ad(...)
        compute_qc_metrics(a, mt_contig=MT_CONTIG)   # metrics, no flags
        add_droplet_qc(a, **params[gsm])             # DropletQC, no flags
        flag_doublets(a)                             # per sample, never concatenated
        qc_embedding(a)                              # throwaway UMAP for inspection
        adatas[gsm] = a

    plot_qc_overview(adatas[gsm])                    # choose thresholds by inspection
    nf_band_report(adatas[gsm], band=0.07)           # what lies in the flat band?
    ambient_top_genes(adatas[gsm])                   # what is the soup made of?

    apply_qc_flags(adatas[gsm], **thresholds[gsm])
    qc_summary(adatas)                               # how much each sample loses

WARNING ON SINGLE-CELL vs SINGLE-NUCLEI
---------------------------------------
Several thresholds here have OPPOSITE meanings depending on the preparation:

    metric           single-cell                  single-nuclei
    ---------------------------------------------------------------------
    nuclear_frac     high = damaged cell          high = normal (it is a nucleus)
                     low = debris/soup            -> the "damaged" branch
                                                     of DropletQC does not apply
    pct_counts_mt    high = damaged cell          high = cytoplasmic
                                                     contamination, poorly
                                                     isolated nucleus
    Malat1           low = anucleate debris       always high (it is nuclear)

This module does NOT decide for you: it requires thresholds to be passed
explicitly and prints a reminder of the issue. The decision is yours, sample by sample.

NOTE ON THE NUCLEAR FRACTION
----------------------------
In single-cell data it is the ONLY criterion in this module that separates real
cells from anucleate debris, because it is the only one that measures something
debris cannot have. Neither CellBender's cell_probability nor total counts do so:
CellBender classifies "above the ambient background", not "is a cell", and a
cytoplasmic fragment with thousands of UMIs looks like a cell to it with p ~ 1.
"""

from __future__ import annotations

import os
import warnings
from typing import Callable, Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


PREP_WARNING = (
    "[qc] Reminder: nuclear_frac, pct_counts_mt and Malat1 are interpreted the "
    "opposite way in single-nuclei and in single-cell data. Review the thresholds for THIS "
    "preparation before applying them (see the module docstring)."
)

# Common mitochondrial contigs (RefSeq accession of the mitochondrial genome)
MT_CONTIG_MOUSE_REFSEQ = "NC_005089.1"
MT_CONTIG_HUMAN_REFSEQ = "NC_012920.1"


# --------------------------------------------------------------------------
# 1. metrics (nothing is flagged)
# --------------------------------------------------------------------------

def compute_qc_metrics(
    adata: ad.AnnData,
    mt_contig: Optional[str] = None,
    mt_prefix: str = "mt-",
    malat_gene: str = "Malat1",
    counts_layer_for_umi: str = "counts_raw",
    verbose: bool = True,
) -> ad.AnnData:
    """
    Add all QC metrics to obs. No flag is created.

    mt_contig
        Accession of the mitochondrial contig, e.g. 'NC_005089.1' in mouse RefSeq.
        This is the preferred route: selects by CONTIG using var['seqname'], which
        comes from annotate_var_from_gtf(). Independent of gene nomenclature.
    mt_prefix
        Name-based fallback when there is no 'seqname' or mt_contig is not given. NOTE:
        the 'mt-' prefix is an Ensembl/GENCODE convention. In NCBI RefSeq the
        mitochondrial genes are named ND1..ND6, CYTB, COX1..3, ATP6/8, TrnX, Rnr1/2
        and this prefix matches NONE of them.
    counts_layer_for_umi
        Layer used for the UMI axis of DropletQC. Defaults to 'counts_raw':
        the thresholds in the DropletQC paper (umi_rescue=3000 and similar) are
        defined on raw counts, not on decontaminated counts.
    """
    # ---- mitochondrial genes: by contig if possible, otherwise by name
    if mt_contig is not None and "seqname" in adata.var.columns:
        adata.var["mt"] = (adata.var["seqname"] == mt_contig).to_numpy()
        mt_how = f"contig {mt_contig}"
    else:
        if mt_contig is not None:
            warnings.warn(
                f"mt_contig='{mt_contig}' was requested but there is no var['seqname']. "
                f"Did you run annotate_var_from_gtf()? Using the prefix '{mt_prefix}'."
            )
        elif "seqname" in adata.var.columns:
            warnings.warn(
                f"var['seqname'] is present but mt_contig was not given: using the prefix "
                f"'{mt_prefix}', which fails with NCBI RefSeq annotations. "
                f"Consider mt_contig='{MT_CONTIG_MOUSE_REFSEQ}' (mouse)."
            )
        adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
        mt_how = f"prefix '{mt_prefix}'"

    n_mt = int(adata.var["mt"].sum())
    if n_mt == 0:
        warnings.warn(
            f"0 mitochondrial genes selected by {mt_how}. pct_counts_mt "
            f"will be 0 for every cell and the criterion will have no effect."
        )
    elif n_mt > 100:
        warnings.warn(f"{n_mt} genes flagged as mitochondrial; this looks like too many")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True
    )

    # counts per layer, for diagnostics and plots
    for layer, name in (("counts_raw", "n_counts_raw"),
                        ("counts_cellbender", "n_counts_cellbender")):
        if layer in adata.layers:
            adata.obs[name] = np.asarray(adata.layers[layer].sum(axis=1)).ravel()

    if "n_counts_raw" in adata.obs and "n_counts_cellbender" in adata.obs:
        raw = adata.obs["n_counts_raw"].to_numpy()
        cb = adata.obs["n_counts_cellbender"].to_numpy()
        with np.errstate(invalid="ignore", divide="ignore"):
            adata.obs["cellbender_removed_frac"] = np.where(raw > 0, 1 - cb / raw, np.nan)

    # ---- velocyto: nuclear fraction
    if all(k in adata.layers for k in ("spliced", "unspliced")):
        for layer in ("spliced", "unspliced", "ambiguous"):
            if layer in adata.layers:
                v = np.asarray(adata.layers[layer].sum(axis=1)).ravel()
                adata.obs[f"n_{layer}"] = v
                adata.obs[f"n_{layer}_log1p"] = np.log1p(v)

        s = adata.obs["n_spliced"].to_numpy()
        u = adata.obs["n_unspliced"].to_numpy()
        denom = s + u
        with np.errstate(invalid="ignore", divide="ignore"):
            adata.obs["nuclear_frac"] = np.where(denom > 0, u / denom, np.nan)

        # Velocyto coverage relative to GeneFull: if it drops sharply in a subpopulation,
        # suspect the alignment before the biology
        if "n_counts_raw" in adata.obs:
            velo = s + u + adata.obs.get("n_ambiguous", 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                adata.obs["velocyto_coverage"] = np.where(
                    adata.obs["n_counts_raw"] > 0, velo / adata.obs["n_counts_raw"], np.nan
                )
    else:
        warnings.warn("No spliced/unspliced layers: no nuclear_frac and no DropletQC")

    # ---- UMIs for the DropletQC axis
    if counts_layer_for_umi in adata.layers:
        umi = np.asarray(adata.layers[counts_layer_for_umi].sum(axis=1)).ravel()
        src = counts_layer_for_umi
    else:
        umi = np.asarray(adata.X.sum(axis=1)).ravel()
        src = "X"
        warnings.warn(
            f"Layer '{counts_layer_for_umi}' does not exist; using X for the UMI "
            f"axis. If X is CellBender output, the DropletQC thresholds do not mean "
            f"the same as in the paper."
        )
    adata.obs["umi_dropletqc"] = umi
    adata.obs["log10_umi"] = np.log10(umi + 1.0)
    adata.uns.setdefault("qc", {})["umi_source"] = src

    # ---- Malat1
    if malat_gene in adata.var_names:
        x = adata[:, malat_gene].X
        x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
        adata.obs["malat1_log1p"] = np.log1p(x)
    else:
        warnings.warn(f"'{malat_gene}' is not in var_names; no Malat1 metric")

    if verbose:
        print(PREP_WARNING)
        print(f"[qc] {adata.n_obs} droplets | UMIs from '{src}' | "
              f"{n_mt} mt genes by {mt_how}")
    return adata


# --------------------------------------------------------------------------
# 2. DropletQC (nothing is flagged: only the classification is stored)
# --------------------------------------------------------------------------

def add_droplet_qc(
    adata: ad.AnnData,
    *,
    nf_rescue: float,
    umi_rescue: float,
    nf_damaged_threshold: float,
    method: str = "kmeans",
    prep: str = "cell",
    classifier: Optional[Callable] = None,
    random_state: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Run classify_empty_and_damaged on (nuclear_frac, umi) and store the
    result in obs['droplet_class']. Does NOT flag cells for removal.

    All thresholds are required and per sample: there are no defaults.

    umi_rescue
        Rescues as a cell anything with low nf but many UMIs. Caution:
        larger debris (e.g. sarcoplasm fragments in muscle)
        has MORE UMIs than a real cell, so this rule may rescue
        precisely the worst droplets. Use nf_band_report() to see what it is rescuing
        before relying on the paper's default value.
    prep
        'cell' or 'nucleus'. Does NOT change the computation — only the warning printed.
    """
    if prep not in ("cell", "nucleus"):
        raise ValueError("prep must be 'cell' or 'nucleus'")
    if "nuclear_frac" not in adata.obs:
        raise ValueError("Missing nuclear_frac; run compute_qc_metrics() first")

    if classifier is None:
        from kranocyte.dropletQC import classify_empty_and_damaged as classifier

    nf_df = pd.DataFrame({
        "nuclear_RNA_frac": adata.obs["nuclear_frac"].to_numpy(),
        "umi": adata.obs["umi_dropletqc"].to_numpy(),
        "logumi": adata.obs["log10_umi"].to_numpy(),
    }, index=adata.obs_names)

    classifier(
        nf_df, method_cells=method, nf_rescue=nf_rescue, umi_rescue=umi_rescue,
        nf_damaged_threshold=nf_damaged_threshold, random_state=random_state,
    )

    col = f"classification_{method}"
    if col not in nf_df.columns:
        raise ValueError(
            f"The classifier did not create '{col}'. Columns: {list(nf_df.columns)}"
        )

    adata.obs["droplet_class"] = pd.Categorical(nf_df[col].values)
    adata.uns.setdefault("qc", {})["dropletqc"] = {
        "method": method, "nf_rescue": nf_rescue, "umi_rescue": umi_rescue,
        "nf_damaged_threshold": nf_damaged_threshold, "prep": prep,
    }

    if verbose:
        print(f"[dropletqc] {dict(adata.obs['droplet_class'].value_counts())}")
        if prep == "nucleus":
            print(
                "[dropletqc] PREPARATION = NUCLEUS: the 'damaged_cell' class has "
                "been computed anyway but is NOT interpretable here — in "
                "single-nuclei data the nuclear fraction is high for everything. Decide "
                "manually whether to ignore it (use high pct_counts_mt as a proxy for "
                "a poorly isolated nucleus) before excluding anything."
            )
    return adata


def nf_band_report(
    adata: ad.AnnData,
    band: float,
    cols: Optional[Sequence[str]] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    What lies inside the flat band of nuclear_frac, and what umi_rescue is
    rescuing.

    Creates obs['nf_group'] with three groups:
        celula      outside the band, classified as a cell
        debris      inside the band, classified as non-cell
        rescatado   inside the band but classified as a CELL
                    (this is what umi_rescue rescues)
        otro        outside the band but classified as non-cell

    How to read the table of medians: if 'rescued' is identical to 'debris' in
    everything except UMIs, umi_rescue is splitting a homogeneous population on
    the only variable that does not discriminate -> raise it or disable it. If 'rescued'
    resembles 'cell' in Malat1 and complexity, it is a real cell type with
    low intronic content and the rescue is doing its job.
    """
    for req in ("nuclear_frac", "droplet_class"):
        if req not in adata.obs:
            raise ValueError(f"Missing obs['{req}']")

    inband = adata.obs["nuclear_frac"] < band
    iscell = adata.obs["droplet_class"] == "cell"

    if verbose:
        print(pd.crosstab(inband, adata.obs["droplet_class"],
                          rownames=[f"nuclear_frac < {band}"],
                          colnames=["droplet_class"]))

    adata.obs["nf_group"] = pd.Categorical(np.select(
        [inband & iscell, inband & ~iscell, ~inband & iscell],
        ["rescued", "debris", "cell"], default="other"))

    if cols is None:
        cols = ["umi_dropletqc", "nuclear_frac", "n_unspliced_log1p",
                "malat1_log1p", "log1p_n_genes_by_counts", "pct_counts_mt",
                "velocyto_coverage", "cell_probability"]
    cols = [c for c in cols if c in adata.obs.columns]

    out = adata.obs.groupby("nf_group", observed=True)[cols].median().round(3)
    out.insert(0, "n", adata.obs["nf_group"].value_counts())

    if verbose:
        print()
        print(out.to_string())
        n_resc = int((adata.obs["nf_group"] == "rescued").sum())
        if n_resc:
            print(f"\n[nf_band] {n_resc} droplets rescued by umi_rescue. "
                  f"Compare their row with 'debris' before accepting them.")
    return out


def ambient_top_genes(adata: ad.AnnData, n: int = 30, verbose: bool = True) -> pd.DataFrame:
    """
    What the soup is made of, according to CellBender.

    This describes the TISSUE, not the sample: in muscle, expect
    fibre transcripts (Ckm, Acta1, troponins, glycolytic genes) because
    fibres are huge syncytia that are torn apart during dissociation. Genes
    that rank high here are suspect in two places:
      - annotation: a cluster with this signature may be soup, not biology
      - velocity: the soup is 100% spliced, so its U/S is biased
        downwards and yields spurious negative velocity -> remove them from velocity_genes
    """
    if "ambient_expression" not in adata.var.columns:
        raise ValueError("No var['ambient_expression']; did you load the CellBender output?")

    cols = ["ambient_expression"]
    for extra in ("ambient_fraction", "cellbender_analyzed", "gene_biotype"):
        if extra in adata.var.columns:
            cols.append(extra)

    out = adata.var[cols].sort_values("ambient_expression", ascending=False).head(n)
    if verbose:
        print(out.to_string())
    return out


# --------------------------------------------------------------------------
# 3. doublets (per sample, never on the concatenated object)
# --------------------------------------------------------------------------

def flag_doublets(
    adata: ad.AnnData,
    use_doubletdetection: bool = True,
    use_scrublet: bool = True,
    seed: int = 0,
    dd_n_iters: int = 10,
    dd_p_thresh: float = 1e-16,
    dd_voter_thresh: float = 0.5,
    dd_n_jobs: int = 1,
    expected_doublet_rate: float = 0.05,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Doublet scores. Runs on ONE sample.

    A doublet forms when two cells enter the same droplet in the same
    10x run: cells from different samples cannot form one. Running
    this on the concatenated object simulates hybrids that do not exist.

    Neither detector needs external clusters: DoubletDetection
    clusters internally at each iteration, and Scrublet builds its own
    KNN graph. Scrublet's 'batch_key' is for SAMPLES; passing clusters to it
    leaves it nearly blind, because the doublets that matter mix cell types.

    Tip: run it AFTER inspecting the debris. A dataset with 40% anucleate
    fragments gives the detectors an expression space that is not
    that of the cells, and the scores become less informative.
    """
    if use_doubletdetection:
        import doubletdetection
        # dd_n_jobs=1 by default: DoubletDetection clusters internally
        # (louvain/phenograph) and with several processes the result can vary
        # between runs even if random_state is passed. If the number of doublets
        # changes from one run to the next, the cells passing QC change and
        # with them EVERYTHING downstream, UMAP included. Increase it only if you have
        # verified that your version is reproducible in parallel.
        np.random.seed(seed)
        clf = doubletdetection.BoostClassifier(
            n_iters=dd_n_iters, clustering_algorithm="louvain",
            standard_scaling=True, pseudocount=0.1, n_jobs=dd_n_jobs,
            random_state=seed,
        )
        calls = clf.fit(adata.X).predict(p_thresh=dd_p_thresh, voter_thresh=dd_voter_thresh)
        score = clf.doublet_score()
        adata.obs["dd_score"] = pd.to_numeric(
            pd.Series(np.asarray(getattr(score, "data", score)).ravel()),
            errors="coerce").to_numpy(dtype=float)
        llamadas = pd.to_numeric(pd.Series(np.asarray(calls).ravel()),
                                 errors="coerce")
        # EXPLICIT bool dtype: if it stays as object, anndata cannot
        # write the h5ad and casting to bool turns NaN into True.
        adata.obs["dd_doublet"] = (llamadas == 1).fillna(False).to_numpy(dtype=bool)
        adata.uns.setdefault("qc", {})["doubletdetection"] = {
            "n_iters": dd_n_iters, "p_thresh": dd_p_thresh,
            "voter_thresh": dd_voter_thresh, "n_jobs": dd_n_jobs, "seed": seed,
        }
        if verbose:
            print(f"[doublets] DoubletDetection: {int(adata.obs['dd_doublet'].sum())} "
                  f"of {adata.n_obs}")

    if use_scrublet:
        import scanpy.external as sce
        res = sce.pp.scrublet(              # no batch_key: already per sample
            adata.copy(), expected_doublet_rate=expected_doublet_rate,
            random_state=seed, knn_dist_metric="cosine", log_transform=True, copy=True,
        )
        adata.obs["scrublet_score"] = pd.to_numeric(
            pd.Series(res.obs["doublet_score"].values),
            errors="coerce").to_numpy(dtype=float)
        pred = pd.Series(res.obs["predicted_doublet"].values)
        n_na = int(pred.isna().sum())
        if n_na:
            warnings.warn(
                f"Scrublet did not classify {n_na} cells (NaN); they are set to False. "
                f"This usually happens when it cannot find a clear bimodal threshold: "
                f"inspect the scrublet_score histogram before relying on it."
            )
        adata.obs["scrublet_doublet"] = pred.fillna(False).astype(bool).to_numpy()
        adata.uns.setdefault("qc", {})["scrublet"] = dict(
            res.uns.get("scrublet", {}).get("parameters", {})
        )
        del res
        if verbose:
            print(f"[doublets] Scrublet: {int(adata.obs['scrublet_doublet'].sum())} "
                  f"of {adata.n_obs}")
    return adata


# --------------------------------------------------------------------------
# 4. throwaway embedding for inspection
# --------------------------------------------------------------------------

def qc_embedding(
    adata: ad.AnnData,
    resolution: float = 0.5,
    n_comps: int = 15,
    n_neighbors: int = 15,
    seed: int = 0,
) -> ad.AnnData:
    """
    Provisional clustering and UMAP, ONLY to see where the flags fall.

    They are stored as 'qc_leiden' and 'X_umap_qc' on purpose, not as 'X_umap'.
    This keeps them out of the final figures after concatenation, a silent
    and common mistake. All of this is recomputed after concatenation.
    """
    tmp = adata.copy()
    if "counts_cellbender" in tmp.layers:
        tmp.X = tmp.layers["counts_cellbender"].copy()
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)
    sc.pp.pca(tmp, n_comps=min(n_comps, min(tmp.shape) - 1), random_state=seed)
    sc.pp.neighbors(tmp, n_neighbors=n_neighbors, random_state=seed)
    sc.tl.leiden(tmp, resolution=resolution, key_added="qc_leiden")
    sc.tl.umap(tmp, min_dist=0.1, random_state=seed)

    adata.obs["qc_leiden"] = tmp.obs["qc_leiden"].values
    adata.obsm["X_umap_qc"] = tmp.obsm["X_umap"]
    del tmp
    return adata


# --------------------------------------------------------------------------
# 5. apply thresholds (required and explicit)
# --------------------------------------------------------------------------

def apply_qc_flags(
    adata: ad.AnnData,
    min_nuclear_frac: float | None = None,
    max_nuclear_frac: Optional[float] | None = None,
    min_pct_mt: Optional[float] | None = None,
    max_pct_mt: Optional[float] | None = None,
    min_log1p_genes: Optional[float] | None = None,
    max_log1p_genes: Optional[float] | None = None,
    min_log1p_counts: Optional[float] | None = None,
    max_log1p_counts: Optional[float] | None = None,
    min_log10_umi: Optional[float] | None = None,
    max_log10_umi: Optional[float] | None = None,
    min_log1p_malat1: Optional[float] | None = None,
    max_log1p_malat1: Optional[float] | None = None,
    max_cellbender_removed: Optional[float] | None = None,
    min_n_unspliced_log1p: Optional[float] | None = None,
    max_n_unspliced_log1p: Optional[float] | None = None,
    min_n_spliced_log1p: Optional[float] | None = None,
    max_n_spliced_log1p: Optional[float] | None = None,
    droplet_exclude: Sequence[str] = (),
    doublet_cols: Sequence[str] = (),
    verbose: bool = True,
) -> ad.AnnData:
    """
    Create the boolean flags and the 'qc_pass' column. Does NOT remove cells.

    All thresholds are keyword-only and required: pass None to
    disable a criterion, but it must be written out. There are no defaults,
    so that no threshold is applied without having been decided for THIS
    sample and THIS preparation.

    min_nuclear_frac
        In single-cell data this is the main criterion: it separates real cells from
        anucleate debris, and it is the only one that does so. A direct cut here is more
        transparent and reproducible than relying on the internal rules of
        DropletQC (which may rescue large debris via umi_rescue).
        In single-nuclei data it does NOT have the same meaning: everything is nuclear.
    max_nuclear_frac
        Damaged cells in single-cell data (they lost cytoplasm). In single-nuclei data
        leave it as None: it would flag the whole sample.
    min_counts
        Applied to 'umi_dropletqc' (raw counts), not to X.
    max_cellbender_removed
        Fraction of counts that CellBender removed from the cell. A cell that
        loses 40% is mostly soup.
    droplet_exclude
        DropletQC classes to flag, e.g. ('empty_droplet',). Empty to not
        use the criterion (reasonable if already cutting on min_nuclear_frac).
    doublet_cols
        Boolean doublet columns to include, e.g. ('dd_doublet',).
    """
    flags = {}

    def _need(col):
        if col not in adata.obs:
            raise ValueError(f"No obs['{col}']; did you run compute_qc_metrics()?")
        return adata.obs[col].to_numpy()

    if min_nuclear_frac is not None:
        flags["flag_low_nf"] = _need("nuclear_frac") < min_nuclear_frac
    if max_nuclear_frac is not None:
        flags["flag_high_nf"] = _need("nuclear_frac") > max_nuclear_frac
    if max_pct_mt is not None:
        flags["flag_high_mt"] = _need("pct_counts_mt") > max_pct_mt
    if min_pct_mt is not None:
        flags["flag_low_mt"] = _need("pct_counts_mt") < min_pct_mt
    if min_log1p_genes is not None:
        flags["flag_low_genes"] = _need("log1p_n_genes_by_counts") < min_log1p_genes
    if max_log1p_genes is not None:
        flags["flag_high_genes"] = _need("log1p_n_genes_by_counts") > max_log1p_genes
    if min_log1p_counts is not None:
        flags["flag_low_counts"] = _need("log1p_total_counts") < min_log1p_counts
    if max_log1p_counts is not None:
        flags["flag_high_counts"] = _need("log1p_total_counts") > max_log1p_counts
    if min_log10_umi is not None:
        flags["flag_low_umi"] = _need("log10_umi") < min_log10_umi
    if max_log10_umi is not None:
        flags["flag_high_umi"] = _need("log10_umi") > max_log10_umi
    if min_log1p_malat1 is not None:
        flags["flag_low_malat1"] = _need("malat1_log1p") < min_log1p_malat1
    if max_log1p_malat1 is not None:
        flags["flag_high_malat1"] = _need("malat1_log1p") > max_log1p_malat1

    if min_n_unspliced_log1p is not None:
        flags["flag_min_n_unspliced_log1p"] = _need("log10_umi") < min_n_unspliced_log1p
    if max_n_unspliced_log1p is not None:
        flags["flag_max_n_unspliced_log1p"] = _need("log10_umi") > max_n_unspliced_log1p
    if min_n_spliced_log1p is not None:
        flags["flag_min_n_spliced_log1p"] = _need("malat1_log1p") < min_n_spliced_log1p
    if max_n_spliced_log1p is not None:
        flags["flag_max_n_spliced_log1p"] = _need("malat1_log1p") > max_n_spliced_log1p

    if max_cellbender_removed is not None:
        flags["flag_high_ambient"] = _need("cellbender_removed_frac") > max_cellbender_removed
    if len(droplet_exclude):
        if "droplet_class" not in adata.obs:
            raise ValueError("No 'droplet_class'; did you run add_droplet_qc?")
        flags["flag_droplet"] = adata.obs["droplet_class"].isin(droplet_exclude).to_numpy()
    for col in doublet_cols:
        if col not in adata.obs:
            raise ValueError(f"No '{col}' in obs")
        flags[f"flag_{col}"] = adata.obs[col].to_numpy().astype(bool)

    checks = [
        ("flag_low_nf",       "nuclear_frac",            min_nuclear_frac,   "lt"),
        ("flag_high_nf",      "nuclear_frac",            max_nuclear_frac,   "gt"),
        ("flag_low_mt",       "pct_counts_mt",           min_pct_mt,         "lt"),
        ("flag_high_mt",      "pct_counts_mt",           max_pct_mt,         "gt"),
        ("flag_low_genes",    "log1p_n_genes_by_counts", min_log1p_genes,    "lt"),
        ("flag_high_genes",   "log1p_n_genes_by_counts", max_log1p_genes,    "gt"),
        ("flag_low_counts",   "log1p_total_counts",      min_log1p_counts,   "lt"),
        ("flag_high_counts",  "log1p_total_counts",      max_log1p_counts,   "gt"),
        ("flag_low_umi",      "log10_umi",               min_log10_umi,      "lt"),
        ("flag_high_umi",     "log10_umi",               max_log10_umi,      "gt"),
        ("flag_low_malat1",   "malat1_log1p",            min_log1p_malat1,   "lt"),
        ("flag_high_malat1",  "malat1_log1p",            max_log1p_malat1,   "gt"),

        ("flag_min_n_unspliced_log1p",   "n_unspliced_log1p",            min_n_unspliced_log1p,   "lt"),
        ("flag_max_n_unspliced_log1p",  "n_unspliced_log1p",            max_n_unspliced_log1p,   "gt"),
        ("flag_min_n_spliced_log1p",   "n_spliced_log1p",            min_n_spliced_log1p,   "lt"),
        ("flag_max_n_spliced_log1p",  "n_spliced_log1p",            max_n_spliced_log1p,   "gt"),

        ("flag_high_ambient", "cellbender_removed_frac", max_cellbender_removed, "gt"),
    ]
    for name, col, thr_val, op in checks:
        if thr_val is None:
            continue
        v = _need(col)
        flags[name] = v < thr_val if op == "lt" else v > thr_val


    if not flags:
        raise ValueError("No criterion has been enabled")

    fail = np.zeros(adata.n_obs, dtype=bool)
    for k, v in flags.items():
        v = np.asarray(v)
        v = np.where(np.isnan(v.astype(float)), False, v) if v.dtype.kind == "f" else v
        adata.obs[k] = v.astype(bool)
        fail |= adata.obs[k].to_numpy()
    adata.obs["qc_pass"] = ~fail

    adata.uns.setdefault("qc", {})["thresholds"] = {
        "min_nuclear_frac": min_nuclear_frac, "max_nuclear_frac": max_nuclear_frac,
        "min_pct_mt": min_pct_mt, "max_pct_mt": max_pct_mt, 
        "min_log1p_genes": min_log1p_genes, "max_log1p_genes": max_log1p_genes,
        "min_log1p_counts": min_log1p_counts, "max_log1p_counts": max_log1p_counts, 
        "min_log1p_malat1": min_log1p_malat1, "max_log1p_malat1": max_log1p_malat1,
        "max_cellbender_removed": max_cellbender_removed,
        "droplet_exclude": list(droplet_exclude), "doublet_cols": list(doublet_cols),
    }

    if verbose:
        print(f"[qc] {int(adata.obs['qc_pass'].sum())}/{adata.n_obs} pass "
              f"({100 * adata.obs['qc_pass'].mean():.1f}%)")
    return adata


# --------------------------------------------------------------------------
# 6. cross-sample summary
# --------------------------------------------------------------------------

def sanitize_for_h5ad(adata: ad.AnnData, verbose: bool = True) -> ad.AnnData:
    """
    Fix the obs/var columns that prevent writing the .h5ad.

    The symptom is this error, which says nothing useful:

        TypeError: Can't implicitly convert non-string objects to strings
        Error raised while writing key 'scrublet_doublet' of ... to /obs

    The cause is a column of dtype 'object'. anndata can only write
    numeric, boolean, categorical and string columns; for an 'object' it assumes strings
    and h5py fails when it finds booleans or NaN inside. This happens
    easily when concatenating samples: if the doublet detector returned NaN in one
    of them, or if a column is missing in one sample and present in another, concatenation
    promotes it to object without warning.

    And it is not only a problem when saving. An 'object' with NaN slips into
    filters: astype(bool) turns NaN into True, because bool(nan) is True.

    Converts in place and returns the same object.
    """
    cambios = []
    for donde, df in (("obs", adata.obs), ("var", adata.var)):
        for col in df.columns:
            if df[col].dtype != object:
                continue
            v = df[col]
            n_na = int(v.isna().sum())
            vals = set(map(type, v.dropna().unique().tolist()))
            if vals <= {bool, np.bool_}:
                df[col] = np.array([bool(x) if pd.notna(x) else False
                                    for x in v], dtype=bool)
                cambios.append((donde, col, "object -> bool", n_na))
            elif vals <= {str}:
                df[col] = v.fillna("").astype(str)
                cambios.append((donde, col, "object -> str", n_na))
            elif all(np.issubdtype(np.dtype(t), np.number) for t in vals) and vals:
                df[col] = pd.to_numeric(v, errors="coerce")
                cambios.append((donde, col, "object -> numeric", n_na))
            else:
                df[col] = v.astype(str)
                cambios.append((donde, col, f"mixed {vals} -> str", n_na))

    if verbose:
        if not cambios:
            print("[sanitize] nothing to fix")
        else:
            for donde, col, que, n_na in cambios:
                extra = f", {n_na} NaN" if n_na else ""
                print(f"[sanitize] {donde}['{col}']: {que}{extra}")
            if any(x[3] for x in cambios):
                print("[sanitize] NaN values in boolean columns were set to False. "
                      "If there were many, check why the detector did not "
                      "classify those cells BEFORE accepting the filtering.")
    return adata


def qc_summary(
    adatas: Mapping[str, ad.AnnData],
    spread_warning: float = 15.0,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Per-sample table: how many cells each criterion flags, how many it flags
    EXCLUSIVELY, and how many survive.

    The '_solo' column is the informative one: if a criterion flags 2000 cells
    but only 30 exclusively, it is redundant with the others and its threshold barely
    matters. If it flags 2000 and 1800 exclusively, that threshold is shaping your
    dataset on its own.

    spread_warning
        Difference in pct_pass, in percentage points, between the best and the worst
        sample above which a warning is issued. In a time course, very uneven
        retention means condition-biased loss: either the sample is
        poor, or the threshold is measuring biology rather than quality.
    """
    rows = []
    for name, a in adatas.items():
        flag_cols = [c for c in a.obs.columns if c.startswith("flag_")]
        row = {"sample": name, "n_cells": a.n_obs}
        for c in flag_cols:
            v = a.obs[c].to_numpy().astype(bool)
            others = np.zeros(a.n_obs, dtype=bool)
            for o in flag_cols:
                if o != c:
                    others |= a.obs[o].to_numpy().astype(bool)
            short = c.replace("flag_", "")
            row[short] = int(v.sum())
            row[short + "_solo"] = int((v & ~others).sum())
        if "qc_pass" in a.obs:
            row["n_pass"] = int(a.obs["qc_pass"].sum())
            row["pct_pass"] = round(100 * float(a.obs["qc_pass"].mean()), 1)
        rows.append(row)

    df = pd.DataFrame(rows).set_index("sample")

    if verbose and "pct_pass" in df.columns and len(df) > 1:
        spread = float(df["pct_pass"].max() - df["pct_pass"].min())
        print(f"[qc_summary] retention between {df['pct_pass'].min()}% "
              f"({df['pct_pass'].idxmin()}) and {df['pct_pass'].max()}% "
              f"({df['pct_pass'].idxmax()}) — spread {spread:.1f} points")
        if spread > spread_warning:
            print(
                f"[qc_summary] WARNING: spread > {spread_warning} points. If the "
                f"samples are different conditions, check whether a threshold is "
                f"measuring biology. Inspect the '_solo' column to see which one."
            )
    return df


# --------------------------------------------------------------------------
# 7. inspection figure
# --------------------------------------------------------------------------

def plot_qc_overview(adata: ad.AnnData, title: str = "", figsize=(13, 12)):
    """
    The nine panels used to choose thresholds for one sample. Draws lines
    for the thresholds already present in uns['qc']['thresholds'].

    Does not modify adata: everything is computed on a copy of obs.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    thr = adata.uns.get("qc", {}).get("thresholds", {})
    obs = adata.obs.copy()                      # <- no side effects
    hue = obs["droplet_class"] if "droplet_class" in obs else None

    fig, axs = plt.subplots(3, 3, figsize=figsize, constrained_layout=True)

    def hline(ax, key, **kw):
        if thr.get(key) is not None:
            ax.axhline(thr[key], c="r", ls="--", **kw)

    # --- row 0: mitochondrial and distributions
    sns.scatterplot(x=obs["log1p_total_counts"], y=obs["pct_counts_mt"],
                    hue=hue, s=3, alpha=.4, ax=axs[0, 0], legend=False)
    hline(axs[0, 0], "max_pct_mt")
    axs[0, 0].set(title="mitochondrial")

    sns.violinplot(obs, x="droplet_class", y="pct_counts_mt", ax=axs[0, 1])
    hline(axs[0, 1], "max_pct_mt")
    axs[0, 1].set(title="mitochondrial by class")

    sns.violinplot(obs, x="droplet_class", y="log1p_n_genes_by_counts", ax=axs[0, 2])
    hline(axs[0, 2], "min_log1p_genes")          # <- was max_pct_mt (bug)
    axs[0, 2].set(title="complexity by class")

    # --- row 1: DropletQC, complexity, Malat1
    sns.scatterplot(x=obs["log10_umi"], y=obs["nuclear_frac"], hue=hue,
                    s=3, alpha=.4, ax=axs[1, 0], legend="brief")
    hline(axs[1, 0], "min_nuclear_frac")
    hline(axs[1, 0], "max_nuclear_frac")
    axs[1, 0].set(xlabel="log10 UMI", ylabel="nuclear fraction", title="DropletQC")

    sns.scatterplot(x=obs["log1p_total_counts"], y=obs["log1p_n_genes_by_counts"],
                    hue=hue, s=3, alpha=.4, ax=axs[1, 1], legend=False)
    hline(axs[1, 1], "min_log1p_genes")
    axs[1, 1].set(title="complexity")

    if "malat1_log1p" in obs:
        sns.scatterplot(x=obs["nuclear_frac"], y=obs["malat1_log1p"], hue=hue,
                        s=3, alpha=.4, ax=axs[1, 2], legend=False)
        hline(axs[1, 2], "min_malat1")
        axs[1, 2].set(title="Malat1")

    # --- row 2: ambient, doublets, spliced/unspliced
    if "cellbender_removed_frac" in obs:
        sns.kdeplot(obs, x="cellbender_removed_frac", hue=hue, ax=axs[2, 0])
        if thr.get("max_cellbender_removed") is not None:
            axs[2, 0].axvline(thr["max_cellbender_removed"], c="r", ls="--")
        axs[2, 0].set(title="fraction removed by CellBender")

    score = next((c for c in ("scrublet_score", "dd_score") if c in obs), None)
    if score:
        sns.scatterplot(x=obs["log1p_n_genes_by_counts"], y=obs[score], hue=hue,
                        s=3, alpha=.4, ax=axs[2, 1], legend=False)
        axs[2, 1].set(title=score)

    if "n_spliced_log1p" in obs and "n_unspliced_log1p" in obs:
        sns.scatterplot(x=obs["n_spliced_log1p"], y=obs["n_unspliced_log1p"],
                        hue=hue, s=3, alpha=.4, ax=axs[2, 2], legend=False)
        axs[2, 2].set(title="spliced vs unspliced")

    fig.suptitle(title or adata.uns.get("qc", {}).get("sample", ""))
    return fig, axs


# --------------------------------------------------------------------------
# 8. debris filtering and assembly
# --------------------------------------------------------------------------

def filter_debris(
    adata: ad.AnnData,
    min_nuclear_frac: float,
    sample: Optional[str] = None,
    report_cols: Optional[Sequence[str]] = None,
    verbose: bool = True,
) -> ad.AnnData:
    """
    The ONLY hard filter in the module, by design: anucleate debris contributes
    nothing downstream and distorts every later step (doublet
    detection, HVGs, PCA, neighbours). Everything else remains a flag.

    Never removes silently: prints how much is removed and the median profile of what
    is removed, so you can check that it is what you think it is. Records
    this in uns['qc']['debris_filter'].

    In single-nuclei data this filter does NOT apply as is: the nuclear fraction is
    high for everything and a low minimum removes nothing. There the proxy is different
    (low counts, high mitochondrial fraction).
    """
    if "nuclear_frac" not in adata.obs:
        raise ValueError("Missing nuclear_frac; run compute_qc_metrics() first")

    nf = adata.obs["nuclear_frac"].to_numpy()
    keep = np.where(np.isnan(nf), False, nf >= min_nuclear_frac)
    n_before, n_keep = adata.n_obs, int(keep.sum())
    n_nan = int(np.isnan(nf).sum())

    if n_keep == 0:
        raise ValueError(
            f"min_nuclear_frac={min_nuclear_frac} removes ALL cells. "
            f"Is this a single-nuclei sample, or is the threshold inverted?"
        )

    if report_cols is None:
        report_cols = ["umi_dropletqc", "nuclear_frac", "malat1_log1p",
                       "log1p_n_genes_by_counts", "pct_counts_mt",
                       "n_unspliced_log1p", "cell_probability"]
    report_cols = [c for c in report_cols if c in adata.obs.columns]

    if verbose:
        tag = f"{sample}: " if sample else ""
        pct = 100 * n_keep / n_before
        print(f"[filter] {tag}{n_before} -> {n_keep} droplets ({pct:.1f}% retained); "
              f"{n_before - n_keep} removed by nuclear_frac < {min_nuclear_frac}"
              + (f" ({n_nan} with undefined nuclear_frac)" if n_nan else ""))
        prof = pd.DataFrame({
            "removed": adata.obs.loc[~keep, report_cols].median(),
            "kept": adata.obs.loc[keep, report_cols].median(),
        }).round(3)
        print(prof.to_string())
        if pct < 40:
            print(f"[filter] WARNING: more than 60% of droplets are removed. Check "
                  f"that the profile of the removed droplets is really debris.")

    out = adata[keep].copy()
    out.uns.setdefault("qc", {})["debris_filter"] = {
        "min_nuclear_frac": min_nuclear_frac, "n_before": n_before,
        "n_after": n_keep, "n_removed": n_before - n_keep,
        "n_nan_nuclear_frac": n_nan,
    }
    return out


def collect_ambient_profiles(
    adatas: Mapping[str, ad.AnnData],
    key: str = "ambient_expression",
) -> pd.DataFrame:
    """
    Per-sample soup profiles in a genes x samples DataFrame.

    They must be extracted BEFORE concatenating: they differ between samples, so
    merge='same' discards them. And they are worth keeping: comparing the soup
    across timepoints tells you how much tissue was destroyed in each, which
    in an injury model is information, not noise.
    """
    cols = {}
    for name, a in adatas.items():
        if key not in a.var.columns:
            warnings.warn(f"{name}: has no var['{key}']")
            continue
        cols[name] = pd.Series(a.var[key].to_numpy(dtype=float), index=a.var_names)
    if not cols:
        raise ValueError(f"No sample has var['{key}']")
    return pd.DataFrame(cols)


def compare_ambient_profiles(
    profiles: pd.DataFrame,
    method: str = "spearman",
    top: int = 10,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Correlation of the soup profile between samples, and the genes that dominate it.

    If all samples correlate highly, the soup is the same destroyed tissue
    in different amounts. If one stands apart, that sample has
    a different contamination composition — which is worth knowing before
    interpreting expression differences between conditions.
    """
    corr = profiles.corr(method=method)
    if verbose:
        print(f"correlation ({method}) of the ambient profile between samples:")
        print(corr.round(3).to_string())
        print(f"\ntop {top} soup genes per sample:")
        print(pd.DataFrame({
            s: profiles[s].sort_values(ascending=False).head(top).index
            for s in profiles.columns
        }).to_string(index=False))
        off = corr.where(~np.eye(len(corr), dtype=bool))
        if off.min().min() < 0.5:
            # the outlier is the one with the lowest MEAN correlation with the rest, not
            # the one in the lowest pair: correlation is symmetric and a pair implicates two
            worst = off.mean().idxmin()
            print(f"\n[ambient] WARNING: '{worst}' is the sample with the most "
                  f"divergent soup (mean correlation {off.mean().min():.3f} with the rest). "
                  f"Its contamination is not of the same kind as that of the others.")
    return corr


def concat_samples(
    adatas: Mapping[str, ad.AnnData],
    sample_key: str = "gsm",
    index_unique: str = "-",
    drop_obsm: Sequence[str] = ("X_umap_qc",),
    prefix_obs: Sequence[str] = ("qc_leiden",),
    verbose: bool = True,
) -> ad.AnnData:
    """
    Concatenate while preserving var and uns, which ad.concat discards by default.

    Three things it does that are easily forgotten when done by hand:

    merge='same' / uns_merge='same'
        Without them the whole of var is lost (entrez_id, seqname, gene_biotype...) as well as
        uns['annotation'] with the GTF provenance. The ad.concat defaults
        are None: discard.
    drop_obsm
        The X_umap_qc of each sample lives in its own space. Concatenating those
        coordinates produces a meaningless embedding that looks valid. Removed.
    prefix_obs
        Cluster '0' in one sample has nothing to do with '0' in another.
        They are prefixed with the sample name so that they are not merged.
    """
    names = list(adatas.keys())
    var_sets = [tuple(adatas[n].var_names) for n in names]
    if len(set(var_sets)) != 1:
        warnings.warn(
            "The samples do not share the same gene set; with join='outer' "
            "missing genes are filled with zeros. Same GTF for all of them?"
        )

    prepared = []
    for n in names:
        a = adatas[n].copy()
        a.obs[sample_key] = n
        for k in drop_obsm:
            a.obsm.pop(k, None)
        for c in prefix_obs:
            if c in a.obs.columns:
                a.obs[c] = pd.Categorical(f"{n}_" + a.obs[c].astype(str))
        prepared.append(a)

    out = ad.concat(
        prepared, label=sample_key, keys=names, index_unique=index_unique,
        join="outer", merge="same", uns_merge="same",
    )

    # soup profiles do not survive merge='same': they are stored separately
    try:
        out.varm["ambient_expression"] = collect_ambient_profiles(adatas).reindex(
            out.var_names).to_numpy(dtype=np.float32)
        out.uns.setdefault("qc", {})["ambient_samples"] = names
    except ValueError:
        pass

    if verbose:
        print(f"[concat] {out.n_obs} cells x {out.n_vars} genes from {len(names)} samples")
        print(out.obs[sample_key].value_counts().to_string())
        kept = [c for c in ("entrez_id", "seqname", "gene_biotype", "gene_symbol")
                if c in out.var.columns]
        print(f"[concat] var kept: {kept if kept else 'NONE — check merge='}")
        if not kept:
            warnings.warn("No annotation column survived in var")
        print(f"[concat] uns kept: {sorted(out.uns)}")
        if "ambient_expression" in out.varm:
            print(f"[concat] soup profiles in varm['ambient_expression'] "
                  f"({out.varm['ambient_expression'].shape})")
    return out


def resolve_cellbender_h5(
    cb_dir: str,
    fpr: Optional[str] = None,
    filtered: bool = False,
    stem: str = "adata_raw_cellbender",
    strict: bool = True,
) -> str:
    """
    Locate the CellBender output regardless of how it was named.

    CellBender changes the file name pattern depending on how many --fpr values it is given:
        --fpr 0.01         ->  <stem>.h5              and  <stem>_filtered.h5
        --fpr 0.01 0.05    ->  <stem>_FPR_0.01.h5     and  <stem>_FPR_0.01_filtered.h5

    So samples run at different times end up with different
    names, and an os.path.exists()-style guard fails silently on the
    mixture. This function accepts both.

    Beware of '_filtered': it does NOT mean "clean cells", it means
    cell_probability > 0.5. In tissues with much anucleate debris that
    probability is ~1 for the fragments too, so the
    filtered file contains them all the same. The actual filter is filter_debris().
    """
    import glob

    suffix = "_filtered" if filtered else ""
    cands = []
    if fpr is not None:
        cands.append(f"{cb_dir}/{stem}_FPR_{fpr}{suffix}.h5")
    cands.append(f"{cb_dir}/{stem}{suffix}.h5")

    for c in cands:
        if os.path.exists(c):
            return c

    # last resort: any FPR present, so that we can report which one exists
    found = sorted(glob.glob(f"{cb_dir}/{stem}_FPR_*{suffix}.h5"))
    if not strict and found:
        warnings.warn(f"FPR={fpr} not found; using {os.path.basename(found[0])}")
        return found[0]

    raise FileNotFoundError(
        f"CellBender output not found in {cb_dir}.\n"
        f"  tried: {[os.path.basename(c) for c in cands]}\n"
        f"  present: {[os.path.basename(f) for f in sorted(glob.glob(f'{cb_dir}/{stem}*.h5'))]}"
    )


# --------------------------------------------------------------------------
# 9. choose the CellBender FPR with a measurable criterion
# --------------------------------------------------------------------------

def fpr_tradeoff(
    adatas_fpr: Mapping[str, ad.AnnData],
    group_key: str,
    negative: Mapping[str, Sequence[str]],
    positive: Mapping[str, Sequence[str]],
    baseline: Optional[str] = None,
    layer: str = "counts_cellbender",
    verbose: bool = True,
):
    """
    The two curves needed to choose the FPR without falling into circularity.

    There is no "correct" FPR recoverable from the data: it is a decision about
    which error you prefer. And the usual internal metrics (silhouette,
    marker specificity, "how clean the clusters look") reward
    homogeneity, so they improve almost monotonically the more you
    decontaminate — they push you to over-correct. They are circular.

    The only thing with known ground truth is the set of genes a population CANNOT
    express: Ckm in a macrophage is 100% soup, and the correct value is zero.
    This yields two measurable, opposing quantities:

        specificity gained    how much the impossible gene drops where it should not be
        signal lost           how much THAT SAME gene drops where it should be

    The criterion is the balance between the two, never one alone. If raising the FPR
    makes the gene drop by 90% in the negative population and by 5% in the positive one, the
    correction works. If it drops by 60% and 40%, you have not gained specificity:
    you have flattened the whole axis.

    Parameters
    ----------
    adatas_fpr
        {fpr: AnnData}, same cells and genes in all of them. e.g.
        {"0.0": a0, "0.01": a1, "0.05": a5}
    group_key
        obs column with the populations. A coarse annotation is enough: you do not
        need the final one to know what a macrophage is.
    negative
        {population: [genes it CANNOT express]}. The correct value is 0.
    positive
        {population: [genes it MUST express]}. Use the SAME genes as in
        'negative': comparing the same gene in both places is what makes
        the balance meaningful.
    baseline
        Reference FPR. Defaults to the lowest, ideally '0.0'. Without an
        uncorrected reference you cannot measure "how much it has dropped".

    Returns
    -------
    (detalle, resumen)
        detalle : one row per (fpr, tipo, poblacion, gen)
        resumen : one row per fpr, with the aggregated balance
    """
    fprs = sorted(adatas_fpr, key=float)
    if baseline is None:
        baseline = fprs[0]
    if baseline not in adatas_fpr:
        raise ValueError(f"baseline '{baseline}' is not in adatas_fpr")
    if float(baseline) > 0 and verbose:
        print(f"[fpr] WARNING: the reference is {baseline}, not 0.0. Drops are "
              f"measured against an already corrected matrix, so they underestimate "
              f"both the specificity gained and the signal lost.")

    ref = adatas_fpr[baseline]
    for f, a in adatas_fpr.items():
        if not a.obs_names.equals(ref.obs_names):
            raise ValueError(f"FPR {f}: cells do not match the reference")
        if not a.var_names.equals(ref.var_names):
            raise ValueError(f"FPR {f}: genes do not match the reference")
    if group_key not in ref.obs:
        raise ValueError(f"No obs['{group_key}'] in the reference")

    groups = ref.obs[group_key].astype(str).to_numpy()   # same cells: valid for all
    gene_pos = {g: i for i, g in enumerate(ref.var_names)}

    rows = []
    for f in fprs:
        a = adatas_fpr[f]
        X = a.layers[layer] if layer and layer in a.layers else a.X
        X = sp.csc_matrix(X)
        tot = np.asarray(X.sum(axis=1)).ravel()
        tot[tot == 0] = np.nan                       # avoid division by zero

        for kind, spec in (("negative", negative), ("positive", positive)):
            for pop, genes in spec.items():
                mask = groups == pop
                if not mask.any():
                    warnings.warn(f"'{pop}' does not appear in obs['{group_key}']")
                    continue
                for g in genes:
                    if g not in gene_pos:
                        warnings.warn(f"'{g}' is not in var_names")
                        continue
                    col = np.asarray(X[:, gene_pos[g]].todense()).ravel()[mask]
                    cpm = np.nanmean(col / tot[mask] * 1e6)
                    rows.append({
                        "fpr": f, "kind": kind, "population": pop, "gene": g,
                        "n_cells_es": int(mask.sum()),
                        "pct_expressing": round(100 * float((col > 0).mean()), 2),
                        "mean_cpm": round(float(cpm), 2),
                        "mean_counts": round(float(col.mean()), 3),
                    })

    detalle = pd.DataFrame(rows)
    if detalle.empty:
        raise ValueError("No population/gene combination produced any data")

    # drop relative to the reference, per (tipo, poblacion, gen)
    idx = pd.MultiIndex.from_frame(detalle[["kind", "population", "gene"]])
    for col_src, col_dst in (("mean_cpm", "frac_removed"),
                             ("mean_counts", "frac_removed_counts")):
        base = (detalle[detalle["fpr"] == baseline]
                .set_index(["kind", "population", "gene"])[col_src])
        ref_vals = base.reindex(idx).to_numpy()
        with np.errstate(invalid="ignore", divide="ignore"):
            detalle[col_dst] = np.where(
                ref_vals > 0, 1 - detalle[col_src].to_numpy() / ref_vals, np.nan)
        detalle[col_dst] = detalle[col_dst].round(3)

    # summary per fpr
    res = []
    for f in fprs:
        d = detalle[detalle["fpr"] == f]
        neg, pos = d[d["kind"] == "negative"], d[d["kind"] == "positive"]
        gain = float(neg["frac_removed"].mean())
        loss = float(pos["frac_removed"].mean())
        loss_counts = float(pos["frac_removed_counts"].mean())
        res.append({
            "fpr": f,
            "n_neg": len(neg), "n_pos": len(pos),
            "specificity_gained": round(gain, 3),
            "signal_lost": round(loss, 3),
            "balance": round(gain - loss, 3),
            # on raw counts: per-cell normalisation can mask
            # uniform removal, so this column is the control
            "signal_lost_counts": round(loss_counts, 3),
            "specificity_gained_counts": round(float(neg["frac_removed_counts"].mean()), 3),
            "neg_pct_expressing": round(float(neg["pct_expressing"].mean()), 2),
            "neg_cpm_residual": round(float(neg["mean_cpm"].mean()), 2),
        })
    resumen = pd.DataFrame(res).set_index("fpr")

    if verbose:
        print(resumen.to_string())
        print()
        cand = resumen.drop(index=baseline, errors="ignore")
        if len(cand):
            best = cand["balance"].idxmax()
            print(f"[fpr] best balance: {best} "
                  f"(specificity +{cand.loc[best, 'specificity_gained']:.2f}, "
                  f"signal -{cand.loc[best, 'signal_lost']:.2f})")
            if cand["balance"].max() <= 0:
                print("[fpr] No FPR gains more specificity than the signal it "
                      "loses: keep the most conservative one.")
            elif cand.loc[best, "signal_lost_counts"] > 0.25:
                print(f"[fpr] WARNING: on raw counts "
                      f"{100*cand.loc[best,'signal_lost_counts']:.0f}% of the real signal "
                      f"of the positive population is lost. In CPM it is barely visible "
                      f"({100*cand.loc[best,'signal_lost']:.0f}%) because "
                      f"per-cell normalisation compensates for uniform "
                      f"removal. For velocity and quantitative kinetics, which "
                      f"work on counts, this loss DOES matter.")
        print("\n[fpr] Reminder: this ranks the FPR values, it does not dictate the value. If your "
              "biological conclusion does not change between them, the choice does not "
              "matter and that should be stated. If it changes, the result depends on a "
              "nuisance parameter and that is itself the finding.")
    return detalle, resumen


def plot_fpr_tradeoff(resumen: pd.DataFrame, figsize=(10, 4)):
    """The two curves and the balance, as a function of the FPR."""
    import matplotlib.pyplot as plt

    x = [float(f) for f in resumen.index]
    fig, axs = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)

    axs[0].plot(x, resumen["specificity_gained"], "o-", label="specificity gained")
    axs[0].plot(x, resumen["signal_lost"], "s-", label="signal lost")
    axs[0].set(xlabel="FPR", ylabel="fraction removed vs reference",
               title="the two curves")
    axs[0].legend(); axs[0].grid(alpha=.3)

    axs[1].plot(x, resumen["balance"], "o-", c="k")
    axs[1].axhline(0, c="r", ls="--", lw=1)
    axs[1].set(xlabel="FPR", ylabel="gained - lost", title="balance")
    axs[1].grid(alpha=.3)
    return fig, axs


# --------------------------------------------------------------------------
# 10. per-sample thresholds from numbers, not from figures
# --------------------------------------------------------------------------
#
# plot_qc_overview is useful for seeing the STRUCTURE of a sample. It is not suitable for
# setting thresholds for several samples at once: densities are scaled
# per panel, axes change between figures, and what looks "a bit further to the
# left" can be 0.05 or 0.4 of nuclear_frac. These functions report the
# same data as numbers so that the threshold is reproducible and recorded.


def _kde_modes(x: np.ndarray, n_grid: int = 512, bw: str = "scott"):
    """
    Modes and valleys of a 1D distribution by KDE.

    Returns (grid, dens, modas, valles) with modes/valleys sorted by
    position. Used to answer two specific questions:
      - are there two nuclear_frac clouds or one? (bimodality)
      - where is the minimum between them? (that is where the threshold goes, not at a
        round number chosen by eye)
    """
    from scipy.stats import gaussian_kde

    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 50 or np.allclose(x, x[0]):
        return None, None, np.array([]), np.array([])

    lo, hi = np.nanmin(x), np.nanmax(x)
    pad = 0.02 * (hi - lo)
    grid = np.linspace(lo - pad, hi + pad, n_grid)
    try:
        dens = gaussian_kde(x, bw_method=bw)(grid)
    except Exception:
        return None, None, np.array([]), np.array([])

    up = np.r_[False, dens[1:] > dens[:-1]]
    down = np.r_[dens[:-1] > dens[1:], False]
    modas = grid[up & down]
    alturas = dens[up & down]
    # discard modes below 10% of the main one: KDE noise
    if modas.size:
        keep = alturas >= 0.10 * alturas.max()
        modas, alturas = modas[keep], alturas[keep]

    valles = []
    for a, b in zip(modas[:-1], modas[1:]):
        m = (grid >= a) & (grid <= b)
        valles.append(grid[m][np.argmin(dens[m])])
    return grid, dens, modas, np.asarray(valles)


def _auc(x_pos: np.ndarray, x_neg: np.ndarray) -> float:
    """
    AUC = P(an element of pos > one of neg). 0.5 = indistinguishable.

    It is the unbiased measure of "does this variable separate the two classes?": it
    is scale-independent and does not assume normality, and unlike comparing
    medians it is not misled by overlapping bimodal distributions.
    """
    from scipy.stats import rankdata

    x_pos = np.asarray(x_pos, float); x_pos = x_pos[np.isfinite(x_pos)]
    x_neg = np.asarray(x_neg, float); x_neg = x_neg[np.isfinite(x_neg)]
    n1, n2 = x_pos.size, x_neg.size
    if n1 == 0 or n2 == 0:
        return np.nan
    r = rankdata(np.concatenate([x_pos, x_neg]))
    return float((r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2))


_TR_COLS = (
    "nuclear_frac", "log1p_total_counts", "log1p_n_genes_by_counts",
    "pct_counts_mt", "malat1_log1p", "cellbender_removed_frac",
    "n_spliced_log1p", "n_unspliced_log1p", "cell_probability",
)


def threshold_report(
    adatas,
    cols: Optional[Sequence[str]] = None,
    class_key: str = "droplet_class",
    debris_nf: float = 0.05,
    verbose: bool = True,
) -> dict:
    """
    The numbers behind plot_qc_overview, for one or several samples.

    adatas
        An AnnData, or a dict {name: AnnData}. If you pass an already
        concatenated object, split it by sample first: thresholds are per sample
        and a pooled table averages them exactly where they differ most.

    Returns a dict with three tables:

    'medians'   median of each metric per sample and droplet_class.
    'separation' AUC(cell vs empty_droplet) per metric. This is the table that
                 decides whether droplet_class can be trusted in that sample:
                   ~0.5  -> the class does not discriminate on that variable
                   >0.9  -> separates well
                   <0.1  -> separates well BUT inverted (the "cells" have
                            LESS of that metric than the "empty"); if this happens for
                            complexity or Malat1, the classification is
                            inverted and you should not exclude by it.
    'modes'      modes and valleys of nuclear_frac and cellbender_removed_frac.
                 Two modes in nuclear_frac above debris_nf = two
                 real cell populations (e.g. mono- vs multinucleated,
                 or nuclei vs cells), NOT debris: there max_nuclear_frac
                 must stay None until you know what they are.
    """
    if isinstance(adatas, ad.AnnData):
        adatas = {"sample": adatas}
    cols = list(cols) if cols is not None else list(_TR_COLS)

    med_rows, sep_rows, mod_rows = [], [], []

    for name, a in adatas.items():
        obs = a.obs
        pres = [c for c in cols if c in obs.columns]

        if class_key in obs.columns:
            g = obs.groupby(class_key, observed=True)[pres].median()
            for cls, row in g.iterrows():
                d = {"sample": name, "class": str(cls),
                     "n": int((obs[class_key] == cls).sum())}
                d.update(row.round(3).to_dict())
                med_rows.append(d)

            es_cell = (obs[class_key] == "cell").to_numpy()
            es_empty = (obs[class_key] == "empty_droplet").to_numpy()
            if es_cell.any() and es_empty.any():
                d = {"sample": name,
                     "n_cell": int(es_cell.sum()),
                     "n_empty": int(es_empty.sum())}
                for c in pres:
                    v = obs[c].to_numpy(dtype=float)
                    d[c] = round(_auc(v[es_cell], v[es_empty]), 3)
                sep_rows.append(d)
        else:
            d = {"sample": name, "class": "(sin droplet_class)", "n": a.n_obs}
            d.update(obs[pres].median().round(3).to_dict())
            med_rows.append(d)

        for c in ("nuclear_frac", "cellbender_removed_frac"):
            if c not in obs.columns:
                continue
            x = obs[c].to_numpy(dtype=float)
            if c == "nuclear_frac":
                x_sub = x[x >= debris_nf]          # without the debris peak
                etiqueta = f"{c} (>= {debris_nf})"
            else:
                x_sub, etiqueta = x, c
            _, _, modas, valles = _kde_modes(x_sub)
            mod_rows.append({
                "sample": name, "variable": etiqueta,
                "n_modes": int(modas.size),
                "modes": ", ".join(f"{m:.3f}" for m in modas),
                "valleys": ", ".join(f"{v:.3f}" for v in valles),
                "frac_below_debris_nf": (round(float((x < debris_nf).mean()), 3)
                                          if c == "nuclear_frac" else np.nan),
            })

    out = {
        "medians": pd.DataFrame(med_rows),
        "separation": pd.DataFrame(sep_rows),
        "modes": pd.DataFrame(mod_rows),
    }

    if verbose:
        for k, df in out.items():
            if df.empty:
                continue
            print(f"--- {k} ---")
            print(df.to_string(index=False))
            print()
        sep = out["separation"]
        if not sep.empty:
            print("[thr] How to read 'separation': AUC of cell versus "
                  "empty_droplet.\n"
                  "      0.5 = droplet_class distinguishes nothing on that metric.\n"
                  "      <0.5 on log1p_n_genes_by_counts or malat1_log1p = the\n"
                  "      classification is INVERTED in that sample; use\n"
                  "      droplet_exclude=() and direct thresholds.")
            for _, r in sep.iterrows():
                malos = [c for c in ("malat1_log1p", "log1p_n_genes_by_counts")
                         if c in r.index and pd.notna(r[c]) and r[c] < 0.60]
                if malos:
                    print(f"      -> {r['sample']}: {malos} with AUC "
                          f"{[r[c] for c in malos]}; do not trust the class here.")
    return out


def suggest_thresholds(
    adata: ad.AnnData,
    class_key: str = "droplet_class",
    debris_nf: float = 0.05,
    mad_k: float = 3.0,
    verbose: bool = True,
) -> dict:
    """
    Proposed thresholds for apply_qc_flags derived from the distribution.

    It is NOT automatic in the sense of "accept it and move on". It is a traceable
    starting point: each value comes from an explicit rule on the data and
    its origin is printed. Change it if you know the biology; but then the
    change is your decision and is recorded, which is different from having
    set a round number because it looked right in the figure.

    Rules:
      min_nuclear_frac      KDE valley between the debris peak and the first
                            cell cloud; if there are not two modes, debris_nf.
      max_nuclear_frac      ALWAYS None. A high nuclear_frac is not a
                            defect: it reflects intronic content, and splitting two
                            cell clouds on it discards a cell type.
      min_log1p_counts      1st percentile of the population above
                            min_nuclear_frac.
      min_log1p_genes       same.
      max_pct_mt            median + mad_k*MAD of that population, clipped to
                            [10, 25]. If the distribution is bimodal with a mode
                            near 100, a warning is given: that is not "high mito", they are
                            droplets containing only mitochondrial transcripts.
      min_log1p_malat1      only if Malat1 separates (AUC >= 0.75 against those with
                            low nuclear_frac); otherwise None.
      max_cellbender_removed  valley if there is a second mode, otherwise 99th percentile.
    """
    obs = adata.obs
    razones = {}

    def _q(col, q):
        return float(np.nanquantile(obs.loc[nucl, col].to_numpy(float), q))

    # --- min_nuclear_frac
    if "nuclear_frac" not in obs:
        raise ValueError("Missing nuclear_frac; run compute_qc_metrics() first")
    x = obs["nuclear_frac"].to_numpy(float)
    _, _, modas, valles = _kde_modes(x)
    if valles.size and modas.size >= 2 and modas[0] < debris_nf * 2:
        min_nf = float(round(valles[0], 3))
        razones["min_nuclear_frac"] = (
            f"KDE valley between the debris mode ({modas[0]:.3f}) and the "
            f"next one ({modas[1]:.3f})")
    else:
        min_nf = debris_nf
        razones["min_nuclear_frac"] = (
            f"no clear bimodality in nuclear_frac; using debris_nf={debris_nf}")

    nucl = obs["nuclear_frac"] >= min_nf
    if nucl.sum() < 100:
        raise ValueError(f"Only {int(nucl.sum())} droplets above "
                         f"min_nuclear_frac={min_nf}; check the threshold")

    out = {"min_nuclear_frac": min_nf, "max_nuclear_frac": None}
    razones["max_nuclear_frac"] = ("None by rule: high nuclear_frac is "
                                   "intronic content, not a defect")

    # --- depth and complexity
    for k, col in (("min_log1p_counts", "log1p_total_counts"),
                   ("min_log1p_genes", "log1p_n_genes_by_counts")):
        if col in obs:
            out[k] = round(_q(col, 0.01), 2)
            razones[k] = f"1st percentile of {col} in droplets with nuclear_frac >= {min_nf}"

    # --- mitochondrial
    if "pct_counts_mt" in obs:
        v = obs.loc[nucl, "pct_counts_mt"].to_numpy(float)
        v = v[np.isfinite(v)]
        med = float(np.median(v))
        mad = float(np.median(np.abs(v - med))) * 1.4826
        prop = med + mad_k * mad
        out["max_pct_mt"] = float(np.clip(round(prop, 1), 10, 25))
        razones["max_pct_mt"] = (f"median {med:.1f} + {mad_k}*MAD {mad:.1f} = "
                                 f"{prop:.1f}, clipped to [10, 25]")
        _, _, m_mt, _ = _kde_modes(obs["pct_counts_mt"].to_numpy(float))
        if m_mt.size >= 2 and m_mt.max() > 80:
            razones["max_pct_mt"] += (
                f"  WARNING: there is a mode at pct_counts_mt={m_mt.max():.0f}. "
                f"These are not stressed cells, they are droplets whose transcriptome "
                f"is almost entirely mitochondrial: lost cytoplasm. Inspect them separately.")

    # --- Malat1
    if "malat1_log1p" in obs:
        v = obs["malat1_log1p"].to_numpy(float)
        auc = _auc(v[nucl.to_numpy()], v[~nucl.to_numpy()])
        if np.isfinite(auc) and auc >= 0.75:
            out["min_log1p_malat1"] = round(_q("malat1_log1p", 0.01), 2)
            razones["min_log1p_malat1"] = (
                f"Malat1 separates (AUC={auc:.2f} against low nuclear_frac); "
                f"1st percentile of the retained population")
        else:
            out["min_log1p_malat1"] = None
            razones["min_log1p_malat1"] = (
                f"Malat1 does NOT separate here (AUC={auc:.2f}); filtering on it would remove "
                f"real cells. None.")

    # --- ambient
    if "cellbender_removed_frac" in obs:
        v = obs["cellbender_removed_frac"].to_numpy(float)
        _, _, m_cb, v_cb = _kde_modes(v)
        if v_cb.size and m_cb.size >= 2:
            out["max_cellbender_removed"] = float(round(v_cb[-1], 2))
            razones["max_cellbender_removed"] = (
                f"valley after the main mode; there are {m_cb.size} modes "
                f"({', '.join(f'{m:.2f}' for m in m_cb)}). A high second mode "
                f"is a group of droplets from which CellBender removed a lot: "
                f"inspect them before discarding them, they may be a low-complexity "
                f"cell type, not junk.")
        else:
            out["max_cellbender_removed"] = float(round(np.nanquantile(v, 0.99), 2))
            razones["max_cellbender_removed"] = "99th percentile (unimodal distribution)"

    # --- can droplet_class be used?
    fiable = None
    if class_key in obs.columns:
        c = (obs[class_key] == "cell").to_numpy()
        e = (obs[class_key] == "empty_droplet").to_numpy()
        if c.any() and e.any():
            aucs = {}
            for col in ("nuclear_frac", "malat1_log1p", "log1p_n_genes_by_counts"):
                if col in obs:
                    v = obs[col].to_numpy(float)
                    aucs[col] = _auc(v[c], v[e])
            fiable = all(a >= 0.75 for a in aucs.values() if np.isfinite(a))
            razones["droplet_exclude"] = (
                f"AUC cell vs empty: "
                + ", ".join(f"{k}={v:.2f}" for k, v in aucs.items())
                + ("  -> the class separates, droplet_exclude can be used"
                   if fiable else
                   "  -> the class does NOT separate; use droplet_exclude=() and filter by "
                     "direct thresholds"))
    out["_droplet_class_fiable"] = fiable

    if verbose:
        print("proposed thresholds:")
        for k, v in out.items():
            if k.startswith("_"):
                continue
            print(f"  {k:24s} = {v!r}")
        print("\nrationale:")
        for k, r in razones.items():
            print(f"  {k}: {r}")
        print("\n[thr] These are a starting point derived from the distribution of "
              "THIS sample. Before accepting them, run apply_qc_flags and check "
              "what each flag removes in qc_summary.")
    out["_razones"] = razones
    return out
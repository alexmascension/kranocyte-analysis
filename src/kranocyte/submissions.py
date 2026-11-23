"""
Preparation of the project deposits.

Three destinations, with different requirements:

* **Zenodo (data)** — the processed objects slimmed according to
  :class:`SlimSpec` below: from 48.65 GB to 32.9 GB. This is the re-entry point
  artefact: with it, notebook 1 onwards can be run.
* **Zenodo (code and results)** — notebooks with outputs, final tables and figures, plus
  the container image that produced the matrices.
* **CELLxGENE Discover** — the in-house object adapted to schema 7.1.0.

No function ever writes over its source object: all of them produce new files under the
deposit directory.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import anndata as ad

__all__ = [
    "SlimSpec", "SPEC_DEFAULT", "slim_adata", "slim_all",
    "write_manifest", "verify_manifest", "cxg_report",
    "CXG_SCHEMA_VERSION", "CXG_MOUSE_ANNOTATION", "CXG_GENCODE_URL",
    "gencode_symbol_map", "cxg_build",
]

#: Version of the CELLxGENE Discover schema the object is prepared against.
#: Checked on 2026-11-20 at
#: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md
CXG_SCHEMA_VERSION = "7.1.0"

#: Gene annotation PINNED by schema 7.1.0 for Mus musculus (NCBITaxon:10090).
#: The project quantified against NCBI RefSeq GRCm39, so a mapping is needed.
CXG_MOUSE_ANNOTATION = "GENCODE vM37 (GRCm39), Ensembl 114"

#: Layers that do not travel to the deposit, with the reason. All three are redundant or
#: never consumed; no information is lost by removing them.
DROP_ALWAYS: Mapping[str, str] = {
    "norm_scran_log1p": "it is X, with unsorted indices; exact duplicate",
    "norm_log1p": "produced and never consumed anywhere",
    "norm_cpm_log1p": "produced and never consumed anywhere",
}

#: RNA velocity layers. Used by `2AA` and `2OP`, but **kept in every object**: there is
#: ample headroom in the Zenodo record, and a third party may want to run velocity on a
#: dataset for which the project did not. Consistency across the eleven objects is worth
#: more than the 15 GB that dropping them would save.
VELOCYTO = ("spliced", "unspliced", "ambiguous")

#: `counts_raw` holds the pre-CellBender counts. Not used downstream of notebook 1, but
#: **kept**: it is what allows auditing what the decontamination removed, and that is
#: exactly the check a reviewer may want to redo.
KEEP_ALWAYS = ("counts_cellbender", "counts_raw", *VELOCYTO)


@dataclass(frozen=True)
class SlimSpec:
    """
    What is removed from an object before depositing it.

    By default only the three layers in `DROP_ALWAYS` are removed. `drop_velocyto_in`
    exists so the policy can be tightened if space is ever needed: it names the files whose
    velocity layers WOULD be dropped. Empty means keeping them in all of them.
    """

    drop_layers: Mapping[str, str] = field(default_factory=lambda: dict(DROP_ALWAYS))
    velocyto_layers: Sequence[str] = VELOCYTO
    drop_velocyto_in: frozenset[str] = frozenset()

    def layers_to_drop(self, filename: str, present: Iterable[str]) -> dict[str, str]:
        """Layers to remove from ``filename``, with the reason for each."""
        present = set(present)
        out = {k: v for k, v in self.drop_layers.items() if k in present}
        if filename in self.drop_velocyto_in:
            for layer in self.velocyto_layers:
                if layer in present:
                    out[layer] = "RNA velocity; dropped by an explicit space decision"
        return out


SPEC_DEFAULT = SlimSpec()


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    """sha256 of a file, read in blocks."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def slim_adata(src: Path, dst: Path, spec: SlimSpec = SPEC_DEFAULT,
               overwrite: bool = False, verbose: bool = True) -> dict:
    """
    Write to ``dst`` a copy of ``src`` without the layers the specification discards.

    ``src`` is not modified. Returns a record of the removed layers and the sizes, meant
    to go into the bundle manifest.

    The object is opened with ``backed='r'`` only to inspect which layers are present;
    writing needs it in memory, because removing a layer from a backed h5ad is not
    supported. The largest of these objects is 11.8 GB on disk.

    Parameters
    ----------
    src : pathlib.Path
        Source object. Never modified.
    dst : pathlib.Path
        Destination inside the deposit bundle.
    spec : SlimSpec
        Which layers to drop.
    overwrite : bool
        Refuse to overwrite an existing destination unless True.

    Returns
    -------
    dict
        file, source_bytes, deposit_bytes, dropped_layers, kept_layers, n_obs, n_vars.

    Raises
    ------
    FileExistsError
        If ``dst`` exists and ``overwrite`` is False.
    """
    src, dst = Path(src), Path(dst)
    if dst.exists() and not overwrite:
        raise FileExistsError(f"{dst} already exists; pass overwrite=True to rebuild it")

    peek = ad.read_h5ad(src, backed="r")
    present = list(peek.layers.keys())
    peek.file.close()

    drop = spec.layers_to_drop(src.name, present)
    if verbose:
        print(f"{src.name}: layers present {present}")
        for layer, reason in drop.items():
            print(f"    drop {layer:18s} - {reason}")

    adata = ad.read_h5ad(src)
    for layer in drop:
        del adata.layers[layer]

    adata.uns["deposit"] = {
        "slimmed_on": str(date.today()),
        "kept_note": ("counts_raw and the velocity layers are kept in all eleven objects, for "
                      "consistency and so a third party can audit the decontamination or run "
                      "RNA velocity on any dataset."),
        "dropped_layers": sorted(drop),
        "source_file": src.name,
        "source_sha256": sha256_of(src),
        "note": ("X is the scran-normalised, log1p matrix; the norm_scran_log1p layer was "
                 "removed because it is identical to X. prepare_cpdb_input detects this "
                 "and uses expm1(X)."),
    }

    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".part")
    adata.write_h5ad(tmp, compression="gzip")
    tmp.replace(dst)

    record = {
        "file": dst.name,
        "source_bytes": src.stat().st_size,
        "deposit_bytes": dst.stat().st_size,
        "dropped_layers": sorted(drop),
        "kept_layers": sorted(adata.layers.keys()),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
    }
    if verbose:
        ratio = 100 * record["deposit_bytes"] / record["source_bytes"]
        print(f"    {record['source_bytes']/1e9:.2f} GB -> {record['deposit_bytes']/1e9:.2f} GB "
              f"({ratio:.0f} %)")
    return record


def slim_all(sources: Sequence[Path], out_dir: Path, spec: SlimSpec = SPEC_DEFAULT,
             overwrite: bool = False) -> list[dict]:
    """Slim a list of objects into the deposit directory."""
    out_dir = Path(out_dir)
    records = []
    for src in sources:
        src = Path(src)
        records.append(slim_adata(src, out_dir / src.name, spec=spec, overwrite=overwrite))
    total_in = sum(r["source_bytes"] for r in records) / 1e9
    total_out = sum(r["deposit_bytes"] for r in records) / 1e9
    print(f"\nTotal: {total_in:.2f} GB -> {total_out:.2f} GB")
    return records


def write_manifest(bundle_dir: Path, manifest_name: str = "MANIFEST.tsv",
                   extra: Mapping[str, str] | None = None) -> Path:
    """
    Write the bundle manifest: relative path, bytes and sha256 of every file.

    It is what lets whoever downloads the deposit check they have what was published, and
    what connects the deposit to the commit that produced it.
    """
    bundle_dir = Path(bundle_dir)
    target = bundle_dir / manifest_name
    rows = ["path\tbytes\tsha256"]
    for path in sorted(p for p in bundle_dir.rglob("*") if p.is_file() and p.name != manifest_name):
        rows.append(f"{path.relative_to(bundle_dir)}\t{path.stat().st_size}\t{sha256_of(path)}")
    if extra:
        rows.append("")
        rows += [f"# {k}: {v}" for k, v in extra.items()]
    target.write_text("\n".join(rows) + "\n")
    print(f"{target}: {len(rows) - 1} entries")
    return target


def verify_manifest(bundle_dir: Path, manifest_name: str = "MANIFEST.tsv") -> list[str]:
    """Check the bundle against its manifest. Returns the list of discrepancies."""
    bundle_dir = Path(bundle_dir)
    problems = []
    for line in (bundle_dir / manifest_name).read_text().splitlines()[1:]:
        if not line or line.startswith("#"):
            continue
        rel, size, digest = line.split("\t")
        path = bundle_dir / rel
        if not path.is_file():
            problems.append(f"missing: {rel}")
        elif path.stat().st_size != int(size):
            problems.append(f"size differs: {rel}")
        elif sha256_of(path) != digest:
            problems.append(f"sha256 differs: {rel}")
    print("manifest matches" if not problems else f"{len(problems)} discrepancies")
    return problems


#: obs fields the curator MUST annotate under schema 7.1.0 (non-spatial).
CXG_OBS_REQUIRED = (
    "assay_ontology_term_id",
    "cell_type_ontology_term_id",
    "development_stage_ontology_term_id",
    "disease_ontology_term_id",
    "donor_id",
    "is_primary_data",
    "self_reported_ethnicity_ontology_term_id",
    "sex_ontology_term_id",
    "suspension_type",
    "tissue_ontology_term_id",
    "tissue_type",
)

#: uns fields the curator is responsible for.
CXG_UNS_REQUIRED = ("organism_ontology_term_id", "title")

#: Values already determined for the in-house dataset.
CXG_KNOWN_VALUES = {
    "organism_ontology_term_id": ("NCBITaxon:10090", "Mus musculus"),
    "sex_ontology_term_id": ("PATO:0000384", "male; determined from Xist and Y-chromosome gene expression"),
    "assay_ontology_term_id": ("EFO:0009922", "10x 3' v3; the SI-TT indices indicate v3.1, curated as v3"),
    "self_reported_ethnicity_ontology_term_id": ("na", "not applicable outside human"),
    "suspension_type": ("cell", "cell suspension, not nuclei"),
    "is_primary_data": ("True", "data generated by this project"),
}


def cxg_report(adata, verbose: bool = True) -> dict:
    """
    Report what the object is missing to meet the CELLxGENE Discover schema.

    Modifies nothing: it says what is there, what is missing and which values are already
    decided. The gene identifier conversion is not done here because it is not a rename:
    the schema pins GENCODE vM37 / Ensembl 114 while the project quantified against NCBI
    RefSeq GRCm39, so the mapping loses genes and has to be reviewed. See ``cxg_build``.
    """
    missing_obs = [c for c in CXG_OBS_REQUIRED if c not in adata.obs.columns]
    missing_uns = [k for k in CXG_UNS_REQUIRED if k not in adata.uns]
    var_looks_ensembl = bool(len(adata.var_names)) and str(adata.var_names[0]).startswith("ENSMUSG")

    report = {
        "schema_version": CXG_SCHEMA_VERSION,
        "required_annotation": CXG_MOUSE_ANNOTATION,
        "missing_obs": missing_obs,
        "missing_uns": missing_uns,
        "var_index_is_ensembl": var_looks_ensembl,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
    }
    if verbose:
        print(f"CELLxGENE schema {CXG_SCHEMA_VERSION} - required annotation: {CXG_MOUSE_ANNOTATION}")
        print(f"  obs: faltan {len(missing_obs)} campos -> {missing_obs}")
        print(f"  uns: faltan {len(missing_uns)} campos -> {missing_uns}")
        print(f"  var indexado por Ensembl: {var_looks_ensembl}"
              f"{'' if var_looks_ensembl else '  <- needs mapping from NCBI RefSeq'}")
        print("\n  Valores ya determinados:")
        for key, (value, why) in CXG_KNOWN_VALUES.items():
            print(f"    {key:42s} = {value:18s} ({why})")
    return report


#: GTF pinned by schema 7.1.0 for Mus musculus. The project quantified against
#: NCBI RefSeq GRCm39, so the mapping goes through gene symbols.
CXG_GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/"
    "gencode.vM37.primary_assembly.annotation.gtf.gz"
)


def gencode_symbol_map(gtf_path: Path, download: bool = True) -> "pd.DataFrame":
    """
    Symbol -> stable Ensembl identifier table, from the GENCODE vM37 GTF.

    The version suffix (``ENSMUSG00000051951.6`` -> ``ENSMUSG00000051951``) is stripped:
    the schema requires stable identifiers without a version.

    Downloads the GTF if absent (about 37 MB).

    Returns
    -------
    pandas.DataFrame
        Columns ``symbol`` (str) and ``ensembl`` (str), one row per unique pair.
    """
    import gzip
    import re as _re
    import urllib.request

    import pandas as pd

    gtf_path = Path(gtf_path)
    if not gtf_path.exists():
        if not download:
            raise FileNotFoundError(f"{gtf_path} is absent and download=False")
        gtf_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"descargando {CXG_GENCODE_URL}")
        urllib.request.urlretrieve(CXG_GENCODE_URL, gtf_path)

    rows = []
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] != "gene":
                continue
            gene_id = _re.search(r'gene_id "([^"]+)"', fields[8])
            gene_name = _re.search(r'gene_name "([^"]+)"', fields[8])
            if gene_id and gene_name:
                rows.append((gene_name.group(1), gene_id.group(1).split(".")[0]))

    table = pd.DataFrame(rows, columns=["symbol", "ensembl"]).drop_duplicates()
    print(f"GENCODE vM37: {len(table)} genes, {table.symbol.nunique()} unique symbols")
    return table


def cxg_build(adata, gene_map, out_path: Path, obs_values: Mapping[str, object],
              title: str, counts_layer: str = "counts_cellbender",
              verbose: bool = True):
    """
    Build the object for CELLxGENE Discover: reindex by Ensembl and drop the genes that
    are not in GENCODE vM37.

    The drop is unavoidable and not neutral, so the function **reports what it loses**:
    how many genes, of which biotype and at what detection level. On the in-house dataset
    18.5 % of genes are lost, but almost all are predicted `Gm*` entries and lncRNAs
    detected in very few cells.

    A symbol mapping to several Ensembl identifiers resolves to the first in the GTF; if
    duplicated identifiers remain after reindexing, the most-detected one is kept and a
    message is printed. No identifier is ever invented.

    ``obs_values`` are the schema fields the curator must annotate: they are written as
    given, without validating them against the ontologies (``cellxgene-schema`` does that).

    Returns
    -------
    anndata.AnnData
        The written object, reindexed by Ensembl identifier.
    """
    import numpy as np
    import pandas as pd

    mapping = gene_map.drop_duplicates("symbol").set_index("symbol")["ensembl"]
    keep = pd.Series(adata.var_names).isin(mapping.index).to_numpy()

    if verbose:
        lost = adata.var[~keep]
        print(f"genes: {keep.sum()} of {adata.n_vars} map to GENCODE vM37 "
              f"({keep.mean():.1%}); dropping {len(lost)}")
        if "gene_biotype" in lost.columns:
            print("  dropped, by biotype:")
            for biotype, n in lost["gene_biotype"].value_counts().head(5).items():
                print(f"    {biotype:28s} {n}")
        if "n_cells" in lost.columns and len(lost):
            print(f"  median cells in which they are detected: {lost['n_cells'].median():.0f} "
                  f"(versus {adata.var.loc[keep, 'n_cells'].median():.0f} among those kept)")

    out = adata[:, keep].copy()
    out.var["gene_symbol_original"] = out.var_names
    out.var_names = pd.Index([mapping[s] for s in out.var_names], name=None)

    duplicated = out.var_names.duplicated(keep=False)
    if duplicated.any():
        order = np.argsort(-out.var["n_cells"].to_numpy()) if "n_cells" in out.var else np.arange(out.n_vars)
        first = ~out.var_names[order].duplicated()
        selected = np.zeros(out.n_vars, bool)
        selected[order[first]] = True
        if verbose:
            print(f"  {int(duplicated.sum())} duplicated identifiers after reindexing; "
                  f"keeping the most-detected of each")
        out = out[:, selected].copy()

    # The schema requires the raw counts in addition to the normalised matrix.
    if counts_layer in out.layers:
        out.raw = out.copy()
        out.raw.X = out.layers[counts_layer]
    out.var["feature_is_filtered"] = False

    for key, value in obs_values.items():
        out.obs[key] = value
    out.uns["title"] = title
    out.uns["organism_ontology_term_id"] = "NCBITaxon:10090"
    out.uns["schema_version"] = CXG_SCHEMA_VERSION
    out.uns["gene_annotation"] = CXG_MOUSE_ANNOTATION
    out.uns["gene_mapping"] = {
        "route": "NCBI RefSeq GRCm39 symbol -> GENCODE vM37 Ensembl identifier",
        "n_in": int(adata.n_vars),
        "n_out": int(out.n_vars),
        "source_gtf": CXG_GENCODE_URL,
    }

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(out_path, compression="gzip")
    if verbose:
        print(f"wrote {out_path} ({out_path.stat().st_size/1e9:.2f} GB, "
              f"{out.n_obs} cells x {out.n_vars} genes)")
        print("  validate with: cellxgene-schema validate " + str(out_path))
    return out

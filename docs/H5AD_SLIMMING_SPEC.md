# Slimming the h5ad files for deposit

From **48.65 GB to 32.88 GB as an upper bound**, with no loss at all: the three layers removed are
redundant or not consumed anywhere.

The Zenodo limit is 50 GB per record, which leaves **17.12 GB of headroom (34 %)**. That is well
above the 10 % below which the bundle should be split. This figure is an **upper bound**: it does
not account for recompression on rewrite. The object that was verified went from 0.77 to
0.15 GB, i.e. 19.5 % and not the 67 % that the bound predicted, so the actual total will be well
below it. I do not extrapolate from a single file: the figure that holds is the bound.

Measured with `h5py` on the **on-disk** size of each element (`get_storage_size()`), without
loading the objects into memory.

## What the 48.65 GB consist of

| Group | GB | % |
|---|---|---|
| `layers` | 41.84 | 86 % |
| `X` | 5.26 | 11 % |
| `obsp` (neighbor graphs) | 0.90 | 2 % |
| `obsm` (embeddings) | 0.33 | 1 % |
| `obs` + `var` + `varm` + `uns` | 0.32 | 1 % |

The problem is `layers` alone: eight layers of about 5 GB each.

## Decision per layer

| Layer | GB | Decision | Evidence |
|---|---|---|---|
| `counts_cellbender` | 5.26 | **keep** | `degs.py:99,183`, `characterization.py:46` and `interactions.py:270` use it as the default `counts_layer` |
| `norm_scran_log1p` | 5.26 | **remove** | **It is `X`.** Identical row by row once the indices are sorted. The layer was saved with unsorted indices, which is the only difference between them |
| `norm_log1p` | 5.26 | **remove** | It does not appear in any notebook or in `pyfuncs`: it is produced and never consumed |
| `norm_cpm_log1p` | 5.26 | **remove** | Same as the previous one |
| `counts_raw` | 5.72 | **keep** | Not used downstream of notebook 1 (`io.py:287,742-749` uses it only in the `0*` notebooks). But it is what allows auditing what CellBender removed, and that is exactly the check a reviewer may want to redo |
| `spliced` / `unspliced` / `ambiguous` | 15.09 | **keep in all eleven** | Required by `normalize_velocyto_layers` + `run_velocity`, which the project only calls in `2AA` and `2OP`. They are kept in all objects anyway: there is plenty of headroom, and a third party may want to run velocyto on a dataset where the project did not |

### Why `X` is enough where `norm_scran_log1p` is removed

`uns/normalization` records `scran`, and the direct comparison confirms it: for rows 0, 1 and 500
of `AA_inhouse_dataset_processed.h5ad`, the sorted indices and the values of `X` and of
`layers/norm_scran_log1p` match exactly. A naive positional comparison says they differ, because
the layer has unsorted indices within each row.

**Adaptation, already applied in the code:** `interactions.prepare_cpdb_input()` had
`log_layer: Optional[str] = "norm_scran_log1p"` and, if the layer was missing, it **silently fell
back** to library-size normalization — a different method, with no warning. Now, when the layer is
absent but `uns['normalization']` records `scran`, it uses `expm1(X)`, which is the same matrix.
If it cannot determine this, it emits a `warning` instead of silently switching method.

### Why `counts_raw` and velocyto are kept in all eleven

They could be removed from the objects that do not use them, saving about 15 GB more, but it is
not worth it:

- **Consistency.** Eleven objects with the same structure are easier to document and use than
  eleven with three variants depending on which notebook opened them.
- **Reuse.** Only `2AA` and `2OP` use the velocyto layers, but someone who downloads Nicoletti or
  Song may want to run RNA velocity there, which the project did not do. Removing them would close
  that door for no reason.
- **Auditability.** `counts_raw` versus `counts_cellbender` is what allows checking what the
  decontamination removed. It is exactly the kind of check a reviewer may ask for.
- **There is room.** 17.12 GB of headroom below the Zenodo limit, and that is before counting
  recompression.

`SlimSpec.drop_velocyto_in` stays in the code, empty, so the spec can be tightened if space is
ever needed without rewriting it.

## Result per file

| File | Current (GB) | Removed (GB) | Remaining (GB) |
|---|---|---|---|
| `OP_oprescu_processed.h5ad` | 11.81 | 3.80 | 8.01 |
| `NI_nicoletti_processed.h5ad` | 8.59 | 2.75 | 5.84 |
| `SU_southerland_processed.h5ad` | 4.60 | 1.48 | 3.12 |
| `OP_oprescu_processed_FAPs.h5ad` | 4.57 | 1.49 | 3.07 |
| `NI_nicoletti_processed_FAPs.h5ad` | 4.51 | 1.47 | 3.04 |
| `SO_song_processed.h5ad` | 3.76 | 1.22 | 2.54 |
| `OP_oprescu_processed_FAPs_cleared.h5ad` | 2.92 | 0.96 | 1.96 |
| `AA_inhouse_dataset_processed.h5ad` | 2.45 | 0.81 | 1.64 |
| `SO_song_processed_FAPs.h5ad` | 2.39 | 0.79 | 1.61 |
| `AA_inhouse_dataset_processed_FAPs.h5ad` | 2.28 | 0.75 | 1.53 |
| `SU_southerland_processed_FAPs.h5ad` | 0.76 | 0.25 | 0.51 |
| **Total** | **48.65** | **15.77** | **32.88** |

## What is kept intact

`X`, the `counts_cellbender`, `counts_raw`, `spliced`, `unspliced` and `ambiguous` layers, `obs`,
`var`, `obsm` (`X_pca`, `X_harmony`, `X_umap`, `X_triku`), `varm`, `obsp` (the three neighbor
graphs: `neighbors_pca`, `neighbors_harmony` and the default one) and the full `uns`
—including `normalization`, `triku_params`, `cellbender`, `qc`, `processing` and `annotation`,
which are the provenance record of the object. Everything other than `X` and the layers — that
is, `obsp` (0.90 GB), `obsm` (0.33 GB) and `obs`/`var`/`varm`/`uns` (0.32 GB) — adds up to
1.55 GB across the eleven objects. There is nothing to gain by trimming that and much to lose;
the 32.88 GB figure above is dominated by `X` and the kept layers.

## What is no longer reproducible from the deposit

Practically nothing. Because `counts_raw` and the velocyto layers are kept, the following remain
possible from the deposit: the per-gene ambient contamination audit (`counts_raw` versus
`counts_cellbender`), the check that `GeneFull` ≈ spliced + unspliced + ambiguous, and running
RNA velocity on any of the eleven objects.

The only thing that cannot be redone is the **comparison between the three normalizations**,
because `norm_log1p` and `norm_cpm_log1p` are not shipped. That decision is made in the `0*`
notebooks, whose entry point is the FASTQ files, and `uns['normalization']` records what was
chosen and with which parameters.

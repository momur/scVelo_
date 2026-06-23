#!/usr/bin/env python3
"""
scvelo_integrated_pipeline.py
------------------------------
RNA velocity analysis pipeline for a Seurat-integrated single-cell dataset.

Workflow:
  1. Convert Seurat MTX/CSV exports -> AnnData (h5ad)
  2. Load per-sample loom files (spliced/unspliced counts)
  3. Merge and pre-process
  4. Stochastic velocity
  5. Dynamical velocity
  6. Visualise and export results

Usage:
    python scvelo_integrated_pipeline.py

Edit the Configuration section below before running.
"""

# ── 1. Imports ──────────────────────────────────────────────────


import os
import numpy as np
import pandas as pd
import anndata
from scipy import io
from scipy.sparse import csr_matrix

import scanpy as sc
import scvelo as scv
import cellrank as cr

# Verbosity / figure settings
scv.settings.verbosity = 3          # 0=errors 1=warnings 2=info 3=hints
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2



# ── 2. Configuration ────────────────────────────────────────────
# Edit the paths below to match your project layout before running the notebook.


# ── Edit these paths ─────────────────────────────────────────────────────────
DATA_DIR = "/path/to/seurat_outputs"        # MTX / CSV exports from Seurat
LOOM_DIR = "/path/to/loom_files"            # .loom files from STARsolo/velocyto
OUT_DIR  = "/path/to/output"                # where h5ad + figures are saved

LOOM_D30 = os.path.join(LOOM_DIR, "gex_possorted_bam_96HZZ_d30.loom")
LOOM_D50 = os.path.join(LOOM_DIR, "d50_gex_possorted_bam_8VM52.loom")

# Cell-type colour palette (12 colours — extend if you have more cell types)
CELLTYPE_COLORS = [
    "#00C19A", "#00A9FF", "#FF61CC", "#F8766D",
    "#00BCD7", "#00BE67", "#E68613", "#ABA300",
    "#C77CFF", "#8494FF", "#0CB702", "#d3d3d3",
]

os.makedirs(OUT_DIR, exist_ok=True)
os.chdir(OUT_DIR)



# ── 3. Convert Seurat output to AnnData ─────────────────────────
# Reads the MTX count matrix, cell metadata, gene names, PCA, and UMAP
# coordinates exported from Seurat and assembles them into an AnnData object.


# Count matrix (transpose: Seurat exports genes × cells; AnnData expects cells × genes)
X = io.mmread(os.path.join(DATA_DIR, "integrated_seurat_rna_counts.mtx"))
adata = anndata.AnnData(X=csr_matrix(X.T))

# Cell metadata
cell_meta = pd.read_csv(os.path.join(DATA_DIR, "integrated_seurat_metadata.csv"))
adata.obs = cell_meta
adata.obs.index = adata.obs["barcode"]

# Gene names
with open(os.path.join(DATA_DIR, "integrated_seurat_rna_counts.csv")) as f:
    gene_names = f.read().splitlines()
adata.var.index = gene_names

# Embeddings
pca = pd.read_csv(os.path.join(DATA_DIR, "integrated_pca.csv"), index_col=0)
pca.index = adata.obs.index
adata.obsm["X_pca"] = pca.to_numpy()
adata.obsm["X_umap"] = np.vstack(
    (adata.obs["UMAP_1"].to_numpy(), adata.obs["UMAP_2"].to_numpy())
).T

# Sanity check
print(adata)
sc.pl.umap(adata, color=["celltype"], frameon=False,
           save=os.path.join(OUT_DIR, "UMAP_seurat_converted.pdf"))

adata.write(os.path.join(OUT_DIR, "integrated_seurat_to_h5ad.h5ad"))



# ── 4. Load loom files and merge with AnnData ───────────────────
# Loom files contain the spliced/unspliced count matrices produced by STARsolo
# or velocyto. Barcodes are renamed to match the Seurat convention
# (`D30_<barcode>-1`) before merging.


adata = sc.read_h5ad(os.path.join(OUT_DIR, "integrated_seurat_to_h5ad.h5ad"))

d30_data = scv.read(LOOM_D30, cache=True)
d50_data = scv.read(LOOM_D50, cache=True)


def rename_barcodes(ldata, prefix):
    """Strip the CB prefix added by STARsolo and apply sample prefix."""
    bcs = [bc.split(":")[1] for bc in ldata.obs.index]
    bcs = [f"{prefix}_{bc[:-1]}-1" for bc in bcs]
    ldata.obs.index = bcs
    return ldata


d30_data = rename_barcodes(d30_data, "D30")
d50_data = rename_barcodes(d50_data, "D50")

d30_data.var_names_make_unique()
d50_data.var_names_make_unique()

ldata = d30_data.concatenate([d50_data])

scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)

adata_merged = scv.utils.merge(adata, ldata)
print(adata_merged)

sc.pl.umap(adata_merged, color=["celltype"], frameon=False,
           save=os.path.join(OUT_DIR, "UMAP_merged.pdf"))

adata_merged.write_h5ad(
    os.path.join(OUT_DIR, "adata_integrated_seurat_and_loom_merged.h5ad")
)



# ── 5. Pre-processing ───────────────────────────────────────────
# Filter, normalise, and compute the neighbourhood graph.
# NOTE: filter_and_normalize already applies log1p — do NOT call scv.pp.log1p separately.
# NOTE: sc.pp.neighbors may crash the kernel unless you are running inside the
# NOTE: dedicated 'scvelo' conda environment (see README).


adata_input = sc.read_h5ad(
    os.path.join(OUT_DIR, "adata_integrated_seurat_and_loom_merged.h5ad")
)

scv.pl.proportions(adata_input, groupby="celltype")

# filter_and_normalize includes log1p — do NOT call scv.pp.log1p separately
scv.pp.filter_and_normalize(adata_input, min_shared_counts=20, n_top_genes=5000)

sc.tl.pca(adata_input)
sc.pp.neighbors(adata_input, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata_input, n_pcs=30, n_neighbors=30)



# ── 6. RNA Velocity — Stochastic model ──────────────────────────


scv.tl.velocity(adata_input, mode="stochastic")
scv.tl.velocity_graph(adata_input)

adata_input.write_h5ad(
    os.path.join(OUT_DIR, "adata_integrated_stochastic_velocity.h5ad")
)



# ── 7. Visualise stochastic velocity ────────────────────────────


adata_input = sc.read_h5ad(
    os.path.join(OUT_DIR, "adata_integrated_stochastic_velocity.h5ad")
)
adata_input.uns["celltype_colors"] = CELLTYPE_COLORS

# UMAP with custom colours
sc.pl.umap(adata_input, frameon=True, color="celltype",
           save=os.path.join(OUT_DIR, "UMAP_stochastic.pdf"),
           legend_loc="on data", legend_fontsize="small")

# Velocity embedding grids at different scales
for scale, tag in [(None, "default"), (0.5, "scale05"), (0.25, "scale025")]:
    kwargs = dict(
        basis="umap", color="celltype",
        title="Integrated Stochastic Velocity",
        save=os.path.join(OUT_DIR, f"stochastic_grid_{tag}.pdf"),
    )
    if scale is not None:
        kwargs["scale"] = scale
    scv.pl.velocity_embedding_grid(adata_input, **kwargs)

# Velocity stream
scv.pl.velocity_embedding_stream(
    adata_input, frameon=True, basis="umap", color="celltype",
    title="Integrated Stochastic Velocity Stream",
    save=os.path.join(OUT_DIR, "stochastic_stream.pdf"),
)



# ── 8. Rank velocity genes ──────────────────────────────────────


scv.tl.rank_velocity_genes(adata_input, groupby="celltype", min_corr=0.3)
df = scv.DataFrame(adata_input.uns["rank_velocity_genes"]["names"])
df.to_csv(os.path.join(OUT_DIR, "stochastic_rank_velocity_genes.csv"), index=False)
df.head()


adata_input.write_h5ad(os.path.join(OUT_DIR, "Integrated_stochastic_anndata.h5ad"))
scv.logging.print_versions()



# ── 9. RNA Velocity — Dynamical model ───────────────────────────
# The dynamical model fits per-gene transcriptional kinetics and is more
# accurate than the stochastic model but significantly slower
# (`recover_dynamics` can take 30–90 min on a large dataset).


dyn_input = sc.read_h5ad(
    os.path.join(OUT_DIR, "adata_integrated_seurat_and_loom_merged.h5ad")
)

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.settings.set_figure_params('scvelo')

scv.pp.filter_and_normalize(dyn_input, min_shared_counts=20, n_top_genes=5000)
scv.pp.moments(dyn_input, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(dyn_input)
scv.tl.velocity(dyn_input, mode='dynamical')
scv.tl.velocity_graph(dyn_input)

dyn_input.write(
    os.path.join(OUT_DIR, "dynamical_velocity_integrated.h5ad"),
    compression='gzip'
)



# ── 10. Visualise dynamical velocity ────────────────────────────


scv.pl.velocity_embedding_stream(
    dyn_input, basis='umap',
    save=os.path.join(OUT_DIR, 'dynamical_stream.pdf')
)
print(f"Number of velocity genes: {dyn_input.var.velocity_genes.sum()}")

# Heatmap — top genes ordered by latent time
top_genes = dyn_input.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(dyn_input, var_names=top_genes, sortby='latent_time',
               col_color='celltype', n_convolve=100,
               save=os.path.join(OUT_DIR, 'dynamical_top_genes_heatmap.pdf'))

# Phase portraits for top 15 genes
top15 = dyn_input.var['fit_likelihood'].sort_values(ascending=False).index[:15]
scv.pl.scatter(dyn_input, basis=top15, ncols=5, frameon=False,
               color='celltype',
               save=os.path.join(OUT_DIR, 'dynamical_top15_scatter.pdf'))

# scVelo Integrated RNA Velocity Analysis

RNA velocity analysis pipeline for a Seurat-integrated single-cell RNA-seq dataset using [scVelo](https://scvelo.readthedocs.io/) and [CellRank](https://cellrank.readthedocs.io/).

## Overview

This notebook takes a Seurat-integrated dataset (exported as MTX + CSV) and per-sample loom files (spliced/unspliced counts), merges them, and computes RNA velocity using both the **stochastic** and **dynamical** models.

```
Seurat MTX/CSV exports
        │
        ▼
  AnnData (h5ad)  ←──  loom files (STARsolo / velocyto)
        │
        ▼
  Pre-processing (filter, normalize, PCA, neighbours, moments)
        │
        ├──▶  Stochastic velocity → embeddings + ranked genes
        │
        └──▶  Dynamical velocity  → latent time + phase portraits
```

## Requirements

### Conda environment

It is strongly recommended to run this notebook inside a dedicated `scvelo` conda environment. The `sc.pp.neighbors` step will crash the kernel in a generic environment.

```bash
conda create -n scvelo python=3.9
conda activate scvelo
pip install scvelo cellrank scanpy anndata session-info
```

### Input files

| File | Source |
|---|---|
| `integrated_seurat_rna_counts.mtx` | Seurat `WriteMMMatrix()` |
| `integrated_seurat_rna_counts.csv` | Gene names (one per line) |
| `integrated_seurat_metadata.csv` | Seurat metadata with `barcode`, `celltype`, `UMAP_1`, `UMAP_2` |
| `integrated_pca.csv` | Seurat PCA embeddings |
| `*.loom` | STARsolo or velocyto per-sample loom files |

## Usage

1. Clone this repo and activate the `scvelo` environment.
2. Open `scVelo_integrated_improved.ipynb` in Jupyter.
3. Edit the **Configuration** cell (Section 2) to set your `DATA_DIR`, `LOOM_DIR`, and `OUT_DIR`.
4. Run all cells top to bottom.

## Output files

| File | Description |
|---|---|
| `integrated_seurat_to_h5ad.h5ad` | Seurat data converted to AnnData |
| `adata_integrated_seurat_and_loom_merged.h5ad` | Merged RNA + spliced/unspliced |
| `adata_integrated_stochastic_velocity.h5ad` | Stochastic velocity computed |
| `Integrated_stochastic_anndata.h5ad` | Final stochastic AnnData |
| `stochastic_rank_velocity_genes.csv` | Per-celltype ranked velocity genes |
| `dynamical_velocity_integrated.h5ad` | Dynamical model output (gzipped) |
| `*.pdf` | Figures (UMAP, grid, stream, heatmap, scatter) |

## Known issues

- **pandas >= 2.0** causes `AttributeError: can't set attribute` in some scVelo plot functions. Downgrade to pandas 1.3.5 if you hit this:
  ```bash
  pip install "pandas==1.3.5"
  ```
- `recover_dynamics` (dynamical model) is slow — allow 30–90 min on a large dataset.

## Citation

If you use scVelo, please cite:

> Bergen et al. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. *Nature Biotechnology*.

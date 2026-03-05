ReDeconv Atlas Runbook
======================

Goal
- Reproduce the current visualization results in `visualization.ipynb` starting from:
  - `/gpfs/gpfs1/scratch/c9881013/felix_granada/granada/atlas/merged_annotated.h5ad`

Working directory
- `cd /gpfs/gpfs1/scratch/c9881013/felix_granada/granada/redeconv/redeconv_atlas`

Environment
- `conda activate redeconv-env`

1) Build demo ReDeconv inputs from atlas
- This uses:
  - `Cell_type <- cell_type`
  - `Sample_ID <- batch`
  - `Sample_ID_original <- sample`
- Command:
  - `python prepare_demo_atlas_from_h5ad.py --h5ad /gpfs/gpfs1/scratch/c9881013/felix_granada/granada/atlas/merged_annotated.h5ad --sample-col batch --sample-original-col sample --max-cells-per-sample-celltype 20 --write-chunk-rows 100`

2) Create ReDeconv-safe 3-column meta (and preserve original sample column)
- Command:
  - `python prepare_redeconv_meta_columns.py --input /gpfs/gpfs1/scratch/c9881013/felix_granada/granada/redeconv/redeconv_atlas/demo_data_4_normalization/Demo_ATLAS_meta_original.tsv --output-3col /gpfs/gpfs1/scratch/c9881013/felix_granada/granada/redeconv/redeconv_atlas/demo_data_4_normalization/Demo_ATLAS_meta.tsv --output-preserved /gpfs/gpfs1/scratch/c9881013/felix_granada/granada/redeconv/redeconv_atlas/demo_data_4_normalization/Demo_ATLAS_meta_original.tsv`

3) Run ReDeconv normalization (all 4 steps)
- `python ReDeconv_Normalization_atlas.py --choice 1`
- `python ReDeconv_Normalization_atlas.py --choice 2`
- `python ReDeconv_Normalization_atlas.py --choice 3`
- `python ReDeconv_Normalization_atlas.py --choice 4`

4) Open visualization notebook
- Notebook:
  - `/gpfs/gpfs1/scratch/c9881013/felix_granada/granada/redeconv/redeconv_atlas/visualization.ipynb`
- It contains:
  - UMAP for normalized data from `Results_4_normalization/scRNA_seq_new_noShift.tsv` + `Results_4_normalization/Meta_data_new.tsv`
  - UMAP for non-normalized data from `demo_data_4_normalization/Demo_ATLAS_scRNAseq.tsv` + `demo_data_4_normalization/Demo_ATLAS_meta_original.tsv`
  - Coloring by `Sample_ID_original`

Key outputs used by visualization
- Normalized:
  - `Results_4_normalization/scRNA_seq_new_noShift.tsv`
  - `Results_4_normalization/Meta_data_new.tsv`
- Non-normalized:
  - `demo_data_4_normalization/Demo_ATLAS_scRNAseq.tsv`
  - `demo_data_4_normalization/Demo_ATLAS_meta_original.tsv`

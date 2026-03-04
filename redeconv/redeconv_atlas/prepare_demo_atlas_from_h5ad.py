#!/usr/bin/env python3
"""
Prepare ReDeconv demo-style input files from a .h5ad atlas.

Outputs under <out_dir>:
  - demo_data_4_normalization/Demo_ATLAS_meta.tsv
  - demo_data_4_normalization/Demo_ATLAS_scRNAseq.tsv
  - demo_data_4_deconvolution/References_Meta_data_subset.tsv
  - demo_data_4_deconvolution/References_scRNA_seq_Nor.tsv
  - demo_data_4_deconvolution/Synthetic_Bulk_RNA_seq_Equal_Fraction_TPM.tsv
"""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path
from typing import Iterable, List

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

try:
    from tqdm import tqdm
except Exception:
    # Fallback keeps script runnable if tqdm is unavailable.
    def tqdm(iterable=None, **kwargs):
        return iterable


def _pick_col(candidates: Iterable[str], cols: Iterable[str]) -> str | None:
    cols = set(cols)
    for c in candidates:
        if c in cols:
            return c
    return None


def _ensure_dirs(base: Path) -> None:
    for d in (
        base / "demo_data_4_normalization",
        base / "demo_data_4_deconvolution",
        base / "Results_4_normalization",
        base / "Results_4_deconvolution",
    ):
        d.mkdir(parents=True, exist_ok=True)


def _make_synthetic_bulk(
    expr_gene_by_cell: pd.DataFrame, meta: pd.DataFrame, out_path: Path, n_bulk: int, seed: int
) -> None:
    # expr_gene_by_cell index=gene, columns=cell ids
    cell_to_type = meta.set_index("Cell_ID")["Cell_type"]
    ct_names: List[str] = sorted(cell_to_type.unique().tolist())

    ct_means = {}
    for ct in tqdm(ct_names, desc="Computing cell-type means", unit="cell_type"):
        cells = cell_to_type[cell_to_type == ct].index
        ct_means[ct] = expr_gene_by_cell[cells].mean(axis=1)
    ct_means = pd.DataFrame(ct_means)

    rng = np.random.default_rng(seed)
    bulk_cols = {}
    for i in tqdm(range(1, n_bulk + 1), desc="Generating synthetic bulk", unit="sample"):
        p = rng.random(len(ct_names))
        p = p / p.sum()
        bulk_cols[f"atlas_bulk_{i}"] = (ct_means.values * p).sum(axis=1)

    bulk = pd.DataFrame(bulk_cols, index=ct_means.index)
    bulk.insert(0, "sample", bulk.index.astype(str))
    bulk.to_csv(out_path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create ReDeconv atlas demo TSVs from merged_annotated.h5ad")
    parser.add_argument(
        "--h5ad",
        default="/gpfs/gpfs1/scratch/c9881013/felix_granada/granada/atlas/merged_annotated.h5ad",
        help="Path to source .h5ad file",
    )
    parser.add_argument(
        "--out-dir",
        default="/gpfs/gpfs1/scratch/c9881013/felix_granada/granada/redeconv/redeconv_atlas",
        help="Output directory (redeconv_atlas root)",
    )
    parser.add_argument("--cell-type-col", default=None, help="obs column for cell type (auto if omitted)")
    parser.add_argument("--sample-col", default=None, help="obs column for sample ID (auto if omitted)")
    parser.add_argument("--min-cells-per-type", type=int, default=300, help="Minimum cells per cell type to keep")
    parser.add_argument("--max-cell-types", type=int, default=10, help="Maximum number of cell types to keep")
    parser.add_argument(
        "--max-cells-per-sample-celltype",
        type=int,
        default=120,
        help="Cap cells sampled per (Sample_ID, Cell_type)",
    )
    parser.add_argument("--n-bulk", type=int, default=20, help="Number of synthetic bulk samples to generate")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    h5ad_path = Path(args.h5ad)
    out_root = Path(args.out_dir)
    _ensure_dirs(out_root)

    norm_dir = out_root / "demo_data_4_normalization"
    deconv_dir = out_root / "demo_data_4_deconvolution"

    print(f"Reading atlas in backed mode: {h5ad_path}")
    atlas = ad.read_h5ad(h5ad_path, backed="r")

    cell_col = args.cell_type_col or _pick_col(
        ["cell_type", "cell_type_manual", "cell_type_coarse", "cellType", "celltype"],
        atlas.obs.columns,
    )
    sample_col = args.sample_col or _pick_col(
        ["sample", "sample_original", "sampleID", "Sample", "batch", "dataset"],
        atlas.obs.columns,
    )
    if not cell_col or not sample_col:
        raise ValueError(f"Could not resolve columns. cell_col={cell_col}, sample_col={sample_col}")

    obs = atlas.obs[[cell_col, sample_col]].copy()
    obs.columns = ["Cell_type", "Sample_ID"]
    obs.index = obs.index.astype(str)
    obs = obs.dropna()

    # Keep most represented cell types above threshold.
    counts = obs["Cell_type"].value_counts()
    keep_ct = counts[counts >= args.min_cells_per_type].index.tolist()[: args.max_cell_types]
    obs = obs[obs["Cell_type"].isin(keep_ct)]

    # Balanced sample: up to N cells per (sample, cell type).
    rng = np.random.default_rng(args.seed)
    selected = []
    groups = list(obs.groupby(["Sample_ID", "Cell_type"]))
    for _, grp in tqdm(groups, desc="Sampling cells", unit="group"):
        n = min(args.max_cells_per_sample_celltype, len(grp))
        # deterministic sampling via numpy indices
        idx = grp.index.to_numpy()
        if len(idx) > n:
            pick = rng.choice(idx, size=n, replace=False)
        else:
            pick = idx
        selected.extend(pick.tolist())

    selected = pd.Index(selected)
    meta = obs.loc[selected].copy()
    meta["Cell_ID"] = meta.index.astype(str)
    meta = meta[["Cell_ID", "Cell_type", "Sample_ID"]]

    print(f"Selected cells: {meta.shape[0]} | genes: {atlas.n_vars}")

    # Materialize only selected slice to memory.
    print("Materializing selected cell/gene block into memory...")
    sub = atlas[selected, :].to_memory()
    X = sub.X
    if sparse.issparse(X):
        print("Converting sparse matrix to dense array...")
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)

    print("Building expression DataFrame (genes x cells)...")
    expr = pd.DataFrame(X.T, index=sub.var_names.astype(str), columns=sub.obs_names.astype(str))
    expr.insert(0, "Gene_sample", expr.index.astype(str))
    expr = expr.reset_index(drop=True)

    # Normalization inputs
    meta_norm = norm_dir / "Demo_ATLAS_meta.tsv"
    expr_norm = norm_dir / "Demo_ATLAS_scRNAseq.tsv"
    print("Writing normalization input files...")
    meta.to_csv(meta_norm, sep="\t", index=False)
    # Compact float formatting reduces file size and write time for large matrices.
    expr.to_csv(expr_norm, sep="\t", index=False, float_format="%.6g")

    # Deconvolution inputs
    meta_ref = deconv_dir / "References_Meta_data_subset.tsv"
    expr_ref = deconv_dir / "References_scRNA_seq_Nor.tsv"
    bulk_out = deconv_dir / "Synthetic_Bulk_RNA_seq_Equal_Fraction_TPM.tsv"
    print("Writing deconvolution reference files...")
    meta.to_csv(meta_ref, sep="\t", index=False)
    # Reuse the large expression TSV instead of writing it twice.
    if expr_ref.exists():
        expr_ref.unlink()
    try:
        os.link(expr_norm, expr_ref)
        print("Linked reference expression TSV to normalization TSV (hardlink).")
    except OSError:
        shutil.copyfile(expr_norm, expr_ref)
        print("Hardlink unavailable; copied expression TSV for reference input.")
    _make_synthetic_bulk(expr.set_index("Gene_sample"), meta, bulk_out, n_bulk=args.n_bulk, seed=args.seed)

    readme = out_root / "README.txt"
    readme.write_text(
        "\n".join(
            [
                f"Prepared from: {h5ad_path}",
                f"Selected obs columns: Cell_type <- {cell_col}, Sample_ID <- {sample_col}",
                f"Cells selected: {meta.shape[0]}",
                f"Genes selected: {sub.n_vars}",
                "Normalization input:",
                f"  - {meta_norm}",
                f"  - {expr_norm}",
                "Deconvolution input:",
                f"  - {meta_ref}",
                f"  - {expr_ref}",
                f"  - {bulk_out}",
            ]
        )
        + "\n"
    )

    try:
        atlas.file.close()
    except Exception:
        pass

    print("Done.")
    print(f"  Meta: {meta_norm}")
    print(f"  Expr: {expr_norm}")
    print(f"  Bulk: {bulk_out}")


if __name__ == "__main__":
    main()

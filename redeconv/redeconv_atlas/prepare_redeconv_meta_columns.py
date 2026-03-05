#!/usr/bin/env python3
"""
Prepare ReDeconv-compatible metadata columns.

Input can contain either:
  - Cell_ID, Cell_type, Sample_ID
  - Cell_ID, Cell_type, Sample_ID, Sample_ID_original

Outputs:
  1) 3-column ReDeconv file: Cell_ID, Cell_type, Sample_ID
  2) optional preserved file including Sample_ID_original when present
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Create ReDeconv-safe 3-column metadata TSV.")
    parser.add_argument(
        "--input",
        default="/gpfs/gpfs1/scratch/c9881013/felix_granada/granada/redeconv/redeconv_atlas/demo_data_4_normalization/Demo_ATLAS_meta_original.tsv",
        help="Input metadata TSV",
    )
    parser.add_argument(
        "--output-3col",
        default="/gpfs/gpfs1/scratch/c9881013/felix_granada/granada/redeconv/redeconv_atlas/demo_data_4_normalization/Demo_ATLAS_meta.tsv",
        help="Output TSV with exactly 3 columns for ReDeconv",
    )
    parser.add_argument(
        "--output-preserved",
        default="",
        help="Optional output TSV preserving Sample_ID_original if present",
    )
    args = parser.parse_args()

    in_path = Path(args.input)
    out_3col = Path(args.output_3col)
    out_preserved = Path(args.output_preserved) if args.output_preserved else None

    df = pd.read_csv(in_path, sep="\t")

    required = ["Cell_ID", "Cell_type", "Sample_ID"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Found: {df.columns.tolist()}")

    # Normalize IDs/labels.
    for c in required:
        df[c] = df[c].astype(str).str.strip()

    # Write ReDeconv-safe metadata (exactly 3 columns).
    meta_3 = df[required].copy()
    out_3col.parent.mkdir(parents=True, exist_ok=True)
    meta_3.to_csv(out_3col, sep="\t", index=False)

    # Optionally preserve the original sample column as a sidecar file.
    if out_preserved is not None:
        keep = required + (["Sample_ID_original"] if "Sample_ID_original" in df.columns else [])
        preserved = df[keep].copy()
        out_preserved.parent.mkdir(parents=True, exist_ok=True)
        preserved.to_csv(out_preserved, sep="\t", index=False)

    print(f"Input rows: {len(df)}")
    print(f"Wrote 3-column metadata: {out_3col}")
    if out_preserved is not None:
        print(f"Wrote preserved metadata: {out_preserved}")


if __name__ == "__main__":
    main()

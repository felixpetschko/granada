#!/usr/bin/env python3
import argparse
import time
from pathlib import Path

from redeconv.__ReDeconv_P import *


def main():
    parser = argparse.ArgumentParser(description="Run ReDeconv deconvolution steps on atlas-prepared inputs.")
    parser.add_argument(
        "--choice",
        type=int,
        required=True,
        choices=[1, 2, 3],
        help="1: initial signature genes, 2: signature matrix, 3: deconvolution",
    )
    args = parser.parse_args()

    base = Path(__file__).resolve().parent
    in_dir = base / "demo_data_4_deconvolution"
    out_dir = base / "Results_4_deconvolution"
    out_dir.mkdir(parents=True, exist_ok=True)

    fn_meta = str(in_dir / "References_Meta_data_subset.tsv")
    fn_exp = str(in_dir / "References_scRNA_seq_Nor.tsv")

    fn_ini_sig = str(out_dir / "Initial_sig_t_test_fd2.0_corr.tsv")
    fn_mean_std = str(out_dir / "Signature_mean_std_fd2.0.tsv")
    fn_heatmap = str(out_dir / "Heatmap_signature_gene_matrix.png")
    fn_extra_info = str(out_dir / "Signature_genes_extra_information.txt")

    fn_bulk_rnaseq_raw = str(in_dir / "Synthetic_Bulk_RNA_seq_Equal_Fraction_TPM.tsv")
    fn_percentage_save = str(out_dir / "ReDeconv_results.tsv")

    st_time = time.mktime(time.gmtime())

    if args.choice == 1:
        l_max_pv = 0.05
        l_min_fold_change = 2.0
        l_celltype_cellno_lb = 30
        l_nosep_sampleno_ub = 2
        l_status_data = check_meta_and_scRNAseq_data(fn_meta, fn_exp)
        if l_status_data > 0:
            get_initial_Signature_Candidates(
                fn_meta,
                fn_exp,
                fn_ini_sig,
                l_max_pv,
                l_min_fold_change,
                l_celltype_cellno_lb,
                l_nosep_sampleno_ub,
            )

    if args.choice == 2:
        l_top_no = 133
        Get_signature_gene_matrix(fn_exp, fn_meta, fn_ini_sig, fn_mean_std, l_top_no, fn_heatmap, fn_extra_info)

    if args.choice == 3:
        ReDeconv(fn_mean_std, fn_bulk_rnaseq_raw, fn_percentage_save)

    end_time = time.mktime(time.gmtime())
    print(f"Total time = {(end_time - st_time) / 60} minutes")


if __name__ == "__main__":
    main()

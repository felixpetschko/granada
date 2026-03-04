#!/usr/bin/env python3
import argparse
import time
from pathlib import Path

from redeconv.__ReDeconv_N import *


def main():
    parser = argparse.ArgumentParser(description="Run ReDeconv normalization steps on atlas-prepared inputs.")
    parser.add_argument(
        "--choice",
        type=int,
        required=True,
        choices=[1, 2, 3, 4],
        help="1: means/counts, 2: heatmap, 3: point-plot, 4: normalization",
    )
    args = parser.parse_args()

    base = Path(__file__).resolve().parent
    in_dir = base / "demo_data_4_normalization"
    out_dir = base / "Results_4_normalization"
    out_dir.mkdir(parents=True, exist_ok=True)

    fn_meta = str(in_dir / "Demo_ATLAS_meta.tsv")
    fn_exp = str(in_dir / "Demo_ATLAS_scRNAseq.tsv")

    fn_ctyp_mean = str(out_dir / "Ctype_size_means.tsv")
    fn_ctyp_count = str(out_dir / "Ctype_cell_counts.tsv")
    fn_cell_transcriptome_size = str(out_dir / "Cell_trans_sizes.tsv")

    st_time = time.mktime(time.gmtime())

    if args.choice == 1:
        l_cell_count_low_bound = 10
        l_status_data = check_meta_and_scRNAseq_data(fn_meta, fn_exp)
        if l_status_data > 0:
            get_sample_cell_type_exp_mean_and_cell_count(
                fn_meta,
                fn_exp,
                fn_ctyp_mean,
                fn_ctyp_count,
                fn_cell_transcriptome_size,
                l_cell_count_low_bound,
            )

    if args.choice == 2:
        fn_heatmap = str(out_dir / "Heatmap_plot.png")
        fn_heatmap_matrix = str(out_dir / "Heatmap_plot_correlation_matrix.csv")
        draw_heatmap_Pearson_all(fn_ctyp_mean, fn_heatmap, fn_heatmap_matrix)

    if args.choice == 3:
        fn_extra_info = str(out_dir / "Extra_information.txt")
        fn_point = str(out_dir / "Points_plot.png")
        l_figure_no_each_row = 2
        l_baseline = 0
        l_pearson_lb = 0.95
        get_sample_cell_type_information_top_Pearson_2(fn_ctyp_count, fn_ctyp_mean, fn_extra_info, l_pearson_lb)
        l_fit_line_with_shift = 1
        l_fit_line_no_shift = 1
        draw_cell_type_size_mean_point_plot(
            fn_ctyp_mean,
            fn_point,
            l_baseline,
            l_pearson_lb,
            l_figure_no_each_row,
            l_fit_line_with_shift,
            l_fit_line_no_shift,
        )

    if args.choice == 4:
        fn_meta_2 = str(out_dir / "Meta_data_new.tsv")
        fn_exp_2 = str(out_dir / "scRNA_seq_temp_file.tsv")
        fn_exp_3 = str(out_dir / "scRNA_seq_new_noShift.tsv")
        fn_exp_4 = str(out_dir / "scRNA_seq_new_withShift.tsv")
        l_pearson_lb = 0.95
        l_baseline = 0
        l_normalization_with_shift = "N"
        if l_normalization_with_shift != "Y":
            get_cell_subset_scRNA_seq_data_normalization_no_shift(
                fn_exp,
                fn_meta,
                fn_ctyp_mean,
                fn_cell_transcriptome_size,
                fn_exp_2,
                fn_meta_2,
                fn_exp_3,
                l_baseline,
                l_pearson_lb,
            )
        else:
            get_cell_subset_scRNA_seq_data_normalization(
                fn_exp,
                fn_meta,
                fn_ctyp_mean,
                fn_cell_transcriptome_size,
                fn_exp_2,
                fn_meta_2,
                fn_exp_4,
                l_baseline,
                l_pearson_lb,
            )

    end_time = time.mktime(time.gmtime())
    print(f"Total time = {(end_time - st_time) / 60} minutes")


if __name__ == "__main__":
    main()

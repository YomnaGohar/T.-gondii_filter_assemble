#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 15:56:25 2024

@author: yomna
"""


import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
def compute_depth_means(bam_file, window_size=100000):
    """
    Compute mean depth for each window in a BAM file.

    Parameters:
        bam_file (str): Path to the BAM file.
        window_size (int): Size of the window for calculating mean depth.

    Returns:
        list: List of mean depths for each window.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    mean_depths = []

    for contig in bam.references:
        contig_length = bam.get_reference_length(contig)
        for start in range(0, contig_length, window_size):
            end = min(start + window_size, contig_length)
            depths = []
            for pileupcolumn in bam.pileup(contig, start, end):
                depths.append(pileupcolumn.nsegments)
            if depths:
                mean_depths.append(np.mean(depths))

    bam.close()
    return mean_depths


def compute_and_plot_depth(bam_file_1, bam_file_2, labels, output_file_plot="output_violin.pdf", output_file_table="output_statistics.tsv", window_size=100000):
    """
    Compute mean depths for each window, plot violin plots for means, and save statistics table.

    Parameters:
        bam_file_1 (str): Path to the first BAM file.
        bam_file_2 (str): Path to the second BAM file.
        labels (list): Labels for the violin plots and statistics table.
        output_file_plot (str): File path to save the violin plot.
        output_file_table (str): File path to save the depth statistics table.
        window_size (int): Size of the window for calculating mean depth.
    """
    # Compute mean depths for each BAM file
    depth_means_1 = compute_depth_means(bam_file_1, window_size)
    depth_means_2 = compute_depth_means(bam_file_2, window_size)

    # Create violin plot
    plt.figure(figsize=(12, 6))
    plt.violinplot([depth_means_1, depth_means_2], showmeans=True, showextrema=True)
    plt.title('Violin Plot of Mean Depths Across Windows')
    plt.xlabel('BAM Files')
    plt.ylabel('Mean Depth')
    plt.xticks([1, 2], labels)
    plt.grid(True)
    plt.savefig(output_file_plot, format="pdf", bbox_inches="tight", dpi=100)
    plt.close()

    # Compute overall statistics for mean depths
    stats = pd.DataFrame([
        {
            "label": labels[0],
            "mean_of_means": np.mean(depth_means_1),
            "median_of_means": np.median(depth_means_1),
            "min_of_means": np.min(depth_means_1),
            "max_of_means": np.max(depth_means_1),
            "Q5_of_means": np.percentile(depth_means_1, 5),
            "Q95_of_means": np.percentile(depth_means_1, 95),
        },
        {
            "label": labels[1],
            "mean_of_means": np.mean(depth_means_2),
            "median_of_means": np.median(depth_means_2),
            "min_of_means": np.min(depth_means_2),
            "max_of_means": np.max(depth_means_2),
            "Q5_of_means": np.percentile(depth_means_2, 5),
            "Q95_of_means": np.percentile(depth_means_2, 95),
        }
    ])

    # Save statistics table
    stats.to_csv(output_file_table, sep="\t", index=False)


# Example usage
bam_path_1 = sys.argv[1]  # First BAM file
bam_path_2 = sys.argv[2]  # Second BAM file
output_file_plot = sys.argv[3]  # Output violin plot file
output_file_table = sys.argv[4]  # Output statistics table file
labels = ["Xie assembly", "Düsseldorf assembly"]

compute_and_plot_depth(bam_path_1, bam_path_2, labels, output_file_plot, output_file_table)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 08:22:56 2024

@author: yomna
"""
import matplotlib.pyplot as plt
import pandas as pd
import sys
import numpy as np
from matplotlib.gridspec import GridSpec
paf_file_my_assembly = sys.argv[1]
paf_file_xie= sys.argv[2]
output = sys.argv[3]
def plot_two_scatter_aligned_read_vs_identity(paf_file1, paf_file2, output_file=output):
    # Helper function to process PAF file
    def process_paf(file):
        columns_to_use = [0, 1, 2, 3, 5, 6, 7, 8, 9]
        paf = pd.read_csv(
            file,
            sep="\t",
            header=None,
            usecols=columns_to_use,
            names=["query_name", "query_length", "query_start", "query_end",
                   "target_name", "target_length", "target_start", "target_end",
                   "match_length"]
        )
        paf["mapped_read_length"] = paf["query_end"] - paf["query_start"]
        paf["read_identity"] = paf["match_length"] / paf["mapped_read_length"] * 100
        paf = paf[(paf["read_identity"] >= 0) & (paf["read_identity"] <= 100)]
        return paf

     # Process both PAF files
    paf1 = process_paf(paf_file1)
    paf2 = process_paf(paf_file2)

    # Create a grid layout for scatter plot and histograms
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(4, 4, fig)

    # Main scatter plot
    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_main.scatter(
        paf1["mapped_read_length"],
        paf1["read_identity"],
        alpha=0.6,
        s=10,
        color="blue",
        label="PAF File 1"
    )
    ax_main.scatter(
        paf2["mapped_read_length"],
        paf2["read_identity"],
        alpha=0.6,
        s=10,
        color="green",
        label="PAF File 2"
    )
    ax_main.set_xscale("log")
    ax_main.set_xlim(10, max(paf1["mapped_read_length"].max(), paf2["mapped_read_length"].max()))
    ax_main.set_ylim(0, 100)
    ax_main.set_xlabel("Aligned Read Length (bp, log scale)", fontsize=12)
    ax_main.set_ylabel("Read Identity (%)", fontsize=12)
    ax_main.legend()
    ax_main.grid(True)

    # Histogram on the x-axis
    ax_xhist = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_xhist.hist(
        [paf1["mapped_read_length"], paf2["mapped_read_length"]],
        bins=30,
        color=["blue", "green"],
        alpha=0.6,
        label=["PAF File 1", "PAF File 2"]
    )
    ax_xhist.set_ylabel("Frequency", fontsize=10)
    ax_xhist.legend()
    ax_xhist.grid(True)

    # Histogram on the y-axis
    ax_yhist = fig.add_subplot(gs[1:4, 3], sharey=ax_main)
    ax_yhist.hist(
        [paf1["read_identity"], paf2["read_identity"]],
        bins=30,
        orientation="horizontal",
        color=["blue", "green"],
        alpha=0.6
    )
    ax_yhist.set_xlabel("Frequency", fontsize=10)
    ax_yhist.grid(True)

    # Tight layout to align subplots
    fig.tight_layout()

    # Save plot to a PDF
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.show()
plot_two_scatter_aligned_read_vs_identity(paf_file_my_assembly, paf_file_xie)



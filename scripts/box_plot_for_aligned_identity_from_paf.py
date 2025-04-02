#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 18:22:49 2024

@author: yomna
"""


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
def read_paf(paf_file):
    """
    Reads a PAF file and calculates the aligned read identity.
    """
    columns_to_use = [0, 1, 2, 3, 5, 6, 7, 8, 9,10]
    paf = pd.read_csv(
        paf_file,
        sep="\t",
        header=None,
        usecols=columns_to_use,
        names=["query_name", "query_length", "query_start", "query_end",
               "target_name", "target_length", "target_start", "target_end",
               "match_length","Gapped_read_alignment"]
    )
    # Calculate mapped read identity
    #paf["mapped_read_length"] = paf["query_end"] - paf["query_start"]
    paf["read_identity"] = paf["match_length"] / paf["Gapped_read_alignment"] * 100
    return paf[["read_identity"]]

def plot_violin(paf_file1, paf_file2, labels, output_file=sys.argv[3]):
    """
    Creates a violin plot comparing the aligned read identity of two PAF files.
    """
    # Read and label data from the two PAF files
    paf1 = read_paf(paf_file1)
    paf1["Dataset"] = labels[0]

    paf2 = read_paf(paf_file2)
    paf2["Dataset"] = labels[1]

    # Combine the data into a single DataFrame
    combined_data = pd.concat([paf1, paf2])

    # Create the violin plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(
        x="Dataset",
        y="read_identity",
        data=combined_data,
        palette="Set2",
        cut=0  # Don't extend violins beyond the data range
    )
    plt.xlabel("Dataset", fontsize=12)
    plt.ylabel("Read Identity (%)", fontsize=12)
    plt.title("Violin Plot of Aligned Read Identity", fontsize=14)
    plt.grid(True, axis="y", linestyle="--", alpha=0.6)

    # Save the plot as a PDF
    plt.savefig(output_file, format="pdf", bbox_inches="tight", dpi=100)
    plt.show()

# Example usage
# Replace 'path/to/paf1.paf' and 'path/to/paf2.paf' with your actual PAF file paths
paf_file1 = sys.argv[1] #"path/to/paf1.paf"
paf_file2 =  sys.argv[2]#"path/to/paf2.paf"
labels = ["PAF File 1", "PAF File 2"]
plot_violin(paf_file1, paf_file2, labels)

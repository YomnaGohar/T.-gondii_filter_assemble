#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 17:52:14 2024

@author: yomna
"""

import pandas as pd
import sys

def process_paf_line(line):
    """
    Process a single PAF line and extract relevant fields.
    """
    fields = line.strip().split("\t")

    if len(fields) < 12:  # PAF lines must have at least 12 columns
        return None

    # Extract required fields
    query_name = fields[0]
    query_length = int(fields[1])
    query_start = int(fields[2])
    query_end = int(fields[3])
    target_name = fields[5]
    target_length = int(fields[6])
    target_start = int(fields[7])
    target_end = int(fields[8])
    residue_matches = int(fields[9])
    alignment_block_length = int(fields[10])
    mapping_quality = int(fields[11])

    # Extract optional CIGAR (cg:Z:) if present
    optional_fields = {field.split(":")[0]: field.split(":")[2] for field in fields[12:] if ":" in field}
    cigar = optional_fields.get("cg", None)

    return {
        "query_name": query_name,
        "query_length": query_length,
        "query_start": query_start,
        "query_end": query_end,
        "target_name": target_name,
        "target_length": target_length,
        "target_start": target_start,
        "target_end": target_end,
        "residue_matches": residue_matches,
        "alignment_block_length": alignment_block_length,
        "mapping_quality": mapping_quality,
        "cigar": cigar,
    }
import re

def parse_cigar(cigar_string):
    """
    Parse a CIGAR string from a PAF cg:Z: tag into a list of tuples.
    
    Parameters:
        cigar_string (str): CIGAR string from a PAF cg:Z: tag (e.g., "cg:Z:5=2I4X3=1D").
    
    Returns:
        list of tuples: Parsed CIGAR operations in the form [(code, length), ...].
    """

    # Use a regular expression to split the CIGAR string into (length, operation) pairs
    parsed_cigar = re.findall(r"(\d+)([=XIDSH])", cigar_string)
    # Convert the parsed operations into tuples of (code, length)
    cigartuples = [(int(length), op) for length, op in parsed_cigar]
    #print(cigartuples)
    return cigartuples
def adjusted_CIGAR(record):
    cigar = parse_cigar(record["cigar"])
    #print(cigar)
    # Add soft clipping if necessary
    if record["query_start"] > 0:
        clipping_length_start = record["query_start"]
        cigar.insert(0, (clipping_length_start, "S"))
    if record["query_end"] < record["query_length"]:
        clipping_length_end = record["query_length"] - record["query_end"]
        cigar.append((clipping_length_end, "S"))
    
    return cigar              

def check_clipping_perc(record):
    adju_cigar_1=record["cigar"]
    clipped_length=0
    for length, op in adju_cigar_1:
        if op == "S":
            clipped_length += length
    perc =  clipped_length / record["query_length"] *100   
    return perc

import numpy as np

def calculate_statistics(clipping_percentages):
    """
    Calculate various statistics for the clipping percentages.
    """
    stats = {
        "min": np.min(clipping_percentages),
        "q5": np.percentile(clipping_percentages, 5),
        "q25": np.percentile(clipping_percentages, 25),
        "mean": np.mean(clipping_percentages),
        "median": np.median(clipping_percentages),
        "q50": np.percentile(clipping_percentages, 50),
        "q75": np.percentile(clipping_percentages, 75),
        "q95": np.percentile(clipping_percentages, 95),
        "max": np.max(clipping_percentages),
    }
    return stats

def classify_positions_from_paf( paf_file_1, paf_file_2, output_file="", suspicious_threshold=0.25, erroneous_threshold=0.5):
    """
    Classify positions in a PAF file as Erroneous, Suspicious, or Fine for each BED position.

    Parameters:
        bed_file (str): Path to the BED file.
        paf_file (str): Path to the PAF file.
        output_file (str): Path to save the classified positions.
        identity_threshold (float): Minimum identity for reads to be considered correct.
        suspicious_threshold (float): Fraction of erroneous reads for suspicious classification.
        erroneous_threshold (float): Fraction of erroneous reads for erroneous classification.
    """
    # Load BED file
    results = {"old": [], "new": []}
    
    for assembly in ["old", "new"]:
        if assembly == "old":
            paf_file = paf_file_1
        else:
            paf_file = paf_file_2
            
        with open(paf_file) as paf:
            paf_lines = paf.readlines()

            # Process PAF records
            for line in paf_lines:
                record = process_paf_line(line)
                if not record:
                    continue
                adj = adjusted_CIGAR(record)
                record["cigar"] = adj
                clipping_perc = check_clipping_perc(record)
                results[assembly].append(clipping_perc)
    
    # Calculate statistics for each assembly
    statistics = {assembly: calculate_statistics(values) for assembly, values in results.items()}
    
    # Optionally save statistics to the output file
    if output_file:
        with open(output_file, "w") as f:
            for assembly, stats in statistics.items():
                f.write(f"{assembly} assembly statistics:\n")
                for stat_name, stat_value in stats.items():
                    f.write(f"{stat_name}: {stat_value:.2f}\n")
                f.write("\n")
    
    return statistics

# Example usage
import sys
paf_file_1 = sys.argv[1]  # First BAM file
paf_file_2 = sys.argv[2]  # Second BAM file
output_file = sys.argv[3]  # Output violin plot file

statistics = classify_positions_from_paf( paf_file_1, paf_file_2, output_file)

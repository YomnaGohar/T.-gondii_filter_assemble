#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:53:07 2024

@author: yomna
"""

import pandas as pd
import sys
import numpy as np
import pysam
def compute_depth_mean_in_window(bam_file, chrom, start, end):
    """
    Compute the mean depth for a specific window in a BAM file.

    Parameters:
        bam_file (str): Path to the BAM file.
        chrom (str): Chromosome name.
        start (int): Start position of the window (0-based).
        end (int): End position of the window (exclusive).

    Returns:
        float: Mean depth for the specified window, or None if no data is available.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    depths = []

    for pileupcolumn in bam.pileup(chrom, start, end):
        depths.append(pileupcolumn.nsegments)

    bam.close()

    if depths:
        return np.mean(depths)
    else:
        return None  # No data available for the specified window

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
def adjusted_CIGAR(record, start, end):
    adju_cigar = []
    cigar = parse_cigar(record["cigar"])
    #print(cigar)
    # Add soft clipping if necessary
    if record["query_start"] > 0:
        clipping_length_start = record["query_start"]
        cigar.insert(0, (clipping_length_start, "S"))
    if record["query_end"] < record["query_length"]:
        clipping_length_end = record["query_length"] - record["query_end"]
        cigar.append((clipping_length_end, "S"))
    
    ref_start = record["target_start"]
    ref_end = record["target_end"]
    window_start = max(start, ref_start)
    window_end = min(end, ref_end)
    adjusted_pos = ref_start  # Reference position
    
    for length, op in cigar:
        #print(adjusted_pos)
        #print(adju_cigar)
        if op in ["=", "X"]:  # Matches or mismatches
            #print( op in ["=", "X"])
            for _ in range(length):
                if window_start <= adjusted_pos < window_end:  # Inclusive of start and end
                    adju_cigar.append((1, op))
                adjusted_pos += 1
        elif op == "I":  # Insertion (does not consume reference)
            #print(op)
            if window_start <= adjusted_pos < window_end:  # Inclusive of start and end
               adju_cigar.append((length, op))
        elif op == "D":  # Deletion (consumes reference)
            #print(op)
            for _ in range(length):
                if window_start <= adjusted_pos < window_end:
                    adju_cigar.append((1, op))
                adjusted_pos += 1
        elif op == "S":  # Soft clipping (query only, no reference position)
            if window_start <= adjusted_pos < window_end:  # Inclusive of start and end
               adju_cigar.append((length, op))
    
    return adju_cigar              

def check_clipping_perc(record,start,end):
    adju_cigar_1=record["cigar"]
    clipped_length=0
    for length, op in adju_cigar_1:
        if op == "S":
            clipped_length += length
    perc =  clipped_length / record["query_length"]    *100   
    return perc
def aligned_read_identiy(record,start,end,chrom):
    #print(chrom,start,end)
    adju_cigar_1=record["cigar"]
    #print(adju_cigar_1)
    #print("adju_cigar_1", adju_cigar_1 )
    #print("")
    aligned_read_length = sum(length for length, op in adju_cigar_1 if op in ["X", "=", "I", "D"])
    matches=sum(length for length, op in adju_cigar_1 if op in [ "="])
    aligned_read_identity =  matches / aligned_read_length * 100     
    return aligned_read_identity

def classify_positions_from_paf(bed_file, paf_file_1,bam_file_1,paf_file_2,bam_file_2, output_file="", suspicious_threshold=0.25, erroneous_threshold=0.5):
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
    bed_positions = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom_xie", "start_xie", "end_xie","chrom_Düs", "start_Düs", "end_Düs"])
    results = {"old":[], "new":[]}
    for assembly in ["old","new"]:
        if assembly == "old":
            paf_file = paf_file_1
            bam_file = bam_file_1
            depth_threshold_upper=25
            depth_threshold_lower=15
            
        else:
            paf_file = paf_file_2
            bam_file = bam_file_2
            depth_threshold_upper=113
            depth_threshold_lower=54
            
        with open(paf_file) as paf:
            paf_lines = paf.readlines()
    
        for _, row in bed_positions.iterrows():
            if assembly == "old":
                chrom, start, end = row["chrom_xie"], row["start_xie"], row["end_xie"]
            else:
                chrom, start, end = row["chrom_Düs"], row["start_Düs"], row["end_Düs"]  
            # Filter PAF records overlapping the BED position
            depth =  compute_depth_mean_in_window(bam_file, chrom, start, end)
                
            relevant_records = []
            spaning_read=[]
            for line in paf_lines:
                record = process_paf_line(line)
                if not record:
                    continue
                # Check overlap
                if (
                    record["target_name"] == chrom and
                    record["target_end"] > start and
                    record["target_start"] < end
                ):
                    relevant_records.append(record)
                    if (
                        record["target_start"] < start and
                        record["target_end"] > end
                    ):  
                     spaning_read.append(record["query_name"])   
    
            if not relevant_records:
                results[assembly].append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "classification": "No Reads"
                })
                continue
    
            # Calculate classifications
            total_alignments = len(relevant_records)
            total_reads =len(list(set(record["query_name"]for record in relevant_records )))
            erroneous_reads_list=[]
            clipped= []
            low_identity= []
            for record in relevant_records:
                adj=adjusted_CIGAR(record, start,end)
                #print(record["cigar"])
                record["cigar"] = adj
                #check clippping of the read 
                clipping_perc= check_clipping_perc(record,start,end)
                if clipping_perc >=10:
                    erroneous_reads_list.append(record["query_name"])
                    clipped.append(record["query_name"])
                #check aligned identity for the given window  
                aligned_read_identity = aligned_read_identiy(record,start,end,chrom)
                if aligned_read_identity < 90:
                   erroneous_reads_list.append(record["query_name"])
                   low_identity .append(record["query_name"])
            erroneous_reads = list(set (erroneous_reads_list))
            erroneous_ratio = len(erroneous_reads) / total_reads
            spaning_read = len(list(set(spaning_read)))
            spaning_ratio = spaning_read / total_reads
            if spaning_ratio > 0.25:
                if erroneous_ratio > erroneous_threshold: #or (depth > depth_threshold_upper or depth < depth_threshold_lower):
                   classification = "Suspicious"   
                else:   
                    classification = "Fine"
            elif erroneous_ratio > erroneous_threshold: #or (depth > depth_threshold_upper or depth < depth_threshold_lower):
                classification = "Erroneous"
            elif suspicious_threshold < erroneous_ratio <= erroneous_threshold:
                classification = "Suspicious"
            else:
                classification = "Fine"
    
            results[assembly].append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "spaning_read":spaning_read,
                "classification": classification,
                 "total_reads": total_reads,
                "total_alignments": total_alignments,
                "erroneous_reads": len(erroneous_reads),
                "reads_that_has_more_that_10%_clippimg":len(list(set(clipped))),
                "less_than_90%_identity": len(list(set(low_identity)))
            })
    
    # Save results
    results_df_1 = pd.DataFrame(results["old"])
    results_df_2= pd.DataFrame(results["new"])
    results_df = pd.concat([results_df_1, results_df_2], axis=1)
    results_df.columns = ["chrom_xie", "start_xie", "end_xie","spaning_read_xie","classification_xie","total_reads_xie","total_alignments_xie","erroneous_reads_xie","clipped","low_identity","chrom_Düs", "start_Düs", "end_Düs","spaning_read_Düs","classification_Düs","total_reads_Düss","total_alignments_Düss","erroneous_reads_Düss","clipped","low_identity"]
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"Results saved to {output_file}")


# Example Usage
if __name__ == "__main__":
    classify_positions_from_paf(
        bed_file=sys.argv[1],
        paf_file_1=sys.argv[2],
        bam_file_1=sys.argv[3],
        paf_file_2=sys.argv[4],
        bam_file_2=sys.argv[5],
        output_file=sys.argv[6]
    )

# #ed_file="/home/yomna/hpc_project/combined_with_pipline_after_rebasecalling_and_herro/after_mit_removal/assembly_to_reference/nano-corr/collective_projected_positions.bed"
# #af_file_1="/home/yomna/Desktop/PhD_Yomna_Gohar/fasta/ME49_genome_from_Third_generation_sequencing_paper/alignment_SRR12363179.paf"
# paf_file_2="/home/yomna/hpc_project/combined_with_pipline_after_rebasecalling_and_herro/after_mit_removal/reads_to_assembly/nano-corr/mapped_reads.paf"
# classify_positions_from_paf( bed_file,paf_file_1,paf_file_2)

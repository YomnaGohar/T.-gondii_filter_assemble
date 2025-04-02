import pysam
import sys
def read_positions(file_path):
    """Read positions from a text file."""
    with open(file_path, 'r') as file:
        positions = [line.strip() for line in file if line.strip()]
    return positions

def check_contig_coverage(bam_path, positions, output_file):
    """Check and order contigs spanning the given positions, write results to a file."""
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    results = {}
    
    for position in positions:
        chr, start_end = position.split(':')[0:2]
        start, end = start_end.split('-')
        start = int(start.replace(',', ''))
        end = int(end.replace(',', ''))
        contigs = bamfile.fetch(chr, start, end)
        if position not in results:
            results[position] = []
        for contig in contigs:
            results[position].append((contig.reference_start, contig.query_name, contig.query_length))
    
    bamfile.close()
    
    # Write results to the specified output file
    with open(output_file, 'w') as file:
        for position in sorted(results):
            file.write(f"Position: {position}\n")
            sorted_contigs = sorted(results[position], key=lambda x: x[0])  # Sort by start position
            for start_pos, name, length in sorted_contigs:
                file.write(f"  Contig {name} starts at {start_pos} with length {length}\n")

    return results

# Usage example
positions = read_positions(sys.argv[1])
bam_path = sys.argv[2]
output_file = sys.argv[3]
contig_order = check_contig_coverage(bam_path, positions, output_file)


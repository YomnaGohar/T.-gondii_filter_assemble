import pysam
import sys

def read_positions(file_path):
    """Read positions from a text file."""
    with open(file_path, 'r') as file:
        positions = [line.strip() for line in file if line.strip()]
    return positions
def check_contig_coverage(bam_path, positions):
    """Check if there is a single contig spanning the given positions."""
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    results = {}
    for position in positions:
        chr, start_end = position.split(':')[0:2]
        start, end = start_end.split('-')
        
        # Remove commas before converting to integers
        start = int(start.replace(',', ''))
        end = int(end.replace(',', ''))

        contigs = bamfile.fetch(chr, start, end)
        contig_names = set()

        for contig in contigs:
            if contig.reference_start <= start and contig.reference_end >= end:
                contig_names.add(contig.query_name)
        
        if len(contig_names) == 1:
            results[position] = f"Covered by single contig: {list(contig_names)[0]}"
        else:
            results[position] = "No single contig covers the range"
    
    bamfile.close()
    return results

# Path to your positions text file
positions_file = sys.argv[1]

# Path to your BAM file
bam_path = sys.argv[2]

# Path to the output file
output_file = sys.argv[3]

# Read positions from the file
positions = read_positions(positions_file)

# Check contig coverage
coverage_results = check_contig_coverage(bam_path, positions)

# Write results to the output file
with open(output_file, 'w') as f:
    for position, result in coverage_results.items():
        f.write(f"{position} : {result}\n")

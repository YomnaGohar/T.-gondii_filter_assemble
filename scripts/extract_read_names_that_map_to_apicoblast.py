import sys

# Check command-line arguments
if len(sys.argv) < 3:
    print("Usage: python extract_read_names.py <input.paf> <output.txt>")
    sys.exit(1)

input_paf = sys.argv[1]  # Input PAF file
output_txt = sys.argv[2]  # Output text file with read names

# Contigs of interest
target_contigs = {"contig_43", "contig_24"}

# Set to store unique read names
read_names = set()

# Process the PAF file
with open(input_paf, 'r') as paf_file:
    for line in paf_file:
        columns = line.strip().split('\t')
        if len(columns) > 5:
            read_name = columns[0]
            contig_name = columns[5]
            # Check if contig matches target
            if contig_name in target_contigs:
                read_names.add(read_name)

# Write the unique read names to the output text file
with open(output_txt, 'w') as output_file:
    for read_name in sorted(read_names):
        output_file.write(f"{read_name}\n")

print(f"Extracted {len(read_names)} read names to {output_txt}")


from Bio import SeqIO
import sys

# Check command-line arguments
if len(sys.argv) < 3:
    print("Usage: python filter_fastq_by_length.py <input.fastq> <output.fastq>")
    sys.exit(1)

input_fastq = sys.argv[1]  # Input FASTQ file
output_fastq = sys.argv[2]  # Output filtered FASTQ file

# Define length threshold
length_threshold = 35000

# Filter and write sequences shorter than the threshold
with open(output_fastq, 'w') as output_handle:
    for record in SeqIO.parse(input_fastq, "fastq"):
        if len(record.seq) < length_threshold:
            SeqIO.write(record, output_handle, "fastq")

print(f"Filtered sequences shorter than {length_threshold} bases written to {output_fastq}")


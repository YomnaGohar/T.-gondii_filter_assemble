import matplotlib.pyplot as plt
from Bio import SeqIO

def plot_read_length_distribution(fastq_file, output_image):
    read_lengths = []

    # Read the FASTQ file and calculate the length of each sequence
    for record in SeqIO.parse(fastq_file, "fastq"):
        read_lengths.append(len(record.seq))

    # Find the maximum read length
    max_read_length = max(read_lengths)

    # Plot the histogram of read lengths
    plt.hist(read_lengths, bins=50, color='blue', edgecolor='black')
    plt.title("Read Length Distribution")
    plt.xlabel("Read Length")
    plt.ylabel("Frequency")
    plt.grid(axis='y', alpha=0.75)

    # Add a vertical line at the maximum read length
    plt.axvline(max_read_length, color='red', linestyle='dashed', linewidth=1, label=f'Max Length: {max_read_length}')

    # Add legend
    plt.legend()

    # Save the plot as an image
    plt.savefig(output_image)
    # plt.show()

if __name__ == "__main__":
    # Specify the input FASTQ file and output image file
    import sys
    fastq_file = sys.argv[1]
    output_image = sys.argv[2]

    # Generate the histogram
    plot_read_length_distribution(fastq_file, output_image)


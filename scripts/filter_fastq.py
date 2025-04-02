
"""                                                                                 
%prog some.fasta wanted-list.txt
for seqIO documentation visit https://biopython.org/wiki/SeqIO                                                    
"""       
from Bio import SeqIO                                                               
import sys                                                                          
                                                                                    
wanted = [line.strip() for line in open(sys.argv[2])]                               
seqiter = SeqIO.parse(open(sys.argv[1]), 'fastq')                                    
SeqIO.write((seq for seq in seqiter if seq.id in wanted), sys.argv[3], "fastq")

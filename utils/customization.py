#!/bin/python3

import sys
import pandas as pd
from Bio import SeqIO
from pysam import FastaFile

options=sys.argv
bedfile=options[1]
fastafile=options[2]
output=options[3]

bed=pd.read_table(bedfile, header=None)
chrom = str(bed.iloc[0,0])
start = int(bed.iloc[0,1])
end = int(bed.iloc[0,2])
print("Extracting reads from ",chrom, start, end)

try:
  sequences_object = FastaFile(fastafile)
  tmp_sequence = sequences_object.fetch(chrom, start, end)
  #print(tmp_sequence)
except:
  print("ERROR : FASTA contig name does not match with BED file")  
  quit()


sequence = ("N"*start) + tmp_sequence + ("N"*(sequences_object.get_reference_length(chrom)-end))
f = open(output, 'w')
f.write(">"+chrom+'\n')
f.write(sequence+'\n')
f.close()
print("Finished creating "+output)  

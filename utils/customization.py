#!/bin/python3
# Customization Script: Region-specific FASTA sequence extraction and padding

# This script is designed to extract a specific genomic region defined by a BED file
# from a reference FASTA file. It then generates a FASTA output file where the
# extracted region is present, and all other parts of the chromosome are padded
# with 'N' (ambiguous nucleotide) characters. This can be useful for creating
# custom reference sequences for specific alignment scenarios, such as aligning
# a gene and its pseudogene within a defined chromosomal context.


import sys
import pandas as pd
from Bio import SeqIO
from pysam import FastaFile
import os 

# --- 1. Argument Parsing ---
# Expects three command-line arguments:
# sys.argv[1]: Path to the BED file (e.g., specifying gene/pseudogene coordinates)
# sys.argv[2]: Path to the comprehensive reference FASTA file
# sys.argv[3]: Desired output FASTA file name
options=sys.argv
bedfile=options[1]
fastafile=options[2]
output=options[3]

# --- 2. Extract Coordinates from BED File ---
# Reads the first line of the BED file (assuming only one region is defined).
# The BED file is expected to be tab-separated with no header.
bed=pd.read_table(bedfile, header=None)
chrom = str(bed.iloc[0,0])
start = int(bed.iloc[0,1])
end = int(bed.iloc[0,2])
print("Extracting reads from ",chrom, start, end)

# --- 3. Fetch Sequence from Reference FASTA ---
try:
  # Opens the reference FASTA file using pysam.FastaFile for efficient random access.
  sequences_object = FastaFile(fastafile)
  # Fetches the sequence for the specified chromosome and coordinates.
  tmp_sequence = sequences_object.fetch(chrom, start, end)
# --- 4. Error Handling: Contig Name Mismatch ---
# Catches exceptions if the chromosome name from the BED file is not found in the FASTA.
except ValueError as e: # Catch specific ValueError for invalid contig
  print(f"ERROR: FASTA contig name '{chrom}' does not match with BED file or coordinates are invalid. Details: {e}")
  quit()
except Exception as e: # Catch any other unexpected errors
  print(f"An unexpected error occurred: {e}")
  quit()

# --- 5. Generate Padded Sequence for Output ---
# This is a key step: it constructs a full-chromosome-length sequence.
# It pads the regions *before* and *after* the extracted segment with 'N's.
# This results in a FASTA entry representing the entire chromosome,
# where only the specified region has its actual sequence, and the rest is 'N's.
chromosome_length = sequences_object.get_reference_length(chrom)
sequence = ("N"*start) + tmp_sequence + ("N"*(chromosome_length-end))


# --- 6. Write Output FASTA File ---
# Writes the generated (and N-padded) sequence to the specified output FASTA file.
f = open(output, 'w')
f.write(">"+chrom+'\n')
f.write(sequence+'\n')
f.close()

print(f"Finished creating {output}")

# --- 7. Generate BWA and Samtools Indexes for the New FASTA File ---
# These index files are necessary for BWA to align reads against this custom reference
# and for samtools/GATK to access specific regions of the FASTA.
print(f"Generating BWA index for {output}...")
# bwa index command: creates .amb, .ann, .bwt, .pac, .sa files.
# 'output' variable holds the path to the newly created FASTA file.
bwa_index_command = f"bwa index {output}"
os.system(bwa_index_command) # Executes the BWA index command

print(f"Generating Samtools FASTA index for {output}...")
# samtools faidx command: creates a .fai index file.
samtools_faidx_command = f"samtools faidx {output}"
os.system(samtools_faidx_command) # Executes the samtools faidx command

print("Indexing complete.")

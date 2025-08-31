#!/bin/python3

# PseudoRA: Core Python Script for Pseudogene ReAlignment
# This script processes NGS alignment data (BAM), variant calls (VCF),
# and performs sophisticated haplotype phasing to reassign reads to their
# true genomic origin (e.g., SBDS vs. SBDSP1), correcting pseudoalignment issues.

import pysam
import sys
import pandas as pd
import numpy as np
from Bio import pairwise2
import warnings, random

# --- 1. Initialization and Input Processing ---
# Retrieves command-line arguments: [VCF/caller_file] [input_bam_file] [output_prefix] [bed_file]
options=sys.argv
f=open(options[1],'r')
counter=0
for line in f.readlines():
  if not line.startswith("#"):
    break
  counter+=1
f.close()
# Loads the VCF data (from the first GATK HaplotypeCaller run) into a pandas DataFrame.
# Filters for single nucleotide variants (SNV) with sufficient depth (>=20) as phasing markers.
caller=pd.read_table(options[1], skiprows=counter-1)
caller[['GT','AD','DP','GQ','PL']] = caller.iloc[:,-1].str.split(":", expand=True)
caller = caller.loc[(caller.REF.str.len()==1)&(caller.ALT.str.len()==1)&(caller.DP.astype(int)>=20)]
caller = caller.reset_index(drop=True)

# Opens the input BAM file (the unified SBDS/SBDSP1 alignment from BWA).
samfile = pysam.AlignmentFile(options[2], "rb" )
output = options[3]
bed = pd.read_table(options[4],header=None)

# Concatenates reference and alternate alleles to create full reference and alternate haplotypes.
# These will be used for scoring during read assignment.
ref_phase = ''.join(caller["REF"])
alt_phase = ''.join(caller["ALT"])

# --- 2. Read Information Extraction and Table Construction ---
# This section directly processes the input BAM file to build a detailed table of read sequences.
# It iterates through genomic positions, extracts nucleotide bases from aligned reads,
# and organizes them into a pandas DataFrame (table).
# This 'table' represents the input for the phasing logic, mapping reads to their observed SNP sequences.
print("=================================================")
print("Phasing all reads")
list_site=[]
tables = dict()

# Iterates through pileup columns for specified genomic region to collect base calls.
caller_index=0
for pileupcolumn in samfile.pileup(str(caller.iloc[0]['#CHROM']), int(caller.iloc[0]['POS'])-10, int(caller.iloc[-1]['POS'])+10, max_depth=100000000, truncate=True, min_base_quality=0, ignore_orphans=False):
  if pileupcolumn.pos+1 == int(caller.loc[caller_index,'POS']):
    for pileupread in pileupcolumn.pileups:
     if not pileupread.is_del and not pileupread.is_refskip: 
      if pileupread.alignment.query_name in tables.keys():
        if len(tables[pileupread.alignment.query_name]) <= len(list_site):
          tables[pileupread.alignment.query_name].append(pileupread.alignment.query_sequence[pileupread.query_position])
      else:
        tables[pileupread.alignment.query_name]=["-"]*len(list_site)
    # Records the current variant site.
    list_site.append(str(pileupcolumn.pos+1)+caller.loc[caller_index,'REF']+caller.loc[caller_index,'ALT'])
    # Fills in gaps for reads that don't cover the current site.
    for read_name in tables.keys():
      if len(tables[read_name]) != len(list_site):
        tables[read_name].append("-")
    caller_index+=1
    if(caller_index>=len(caller)): break

# Converts the collected read data into a pandas DataFrame (the 'table' for phasing).
table = pd.DataFrame.from_dict(tables, orient='index', columns=list_site)
table = table.sort_values(list(table.columns), ascending=False)
# Concatenates all SNP bases for each read into a single string, representing the read's observed haplotype sequence.
table['concat']=table[table.columns].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
# Saves the constructed table to a TSV file for potential debugging or further analysis.
table.to_csv(output+".tsv",sep='\t')

# --- 3. Sequence Alignment Scoring Functions ---
# These functions define custom scoring matrices for Biopython's pairwise alignment.
# Scores matches/mismatches with gaps receiving zero. Used for general sequence similarity.
def matrix_get(nuc1,nuc2):
   if nuc1=="-" or nuc2=="-":
     return 0
   elif nuc1==nuc2:
     return 1
   else:
     return -1

# Scores only non-gap matches. Used to count comparable bases for normalization.
def matrix_get2(nuc1,nuc2):
   if nuc1=="-" or nuc2=="-":
     return 0
   else:
     return 1

# Combines two sequences, prioritizing non-gap bases from seq1, or taking from seq2 if seq1 has a gap.
# Used to build consensus sequences during haplotype extension.
def combine_two_seqs(seq1, seq2):
  seq3=''
  for i in range(0,len(seq1)):
    if seq1[i]=="-":
      seq3+=seq2[i]
    else:
      seq3+=seq1[i]
  return seq3

# --- 4. Exon Segmentation and Haplotype Phasing Logic ---
# Divides the main table into 'exons' based on non-gap patterns.
# Phasing is then performed independently within each exon to build haplotype groups.
# This handles potential structural variations or discontinuities between regions.
columns_to_drop_all_gap = [col for col in table.columns if (table[col] == '-').all()]
if len(columns_to_drop_all_gap) > 0:
    print(f"Dropping {len(columns_to_drop_all_gap)} columns with all '-' values.")
    table = table.drop(columns=columns_to_drop_all_gap, axis=1)
else:
    print("No columns with all '-' values found to drop.")
table2 = table.reset_index()
exons = []
before = 0
# Iterates through columns (variant sites) to identify boundaries of contiguous non-gapped regions (exons).
for i in range(2,len(table.columns)-1):
  # Finds last non-gap position in previous column and first non-gap in current column.
  p1 = table2.loc[table2.iloc[:,i-1]!="-"].index[-1]
  p2 = table2.loc[table2.iloc[:,i]!="-"].index[0] 
  if p2 > p1:
    exons.append(table.iloc[before:p2])
    before=p2
exons.append(table.iloc[before:len(table)])

# --- 5. Haplotype Inference and Read Assignment ---
# This is the core logic for probabilistic haplotype phasing using pairwise alignment.
# It assigns reads to existing haplotypes or creates new ones based on sequence similarity.
phases=[]
for exon in exons:
 if not exon.empty:
  # Initialize phasing for the exon with the first read's concatenated sequence.
  phase={exon.iloc[0]['concat']:[exon.index[0]]}
  # Iterate through remaining reads in the exon.
  for i in range(1, len(exon)):
    found=[]
    cand=[]
    seq2=exon.iloc[i]['concat']
    # Compare current read (seq2) to all established haplotypes (seq1) in the 'phase' dictionary.
    for j in list(phase.keys()):
      seq1=j
      # Perform global alignments using custom scoring functions.
      alignments = pairwise2.align.globalcs(seq1, seq2,matrix_get,-10,-10) # Measures overall similarity.
      compare1 = pairwise2.align.globalcs(seq1, seq2,matrix_get2,-10,-10) # Counts comparable bases.
      if compare1[0].score == 0: # If no common comparable bases, skip comparison.
        continue
      # Determine if it's an exact match or a partial candidate.
      if alignments[0].score == compare1[0].score: # Exact match implies all comparable bases match perfectly.
        found.append(seq1)
      elif alignments[0].score > -1*(compare1[0].score): # Partial match: enough similarity, but not perfect.
        cand.append(seq1)
    # Decision logic for assigning current read to a haplotype.
    if len(found) == 0 and len(cand) == 0: # No match, create a new haplotype.
      phase[seq2]=[exon.index[i]]
    elif len(found) == 0 and len(cand) == 1: # Only one partial candidate, extend it.
       ind_blank=0 # Find the first non-gap position
       for char1 in range(0,len(seq2)):
         if seq2[char1]!="-":
           ind_blank=char1
           break
       # Combine: take prefix from candidate, suffix from current read.
       combined=cand[0][:char1]+seq2[char1:]
       phase[combined]=[exon.index[i]] ################################################
    elif len(found) == 0 and len(cand) > 1: # Multiple partial candidates, choose the best one (greedy).
      cand1=dict()
      for cand_curr in cand: 
       for char1 in range(0,len(seq2)):
         if seq2[char1]=="-":
           continue
         elif seq2[char1]!=cand_curr[char1]:
           cand1[cand_curr]=char1
           break
       # Select candidate with longest identical prefix (i.e., highest mismatch position).
       max_seq = max(cand1, key=cand1.get)
       max_index = cand1[max_seq]
       combined=max_seq[:max_index]+seq2[max_index:]
       # Replace old haplotype if sequence changed, assign read.
       phase[combined]=[exon.index[i]] ########################################      
    else: # If exact matches are found, randomly assign to one and update.
      seq1 = random.choice(found) # Randomly choose one exact match.
      combined = combine_two_seqs(seq1, seq2) # Merge sequences.
      if combined != seq1: # If the sequence actually changed (e.g., added non-gap from seq2)
        res = dict()
        for key in phase:
          if key == seq1:
            res[combined] = phase[key]+[exon.index[i]]
          else: 
            res[key] = phase[key]
        phase=res.copy()
      else: # If sequence did not change, just add read to existing haplotype.
        phase[combined]+=[exon.index[i]]
  phases.append(phase) # Store phasing results for the current exon.

# --- 6. Read Selection for Final Output BAM ---
# Based on the derived haplotypes, reads are selected for the final corrected BAM file.
# This selection prioritizes reads aligning best to the reference/alternate phases and total read coverage.
print("Writing corrected BAM file")
#warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 
#warnings.filterwarnings("ignore", category=DeprecationWarning)
df_interprets=[]
final_list=[]
for index in range(0, len(phases)):
   # Create a DataFrame to interpret phasing results (scores, read counts).
   df_interpret=pd.DataFrame(index=phases[index].keys(), columns=['read_names','ref_score','alt_score','score','reads','checked'])
   # Calculate scores for each derived haplotype against the reference and alternate phases.
   for phase in phases[index].keys():
     ref_score=float(pairwise2.align.globalcs(ref_phase, phase,matrix_get,-1,-1)[0].score)
     alt_score=float(pairwise2.align.globalcs(alt_phase, phase,matrix_get,-1,-1)[0].score)
     df_interpret.loc[phase] = [phases[index][phase], ref_score, alt_score, ref_score-alt_score, len(phases[index][phase]),'']#score, depth
     #df_interpret = df_interpret.sort_values(by=['ref_score', 'alt_score'], ascending=[False,True])
   # Sort haplotypes by 'score' (preference for ref_phase over alt_phase).
   df_interpret = df_interpret.sort_values(by=['score'], ascending=[False])
   # Calculate total reads and normalized total reads for the gene.
   sum_total = df_interpret.reads.sum()
   sum_total_by_genes = sum_total/len(bed)
   sum_now=0
   # Iterate through sorted haplotypes to select reads for final list.
   # This logic aims to select reads that best represent the dominant haplotype(s) while covering sufficient depth.
   for i in range(0, len(df_interpret)):
     if i==0:
       final_list += df_interpret.iloc[i]['read_names']
       df_interpret.loc[df_interpret.index[i], 'checked']="O"
       continue
     # Calculate cumulative sum of reads.
     sum_now += int(df_interpret.iloc[i-1]['reads'])
     sum_later = sum_now + int(df_interpret.iloc[i]['reads'])
     # Logic to select additional haplotypes based on coverage threshold (sum_total_by_genes).
     if (sum_now-sum_total_by_genes)<0 and (sum_later-sum_total_by_genes)<0:
       # If both cumulative sums are below target, include current haplotype.
       final_list += df_interpret.iloc[i]['read_names']
       df_interpret.loc[df_interpret.index[i], 'checked']="O"
     elif (sum_now-sum_total_by_genes)<=0 and (sum_later-sum_total_by_genes)>0 and abs(sum_now-sum_total_by_genes)>abs(sum_later-sum_total_by_genes):
       # If crossing threshold, choose based on which side is closer.
       final_list += df_interpret.iloc[i]['read_names']
       df_interpret.loc[df_interpret.index[i], 'checked']="O"
     else: # Both positive, implies threshold met/exceeded in previous steps.
       break
   df_interprets.append(df_interpret) # Store interpretation for current exon.

# --- 7. Output Result Files and Corrected BAM ---
# Writes interpretation results to TSV files and generates the final corrected BAM file.
for i in range(0,len(df_interprets)):
  if i==0:
    df_interprets[i].to_csv(output+".result.tsv", sep='\t')
  else:
    df_interprets[i].to_csv(output+".result.tsv", mode="a", sep='\t', header=False)

samfile.close()
# Re-open original BAM file to filter and write selected reads (more efficient).
samfile = pysam.AlignmentFile(options[2], "rb" )
outfile = pysam.AlignmentFile(output+".correct.bam", "wb", template=samfile)
for s in samfile:
  if s.query_name in final_list:
    outfile.write(s)

print("=================================================")
outfile.close()

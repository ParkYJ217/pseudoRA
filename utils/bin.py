#!/bin/python3

import pysam
import sys
import pandas as pd
import numpy as np
from Bio import pairwise2
import warnings, random

options=sys.argv
print("=================================================")
print("Phasing all reads")
f=open(options[1],'r')
counter=0
for line in f.readlines():
  if not line.startswith("#"):
    break
  counter+=1
f.close()
caller=pd.read_table(options[1], skiprows=counter-1)
samfile = pysam.AlignmentFile(options[2], "rb" )
output = options[3]
bed = pd.read_table(options[4],header=None)

caller[['GT','AD','DP','GQ','PL']] = caller.iloc[:,-1].str.split(":", expand=True)
caller = caller.loc[(caller.REF.str.len()==1)&(caller.ALT.str.len()==1)&(caller.DP.astype(int)>=20)]
caller = caller.reset_index(drop=True)

ref_phase = ''.join(caller["REF"])
alt_phase = ''.join(caller["ALT"])

list_site=[]
tables = dict()
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
    list_site.append(str(pileupcolumn.pos+1)+caller.loc[caller_index,'REF']+caller.loc[caller_index,'ALT'])
    for read_name in tables.keys():
      if len(tables[read_name]) != len(list_site):
        tables[read_name].append("-")
    caller_index+=1
    if(caller_index>=len(caller)): break

table = pd.DataFrame.from_dict(tables, orient='index', columns=list_site)
table = table.sort_values(list(table.columns), ascending=False)
table['concat']=table[table.columns].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
#print(table)
table.to_csv(output+".tsv",sep='\t')

##################################
def matrix_get2(nuc1,nuc2):
   if nuc1=="-" and nuc2=="-":
     return 0
   elif nuc1=="-" and nuc2!="-":
     return 0
   elif nuc1!="-" and nuc2=="-":
     return 0
   elif nuc1==nuc2:
     return 1
   else:
     return 1

table2 = table.reset_index()
exons = []
before = 0
for i in range(2,len(table.columns)-1):
  p1 = table2.loc[table2.iloc[:,i-1]!="-"].index[-1]
  p2 = table2.loc[table2.iloc[:,i]!="-"].index[0] 
  if p2 > p1:
    exons.append(table.iloc[before:p2])
    before=p2

exons.append(table.iloc[before:p2])

#################################
def matrix_get(nuc1,nuc2):
   if nuc1=="-" or nuc2=="-":
     return 0
   elif nuc1==nuc2:
     return 1
   else:
     return -1

def matrix_get2(nuc1,nuc2):
   if nuc1=="-" or nuc2=="-":
     return 0
   else:
     return 1

def combine_two_seqs(seq1, seq2):
  seq3=''
  for i in range(0,len(seq1)):
    if seq1[i]=="-":
      seq3+=seq2[i]
    else:
      seq3+=seq1[i]
  return seq3


phases=[]
for exon in exons:
  phase={exon.iloc[0]['concat']:[exon.index[0]]}
  for i in range(1, len(exon)):
    found=[]
    cand=[]
    seq2=exon.iloc[i]['concat']
    for j in list(phase.keys()):
      seq1=j
      alignments = pairwise2.align.globalcs(seq1, seq2,matrix_get,-10,-10) #check if bps are same
      compare1 = pairwise2.align.globalcs(seq1, seq2,matrix_get2,-10,-10) #check how many bps to compare
      if compare1[0].score == 0: #if no bp to compare
        continue
      if alignments[0].score == compare1[0].score: #if exact match, add to found
        found.append(seq1)
      elif alignments[0].score > -1*(compare1[0].score): # if partial match, add to cand
        cand.append(seq1)
    #
    if len(found) == 0 and len(cand) == 0: # if not found, add new phase
      phase[seq2]=[exon.index[i]]
    elif len(found) == 0 and len(cand) == 1: # if only 1 cand found, combine line
       ind_blank=0
       for char1 in range(0,len(seq2)):
         if seq2[char1]!="-":
           ind_blank=char1
           break
       combined=cand[0][:char1]+seq2[char1:]
       phase[combined]=[exon.index[i]]
    elif len(found) == 0 and len(cand) > 1: # if more than 1 cand found, find candidate with maximum match, then combine line
      cand1=dict()
      for cand_curr in cand: 
       for char1 in range(0,len(seq2)):
         if seq2[char1]=="-":
           continue
         elif seq2[char1]!=cand_curr[char1]:
           cand1[cand_curr]=char1
           break
       max_seq = max(cand1, key=cand1.get)
       max_index = cand1[max_seq]
       combined=max_seq[:max_index]+seq2[max_index:]
       phase[combined]=[exon.index[i]]       
    else: #if found is more than 1, randomly put in one of the founds
      seq1 = random.choice(found)
      combined = combine_two_seqs(seq1, seq2)
      if combined != seq1:
        res = dict()
        for key in phase:
          if key == seq1:
            res[combined] = phase[key]+[exon.index[i]]
          else:
            res[key] = phase[key]
        phase=res.copy()
      else:
        phase[combined]+=[exon.index[i]]
  phases.append(phase)

##################################

print("Writing corrected BAM file")
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 
df_interprets=[]
final_list=[]
for index in range(0, len(phases)):
   df_interpret=pd.DataFrame(index=phases[index].keys(), columns=['read_names','ref_score','alt_score','score','reads','checked'])
   for phase in phases[index].keys():
     ref_score=float(pairwise2.align.globalcs(ref_phase, phase,matrix_get,-1,-1)[0].score)
     alt_score=float(pairwise2.align.globalcs(alt_phase, phase,matrix_get,-1,-1)[0].score)
     df_interpret.loc[phase] = [phases[index][phase], ref_score, alt_score, ref_score-alt_score, len(phases[index][phase]),'']#score, depth
     #df_interpret = df_interpret.sort_values(by=['ref_score', 'alt_score'], ascending=[False,True])
     df_interpret = df_interpret.sort_values(by=['score'], ascending=[False])
   sum_total = df_interpret.reads.sum()
   sum_total_by_genes = sum_total/len(bed)
#   sum_now = int(df_interpret.iloc[0]['reads'])
#   final_list += df_interpret.iloc[0]['read_names']
#   df_interpret.iloc[0]['checked']="O"
   sum_now=0
   for i in range(0, len(df_interpret)):
     if i==0:
       final_list += df_interpret.iloc[i]['read_names']
       df_interpret.iloc[i]['checked']="O"
       continue
     sum_now += int(df_interpret.iloc[i-1]['reads'])
     sum_later = sum_now + int(df_interpret.iloc[i]['reads'])
     if (sum_now-sum_total_by_genes)<0 and (sum_later-sum_total_by_genes)<0:
       final_list += df_interpret.iloc[i]['read_names']
       df_interpret.iloc[i]['checked']="O"
     elif (sum_now-sum_total_by_genes)<=0 and (sum_later-sum_total_by_genes)>0 and abs(sum_now-sum_total_by_genes)>abs(sum_later-sum_total_by_genes):
       final_list += df_interpret.iloc[i]['read_names']
       df_interpret.iloc[i]['checked']="O"
     else: #BOTH positives
       break
#   df_real_gene = df_interpret.loc[df_interpret['ref_score']>=0]
#   for i in df_real_gene.index:
#     final_list += df_real_gene.loc[i,'read_names']
   df_interprets.append(df_interpret)

for i in range(0,len(df_interprets)):
  if i==0:
    df_interprets[i].to_csv(output+".result.tsv", sep='\t')
  else:
    df_interprets[i].to_csv(output+".result.tsv", mode="a", sep='\t', header=False)

samfile.close()
samfile = pysam.AlignmentFile(options[2], "rb" )
outfile = pysam.AlignmentFile(output+".correct.bam", "wb", template=samfile)
for s in samfile:
  if s.query_name in final_list:
    outfile.write(s)
print("=================================================")
outfile.close()
#pysam.index(output+".correct.bam")

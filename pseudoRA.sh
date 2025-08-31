#!/bin/bash

# PseudoRA: A tool for realigning reads in gene-pseudogene pairs to improve variant calling accuracy.
# Version: 1.0

# --- Configuration and Argument Parsing ---
# Set script-related paths
script=$(realpath "${BASH_SOURCE:-$0}")
scriptDIR=$(dirname $script)

# Default paths for reference files. These can be overridden by command-line arguments.
RefFasta=$scriptDIR/reference/chr7_b37_SBDS.fa
roi=$scriptDIR/reference/SBDS.bed

# Define usage instructions for the script
usage="
Program: pseudoRA (Tool for realigning reads in gene-pseudogene pair)
Version: 1.0

$(basename "$0") [options]

where:
    -h  show this help text
    -t  test with demo/demo.bam
    -i  destination of original BAM file [REQUIRED]
    -r  reference FASTA [OPTIONAL]
    -b  region of interest BED [OPTIONAL]
"

# Parse command-line arguments
while getopts ':hti:r:b:' option; do
  echo $option
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    t) ID=${scriptDIR}/demo/demo
       ;;
    i) ID_candidate=$OPTARG
       if [[ $ID_candidate != *".bam" ]]; then
         echo "Incorrect BAM file.."
         echo "Exitting.."
         exit 1
       fi
       ID=${ID_candidate/.bam/}
       ;;
    r) RefFasta=$OPTARG
       if [[ $RefFasta != *".fa" ]]; then
         echo "Incorrect reference file.."
         echo "Exitting.."
         exit 1
       fi       
       ;;
    b) roi=$OPTARG
       if [[ $roi != *".bed" ]]; then
         echo "Incorrect region of interest BED file.."
         echo "Exitting.."
         exit 1
       fi       
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
if [ $OPTIND -eq 1 ]; then echo "$usage"; exit; fi
shift $((OPTIND - 1))

# --- PseudoRA Pipeline Execution ---

# 1. Sample Filter: Extract reads mapping to the defined region of interest (ROI).
# This focuses the analysis on SBDS and SBDSP1 regions, reducing processing time.
echo "Step 1: Filtering reads to ROI..."
samtools view -b -L ${roi} ${ID}.bam > ${ID}_filtered.bam

# 2. Convert Filtered BAM to FASTQ and Sort by Read Name.
# This prepares reads for re-alignment by converting them to FASTQ format and ensures
# proper pairing for downstream BWA processing by sorting based on read names.
echo "Step 2: Converting BAM to FASTQ and sorting..."
# Convert BAM to FASTQ
samtools fastq ${ID}_filtered.bam -1 ${ID}.R1.pre.fastq.gz -2 ${ID}.R2.pre.fastq.gz 
zcat "${ID}.R1.pre.fastq.gz" "${ID}.R2.pre.fastq.gz" | paste - - - - | LC_ALL=C sort -k1,1 -S 3G | tr '\t' '\n' | gzip > "${ID}.mixed.fastq.gz"

# 3. Align to Combined Reference (SBDS).
# All filtered reads (from both SBDS and SBDSP1) are re-aligned to the primary SBDS reference.
# This all-inclusive alignment is critical for initial VAF estimation and haplotype inference,
# by pooling all homologous reads into one context.
echo "Step 3: Re-aligning all reads to SBDS reference using BWA-MEM..."
bwa mem -p "${RefFasta}" "${ID}.mixed.fastq.gz" > "${ID}.sam"
samtools view -bhS ${ID}.sam > ${ID}.pre.bam
samtools sort ${ID}.pre.bam > ${ID}.rg.bam
java -XX:-DoEscapeAnalysis -jar ${scriptDIR}/jar/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=SILENT I=${ID}.rg.bam O=${ID}.out.bam SO=coordinate RGID=${ID} RGLB=${ID} RGPL=Illumina,paired-end RGPU=${ID} RGSM=${ID} CREATE_INDEX=true;

# 4. Initial Variant Calling.
# Perform the first pass of variant calling on the unified BAM. The VCF generated here
# will contain variants from both SBDS and SBDSP1, with VAFs that reflect the combined pool.
# This VCF is crucial input for the PseudoRA Python script for detailed phasing.
echo "Step 4: Adding read groups and indexing combined BAM..."
java -XX:-DoEscapeAnalysis -jar ${scriptDIR}/jar/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${RefFasta} -I ${ID}.out.bam -L ${roi} -o ${ID}.hc.vcf -stand_call_conf 20;

# 5. Run PseudoRA's Reanalysis Script (bin.py).
# This is the core of PseudoRA, where the custom Python script performs
# probabilistic haplotype phasing and reassigns reads to their true genomic origin (SBDS or SBDSP1).
echo "Step 5: Running PseudoRA's core reanalysis script (bin.py)..."
python3 ${scriptDIR}/utils/bin.py ${ID}.hc.vcf ${ID}.out.bam ${ID} ${roi}
samtools index ${ID}.correct.bam

# 6. Final Variant Calling.
# Perform the final, refined variant calling on the PseudoRA-processed BAM.
# This VCF contains accurate variant calls with corrected VAFs for the SBDS gene,
# as pseudoalignment issues have been resolved.
echo "Step 6: Performing final variant calling on corrected BAM..."
java -XX:-DoEscapeAnalysis -jar ${scriptDIR}/jar/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${RefFasta} -I ${ID}.correct.bam -L ${roi} -o ${ID}.correct.hc.vcf -stand_call_conf 20;

# 7. Remove all temporary and intermediate files to keep the directory clean.
echo "Step 7: Cleaning up temporary files..."

rm ${ID}_filtered.bam
rm ${ID}.R1.pre.fastq.gz ${ID}.R2.pre.fastq.gz ${ID}.mixed.fastq.gz
rm ${ID}.pre.bam ${ID}.sam ${ID}.rg.bam
rm ${ID}.hc.vcf* ${ID}.out.ba*

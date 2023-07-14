#!/bin/bash

script=$(realpath "${BASH_SOURCE:-$0}")
scriptDIR=$(dirname $script)
RefFasta=$scriptDIR/reference/chr7_b37_SBDS.fa
roi=$scriptDIR/reference/SBDS.bed

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

#sample filter
samtools view -b -L ${roi} ${ID}.bam > ${ID}_filtered.bam

#make to fastq
samtools fastq ${ID}_filtered.bam -1 ${ID}.R1.pre.fastq.gz -2 ${ID}.R2.pre.fastq.gz
#samtools fastq ${ID}_filtered.bam -1 ${ID}.R1.pre.fastq.gz -2 ${ID}.R2.pre.fastq.gz
zcat ${ID}.R1.pre.fastq.gz | paste - - - - | LC_ALL=C sort -k1,1 -S 3G | tr '\t' '\n' | gzip > ${ID}.R1.fastq.gz
zcat ${ID}.R2.pre.fastq.gz | paste - - - - | LC_ALL=C sort -k1,1 -S 3G | tr '\t' '\n' | gzip > ${ID}.R2.fastq.gz

#align to combined
bwa mem ${RefFasta} ${ID}.R1.fastq.gz ${ID}.R1.fastq.gz > ${ID}.sam
samtools view -bhS ${ID}.sam > ${ID}.pre.bam
samtools sort ${ID}.pre.bam > ${ID}.rg.bam
java -XX:-DoEscapeAnalysis -jar ${scriptDIR}/jar/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=SILENT I=${ID}.rg.bam O=${ID}.out.bam SO=coordinate RGID=${ID} RGLB=${ID} RGPL=Illumina,paired-end RGPU=${ID} RGSM=${ID} CREATE_INDEX=true;

#VCF
java -XX:-DoEscapeAnalysis -jar ${scriptDIR}/jar/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${RefFasta} -I ${ID}.out.bam -L ${roi} -o ${ID}.hc.vcf -stand_call_conf 20;

#reanalyze
python3 ${scriptDIR}/utils/bin.py ${ID}.hc.vcf ${ID}.out.bam ${ID} ${roi}
samtools index ${ID}.correct.bam
java -XX:-DoEscapeAnalysis -jar ${scriptDIR}/jar/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${RefFasta} -I ${ID}.correct.bam -L ${roi} -o ${ID}.correct.hc.vcf -stand_call_conf 20;

#remove files
rm ${ID}_filtered.bam
rm ${ID}.R1.pre.fastq.gz ${ID}.R2.pre.fastq.gz
rm ${ID}.pre.bam ${ID}.sam ${ID}.rg.bam
rm ${ID}.R1.fastq.gz ${ID}.R2.fastq.gz ${ID}.hc.vcf* ${ID}.out.ba*

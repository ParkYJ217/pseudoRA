# PseudoRA : Pseudogene Re-Aligner for short-read next-generation sequencing data

## 1. What is PseudoRA?
**PseudoRA** aims to realign mismapped reads from the general short-read NGS BAM files due to highly-homologous sequences between gene-pseudogene pairs. Briefly, the software works in following steps: (1) Take all the reads from both functional and pseudogene, (2) Phase each read independently, (3) Rank each read by haplotype, (4) Remake BAM and VCF. The software outputs (1) correctly aligned BAM file with only reads belonging to the functional gene and (2) vcf file made by HaplotypeCaller(gatk). These files can be merged into your existing BAM and VCF for further analysis. Currently, the program's default is set for SBDS and SBDSP1 region, user can change these settings.

## 2. Installation
Just download the latest release from the GitHub repository and uncompress the tarball in a suitable directory. The tarball includes the PseudoRA script as well as the third-party software redistributed with PseudoRA (see section 5). The INSTALL files contain detailed installation instructions, including all the external libraries required to make PseudoRA run in Ubuntu.

The test_install.pl script can be run in order to check whether the required dependencies are available in your environment.

    </path/to/PseudoRA/>utils/test_install.pl

## 3. Testing PseudoRA (DEMO)

    bash </path/to/PseudoRA/>utils/pseudoRA.sh -t

## 4. Running scripts

The command for running PseudoRA has the following syntax:

    bash </path/to/PseudoRA/>utils/pseudoRA.sh -i <input BAM>

>Arguments Mandatory parameters
>
>-i <string>: Destination of original BAM file (REQUIRED)
>-r <string>: Reference FASTA (OPTIONAL) [Default:reference/chr7_b37_SBDS.fa]
>-b <string>: Region of interest BED (OPTIONAL)  [Default:reference/SBDS.bed]

## 5. Customization for other genes

    python3 </path/to/PseudoRA/>utils/customization.py -r <reference FASTA> -b <region-of-interest BED> -o <output FASTA>
    bwa index <output FASTA>
    samtools faidx <output FASTA>
    java -jar </path/to/PseudoRA/>jar/picard.jar CreateSequenceDictionary R=<output FASTA> O=<output DICT>

## 6. License and third-party software

PseudoRA is distributed under a GPL-3 license. Additionally, SqueezeMeta redistributes the following third-party software:
[bwa](https://github.com/lh3/bwa)
[gatk](https://gatk.broadinstitute.org/hc/en-us)
[picard](https://broadinstitute.github.io/picard/)
[samtools](http://www.htslib.org/)

## 7. Reference

Acknowledgements
Author of pipeline: Yu Jin Park (parkyj217@gmail.com)
Principal Investigators: Saeam Shin and Seung-Tae Lee
Institution: Yonsei University, College of Medicine, Department of Laboratory Medicine

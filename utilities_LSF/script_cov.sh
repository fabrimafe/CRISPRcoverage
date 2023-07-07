#!/bin/bash

#script to compute genome-wide coverage using GATK
#you can run this by:
#./script_cov.sh myinput.bam output.tab chrX
#ARGUMENTS:
#myinput.bam: a bam file
#output.tab: the output file
#chrX: the chromosome for which to compute coverage
module load jdk/8.111
module load GATK/3.7
module load picard-tools/1.119
module load samtools/1.9
module load bwa
inputbam=$1
outputfile=$2
chromosome=$3
#samtools index $inputbam
java -jar /apps/RH7U2/general/GATK/3.7/GenomeAnalysisTK.jar -T DepthOfCoverage -R /home/labs/alevy/Collaboration/Tomato_WGS/Tomato_plasmid_gRNA_Genome_Or_Aviva/reference_genome.fa -I ${inputbam} -o ${outputfile} -L ${chromosome}

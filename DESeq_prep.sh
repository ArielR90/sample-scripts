#!/bin/bash

#Shell commands to get and process sequencing files for DESeq on ACISS, these are just the 
#files that we were using in comparing wt and myd88. Once you're logged on and have an interactive 
# node running("qsub -I" on the head node), you can make a directory for the project and then 
# run these commands, which should give you a .tab file at the end to send to DESeq. I've used 
# some shortcuts, which should all work, but if not, let me know, and I'll fix it!

wget http://gc3fstorage.uoregon.edu/HHHK2BGXX/131/vibrio-mono-24hr-1_S6_R1_001.fastq.gz
wget http://gc3fstorage.uoregon.edu/HHHK2BGXX/131/vibrio-mono-24hr-2_S7_R1_001.fastq.gz
wget http://gc3fstorage.uoregon.edu/HHHK2BGXX/131/vibrio-mono-24hr-3_S9_R1_001.fastq.gz
wget http://gc3fstorage.uoregon.edu/HHHK2BGXX/131/vibrio-myd88-24hr-1_S12_R1_001.fastq.gz
wget http://gc3fstorage.uoregon.edu/HHHK2BGXX/131/vibrio-myd88-24hr-2_S13_R1_001.fastq.gz
wget http://gc3fstorage.uoregon.edu/HHHK2BGXX/131/vibrio-myd88-24hr-3_S14_R1_001.fastq.gz
wget http://gc3fstorage.uoregon.edu/HHHK2BGXX/131/vibrio-myd88-72hr-1_S15_R1_001.fastq.gz
wget http://gc3fstorage.uoregon.edu/HHHK2BGXX/131/vibrio-myd88-72hr-2_S16_R1_001.fastq.gz
wget http://gc3fstorage.uoregon.edu/HHHK2BGXX/131/vibrio-myd88-72hr-3_S17_R1_001.fastq.gz
wget http://meta.uoregon.edu/genomes/ZWU0020/ZWU0020.scaffolds.fasta
wget http://meta.uoregon.edu/genomes/ZWU0020/ZWU0020.gff

#these load the necessary programs 
module load bowtie
module load samtools
module load bedtools

# this builds a reference index for the vibrio genome for bowtie2, then aligns all the files
bowtie2-build ZWU0020.scaffolds.fasta ZWU0020
ls -1 *gz | while read file; do bowtie2 -x ./ZWU0020 -U $file -S ${file/_S*fastq.gz/.sam}; done

# this reformats the sam files to bam files, and indexes them for bedtools to use
ls -1 *sam | while read file; do samtools view -Sb $file > ${file/sam/bam}; done
ls -1 *bam | while read file; do samtools sort $file ${file%.bam}; done
ls -1 *bam | while read file; do samtools index $file; done

#this counts the reads at each gene in each bam file, giving you a counts table as output.
#the order you give it the bam files is the order it places the sampels columns, even though 
#there isn't a header

head -n 3804 ZWU0020.gff > ZWU0020.nofa.gff

bedtools multicov -bams vibrio-mono-24hr-1.bam vibrio-mono-24hr-2.bam vibrio-mono-24hr-3.bam vibrio-myd88-24hr-1.bam vibrio-myd88-24hr-2.bam vibrio-myd88-24hr-3.bam vibrio-myd88-72hr-1.bam vibrio-myd88-72hr-2.bam vibrio-myd88-72hr-3.bam -bed ZWU0020.nofa.gff > gene_counts.gff

#this reformats the file so you can give it directly to DESeq
sed 's/^.*locus_tag=//' gene_counts.gff > gene_counts.tab
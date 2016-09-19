#!/bin/bash 

module load bbmap/34.41
module load usearch
module load bwa
module load bedtools
module load samtools/1.1
module load lofreq
module load java/1.7
module load vcftools/0.1.11

src=$HOME/ddRAD/tests
ref=$HOME/ref/eukaryotic/Danio_rerio.GRCz10.dna.toplevel.fa
bed=$HOME/ddRAD/bwa_tis/newref/hc_tis.bed

# Merge reads, trim overhang
gunzip $src/*gz
ls -1 $src/*1.fq | while read file; do usearch -fastq_mergepairs $file -reverse ${file/1.fq/2.fq} -fastqout ${file/1.fq/mer.fq} -fastq_merge_maxee 1.0 -fastq_qmaxout 92 -fastq_qmax 92; done
ls -1 $src/*mer.fq | while read file; do reformat.sh in=$file qin=33 qtrim=t trimq=40 minlength=40 out=${file/mer/tr}; done

# Align files to reference
ls -1 $src/*tr.fq | while read file; do bwa mem $HOME/ref/eukaryotic/bwa/dan_rer_10 -t 12 $file > ${file/tr.fq/sam}; done

# Create bam, index
ls -1 $src/*1.sam | while read file; do java -jar /usr/local/packages/picard/1.92/MergeSamFiles.jar INPUT=$file INPUT=${file/1.s/2.s} OUTPUT=${file/1.s/mer.s}; done
ls -1 $src/*.sam  | while read file; do samtools view -Sb $file > ${file/sam/bam}; done
ls -1 $src/*.bam  | while read file; do samtools sort $file  ${file%.bam}; done
ls -1 $src/*.bam  | while read file; do samtools index $file; done

# Restrict to reads in highcov regions
ls -1 $src/*.bam  | while read file; do samtools view -bL $bed -o ${file/.bam/_sites.bam} $file; done

# Downsample to 250,000 reads
ls -1 $src/*sites.bam  | while read file; do reformat.sh in=$file out=${file/sites.bam/ds250.sam} qin=33 samplereads=250000; done

# Create Bam, sort
ls -1 $src/*0.sam | while read file; do samtools view -Sb -o ${file/sam/bam} $file; done
ls -1 $src/*0.bam | while read file; do samtools sort $file ${file%.bam}; done

# Realign, assign indel qualities, call SNPs
ls -1 $src/*0.bam | while read file; do lofreq viterbi -f $ref $file | samtools sort - ${file/ds250.bam/ra}; done
ls -1 $src/*ra.bam | while read file; do lofreq indelqual --dindel -f $ref -o ${file/ra/iq} $file ; done
ls -1 $src/*iq.bam  | while read file; do samtools index $file; done
ls -1 $src/*iq.bam | while read file; do lofreq call-parallel --pp-threads 12 -f $ref -l $bed --no-default-filter --call-indels -o ${file/iq.bam/all.vcf} $file; done
ls -1 $src/*all.vcf | while read file; do lofreq filter -i $file -v 100 -o ${file/_all/_flt}; done

# Filter SNPs by barcodes
ls -1 $src/*flt.vcf | while read file; do bgzip -c $file > $file.gz; done
ls -1 $src/*vcf.gz | while read file; do tabix -p vcf $file; done
ls -1 $src/*1_flt.vcf.gz | while read file; do lofreq vcfset -a intersect -1 $file -2 ${file/1_flt/2_flt} -o ${file/1_flt/shd}; done
ls -1 $src/*shd.vcf.gz | while read file; do tabix -p vcf $file; done
ls -1 $src/*mer_flt.vcf.gz | while read file; do lofreq vcfset -a intersect -1 $file -2 ${file/mer_flt/shd} -o ${file/mer_flt.vcf.gz/PELE.vcf}; done







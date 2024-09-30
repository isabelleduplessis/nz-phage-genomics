#!/bin/bash

# Commands to map raw reads to viral seqs to get depth and coverage info for relative abundance analysis

bwa index -p nztotalseqs nztotalseqs.fna

gunzip -c /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/NZ/rawreads/NZ${samplenumber}_*R1*.fastq.gz > R1_${samplenumber}.fastq
gunzip -c /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/NZ/rawreads/NZ${samplenumber}_*R2*.fastq.gz > R2_${samplenumber}.fastq
bwa mem -t 4 nztotalseqs R1_${samplenumber}.fastq R2_${samplenumber}.fastq > sam/sample_${samplenumber}.sam 
rm R1_${samplenumber}.fastq
rm R2_${samplenumber}.fastq

samtools view -bS sample_${samplenumber}.sam > sample_${samplenumber}.bam
samtools sort sample_${samplenumber}.bam > sample_${samplenumber}_sorted.bam
samtools depth sample_${samplenumber}_sorted.bam > sample_${samplenumber}_depth.txt
samtools coverage sample_${samplenumber}_sorted.bam > sample_${samplenumber}_coverage.txt


#Use the coverage information to get the viral-sample pairs where the virus is covered at a least 75% sequence length
#The following lines output a file that contains rpkm and meandepth information
cat file.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'

#get the total number of reads in each sample 
cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance/bam
for i in {1..22}; do totreads=$(samtools view -c NZ${i}_sorted.bam); echo "NZ${i},$totreads"; done > totalreads.csv

#Get the viral mappings that pass the 75% threshold with their meandepth and  RPKM
cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance/coverage
for ARG in {1..22}; do totreads=$(grep "NZ$ARG," totalreads.csv | cut -f2 -d ","); awk -F'\t' 'NR==1 || $6>=75' "NZ${ARG}_coverage.txt" | grep -v "^#" | sed -e "s/^/NZ$ARG\t/" | awk -v t=$totreads '{print $1"\t"$2"\t"$7"\t"$8"\t"($5*10**9)/($4*t)}'; done > votus_cov75thres.txt;
#add header
sed -i '1 i\#sample\tvotus\tcoverage\tmeandepth\trpkm' votus_cov75thres.txt 

#Step 5: get the number of reads that mapped to viral sequences, for plotting purposes.

#Get the reads that mapped to something viral
cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance/bam
for ARG in {1..22}; do samtools view -b -F 4 NZ${ARG}_sorted.bam > NZ${ARG}_mapped.bam; done
#Get number of mapped reads 
for i in {1..22}; do totreads=$(samtools view -c NZ${i}_mapped.bam); echo -e "NZ${i}\t$totreads"; done > totalreads_mapped.csv  

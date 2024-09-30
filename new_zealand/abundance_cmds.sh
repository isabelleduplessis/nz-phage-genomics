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
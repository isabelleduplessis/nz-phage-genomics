#!/bin/bash

# Commands used to increase viral seqs by:
# Mapping rawreads to viral predictions
# Taking all the reads that mapped from each sample at all and reassemble with megahit
# Running virsorter on those
# Filtering for 5000 bp
# Select non not detrmined checkv sequences
# CDhit with those and original 8000 viral seqs


for i in {1..22}; do 
	samtools sort -n -o NZ$i\sorted NZ$i\_mapped.bam
	bedtools bamtofastq -i NZ$i\sorted -fq NZ$i\_R1 -fq2 NZ$i\_R2
	done

for f in sample*; do mv "$f" "$f.fq"; done

for i in {2..22}; do sbatch --export=samplenumber=$i megahit.sbatch; done

sbatch --export=samplenumber=1 megahit.sbatch

megahit -1 sample_${samplenumber}_R1.fq -2 sample_${samplenumber}_R2.fq -o megahit_${samplenumber} -t 4

seqkit split2 -s 20000 /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/salt/PIG_GOV_db/bwaPIGEONindex/PIGEONv1.0.fa

virsorter run -w viral/sample_${samplenumber}_viral.fa -i megahit_${samplenumber}/final.contigs.fa --min-length 1500 -j 4 all

sbatch --export=samplenumber=1 virsort.sbatch

for i in {2..22}; do sbatch --export=samplenumber=$i virsort.sbatch; done

for f in *.faa; do sed -i "s/^>/>${f}_/" "$f"; done

for f in *.fa/final-viral-combined.fa; do sed -i "s/^>/>${f%.fa}/g" "${f}"; done > viralseqs.fa

# add sample number to each header

for i in {1..22}; do sed "s%^>\(.*\)%>NZ${i}/\1/%I" sample_${i}_viral.fa/final-viral-combined.fa; done > viralseqs.fa










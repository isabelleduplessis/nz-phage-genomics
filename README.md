# Introduction

These are the commands used to characterize viral sequences from marine samples near New Zealand. The goal of this project is to understand how different temperatures, salinity, nutrient levels, and water depths impact viral genomes. This is part of the initial stages of a research project. 

# Methods

## Packages Used:
need to update

## Quality Checking Viral Sequences

Started with 58643 viral sequences predicted with VirSorter: all_virsorter_contigs_Daniel.fa

Samples come from https://www.nature.com/articles/s41467-022-35204-6#Sec10

### Filter Sequence Lengths - 5KB

Filtered to sequences above length 5000 using seqkit v2.4.0 with conda:

    seqkit seq -m 5000 -g all_virsorter_contigs_Daniel.fa > nzvirsorterseqs5000.fna

This resulted in 10494 sequences.

### CDHIT clustering to get unique viral populations

Next, clustered sequences with CDHIT v4.8.1:

    cd-hit-est -i nzvirsorterseqs5000.fna -o nzvirsorterseqscdhit.fna -c 0.95 -s 0.8

This resulted in 8229 sequences.

### CheckV to keep quality viral sequences

Next, CheckV v1.0.1 was run on local computer with the following
commands:

    conda activate checkv
    export CHECKVDB=checkv-db-v1.5
    checkv end_to_end nzvirsorterseqscdhit.fna checkvresults -t 2

Get non \"Not-determined\" sequences:

    grep -v "Not-determined" checkvresults/quality_summary.tsv > quality_summary_filtered.tsv
    grep -A 1 -f <(awk '(NR>1) {print $1}' checkv/quality_summary_filtered.tsv) nzvirsorterseqscdhit.fna | sed '/^-/d' > nzseqs7890.fna


This resulted in 7890 viral sequences, which can be found at:

    nzseqs7890.fna

## Increasing Viral Sequences

### Read Mapping

First, the original viral sequences were indexed with BWA.

    bwa index -p all_virsorter_contigs_Daniel all_virsorter_contigs_Daniel.fa

Then, the raw reads were mapped to the 7890 sequences.

    bwa mem -t 4 all_virsorter_contigs_Daniel R1_${samplenumber}.fastq R2_${samplenumber}.fastq > sam/NZ${samplenumber}.sam 

Sam files were converted to sorted bam files.

    samtools view -bS NZ${samplenumber}.sam > NZ${samplenumber}.bam
    samtools sort NZ${samplenumber}.bam > NZ${samplenumber}_sorted.bam

Bam files were generated for the reads that mapped.

    for ARG in {1..22}; do samtools view -b -F 4 NZ${ARG}_sorted.bam > NZ${ARG}_mapped.bam; done

Next, those bam files were used to create paired fastq files for
reassembly.

    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/nz/increasingseqs/sam/mappedbam
    module load bedtools2/2.30.0-mlxpvg
    module load samtools/1.14-3alo66
    for i in {1..22}; do 
        samtools sort -n -o NZ$i\sorted NZ$i\_mapped.bam
        bedtools bamtofastq -i NZ$i\sorted -fq NZ$i\_R1 -fq2 NZ$i\_R2
        done

### Megahit

The fq files were reassembled with megahit.

    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/nz/increasingseqs/sam/mappedbam/fq
    module load anaconda3
    conda activate megahit
    megahit -1 NZ${samplenumber}_R1.fq -2 NZ${samplenumber}_R2.fq -o megahit_${samplenumber} -t 4

### VirSorter

Run virsorter on each of the 22 outputs.

    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/nz/increasingseqs/sam/mappedbam/fq
    module load anaconda3
    conda activate vs2
    virsorter run -w viral/NZ${samplenumber}_viral.fa -i megahit_${samplenumber}/final.contigs.fa --min-length 1500 -j 4 all

Add sample number to headers and concatenate viral sequences:

    for i in {1..22}; do sed "s%^>\(.*\)%>NZ${i}/\1/%I" sample_${i}_viral.fa/final-viral-combined.fa; done > viralseqs.fa

This resulted in 112825 sequences.

### 5KB

    seqkit seq -m 5000 -g viralseqs.fa > viralseqs5000.fa

This resulted in 21587 sequences.

### CDHIT

    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/nz/increasingseqs/sam/mappedbam/fq/viral
    module load cdhit/4.8.1-n4smhr
    cd-hit-est -i viralseqs5000.fa -o increasedviralseqs.fa -c 0.95 -s 0.8

This resulted in 14092 sequences

### CheckV

    conda activate checkv
    export CHECKVDB=/Users/isabelleduplessis/phage/NZ/checkv-db-v1.5
    checkv end_to_end increasedviralseqs.fa checknewseqs -t 2

    awk '(NR>1) {print $1}' checknewseqs/quality_summary_filtered.tsv > increasedseqids.txt
    seqkit grep -f increasedseqids.txt increasedviralseqs.fa > nzincreasedseqs.fna
    cat nzincreasedseqs.fna nzseqs7890.fna > nzcombinedseqs.fna

This resulted in 13586 sequences. After combining the files and before
clustering, there were 21476 sequences.

### CDHIT

    cd-hit-est -i nzcombinedseqs.fna -o nztotalseqs.fna -c 0.95 -s 0.8

This resulted in 18490 sequences.

Get sequence lengths:

    awk '/^>/ { print (NR==1 ? "" : RS) $0; next } \
    { printf "%s", $0 } END { printf RS }' nztotalseqs.fna | \
    awk '/^>/{sub(/^>/,"");val=$0;next} \ 
    {print val,length($0);val=""} \
    END{if(val!=""){print val}}' | \
    awk '{print $2}'

## Abundance and Diversity

### Read Mapping

The nztotalseqs.fna file was indexed using BWA:

    module load bwa/0.7.17-r23kej
    bwa index -p nztotalseqs nztotalseqs.fna

Then, the raw reads were mapped to this index:

    bwa mem -t 4 nztotalseqs R1_${samplenumber}.fastq R2_${samplenumber}.fastq \
    > sam/sample_${samplenumber}.sam 

### Coverage

    samtools view -bS sample_${samplenumber}.sam > sample_${samplenumber}.bam
    samtools sort sample_${samplenumber}.bam > sample_${samplenumber}_sorted.bam
    samtools depth sample_${samplenumber}_sorted.bam > sample_${samplenumber}_depth.txt
    samtools coverage sample_${samplenumber}_sorted.bam > sample_${samplenumber}_coverage.txt

After running coverage.sbatch:

    for i in {1..22}; do totreads=$(samtools view -c sample_${i}_sorted.bam); \
    echo "sample_${i},$totreads"; done > totalreads.csv

    for ARG in {1..22}; do totreads=$(grep "sample_$ARG," totalreads.csv | \
    cut -f2 -d ","); awk -F'\t' 'NR==1 || $6>=75' "sample_${ARG}_coverage.txt" | grep -v "^#" | \
    sed -e "s/^/sample_$ARG\t/" | \
    awk -v t=$totreads '{print $1"\t"$2"\t"$7"\t"$8"\t"($5*10**9)/($4*t)}'; \
    done > votus_cov75thres.txt;

    sed -i '1 i\#sample\tvotus\tcoverage\tmeandepth\trpkm' votus_cov75thres.txt

    for ARG in {1..22}; do samtools view -b -F 4 sample_${ARG}_sorted.bam \
    > sample_${ARG}_mapped.bam; done

    for i in {1..22}; do totreads=$(samtools view -c sample_${i}_mapped.bam); \
    echo -e "sample_${i}\t$totreads"; done > totalreads_mapped.csv 

# Visualizing Results
add information for running R plot


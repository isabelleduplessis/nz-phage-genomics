#!/bin/bash

# start with 58643 sequences in all_virsorter_contigs_Daniel.fa

# get sequences above length 5000 (using seqkit v2.4.0 with conda)
seqkit seq -m 5000 -g /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/NZ/all_virsorter_contigs_Daniel.fa > nzvirsorterseqs5000.fna
# 10494 sequences

# cdhit v4.8.1
cd-hit-est -i nzvirsorterseqs5000.fna -o nzvirsorterseqscdhit.fna -c 0.95 -s 0.8
# 8229 sequences

# Run checkv v1.0.1 on local computer
scp iplessis3@login-phoenix-slurm.pace.gatech.edu:/storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/nz/nzvirsorterseqscdhit.fna /Users/isabelleduplessis/phage/NZ
cd /Users/isabelleduplessis/phage/NZ
conda activate checkv
export CHECKVDB=/Users/isabelleduplessis/phage/NZ/checkv-db-v1.5
checkv end_to_end nzvirsorterseqscdhit.fna checkvresults -t 2

# get non "Not-determined" sequences. These are final viral sequences to be used for analysis
grep -v "Not-determined" checkvresults/quality_summary.tsv > quality_summary_filtered.tsv
grep -A 1 -f <(awk '(NR>1) {print $1}' checkv/quality_summary_filtered.tsv) nzvirsorterseqscdhit.fna | sed '/^-/d' > nzfinalseqs.fna
# 7890 final sequences

# get sequence lengths
cat nzfinalseqs.fna | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | awk '{print $2}' > nzseqlens.txt

# final sequences and sequence lengths can be found at /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/nz

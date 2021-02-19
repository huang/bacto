# bacto
A pipeline for clinical bacteria NGS data

#1. Run FastQC to allow manual inspection of the quality of sequences
mkdir fastqc_out
fastqc -t 4 raw_data/* -o fastqc_out/

#2. Rename the files
cd raw_data
for file in *.fastq.gz; do mv $file $(echo $file | cut -d'_' -f1 | cut -d'-' -f1-2)_$(echo $file | cut -d'_' -f4).fastq.gz; done
cd ..

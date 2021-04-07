# bacto
A pipeline for clinical bacteria NGS data

## Run FastQC to allow manual inspection of the quality of sequences
```sh
mkdir fastqc_out
fastqc -t 4 raw_data/* -o fastqc_out/
```

## Rename the files
```sh
cd raw_data
for file in *.fastq.gz; do mv $file $(echo $file | cut -d'_' -f1 | cut -d'-' -f1-2)_$(echo $file | cut -d'_' -f4).fastq.gz; done
cd ..
```

## Run snakemake
```sh
git clone https://github.com/huang/bacto
conda activate bengal3_ac3
snakemake --printshellcmds
#Note that we need delete/add files in fastq if we want to change the results in roary, variants, fasttree, raxml-ng, and gubbins.
```

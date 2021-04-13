# bacto
A pipeline for clinical bacteria NGS data

## 1, run FastQC to allow manual inspection of the quality of sequences
```sh
mkdir fastqc_out
fastqc -t 4 raw_data/* -o fastqc_out/
```

## 2, rename the files
```sh
cd raw_data
for file in *.fastq.gz; do mv $file $(echo $file | cut -d'_' -f1 | cut -d'-' -f1-2)_$(echo $file | cut -d'_' -f4).fastq.gz; done
cd ..
```

## 3, run snakemake
```sh
git clone https://github.com/huang/bacto
conda activate bengal3_ac3
snakemake --printshellcmds
#Note that we need delete/add files in fastq if we want to change the results in roary, variants, fasttree, raxml-ng, and gubbins.
#https://comparative-genomics.readthedocs.io/en/latest/day2_afternoon.html
```

## 4, Identify species
- StrainSeeker: fast identification of bacterial strains from unassembled sequencing reads using user-provided guide trees
- Reads2Type: a web application for rapid microbial taxonomy identification
- https://github.com/tseemann/sixess (sixess Hentschke_22927_R1.fastq.gz Hentschke_22927_R2.fastq.gz)
- https://pubmlst.org/bigsdb?db=pubmlst_rmlst_seqdef_kiosk

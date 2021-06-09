# bacto
A pipeline for clinical bacteria NGS data

```sh
#…or create a new repository on the command line
#echo "# damian_extended" >> README.md
#git init
#git add README.md
#git commit -m "first commit"
#git remote add origin https://github.com/huang/damian_extended.git
#git push -u origin master
#.......or push an existing repository from the command line
#git remote add origin https://github.com/huang/damian_extended.git or git remote set-url origin https://github.com/huang/damian_extended.git
#git push -u origin master

#git add --all
#git commit -am "<commit message>"
#git push

#https://www.heise.de/tipps-tricks/SSH-Key-fuer-GitHub-Repositories-einrichten-4627459.html
```

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
- ReferenceSeeker: rapid determination of appropriate reference genomes    https://github.com/oschwengers/referenceseeker
```sh
Installation:

$ conda install -c bioconda referenceseeker
$ wget https://zenodo.org/record/4415843/files/bacteria-refseq.tar.gz
$ tar -xzf bacteria-refseq.tar.gz
$ rm bacteria-refseq.tar.gz

Simple:

$ # referenceseeker <REFERENCE_SEEKER_DB> <GENOME>
$ referenceseeker bacteria-refseq/ genome.fasta

Expert: verbose output and increased output of candidate reference genomes using a defined number of threads:

$ # referenceseeker --crg 500 --verbose --threads 8 <REFERENCE_SEEKER_DB> <GENOME>
$ referenceseeker --crg 500 --verbose --threads 8 bacteria-refseq/ genome.fasta

(bengal3_ac3) jhuang@hamburg:~/Tools$ referenceseeker ~/Tools/referenceseeker/test/db ~/Tools/referenceseeker/test/data/Salmonella_enterica_CFSAN000189.fasta
ERROR: failed to execute nucmer!
exit=1
cmd=['nucmer', '--threads=1', '/home/jhuang/Tools/referenceseeker/test/db/GCF_000439415.1.fna', '/tmp/tmp0ljafd1l/dna-fragments.fasta']
```
- Reads2Type: a web application for rapid microbial taxonomy identification
- https://github.com/tseemann/sixess (sixess Hentschke_22927_R1.fastq.gz Hentschke_22927_R2.fastq.gz)
- https://pubmlst.org/bigsdb?db=pubmlst_rmlst_seqdef_kiosk
```sh
spades.py --meta -t 16 -1 trimmed/Hentschke_22927_trimmed_P_1.fastq -2 trimmed/Hentschke_22927_trimmed_P_2.fastq -o spades

bowtie -p 14 -S -f -v 3 -m 1 --best --strata extracted_asso_genes Pantoea_alhagi_LTYR-11Z.fasta > kmer_on_assoc_genes_3maxMismatches.sam
bowtie -p 14 -S -f -v 3 -m 1 --best --strata extracted_asso_genes 96_training_kmers.fasta > 575genes_96kmers_3maxMismatches.sam
bowtie -p 14 -S -f -v 3 -m 1 --best --strata extracted_asso_genes_BH 96_training_kmers.fasta > 205genesBH_96kmers_3maxMismatches.sam

bowtie2-build Pantoea_alhagi_LTYR-11Z.fasta Pantoea_alhagi_LTYR-11Z_index
bowtie2-build Mixta_theicola_SRCM103227.fasta Mixta_theicola_SRCM103227_index
bowtie2-build both.fasta both_index
bowtie2 -x Pantoea_alhagi_LTYR-11Z_index -1 trimmed/Hentschke_22927_trimmed_P_1.fastq -2 trimmed/Hentschke_22927_trimmed_P_2.fastq --threads 15 --very-sensitive --al-conc-gz Pantoea_alhagi_mapped.fastq.gz --un-conc-gz Pantoea_alhagi_unmapped.fastq.gz > Pantoea_alhagi_LTYR-11Z.sam
bowtie2 -x Mixta_theicola_SRCM103227_index -1 trimmed/Hentschke_22927_trimmed_P_1.fastq -2 trimmed/Hentschke_22927_trimmed_P_2.fastq --threads 15 --very-sensitive --al-conc-gz Mixta_theicola_mapped.fastq.gz --un-conc-gz Mixta_theicola_unmapped.fastq.gz > Mixta_theicola_SRCM103227.sam
bowtie2 -x both_index -1 trimmed/Hentschke_22927_trimmed_P_1.fastq -2 trimmed/Hentschke_22927_trimmed_P_2.fastq --threads 15 --very-sensitive --al-conc-gz both_mapped.fastq.gz --un-conc-gz both_unmapped.fastq.gz > both.sam
samtools view -bS -h Pantoea_alhagi_LTYR-11Z.sam | samtools sort - Pantoea_alhagi_LTYR-11Z_sorted
samtools view -bS -h Mixta_theicola_SRCM103227.sam | samtools sort - Mixta_theicola_SRCM103227_sorted
samtools view -bS -h both.sam | samtools sort - both_sorted

makeblastdb -in Pantoea_alhagi_LTYR-11Z.fasta -dbtype 'nucl' -out Pantoea_alhagi_LTYR-11Z.fa.db 
blastn -db Pantoea_alhagi_LTYR-11Z.fa.db -query contigs_Pantoea_alhagi.fasta -out newcontigs_Pantoea_alhagi.blastn -evalue 1e-10  -num_threads 15 -outfmt 6 
=(SUM(E1:E1054)/SUM(D1:D1054))=330471920,87/3721906=88,791044392=88.8%
makeblastdb -in Mixta_theicola_SRCM103227.fasta -dbtype 'nucl' -out Mixta_theicola_SRCM103227.fa.db 
blastn -db Mixta_theicola_SRCM103227.fa.db -query contigs_Mixta_theicola.fasta -out newcontigs_Mixta_theicola.blastn -evalue 1e-10  -num_threads 15 -outfmt 6 
=84.3%
```

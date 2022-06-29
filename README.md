# bacto
A pipeline for clinical bacteria NGS data

WARNING: used the absolute path '-hmm /media/jhuang/Titisee/GAMOLA2/TIGRfam_db/TIGRFAMs_15.0_HMM.LIB' in Snakefile

```sh
#https://www.cyberciti.biz/faq/installing-latest-stable-mainline-linux-kernel-on-ubuntu-with-apt-get/
#https://askubuntu.com/questions/89710/how-do-i-free-up-more-space-in-boot
#dpkg -l 'linux-*' | sed '/^ii/!d;/'"$(uname -r | sed "s/\(.*\)-\([^0-9]\+\)/\1/")"'/d;s/^[^ ]* [^ ]* \([^ ]*\).*/\1/;/[0-9]/!d' | xargs sudo apt-get -y purge
genbankdownload.py -t fasta CP000851
(optional) download genomes https://github.com/kblin/ncbi-genome-download
ncbi-genome-download --species-taxids 562 --formats fasta --refseq-categories reference --assembly-levels complete --parallel 4  bacteria
ncbi-genome-download --species-taxids 562 --formats fasta  --assembly-levels complete --parallel 10 bacteria
ncbi-genome-download -T 1423 -F fasta -l complete --parallel 10 bacteria

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
mv bacto/* ./
rm -rf bacto
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

## 5, variant calling
```sh
[annotation of contigs]
python ~/Tools/VAPiD/vapid3.py --db ~/REFs/all_virus/all_virus.fasta contigs.fa ~/REFs/template_Holger.sbt
#EF9117
https://www.nature.com/articles/srep34963
#GFinisher: a new strategy to refine and finish bacterial genome assemblies
java -Xms2G -Xmx4G -jar GenomeFinisher_1.4/GenomeFinisher.jar  -i shovill/noAB_wildtype/contigs.fa  -ref CP040849.fasta  -o outputDirectory
http://lu168.cs.nthu.edu.tw/CAR/index.php
~/Tools/CONTIGuator_v2.7/CONTIGuator.py -r CP040849.fasta -a /home/jhuang/Tools/act.jar -c shovill/noAB_wildtype/contigs.fa 
~/Tools/CONTIGuator_v2.7/CONTIGuator.py -r CP059793.fasta -a /home/jhuang/Tools/act.jar -c shovill/noAB_wildtype/contigs.fa 
E50862
(referenceseeker) ~/Tools/referenceseeker/bin/referenceseeker -v ~/REFs/bacteria-refseq/ shovill/E50862/contigs.fa
#ID     Mash Distance   ANI     Con. DNA        Taxonomy ID     Assembly Status Organism
GCF_001281065.1 0.00712 98.70   91.86   1747    complete        Cutibacterium acnes KCOM 1861 (= ChDC B594)
GCF_000231215.1 0.00710 98.51   91.85   1091045 complete        Cutibacterium acnes subsp. defendens ATCC 11828
GCF_000240055.1 0.01125 98.42   91.54   1114969 complete        Cutibacterium acnes TypeIA2 P.acn31
GCF_006739385.1 0.01108 98.42   91.48   1734925 complete        Cutibacterium acnes subsp. acnes NBRC 107605
GCF_000008345.1 0.01085 98.43   91.41   267747  complete        Cutibacterium acnes KPA171202
~/Tools/CONTIGuator_v2.7/CONTIGuator.py -r ATCC11828.fasta -a /home/jhuang/Tools/act.jar -c shovill/E50862/scaffolds.fasta
merge_seq.py Excluded.fsa > ../Excluded_sequence.fasta

[genbank copying]
#rm -rf ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/data/CP040849.1
mkdir ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/data/noAB_wildtype
cp PROKKA_01242022/PROKKA_01242022.gbk ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/data/noAB_wildtype/genes.gbk
vim ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/snpEff.config
/home/jhuang/anaconda3/envs/spandx/bin/snpEff build -genbank noAB_wildtype      -d

[spandx calling]
conda activate spandx
#NOTE that Realignment disappeared in gatk4.
#-In GATK4 the indel realignment step will no longer be part of the pipeline. With my great surprise the tool has disappeared. My question, the algorithm is no more there because it doesn't improve the alignment? Do I need to search for a replace or should I just remove it from my pipeline since it is not useful? Any suggestions?
#-Basically, the realignment step has been integrated into the GATK variant callers and they can output their own realigned BAMs. However, regardless of how you perform realignment, it is still preferred.

#NOTE that it requires two samples for merging! --> VERY IMPORTANT!!!!
nextflow run spandx/main.nf --fastq "raw_data/*_R{1,2}_001.fastq.gz" --ref NC_045512.fasta --annotation --database NC_045512 -resume
```

## 6,
```sh
https://cge.cbs.dtu.dk/services/
https://cge.cbs.dtu.dk/services/MLST/
https://pubmlst.org/bigsdb?db=pubmlst_pacnes_seqdef&page=batchProfiles&scheme_id=3
```

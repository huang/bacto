# bacto
A pipeline for clinical bacteria NGS data

WARNING: used the absolute path '-hmm /media/jhuang/Titisee/GAMOLA2/TIGRfam_db/TIGRFAMs_15.0_HMM.LIB' in Snakefile

```sh
#bioconda tutorial; TODO: prepare damian package on bioconda
https://bioconda.github.io/tutorials/gcb2020.html
#extract features from genbank-file
https://services.healthtech.dtu.dk/service.php?FeatureExtract-1.2
#extract gff3 from Genbank-file
bp_genbank2gff3 roary_.fa.gb
#https://stackoverflow.com/questions/63990487/package-conflict-with-anaconda
#IMPORTANT: By telling conda that it can visit multiple channels (conda-forge, bioconda and defaults), it has some extra options to resolve the dependency conflicts.

cd raw_data
for file in *.fastq.gz; do mv $file $(echo $file | cut -d'_' -f1)_$(echo $file | cut -d'_' -f4).fastq.gz; done
cd ..
git clone https://github.com/huang/bacto
mv bacto/* ./
rm -rf bacto
rm Comparative_Genomics_Tutorial_v2.pdf dN_dS.pdf call_shell_from_Ruby.png install_linux-generic_LOG bacterial_genomics_desc.txt

#-- we should create two environments: prokka and bacto2 --> SUCCESSFUL
# ---------------------------------------------------------------
# ---------------------- install env prokka ---------------------
conda create -n prokka -c conda-forge -c bioconda -c defaults python=3.7 prokka
(prokka) conda install -c conda-forge -c bioconda -c defaults shovill mlst trimmomatic
#(prokka) conda install -c conda-forge -c bioconda -c defaults snakemake

## ---- install env snippy on conda or host-env (NOT NEED ANYMORE, Since bacto2 get repaired by 'export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0/') ----
#-- install snippy using conda failed --
#conda activate base
#conda update --all --yes
#conda create -n snippy
#conda activate snippy
#(snippy) conda install -c conda-forge -c bioconda -c defaults snippy=4.6 roary
#(snippy) conda list (To check version is above 4)
#(snippy) snippy --check

##-- install brew and snippy to ubuntu host-env --
##https://www.how2shout.com/linux/how-to-install-brew-ubuntu-20-04-lts-linux/#google_vignette
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
##add 'eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"' to .bashrc
#source ~/.bashrc
#brew install brewsci/bio/snippy

# ---------------------------------------------------------------
# ---------------------- install env bacto2 ---------------------
conda create -n bacto2 -c conda-forge -c bioconda -c defaults python=3.7
conda activate bacto2
(bacto2) conda install -c conda-forge -c bioconda -c defaults trimmomatic fastqc shovill 
(bacto2) conda install -c conda-forge -c bioconda -c defaults mlst roary snippy fasttree 
#(bacto2) cpanm Bio::Roary  #DEBUG NOT WORKING
(bacto2) export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0/  #DEBUG SUCCESSFUL, repaired also snippy and snippy-core
(bacto2) conda install -c conda-forge -c bioconda -c defaults raxml-ng gubbins snp-sites
#(bacto2) conda install -c conda-forge -c bioconda -c defaults perl-bioperl
#(bacto2) cpanm Bio::SeqIO
#(bacto2) cpanm Bio::Perl
#(bacto2) conda install -c conda-forge -c bioconda -c defaults biopython
#(bacto2) conda install -c conda-forge -c bioconda -c defaults snakemak

# ------------------------------------------------------------------------------------------
# ------------ using prokka for first three steps and bacto2 for last 5 steps --------------
(prokka) for the first 3 steps (f,f,t,t,t,f,f,f,f,f in bacto-0.1.json)
(prokka) run prokka manually if needed: 
for sample in 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87; do
    #under the environment 'bincs'
    #bakta --db ~/REFs/bakta/db_1_1 --output bakta/${sample} --genus Pseudomonas --species aeruginosa --prefix ${sample} --locus-tag ${sample} shovill/${sample}/contigs.fa
    prokka --force --outdir prokka/${sample} --cpus 2 --usegenus --genus Staphylococcus --kingdom Bacteria --species epidermidis --addgenes --addmrna --prefix ${sample} --locustag ${sample} shovill/${sample}/contigs.fa #-hmm /media/jhuang/Titisee/GAMOLA2/TIGRfam_db/TIGRFAMs_15.0_HMM.LIB
done

#1. if the reference is given in bacto-0.1.json, automatically run all the remaining steps.
(bacto2 on titisee) for the remaining steps (f,f,f,f,f,t,f,t,t,t in bacto-0.1.json)

#2. if the reference is not given, manually run until raxml, then call plotTreeHeatmap and scoary
(bacto2 on titisee) for the step roary (f,f,f,f,f,t,f,f,f,f in bacto-0.1.json) generating roary, gathering samples only from raw_data/*.fastq.gz, be careful need to delete all roary* directories before, otherwise it will occur errors during running!
(bacto2 on titisee) export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0/
(bacto2 on titisee) roary -p 15 -f ./roary -i 95 -cd 99 -s -e -n -v  prokka/E47727/E47727.gff prokka/BK16922/BK16922.gff ...
./roary$ samtools faidx core_gene_alignment.aln
(bacto2 on titisee) ./roary$ fasttree -gtr -nt core_gene_alignment.aln > core_gene_alignment.tree
(bacto2 on titisee) ./roary$ snp-sites core_gene_alignment.aln > core_gene_alignment_.aln     #372992 --> 16848
(bacto2 on titisee) ./roary$ raxml-ng --all --model GTR+G+ASC_LEWIS --prefix core_gene_raxml --threads 6 --msa core_gene_alignment_.aln --bs-trees 1000
---- plotTreeHeatmap (optional) ----
#using tree generated in last step or generated by spandx 
cp Outputs/Phylogeny_and_annotation/ML_phylogeny.tre plotTreeHeatmap/ML_phylogeny_tree.txt

---- scoary (optional) ----
#cp ../plotTreeHeatmap/typing_188.csv ./
#,If_inf
Command: /home/jhuang/anaconda3/envs/bengal3_ac3/bin/scoary -g ../roary/gene_presence_absence.csv -t f1_13 -r restrict_to_186.csv
(bengal3_ac3) scoary -g ../roary/gene_presence_absence.csv -t f1_13 -r restrict_to_186.csv

#3. seek the next closely reference (optional)
#download referenceseeker reference databases under https://zenodo.org/record/3725706#.YuPMtGFBzCI
(referenceseeker) ~/Tools/referenceseeker/bin/referenceseeker -v ~/REFs/bacteria_refseeker_db/ shovill/E50862/contigs.fa
#refseeker_results  #ID     Mash Distance   ANI     Con. DNA        Taxonomy ID     Assembly Status Organism
#refseeker_results  GCF_001281065.1 0.00712 98.70   91.86   1747    complete        Cutibacterium acnes KCOM 1861 (= ChDC B594)
#refseeker_results  GCF_000231215.1 0.00710 98.51   91.85   1091045 complete        Cutibacterium acnes subsp. defendens ATCC 11828
#refseeker_results  GCF_000240055.1 0.01125 98.42   91.54   1114969 complete        Cutibacterium acnes TypeIA2 P.acn31
#refseeker_results  GCF_015697645.1	0.00138	99.92	94.97	287	complete	Pseudomonas aeruginosa B-I-1
#refseeker_results  GCF_003812885.1	0.00137	99.86	94.80	287	complete	Pseudomonas aeruginosa FDAARGOS_571


#4. varicant calling using spandx
#-4.1. generate PseudoContig_wildtype.fasta --
#./shovill/noAB_wildtype/contigs.fasta
~/Tools/CONTIGuator_v2.7/CONTIGuator.py -r CP040849.fasta -a /home/jhuang/Tools/act.jar -c shovill/noAB_wildtype/contigs.fa 
~/Tools/CONTIGuator_v2.7/CONTIGuator.py -r CP059793.fasta -a /home/jhuang/Tools/act.jar -c shovill/noAB_wildtype/contigs.fa 
#or using 'java -jar r2cat.jar' or CAR: http://lu168.cs.nthu.edu.tw/CAR/index.php 
merge_seq.py Excluded.fsa > ../Excluded_Sequence.fasta
#add Excluded_sequence.fasta to PseudoContig.fsa or only use PseudoContig.fsa
cat PseudoContig.fsa Excluded_Sequence.fasta > PseudoContig_Wildtype_.fasta
seqkit seq -w 80 PseudoContig_Wildtype_.fasta > PseudoContig_Wildtype.fasta

#prokka PseudoContig_Wildtype.fasta #result in PROKKA_01242022
for sample in PseudoContig_Wildtype; do
    #under the environment 'bincs'
    bakta --db ~/REFs/bakta/db_1_1 --output ${sample} --genus Pseudomonas --species aeruginosa --prefix ${sample} --locus-tag ${sample} ${sample}.fasta
done
#../../vrap/external_tools/blast/makeblastdb -in nucleotide.fa -dbtype nucl -parse_seqids -out nucleotide.fasta
python ~/Tools/VAPiD/vapid3.py --db ~/REFs/all_virus/all_virus.fasta contigs.fa ~/REFs/template_Fischer.sbt
python ~/Tools/VAPiD/vapid3.py --db ~/REFs/all_virus/all_virus.fasta contigs.fa ~/REFs/template_Holger.sbt

#-4.2-1. generate genebank in snpEff under the example noAB_wildtype --
mkdir ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/data/noAB_wildtype
cp PROKKA_01242022/PROKKA_01242022.gbk ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/data/noAB_wildtype/genes.gbk
vim ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/snpEff.config
/home/jhuang/anaconda3/envs/spandx/bin/snpEff build -genbank noAB_wildtype      -d

#-4.2-2. generate genbank in snpEff under the example Holger_Kpneumoniae_SNP --
cat reordered.fasta unmatched.fasta > PseudoContig_wildtype.fasta
echo ">wildtype_150" > wildtype_150.fasta;
~/Scripts/merge_seq.py PseudoContig_wildtype.fasta >> wildtype_150.fasta;
prokka wildtype_150.fasta
conda activate spandx
mkdir ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/data/wildtype_150
cp PROKKA_09062022/PROKKA_09062022.gbk ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/data/wildtype_150/genes.gbk
vim ~/anaconda3/envs/spandx/share/snpeff-4.3.1t-5/snpEff.config
/home/jhuang/anaconda3/envs/spandx/bin/snpEff build -genbank wildtype_150      -d
gzip 148_trimmed_P_1.fastq 148_trimmed_P_2.fastq 149_trimmed_P_1.fastq 149_trimmed_P_2.fastq

#-4.3. using spandx calling variants --
#ln -s /home/jhuang/Tools/spandx/ spandx
#-t
snpEff eff -nodownload -no-downstream -no-intergenic -ud 100 -v CP040849.1 noAB_wildtype_trimmed.PASS.snps.vcf > noAB_wildtype_trimmed.PASS.snps.annotated.vcf
#
#               'CP040849'      2659111 Standard
(spandx) nextflow run spandx/main.nf --fastq "trimmed/*_P_{1,2}.fastq.gz" --ref PseudoContig_wildtype.fasta --annotation --database noAB_wildtype -resume
(spandx) nextflow run spandx/main.nf --fastq "trimmed/*_P_{1,2}.fastq.gz" --ref wildtype_150.fasta --annotation --database wildtype_150 -resume

#-4.4. using snippy call variants --
#in bacto-0.1.json
##   "variants_calling": true,
#    "genus": "Klebsiella",
#    "kingdom": "Bacteria",
#    "species": "pneumoniae",
#    "reference": "PROKKA_09062022/PROKKA_09062022.gbk",
snakemake --printshellcmds
#-->extract common SNPs in Outputs/Phylogeny_and_annotation/All_SNPs_indels_annotated.txt, save as All_SNPs_annotated.xlsx

#-4.5. post-processing --
awk '{if($6!=$7) print}' < All_SNPs_indels_annotated.txt > All_SNPs_indels_annotated_.txt
cut -d$'\t' -f1-7 All_SNPs_indels_annotated_.txt > f1_7
grep -v "/" f1_7 > f1_7_
grep -v "\." f1_7_ > f1_7__
grep -v "*" f1_7__ > f1_7___

grep -v "INDEL" f1_7___ > f1_7_SNPs
cut -d$'\t' -f2-2 f1_7_SNPs > f1_7_SNPs_f2
\n --> " All_SNPs_indels_annotated_.txt >> All_SNPs_annotated.txt\ngrep "
#grep "CHROM" All_SNPs_indels_annotated_.txt > All_SNPs_annotated.txt
#grep "46882" All_SNPs_indels_annotated_.txt >> All_SNPs_annotated.txt
#grep "46892" All_SNPs_indels_annotated_.txt >> All_SNPs_annotated.txt

grep -v "SNP" f1_7___ > f1_7_indels
cut -d$'\t' -f2-2 f1_7_indels > f1_7_indels_f2
\n --> " All_SNPs_indels_annotated_.txt >> All_indels_annotated.txt\ngrep "
#grep "CHROM" All_SNPs_indels_annotated_.txt > All_indels_annotated.txt
#grep "49848" All_SNPs_indels_annotated_.txt >> All_indels_annotated.txt
#grep "53872" All_SNPs_indels_annotated_.txt >> All_indels_annotated.txt
```






```sh
#------------------------------------------
#--Always install essential after upgrade--
sudo apt install build-essential
sudo apt install gfortran

#------------------------
#--MultiDisplay Lubuntu--
https://askubuntu.com/questions/657055/dual-monitor-extended-desktop-in-lubuntu
https://help.ubuntu.com/community/Lubuntu/MultiDisplay
https://askubuntu.com/questions/1230924/ubuntu-20-04-does-not-recognize-second-monitor

xrandr --verbose
lspci | egrep 'VGA|3D'
lshw -c video
xrandr
   Screen 0: minimum 8 x 8, current 1920 x 1080, maximum 16384 x 16384
   VGA-0 disconnected (normal left inverted right x axis y axis)
-->DVI-D-0 connected primary 1920x1080+0+0 (normal left inverted right x axis y axis) 521mm x 293mm
   1920x1080     60.00*+  50.00  
   HDMI-0 disconnected (normal left inverted right x axis y axis)
-->VGA-1-0 connected (normal left inverted right x axis y axis)
   1920x1080     60.00 +

sudo apt install arandr
arandr
#try to save the configuration file under ~/.screenlayout$ vim xx.sh
##!/bin/sh
#xrandr --output VGA-0 --off --output DVI-D-0 --primary --mode 1920x1080 --pos 0x0 --rotate normal --output HDMI-0 --off --output VGA-1-0 --mode 1920x1080 --pos 1920x0 --rotate normal --output DVI-D-1-0 --off --output HDMI-1-0 --off
xrandr --output DVI-D-0 --primary --mode 1920x1080 --pos 0x0 --rotate normal          --output VGA-1-0 --mode 1920x1080 --pos 1920x0 --rotate normal     
#-->xrandr: Configure crtc 4 failed
#NOTE: can heat-detected the cables as follows.

sudo apt update
sudo apt remove '^nvidia'
sudo apt autoremove
sudo apt install nvidia-driver-460


#------------------
#--install kernel--
# VERY IMPORTANT: firstly ensuring suitable additional driver for nvidia installed for example 340.108.
# nvidia-driver-470 (proprietary) --> not good, always frozen after few minutes!          
# driver installation via GUI "Additional Drivers" or via command lines as follows.
sudo ubuntu-drivers devices
#and either you install the one you want by hand:
sudo apt install (your driver)
#or you let the system install the recommended ones:
sudo ubuntu-drivers autoinstall
#If you installed any new drivers, do reboot afterwards!

#cd /tmp/
#wget -c https://kernel.ubuntu.com/~kernel-ppa/mainline/v5.7/amd64/linux-headers-5.7.0-050700_5.7.0-050700.202006082127_all.deb
#wget -c https://kernel.ubuntu.com/~kernel-ppa/mainline/v5.7/amd64/linux-headers-5.7.0-050700-generic_5.7.0-050700.202006082127_amd64.deb   
#wget -c https://kernel.ubuntu.com/~kernel-ppa/mainline/v5.7/amd64/linux-image-unsigned-5.7.0-050700-generic_5.7.0-050700.202006082127_amd64.deb   
#wget -c https://kernel.ubuntu.com/~kernel-ppa/mainline/v5.7/amd64/linux-modules-5.7.0-050700-generic_5.7.0-050700.202006082127_amd64.deb

#-before install 5.8.0
total 115M
drwx------ 2 root root  12K Nov 10  2017 lost+found
-rw-r--r-- 1 root root 181K Aug 18  2020 memtest86+_multiboot.bin
-rw-r--r-- 1 root root 181K Aug 18  2020 memtest86+.elf
-rw-r--r-- 1 root root 179K Aug 18  2020 memtest86+.bin
-rw------- 1 root root 4,6M Nov 26  2021 System.map-5.4.0-92-generic
-rw-r--r-- 1 root root 233K Nov 26  2021 config-5.4.0-92-generic
-rw------- 1 root root  14M Nov 26  2021 vmlinuz-5.4.0-92-generic
-rw-r--r-- 1 root root 8,3M Jun 27 17:06 initrd.img-4.10.0-38-generic
-rw-r--r-- 1 root root 8,3M Jun 27 17:06 initrd.img-4.10.0-28-generic
lrwxrwxrwx 1 root root   24 Jun 30 16:18 vmlinuz.old -> vmlinuz-5.4.0-92-generic
-rw-r--r-- 1 root root  80M Jun 30 16:19 initrd.img-5.4.0-92-generic
lrwxrwxrwx 1 root root   24 Jun 30 16:56 vmlinuz -> vmlinuz-5.4.0-92-generic
lrwxrwxrwx 1 root root   27 Jun 30 16:56 initrd.img.old -> initrd.img-5.4.0-92-generic
lrwxrwxrwx 1 root root   27 Jun 30 16:56 initrd.img -> initrd.img-5.4.0-92-generic
drwxr-xr-x 5 root root 1,0K Jun 30 16:56 grub

sudo mv initrd.img-5.4.0-92-generic ../ref/initrd.img-5.4.0-92-generic
sudo ln -s ../ref/initrd.img-5.4.0-92-generic initrd.img-5.4.0-92-generic
# --> total 35M under /boot
#-delete other versions
#dpkg -l 'linux-*' | sed '/^ii/!d;/'"$(uname -r | sed "s/\(.*\)-\([^0-9]\+\)/\1/")"'/d;s/^[^ ]* [^ ]* \([^ ]*\).*/\1/;/[0-9]/!d' | xargs sudo apt-get -y purge
sudo apt -y purge linux-hwe-5.13-headers-5.13.0-52 linux-{headers,image,modules,modules-extra}-5.13.0-52-generic

#sudo dpkg -i *.deb
sudo apt install linux-headers-5.8.0-55-generic
sudo apt install linux-image-5.8.0-55-generic
sudo apt install linux-hwe-5.13-headers-5.13.0-52 linux-{headers,image,modules,modules-extra}-5.13.0-52-generic

total 131M
drwx------ 2 root root  12K Nov 10  2017 lost+found
-rw-r--r-- 1 root root 181K Aug 18  2020 memtest86+_multiboot.bin
-rw-r--r-- 1 root root 181K Aug 18  2020 memtest86+.elf
-rw-r--r-- 1 root root 179K Aug 18  2020 memtest86+.bin
-rw------- 1 root root 5,3M Jun  1  2021 System.map-5.8.0-55-generic
-rw-r--r-- 1 root root 243K Jun  1  2021 config-5.8.0-55-generic
-rw------- 1 root root 9,4M Jun  1  2021 vmlinuz-5.8.0-55-generic
-rw------- 1 root root 4,6M Nov 26  2021 System.map-5.4.0-92-generic
-rw-r--r-- 1 root root 233K Nov 26  2021 config-5.4.0-92-generic
-rw------- 1 root root  14M Nov 26  2021 vmlinuz-5.4.0-92-generic
-rw-r--r-- 1 root root 8,3M Jun 27 17:06 initrd.img-4.10.0-38-generic
-rw-r--r-- 1 root root 8,3M Jun 27 17:06 initrd.img-4.10.0-28-generic
lrwxrwxrwx 1 root root   24 Jun 30 16:18 vmlinuz.old -> vmlinuz-5.4.0-92-generic
lrwxrwxrwx 1 root root   27 Jun 30 16:56 initrd.img.old -> initrd.img-5.4.0-92-generic
lrwxrwxrwx 1 root root   34 Jun 30 16:59 initrd.img-5.4.0-92-generic -> ../ref/initrd.img-5.4.0-92-generic
lrwxrwxrwx 1 root root   24 Jun 30 17:10 vmlinuz -> vmlinuz-5.8.0-55-generic
lrwxrwxrwx 1 root root   27 Jun 30 17:10 initrd.img -> initrd.img-5.8.0-55-generic
-rw-r--r-- 1 root root  82M Jun 30 17:11 initrd.img-5.8.0-55-generic
drwxr-xr-x 5 root root 1,0K Jun 30 17:11 grub
#sudo apt install nvidia-cuda-toolkit


#-----------------------------------------------------------------
#--boot from a live-usb or lubuntu 20.04 (black small USB-stick)--
#http://archive.ubuntu.com/ubuntu/dists/bionic-updates/main/installer-amd64/current/images/netboot/
#https://help.ubuntu.com/community/Installation/MinimalCD#mini_system_in_UEFI_mode
#http://cdimages.ubuntu.com/netboot/focal/
#https://forums.linuxmint.com/viewtopic.php?p=1570850#p1570850
I don't know if there's a simple guide, but I pieced this together from a lot of Googling and testing. If anyone has any corrections to this unprofessional and informal mini-guide, please don't hesitate to correct me, because I'm likely missing something important. :D

For most systems, the basic process should be really simple.
1. Boot into your live USB of Mint
2. Figure out which partition is your root partition. On a default installation of Mint where you just let the installer use the whole drive, it will probably be /dev/sda1. It's easy to check if this is the case. Just open up gparted (this app is already part of the live USB), then look for the partition with the /boot flag. It should be the same partition as your root partition unless you did a special install with a separate /boot partition.
3. Once you've figured out where your root partition is, run this:
sudo mount /dev/sda1 /mnt
Replace /dev/sda1 with whatever the location of your root partition is, if you have a nonstandard installation.
3a. If you have a separate /boot partition somewhere else (let's say /dev/sda4), you can mount that as well:
sudo mount /dev/sda4 /mnt/boot
4. If you need Internet access while doing things to your unbootable installation, such as to download packages, run this:
sudo cp /etc/resolv.conf /mnt/etc/resolv.conf
5. sudo chroot /mnt
6. If everything worked out, you will now see a # symbol instead of a $ symbol in front of your next line. This # symbol means you have chrooted into your installed system, and the commands you run will apply to that system instead of to your live environment. You already have root access when you are chrooted, so you don't need to use sudo. So for example, to remove a package called gtkhash, you'd just run apt-get remove gtkhash, as if you were booted into the system instead of into your live USB. Note that I didn't use sudo.
If you do something special like edit /etc/default/grub, you will also need to do the same things you normally would in the system, but without sudo. So for example, after making a change, I would run update-grub while in chroot.
7. If you forgot something, you can exit the chroot environment (change the # back into a $) by just typing exit. For example, you might do this if you forgot step 4 to get Internet access, and you find that you can't download packages into the chrooted environment. After doing what you forgot, just chroot back in with sudo chroot /mnt
8. That's it! Do what you need to do: remove driver packages, install driver packages, or whatever, then reboot into your fixed system.
sudo mount /dev/sdb1 /mnt
sudo mount /dev/sda2 /mnt/boot
sudo mount /dev/sda1 /mnt/home
sudo mount /dev/sdb5 /mnt/ref
#sudo mount /dev/sdb6 /mnt/unused
sudo mount /dev/sdb7 /mnt/tmp

mount -o bind /dev /mnt/dev
mount -o bind /proc /mnt/proc
mount -o bind /sys /mnt/sys

sudo chroot /mnt
#su jhuang
#top -c --> hinits "mount -t proc proc /proc"
mkinitrd
check, exit, reboot.

#--install graphics driver--
#https://linuxconfig.org/how-to-install-the-nvidia-drivers-on-ubuntu-20-04-focal-fossa-linux
sudo add-apt-repository ppa:graphics-drivers/ppa
ubuntu-drivers devices
sudo ubuntu-drivers autoinstall
or sudo apt install nvidia-driver-390
#https://www.cyberciti.biz/faq/ubuntu-linux-install-nvidia-driver-latest-proprietary-driver/
sudo apt install nvidia-driver-510 nvidia-dkms-510
sudo reboot
#https://www.reddit.com/r/Ubuntu/comments/nuhurp/downgrading_linux_kernel_ubuntu_2004/
#sudo apt install linux-generic
#sudo apt install linux-headers-5.4.0-92-generic
#sudo apt install linux-image-5.4.0-92-generic
#sudo apt -y purge linux-{headers,image,modules}-5.8*
dpkg -l linux-\* | grep ^ii
sudo apt remove linux-hwe-5.8-headers-5.8.0-55
#https://askubuntu.com/questions/75709/how-do-i-install-kernel-header-files
sudo apt remove linux-{headers,image,modules}-5.4-0-101-*
sudo apt install linux-generic

#https://www.cyberciti.biz/faq/installing-latest-stable-mainline-linux-kernel-on-ubuntu-with-apt-get/
#https://askubuntu.com/questions/89710/how-do-i-free-up-more-space-in-boot
```









```sh
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

## 6, PubMLST, ResFinder or RGI calling
```sh
https://cge.cbs.dtu.dk/services/
https://cge.cbs.dtu.dk/services/MLST/
https://pubmlst.org/bigsdb?db=pubmlst_pacnes_seqdef&page=batchProfiles&scheme_id=3
https://cge.food.dtu.dk/services/ResFinder/
https://card.mcmaster.ca/home
https://www.pseudomonas.com/
https://www.pseudomonas.com/amr/list
https://card.mcmaster.ca/analyze/rgi
The antimicrobial resistance genes were predicted by the tool RGI 6.0.0  which utilizes the database CARD 3.2.5 (PMID: 31665441).
!!!!#check the difference of SNPs in the two files, and then look for the SNPs in my tables comparing the two isolates!!!!

cut -d$'\t' -f9-21 BK20399_contigs.fa.txt > f9_21
cut -d$'\t' -f2-6 BK20399_contigs.fa.txt > f2_6
#Best_Hit_ARO --> Antibiotic Resistance Ontology (ARO) Term
paste f9_21 f2_6 > BK20399_.txt
#grep "CARD_Protein_Sequence" BK20399_.txt > header
grep -v "CARD_Protein_Sequence" BK20399_.txt > BK20399__.txt
sort BK20399__.txt > BK20399_sorted.txt

cut -d$'\t' -f9-21 GE3138_contigs.fa.txt > f9_21
cut -d$'\t' -f2-6 GE3138_contigs.fa.txt > f2_6
paste f9_21 f2_6 > GE3138_.txt
grep -v "CARD_Protein_Sequence" GE3138_.txt > GE3138__.txt
sort GE3138__.txt > GE3138_sorted.txt

cut -d$'\t' -f1-6 BK20399_sorted.txt > BK20399_f1_6.txt  
cut -d$'\t' -f1-6 GE3138_sorted.txt > GE3138_f1_6.txt
diff BK20399_f1_6.txt GE3138_f1_6.txt
#-->< APH(3'')-Ib	99.63	3002639	protein homolog model	n/a	n/a

cat header BK20399_sorted.txt > BK20399.txt
cat header GE3138_sorted.txt > GE3138.txt

~/Tools/csv2xls-0.4/csv_to_xls.py BK20399.txt GE3138.txt -d$'\t' -o ARG_calling.xls
```

## 7, calculate mapping table between GWAS and reference Genbank
```sh
echo "##FASTA" gene_presence_absence.csv >> CP023676.gff3
cat CP023676.1.fasta gene_presence_absence.csv >> CP023676.gff3

echo "##FASTA" gene_presence_absence.csv >> CP012351.gff3
cat CP012351.1.fasta gene_presence_absence.csv >> CP012351.gff3

echo "##FASTA" gene_presence_absence.csv >> CP003084.gff3
cat CP003084.1.fasta gene_presence_absence.csv >> CP003084.gff3

#modify unnamed to roary_189_selected_genes

roary -p 15 -f ./roary -i 95 -cd 99 -s -e -n -v  roary_.fa.gb.gff CP012351.gff3 CP003084.gff3 CP023676.gff3 
```


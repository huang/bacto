import os

#######################################################
############### Snakefile Configuration ###############
#######################################################

configfile: "bacto-0.1.json"

SAMPLES, PAIRS, = glob_wildcards("raw_data/{sample}_{pair}.fastq.gz")

SAMPLES = list(set(SAMPLES))

#https://github.com/sanger-pathogens/ariba/wiki/Task%3A-getref
#reference_name is one of: argannot, card, megares, plasmidfinder, resfinder, srst2_argannot, vfdb_core, vfdb_full, virulencefinder
#VFDB: a reference database for bacterial virulence factors.
#VFDB 2016: hierarchical and refined dataset for big data analysis-10 years on"
#ariba getref vfdb_core vfdb_core
#mv vfdb_core.tsv db/vfdb_core/vfdb_core.tsv
#mv vfdb_core.fa db/vfdb_core/vfdb_core.fa
#MODIFIED: trimmomatic threads=6 doesn't work, add hard-coding {threads} as 12 in line 136
#DB = ["argannot", "card", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_core", "vfdb_full", "virulencefinder"]
DB = ["megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_core", "virulencefinder"]

#########################################################################
### Helper functions to access and construct commands from parameters ###
### in the configuration file for this Snakemake Pipeline             ###
#########################################################################

def get_params_trailing():

	return "TRAILING:" + str(config["trimmomatic"]["trailing"])

def get_params_leading():

	return "LEADING:" + str(config["trimmomatic"]["leading"])

def get_params_minlen():

	return "MINLEN:" + str(config["trimmomatic"]["min_len"])

def get_params_illuminaclip():

    params = config["trimmomatic"]

    return "ILLUMINACLIP:" + params["adapter_path"] + ":" \
           + str(params["seed_mismatch"]) + ":" + str(params["palindrome_threshold"]) \
           + ":" + str(params["clip_threshold"]) + ":" + str(params["min_adapter_length"]) \
		   + ":" + str(params["keep_both_reads"])

def get_params_slidingwindow():

    params = config["trimmomatic"]

    return "SLIDINGWINDOW" + ":" + str(params["window_size"]) + ":" + str(params["window_quality"])


def get_genus_options(wildcards):

    if config["prokka"]["genus"]:
        cmd = "--usegenus --genus " + config["prokka"]["genus"]
    else:
        cmd = ""

    return cmd

def get_reference_options(wildcards):

    if config["quast"]["reference"]:
        cmd = "-R " + config["quast"]["reference"]
    else:
        cmd = ""

    return cmd

def get_forward_files(wildcards):

    return "raw_data/{id}_R1.fastq.gz".format(id=wildcards.sample)

def get_reverse_files(wildcards):

    return "raw_data/{id}_R2.fastq.gz".format(id=wildcards.sample)


##############################################
############### Consuming Rule ###############
##############################################

rule all:
    input:
        expand("kraken/{sample}_report.txt", sample=SAMPLES) if config["taxonomic_classifier"] else [],
        #expand("trimmed/{sample}_trimmed_P_1.fastq.gz", sample=SAMPLES) if config["trim"] else [],
        expand("fastqc/{sample}_fastqc.html", sample=SAMPLES) if config["fastqc"] else [],
        expand("shovill/{sample}/contigs.fa", sample=SAMPLES) if config["typing_ariba"] else [],
        expand("prokka/{sample}/{sample}.gbk", sample=SAMPLES) if config["typing_ariba"] else [],
        
        #expand("mykrobe/{sample}/{sample}.json", sample=SAMPLES) if config["typing"] else [],
        expand("ariba/{db}/{sample}.report.txt", db=DB, sample=SAMPLES) if config["typing_ariba"] else [],
        expand("mlst/{sample}.mlst.txt", sample=SAMPLES) if config["assembly"] and config["typing_mlst"] else [],
        
	"roary/gene_presence_absence.csv" if config["pangenome"] else [],
	#"core_alignment/core_alignment_roary.fasta" if config["pangenome"] else [],
	
        "variants/snippy.core.full.aln" if config["variants_calling"] else [],
        "variants/snippy.core.aln" if config["variants_calling"] else [],
        
        "fasttree/snippy.core.tree" if config["phylogeny_fasttree"] else [],
        "raxml-ng/snippy.core.aln.raxml.bestTree" if config["phylogeny_raxml"] else [],
        "gubbins/recomb.final_tree.tre"  if config["recombination"] else [],

        #conda install -c anaconda seaborn
        #conda install libgcc
        #conda install matplotlib biopython numpy pandas
        #"fasttree_matrix.png"

########################################
############### QC Rules ###############
########################################

#if config["trim"]:
rule trimmomatic:
    input:
        fwd = get_forward_files,
        rev = get_reverse_files
    params:
        illumina_clip=get_params_illuminaclip(),
        sliding_window=get_params_slidingwindow(),
        minlen=get_params_minlen(),
        trailing=get_params_trailing(),
        leading=get_params_leading()
    output:
        fwd_p="trimmed/{sample}_trimmed_P_1.fastq",
        fwd_u="trimmed/{sample}_trimmed_U_1.fastq",
        rev_p="trimmed/{sample}_trimmed_P_2.fastq",
        rev_u="trimmed/{sample}_trimmed_U_2.fastq",
        fwd_fq="fastq/{sample}_1.fastq",
        rev_fq="fastq/{sample}_2.fastq"
    threads:
        config["trimmomatic"]["cpu"]
    shell:
        "trimmomatic PE -threads 12 {input.fwd} {input.rev} {output.fwd_p} {output.fwd_u}"
        " {output.rev_p} {output.rev_u} {params.illumina_clip} {params.sliding_window} {params.leading}"
        " {params.trailing} {params.minlen} && "
        "ln -s $(pwd)/{output.fwd_p} $(pwd)/{output.fwd_fq} && ln -s $(pwd)/{output.rev_p} $(pwd)/{output.rev_fq} && "
        "touch -h {output.fwd_fq} && touch -h {output.rev_fq}"

if config["fastqc"]:
    rule fastqc:
        input:
            fwd_after="fastq/{sample}_1.fastq",
            rev_after="fastq/{sample}_2.fastq"
        output:
            "fastqc/{sample}_1_fastqc.html",
            "fastqc/{sample}_2_fastqc.html",
        shell:
            "fastqc --outdir fastqc {input.fwd_after}"
            " && fastqc --outdir fastqc {input.rev_after}"

if config["taxonomic_classifier"]:
    rule kraken:
        input:
            fwd_kra="fastq/{sample}_1.fastq",
            rev_kra="fastq/{sample}_2.fastq"
        params:
            db=config["kraken"]["db_path"]
        output:
            tax="kraken/{sample}_tax.out",
            report="kraken/{sample}_report.txt"
        threads:
            config["kraken"]["cpu"]
        shell:
            "cat {params.db}/database.* > /dev/null"
            " && kraken --db {params.db} --threads {threads} --output {output.tax}"
            " --fastq-input --paired {input.fwd_kra} {input.rev_kra}"
            " && kraken-report --db {params.db} {output.tax} > {output.report}"

##############################################
########### Assembly and Annotation ##########
##############################################

if config["assembly"]:

    rule shovill:
        input:
            forward = "fastq/{sample}_1.fastq",
            reverse = "fastq/{sample}_2.fastq"
        params:
            spades = config["shovill"]["spades"],
            cpu = config["shovill"]["cpu"],
            depth = config["shovill"]["depth"],
            other = config["shovill"]["other"]
        output:
            "shovill/{sample}/contigs.fa"
        shell:
            "shovill --outdir shovill/{wildcards.sample} --depth {params.depth} --force --R1 {input.forward} "
            "--R2 {input.reverse} --mincov 1 --cpus {params.cpu} {params.other}"

    #DEBUG: replace the tbl2asn with the newest version: https://github.com/tseemann/prokka/issues/139 --> SOLVED!
    rule prokka:
        input:
            "shovill/{sample}/contigs.fa"
        params:
            cpu = config["prokka"]["cpu"],
            evalue = config["prokka"]["evalue"],
            genus_options = get_genus_options,
            kingdom = config["prokka"]["kingdom"],
            species = config["prokka"]["species"],
            other = config["prokka"]["other"]
        output:
            "prokka/{sample}/{sample}.gbk",
            "prokka/{sample}/{sample}.gff"
        shell:
             """ 
             prokka --force --outdir prokka/{wildcards.sample} --cpus {params.cpu} {params.genus_options} --kingdom {params.kingdom} --species {params.species} --addgenes --addmrna --prefix {wildcards.sample} --locustag {wildcards.sample} {params.other} {input} -hmm /media/jhuang/Titisee/GAMOLA2/TIGRfam_db/TIGRFAMs_15.0_HMM.LIB
             """
#tbl2asn -V b -a r10k -l paired-ends -M n -N 1 -y 'Annotated using prokka 1.12 from https://github.com/tseemann/prokka' -Z prokka\/HD31_2\/HD31_2\.err -i prokka\/HD31_2\/HD31_2\.fsa
#wget -O tbl2asn.gz ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz
#ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/DOCUMENTATION/VERSIONS
#ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/
#gunzip tbl2asn.gz
#chmod +x tbl2asn
#mv linux64.tbl2asn /home/jhuang/anaconda3/envs/bengal/bin/tbl2asn 
#prokka 1.12 has problem!
#conda install prokka

    if config["typing_mlst"]:
        rule mlst:
            input:
                "shovill/{sample}/contigs.fa"
            output:
                "mlst/{sample}.mlst.txt"
            shell:
                "mlst {input} > {output}"

########################################
############### Typing #################
########################################

if config["typing_ariba"]:
    # Warning. This is some black magic with the LD_LIBRARY_PATH
    # to make MykrobePredictor work in Conda. MCCORTEX31 fails for
    # ZLIB library dependency, which is present in CONDA_PREFIX/lib
    # mykrobe predict --skeleton_dir ./mykrobe/HD31N3_S93 HD31N3_S93 staph --mccortex31_path /home/jhuang/anaconda3/envs/bengal/bin/mccortex31 -1 fastq/HD31N3_S93_1.fastq fastq/HD31N3_S93_2.fastq > mykrobe/HD31N3_S93/HD31N3_S93.json
    #rule mykrobe_predictor:
    #    input:
    #        forward = "fastq/{sample}_1.fastq",
    #        reverse = "fastq/{sample}_2.fastq"
    #    params:
    #        species=config["mykrobe"]["species"]
    #    output:
    #        "mykrobe/{sample}/{sample}.json"
    #    conda:
    #        "envs/mykrobe.yaml"
    #    shell:
    #        "LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH && "
    #        "echo $LD_LIBRARY_PATH && "
    #        "mykrobe predict --skeleton_dir ./mykrobe/{wildcards.sample} {wildcards.sample} {params.species} "
    #        "-1 {input.forward} {input.reverse} > {output}"


    #rule ariba_prepareref:
#        input:
#                gene_dbs="db/{db}/{db}.fa",
#                gene_tsv="db/{db}/{db}.tsv"
#        output:
#                gene_preps="db/{db}/prepareref.{db}"
#        shell:
#                "ariba prepareref -f {input.gene_dbs} -m {input.gene_tsv} {output.gene_preps}"

#argannot, card, megares, plasmidfinder, resfinder, srst2_argannot, vfdb_core, vfdb_full, virulencefinder
##ariba getref argannot argannot  
##mkdir db/argannot
##mv argannot* db/argannot
##ariba prepareref -f db/argannot/argannot.fa -m db/argannot/argannot.tsv db/argannot/prepareref.argannot
##-------
##ariba getref card card  
##mkdir db/card
##mv card* db/card
##ariba prepareref -f db/card/card.fa -m db/card/card.tsv db/card/prepareref.card
##-------
##ariba getref megares megares  
##mkdir db/megares
##mv megares* db/megares
##ariba prepareref -f db/megares/megares.fa -m db/megares/megares.tsv db/megares/prepareref.megares
##-------
##ariba getref plasmidfinder plasmidfinder  
##mkdir db/plasmidfinder
##mv plasmidfinder* db/plasmidfinder
##ariba prepareref -f db/plasmidfinder/plasmidfinder.fa -m db/plasmidfinder/plasmidfinder.tsv db/plasmidfinder/prepareref.plasmidfinder
##-------
##ariba prepareref -f db/resfinder/resfinder.fa -m db/resfinder/resfinder.tsv db/resfinder/prepareref.resfinder
##-------
##ariba getref srst2_argannot srst2_argannot  
##mkdir db/srst2_argannot
##mv srst2_argannot* db/srst2_argannot
##ariba prepareref -f db/srst2_argannot/srst2_argannot.fa -m db/srst2_argannot/srst2_argannot.tsv db/srst2_argannot/prepareref.srst2_argannot
##-------
##ariba prepareref -f db/vfdb_core/vfdb_core.fa -m db/vfdb_core/vfdb_core.tsv db/vfdb_core/prepareref.vfdb_core
##-------
##ariba getref vfdb_full vfdb_full  
##mkdir db/vfdb_full
##mv vfdb_full* db/vfdb_full
##ariba prepareref -f db/vfdb_full/vfdb_full.fa -m db/vfdb_full/vfdb_full.tsv db/vfdb_full/prepareref.vfdb_full
##-------
##ariba getref virulencefinder virulencefinder  
##mkdir db/virulencefinder
##mv virulencefinder* db/virulencefinder
##ariba prepareref -f db/virulencefinder/virulencefinder.fa -m db/virulencefinder/virulencefinder.tsv db/virulencefinder/prepareref.virulencefinder
#

    rule ariba:
        input:
                prepref = "db/{db}/prepareref.{db}",
                forward = "fastq/{sample}_1.fastq",
                reverse = "fastq/{sample}_2.fastq"
        params:
                outdir="ariba/{db}/{sample}",
                report="ariba/{db}/{sample}/report.tsv"
        output:
                "ariba/{db}/{sample}.report.txt"
        shell:
                "ariba run --force {input.prepref} {input.forward} {input.reverse} {params.outdir} && "
                "mv {params.report} {output}"
        

#####################################################
############### Pangenome using roary ###############
#####################################################

#http://sepsis-omics.github.io/tutorials/modules/roary/
if config["pangenome"]:
    rule roary:
        input:
            expand("prokka/{sample}/{sample}.gff", sample=SAMPLES)
        params:
            core=config["roary"]["core"],
            identity=config["roary"]["identity"],
            other=config["roary"]["other"]
        output:
            "roary/gene_presence_absence.csv"
            #"core_alignment/core_alignment_roary.fasta"
        threads:
            config["roary"]["cpu"]
        shell:
            #{threads} doesn't work, using hard-coding 15
            "roary -p 15 -f ./roary -i {params.identity} -cd {params.core} -s -e -n -v {params.other} {input} && "
            "mv roary_*/* ./roary && rm -rf ./roary_*"
            #"cp roary/core_gene_alignment.aln core_alignment/core_alignment_roary.fasta"
            #roary -e --mafft -p 8 *.gff
            #roary -p 4 -f {params.outdir} -s -e -n -v {input} && touch {output}	
	
########################################
############### Variants ###############
########################################

# Snippy in the environment requires the following:
# Install latest Python 3.5 freebayes=1.1.0=3 from BioConda
# Install bzip2 from Conda-Forge (add to channels in YAML)
# Then install standard Snippy from BioConda - for some reason the
# normal install does not work unless using this sequence of installations.

#TODO: using own refs from this STEP, using panphlan
if config["variants_calling"]:
    rule snippy:
        input:
            forward="fastq/{sample}_1.fastq",
            reverse="fastq/{sample}_2.fastq"
        params:
            outdir = "snippy/{sample}",
            reference = config["snippy"]["reference"],
            mincov = config["snippy"]["mincov"],
            minfrac = config["snippy"]["minfrac"],
            mapqual = config["snippy"]["mapqual"],
            other = config["snippy"]["other"],
            cpu = config["snippy"]["cpu"]
        output:
            "snippy/{sample}/{sample}.txt",
            #"snippy/{sample}/{sample}.depth"
        shell:
            "LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH && "
            "snippy --force --outdir {params.outdir} --ref {params.reference} "
            "--R1 {input.forward} --R2 {input.reverse} --cpus {params.cpu} "
            "--mincov {params.mincov} --minfrac {params.minfrac} "
            "--mapqual {params.mapqual} --prefix {wildcards.sample} {params.other} "
            #"&& gzip -d snippy/{wildcards.sample}/{wildcards.sample}.depth.gz"
            
    rule snippy_core:
        input:
            expand("snippy/{sample}/{sample}.txt", sample=SAMPLES)
        params:
            reference = config["snippy"]["reference"]
        output:
            "variants/snippy.core.full.aln",
            "variants/snippy.core.aln"
        shell:
            "snippy-core --ref {params.reference} --prefix snippy.core snippy/* && "
            "mv snippy.core.* variants/"

##########################################################################
########### Set Analysis based on variants/snippy.core.aln ###############
##########################################################################
if config["phylogeny_fasttree"]:
    rule fasttree_snp:
        input:
            "variants/snippy.core.aln"
        output:
            "fasttree/snippy.core.tree"
        shell:
            "FastTree -gtr -nt {input} > {output}"

#DEBUG: threads of raxml_ng doesn't work!
if config["phylogeny_raxml"]:
    rule raxml_ng:
        input:
            "variants/snippy.core.aln"
        params:
            model = config["raxml_ng"]["model"],
            correction = config["raxml_ng"]["correction"],
            bootstrap = config["raxml_ng"]["bootstrap"],
            other = config["raxml_ng"]["other"]
        output:
            "raxml-ng/snippy.core.aln.raxml.bestTree"
        threads:
            config["raxml_ng"]["cpu"]
        shell:
	    #{threads} doesn't work, using hard-coding 8
            "raxml-ng --all --model {params.model}{params.correction} --prefix raxml-ng/snippy.core.aln --threads 8 "
            "--msa {input} --bs-trees {params.bootstrap} {params.other}"

if config["recombination"]:
    rule gubbins:
        input:
            "variants/snippy.core.aln"
        params:
            model = config["gubbins"]["model"],
            prefix = "recomb",
            tree_builder = config["gubbins"]["tree_builder"],
            iterations = config["gubbins"]["iterations"],
            min_snps = config["gubbins"]["min_snps"],
            min_window_size = config["gubbins"]["min_window_size"],
            max_window_size = config["gubbins"]["max_window_size"],
            filter_percentage = config["gubbins"]["filter_percentage"],
            other = config["gubbins"]["other"]
        output:
            "gubbins/recomb.final_tree.tre"
        shell:
            "run_gubbins.py --tree_builder {params.tree_builder} --iterations {params.iterations} --raxml_model "
            "{params.model} --min_snps {params.min_snps} --min_window_size {params.min_window_size} "
            "--max_window_size {params.max_window_size} --filter_percentage {params.filter_percentage} "
            "--prefix {params.prefix} {params.other} {input} && mv {params.prefix}* gubbins"

            

#cp variants/snippy.core.full.aln variants/snippy.core.full.with_ref.aln
#cp variants/snippy.core.aln variants/snippy.core.with_ref.aln
#rule visualize_roary:
#    input:
#        tree = "fasttree/snippy.core.tree",
#        csv = "roary/gene_presence_absence.csv",
#        plot_script_fp = "local/roary_plots"
#    output:
#        "fasttree_matrix.png"
#    shell:
#        """
#        python {input.plot_script_fp}/roary_plots.py {input.tree} {input.csv}
#        mv pangenome_matrix.png matrix_fasttree.png
#        """
rule visualize_roary_fasttree:
    input:
        tree = "fasttree/snippy.core.tree",
        csv = "roary/gene_presence_absence.csv",
        plot_script_fp = "local/roary_plots"
    output:
        "visualize.fasttree.done"
    shell:
        """
        python {input.plot_script_fp}/roary_plots.py {input.tree} {input.csv} && touch {output}
        mv pangenome_matrix.png pangenome_matrix_fasttree.png
        """

rule visualize_roary_raxml:
    input:
        tree = "raxml-ng/snippy.core.aln.raxml.bestTree",
        csv = "roary/gene_presence_absence.csv",
        plot_script_fp = "local/roary_plots"
    output:
        "visualize.raxml.done"
    shell:
        """
        python {input.plot_script_fp}/roary_plots.py {input.tree} {input.csv} && touch {output}
        mv pangenome_matrix.png pangenome_matrix_raxml.png
        """
  
rule visualize_roary_gubbins:
    input:
        tree = "gubbins/recomb.final_tree.tre",
        csv = "roary/gene_presence_absence.csv",
        plot_script_fp = "local/roary_plots"
    output:
        "visualize.gubbins.done"
    shell:
        """
        python {input.plot_script_fp}/roary_plots.py {input.tree} {input.csv} && touch {output}
        mv pangenome_matrix.png pangenome_matrix_gubbins.png
        """

rule run_piggy:
    input:
        gff="prokka_gffs/",
        roa="roary/"
    output:
        "piggy.done"
    params:
        outdir = "piggy"
    shell:
        """
        piggy -i {input.gff} -r {input.roa} -o {params.outdir} && touch {output}
        """

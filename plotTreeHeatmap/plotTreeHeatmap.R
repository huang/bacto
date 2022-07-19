# https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html#fig:curatedgene
# 
# This example introduces annotating a tree with various sources of data (e.g., location, sampling year, curated genotype information, etc.).
# 
# ## use the following command to get the data from Holt Lab
# ## git clong https://github.com/katholt/plotTree.git
# .libPaths()
# 
# TODO:
# - turn the tree of IA2 to another directory. Change the tree information directly!
# 
# (((BK27703:0.000644,BK23521:0.000846):0.000024,((BK003180:0.010535,(BK14772:0.010045,((((BK16920:0.000001,BK16922:0.000001):0.000001,BK16921:0.000001):0.000202,(BK32250:0.000120,E53234:0.000099):0.000001):0.012081,((E50862:0.001032,(((BK3089:0.000446,(X33027:0.000222,(BK13142-18:0.000121,BK15782:0.000303):0.000020):0.000101):0.000764,(BK23837:0.000424,BK2371:0.000404):0.000735):0.187659,( (((((((GE443-21:0.000682,BK33325:0.000141):0.000001,(BK12403:0.000020,BK34618:0.000020):0.000162):0.000183,((BK13171:0.000162,((BK21716:0.000040,(E41351:0.000162,BK020067-19:0.000222):0.000020):0.000001,E54403:0.000001):0.000020):0.000001,D4102:0.000323):0.000101):0.000305,((E56292:0.000456,E50388:0.000543):0.000001,BK002172:0.000769):0.000023):0.000083,(E47914:0.000764,(BK020672-18:0.000141,BK13406-18:0.000202):0.000324):0.000021):0.000231,((((A3794:0.000247,E47727:0.000182):0.000001,D6288:0.000343):0.000123,X32299:0.000243):0.000826,(((((E38935:0.000375,(BK004881:0.000001,BK004880:0.000020):0.000263):0.000395,VA40274-18:0.000566):0.000020,((((A2034:0.000283,(BK010583-18:0.000267,BK22569:0.000502):0.000001):0.000023,BK21333:0.000283):0.000143,E44070:0.000525):0.000001,(BK32630:0.003713,BK008352-18:0.000458):0.000068):0.000062):0.000041,(BK16937:0.000222,((E50389:0.002053,(((((BK24000:0.000358,BK002557:0.000182):0.000001,BK14321:0.000323):0.000126,(((BK008238:0.000525,BK29927:0.000445):0.000060):0.000001,((BK16375:0.000222,BK005978:0.001699):0.000061,BK27248:0.000424):0.000040):0.000040):0.000010,BK26765:0.000344):0.000010,((X28484:0.000384,(E40836:0.000222,(X34222:0.000121,E37403:0.000280):0.000022):0.000021):0.000121,A010781:0.000384):0.000001):0.000020):0.000040,(BK002600:0.000101,BK14441:0.000182):0.000162):0.000020):0.000040):0.000021,(BK1494:0.000603,((((((X32672:0.000182,E37240:0.000340):0.000001,BK007118:0.000283):0.000001,E47153:0.000141):0.000001,(X33762:0.000143,BK26618:0.000182):0.000061):0.000020,BK012500:0.000182):0.000061,((E54438:0.000001,E54438_20:0.000001):0.000222,(X35089:0.000776,X35718:0.000283):0.000222):0.000101):0.000040):0.000043):0.000553):0.000229):0.001906,E50391:0.002240):0.004950, (((((BK33082:0.000325,(BK19775:0.000101,D515:0.000445):0.000019):0.002025,BK32277:0.001653):0.004822,((((((((((((BK9651:0.000203,E40763:0.000344):0.000020,C6573:0.000122):0.001359,(((BK14636:0.000445,(BK27232:0.000390,E39190:0.007252):0.000056):0.000061,(BK002523:0.000263,(D8004:0.000373,E47207:0.000344):0.000093):0.000020):0.000001,X37677:0.000344):0.000020):0.000081,((E47475:0.000729,(E40618:0.000364,((BK14042:0.000271,BK012799:0.001223):0.000357,(BK14738-18:0.001464,E41797:0.000304):0.000082):0.000083):0.000001):0.000001,(E44849:0.000020,E43256:0.000041):0.000527):0.000021):0.000001,BK003184:0.000770):0.000001,(E41484:0.000568,((E52542:0.000001,E52430:0.000020):0.000203,X28149:0.000182):0.000466):0.000101):0.000001,E42226:0.000527):0.000001,(D6657:0.000487,C009283:0.000790):0.000081):0.000001,((((((E40343:0.000263,BK28726:0.000284):0.000101,E52304:0.000527):0.000101,E41905:0.000446):0.000001,E42498:0.000689):0.000041,((E44126:0.000344,E45124:0.000791):0.000001,A010961:0.000466):0.000020):0.000001,(E50396:0.001175,BK011007:0.000628):0.000001):0.000041):0.000383,(X33874:0.000710,BK12722:0.000203):0.000196):0.000562,BK4738:0.000905):0.003539,(((X37824:0.000407,(BK20140:0.000668,BK35047:0.000376):0.000167):0.000040,((((BK24378:0.000202,BK34454:0.000960):0.000235,X36694:0.000264):0.000061,BK31654:0.000506):0.000020,((E44544:0.000263,D4792:0.000406):0.000001,X32508:0.000222):0.000122):0.000042):0.000186,X37153:0.000886):0.007939):0.002202):0.001500):0.006725,(((((E38448:0.002073,D3376:0.000141):0.000001,(((BK32875:0.000020,BK31875:0.000001):0.000182,(E37903:0.008018,E45005:0.000222):0.000001):0.000020,(E50859:0.000001,E50861:0.000001):0.000243):0.000001):0.000001,(GE506-20:0.000222,(BK30307:0.000426,(BK002382:0.000283,((X37209:0.000222,E37795:0.000344):0.000040,(E43122:0.000001,E43120:0.000122):0.000509):0.000122):0.000001):0.000001):0.000081):0.000837,BK003501:0.000981):0.002777,(D3286:0.000725,BK28348:0.001080):0.002133):0.008249):0.006791 ):0.041338):0.038080):0.041448,((((((E43431:0.000080,E56457:0.000161):0.000001,(((E46791:0.000001,E46789:0.000001):0.000060,E54491:0.000100):0.000181,(P35022:0.000040,BK14141-2-18:0.000743):0.000031):0.000029):0.000001,(E40247:0.000100,Ask:0.000080):0.000020):0.000021,E38240:0.000851):0.000464,((BK17415-18:0.000121,E40879:0.000583):0.000241,X35753:0.000325):0.001464):0.001000,((((X36290:0.000241,E47627:0.000140):0.000120,D770:0.000201):0.000001,E55198:0.000241):0.000805,(E42985:0.000001,E42986:0.000059):0.000805):0.000334):0.000913):0.003425):0.004137):0.001462):0.001653,(BK30006:0.000463,E38486:0.000585):0.000951):0.000271):0.001185,GE1700-19:0.004409,X31274:0.000916):0.0;
# 
# ((((BK33151:0.000384,(X34110:0.000262,(A5284:0.000142,GE738-20:0.000081):0.000081):0.000041):0.001225,BK12078:0.001561):0.000365,BK29232:0.001711):0.006892,E39825:0.004146):0.003604,
# 
# 
# 
# 
# ( 
# -1- ((((BK33082:0.000325,(BK19775:0.000101,D515:0.000445):0.000019):0.002025,BK32277:0.001653):0.004822,((((((((((((BK9651:0.000203,E40763:0.000344):0.000020,C6573:0.000122):0.001359,(((BK14636:0.000445,(BK27232:0.000390,E39190:0.007252):0.000056):0.000061,(BK002523:0.000263,(D8004:0.000373,E47207:0.000344):0.000093):0.000020):0.000001,X37677:0.000344):0.000020):0.000081,((E47475:0.000729,(E40618:0.000364,((BK14042:0.000271,BK012799:0.001223):0.000357,(BK14738-18:0.001464,E41797:0.000304):0.000082):0.000083):0.000001):0.000001,(E44849:0.000020,E43256:0.000041):0.000527):0.000021):0.000001,BK003184:0.000770):0.000001,(E41484:0.000568,((E52542:0.000001,E52430:0.000020):0.000203,X28149:0.000182):0.000466):0.000101):0.000001,E42226:0.000527):0.000001,(D6657:0.000487,C009283:0.000790):0.000081):0.000001,((((((E40343:0.000263,BK28726:0.000284):0.000101,E52304:0.000527):0.000101,E41905:0.000446):0.000001,E42498:0.000689):0.000041,((E44126:0.000344,E45124:0.000791):0.000001,A010961:0.000466):0.000020):0.000001,(E50396:0.001175,BK011007:0.000628):0.000001):0.000041):0.000383,(X33874:0.000710,BK12722:0.000203):0.000196):0.000562,BK4738:0.000905):0.003539,(((X37824:0.000407,(BK20140:0.000668,BK35047:0.000376):0.000167):0.000040,((((BK24378:0.000202,BK34454:0.000960):0.000235,X36694:0.000264):0.000061,BK31654:0.000506):0.000020,((E44544:0.000263,D4792:0.000406):0.000001,X32508:0.000222):0.000122):0.000042):0.000186,X37153:0.000886):0.007939):0.002202):0.001500):0.006725, 
# -2- (((((E38448:0.002073,D3376:0.000141):0.000001,(((BK32875:0.000020,BK31875:0.000001):0.000182,(E37903:0.008018,E45005:0.000222):0.000001):0.000020,(E50859:0.000001,E50861:0.000001):0.000243):0.000001):0.000001,(GE506-20:0.000222,(BK30307:0.000426,(BK002382:0.000283,((X37209:0.000222,E37795:0.000344):0.000040,(E43122:0.000001,E43120:0.000122):0.000509):0.000122):0.000001):0.000001):0.000081):0.000837,BK003501:0.000981):0.002777,(D3286:0.000725,BK28348:0.001080):0.002133):0.008249, 
# -3- ((((BK33151:0.000384,(X34110:0.000262,(A5284:0.000142,GE738-20:0.000081):0.000081):0.000041):0.001225,BK12078:0.001561):0.000365,BK29232:0.001711):0.006892,E39825:0.004146):0.003604  
# ):0.006791
# 
# ---->
# 
# (((((BK33082:0.000325,(BK19775:0.000101,D515:0.000445):0.000019):0.002025,BK32277:0.001653):0.004822,((((((((((((BK9651:0.000203,E40763:0.000344):0.000020,C6573:0.000122):0.001359,(((BK14636:0.000445,(BK27232:0.000390,E39190:0.007252):0.000056):0.000061,(BK002523:0.000263,(D8004:0.000373,E47207:0.000344):0.000093):0.000020):0.000001,X37677:0.000344):0.000020):0.000081,((E47475:0.000729,(E40618:0.000364,((BK14042:0.000271,BK012799:0.001223):0.000357,(BK14738-18:0.001464,E41797:0.000304):0.000082):0.000083):0.000001):0.000001,(E44849:0.000020,E43256:0.000041):0.000527):0.000021):0.000001,BK003184:0.000770):0.000001,(E41484:0.000568,((E52542:0.000001,E52430:0.000020):0.000203,X28149:0.000182):0.000466):0.000101):0.000001,E42226:0.000527):0.000001,(D6657:0.000487,C009283:0.000790):0.000081):0.000001,((((((E40343:0.000263,BK28726:0.000284):0.000101,E52304:0.000527):0.000101,E41905:0.000446):0.000001,E42498:0.000689):0.000041,((E44126:0.000344,E45124:0.000791):0.000001,A010961:0.000466):0.000020):0.000001,(E50396:0.001175,BK011007:0.000628):0.000001):0.000041):0.000383,(X33874:0.000710,BK12722:0.000203):0.000196):0.000562,BK4738:0.000905):0.003539,(((X37824:0.000407,(BK20140:0.000668,BK35047:0.000376):0.000167):0.000040,((((BK24378:0.000202,BK34454:0.000960):0.000235,X36694:0.000264):0.000061,BK31654:0.000506):0.000020,((E44544:0.000263,D4792:0.000406):0.000001,X32508:0.000222):0.000122):0.000042):0.000186,X37153:0.000886):0.007939):0.002202):0.001500):0.005725,(((((E38448:0.002073,D3376:0.000141):0.000001,(((BK32875:0.000020,BK31875:0.000001):0.000182,(E37903:0.008018,E45005:0.000222):0.000001):0.000020,(E50859:0.000001,E50861:0.000001):0.000243):0.000001):0.000001,(GE506-20:0.000222,(BK30307:0.000426,(BK002382:0.000283,((X37209:0.000222,E37795:0.000344):0.000040,(E43122:0.000001,E43120:0.000122):0.000509):0.000122):0.000001):0.000001):0.000081):0.000837,BK003501:0.000981):0.002777,(D3286:0.000725,BK28348:0.001080):0.002133):0.007249):0.001000,((((BK33151:0.000384,(X34110:0.000262,(A5284:0.000142,GE738-20:0.000081):0.000081):0.000041):0.001225,BK12078:0.001561):0.000365,BK29232:0.001711):0.006892,E39825:0.004146):0.003604  
# 
# 
# 
# Hi Jiabin!
# - How are you? I am now finally writing the C. acnes paper before I leave for my mission. 
# - I was wondering if you could send me a list of all the strains that went into the GWAS analysis from table  If_inf_02_09_2021_1151.results.xlsx and their clonal complexes and phylotype and if you sorted them into a certain group (I am not quite sure what the GROUP 1 to 3 is. Do they correspond to the phylotypes?).
# - I would also like to make a circle phylogenetic tree with some colour code for the infection isolates (all that are not BK or X, 1 as in column 2 of 2021_Metadata_C. acnes.xlsx) and the commensal isolates (BK and X, 0 as in column 2 of 2021_Metadata_C. acnes.xlsx).
# 
# phylotypes
# 
# ST
# Clonal complex
# Phylotype
# Type
# 
# #https://www.nature.com/articles/s41598-017-04081-1
# Commensal-to-pathogen transition: One-single transposon insertion results in two pathoadaptive traits in Escherichia coli -macrophage interaction
# commensal pathogen typing
# 
# #Commensals are those type of microbes that reside on either surface of the ...
# Type
# Commensals
# pathogens
# 
# infection
# commensal
# 
# 
# #E39883 is a CC5 Type 1B but contains 2 sodA alleles
# BK012348 has only partial sequences for gmk and guaA (possible ST49).
# 
# 
# E39883 is a CC5 Type 1B but contains 2 sodA alleles
# BK012348 has only partial sequences for gmk and guaA (possible ST49).




#BiocManager::install("ggtree")
library(ggtree)
#library(ggimage)

setwd("~/DATA/Data_Anna_C.acnes/plotTree/")
#cp ~/DATA/Data_Anna_C.acnes/182samples_roaries/roary_182s_95/core_gene_alignment.tree ./
#cp ~/DATA/Data_Anna_C.acnes/BIGSdb_040564_9438765207_96860_complete2.csv ./


##Interpolation
#E50396 --> 1,CC1,IA1
#X35718 --> 5,CC5,IB
#X35089 --> 5,CC5,IB
#BK32630 --> 5,CC5,IB

#info2 <- read.csv("info.csv")
info <- read.csv("completed2_interpolated_kick2.csv")
info$name <- info$Isolate
#tree <- read.tree("tree.nwk")
tree <- read.tree("core_gene_raxml.kick2.tree")
cols <- c(infection='purple2', commensal='skyblue2')

##offset=4, hjust=0.5, , yscale
#ggtree(tree, layout='circular', branch.length='none') %<+% info + 
#  geom_tippoint(aes(color=Type)) + 
#  scale_color_manual(values=cols) + geom_tiplab2(aes(label=name), align=TRUE, linetype=NA, hjust=1.5, size=2.5, offset=10)  #
##  + geom_tiplab2(aes(label=year), align=T, linetype=NA, size=2, offset=8, hjust=0.5)


#- The tree was visualized in circular layout and attached with the strain sampling location information. 
#- A geom_tippoint layer added circular symbolic points to tree tips and colored them by their locations. 
#- Two geom_tiplab2 were added to display taxon names and sampling years.


#heatmapData=read.csv("res_genes.csv", row.names=1)
#rn <- rownames(heatmapData)
#heatmapData <- as.data.frame(sapply(heatmapData, as.character))
#rownames(heatmapData) <- rn

library(dplyr)
heatmapData2 <- info %>% select(Isolate, ST, Clonal.Complex, Phylotype)
rn <- heatmapData2$Isolate
heatmapData2$Isolate <- NULL
heatmapData2 <- as.data.frame(sapply(heatmapData2, as.character))
rownames(heatmapData2) <- rn

#heatmapData3 <- info %>% select(Isolate, ST)
#rn <- heatmapData3$Isolate
#heatmapData3$Isolate <- NULL
#heatmapData3 <- as.data.frame(sapply(heatmapData3, as.character))
#rownames(heatmapData3) <- rn

#heatmapData4 <- info %>% select(Isolate, Clonal.Complex)
#rn <- heatmapData4$Isolate
#heatmapData4$Isolate <- NULL
#heatmapData4 <- as.data.frame(sapply(heatmapData4, as.character))
#rownames(heatmapData4) <- rn

#heatmap.colours <- c("white","grey","seagreen3","darkgreen")
#heatmap.colours2 <- c("green","brown","tan","red","orange",
#                    "pink","magenta","purple","blue","skyblue3",
#                    "blue","skyblue2")
#names(heatmap.colours) <- 0:3
#names(heatmap.colours2) <- 4:15

#-"purple",
heatmap.colours <- c("cornflowerblue","darkgreen","seagreen3","tan","red","green","orange","pink","purple","magenta","brown", "darksalmon","chocolate4","darkkhaki",   "blue","cyan", "skyblue2", "azure3","blueviolet","darkgoldenrod",       "cornflowerblue","darkgreen","seagreen3","tan","red","green","orange","pink","brown","magenta",     "cornflowerblue","darkgreen","red","tan","brown",      "darkgrey")

#"CC1 (type IA1)", "CC2 (type IA2)", "CC3 (type IA1)","CC4 (type IA1)","CC5 (type IB)","CC6 (type II)","CC30 (previously CC72) (type II)","CC72 (type II)","CC77 (type III)", "Singleton (type IA2)", "ST linking CC1 & CC4",
                    #"coral1","chocolate4","chartreuse4","cadetblue4","burlywood4","cornflowerblue","cornsilk4","antiquewhite3","bisque4","darkkhaki","coral4",   
                    #"seagreen3","darkgreen","green","brown", "tan","red","orange","pink", "magenta","purple","blue","darkgrey")
#-"ST linking CC1 & CC4",
names(heatmap.colours) <- c("1","2","3","4","5", "20","21","22",   "28","30","33","42","52","53",  "100","105","124","133","134","135",   "CC1","CC2","CC3","CC4","CC5","CC6","CC30","CC72","CC77","Singleton",    "IA1","IA2","IB","II","III",    "NA")



#https://deep-dv.org/wp/

#heatmapData4, 
#, branch.length='none'

#svg("xxx.svg")
#%<+%
#geom_tiplab!=geom_tiplab2(size=3.3, ... colour='pink',  
p <- ggtree(tree, layout='circular', branch.length='none') %<+% info + 
  geom_tippoint(aes(color=Type)) + 
  scale_color_manual(values=cols) + geom_tiplab(aes(label=name), geom='text', align=TRUE,  linetype=NA, hjust=1.8, offset=10,check.overlap=TRUE, size=3.3)
  #difference between geom_tiplab and geom_tiplab2?
  #+ theme(axis.text.x = element_text(angle = 30, vjust = 0.5)) + theme(axis.text = element_text(size = 20))  + scale_size(range = c(1, 20))
  #font.size=10, 
png("ggtree.png", width=1200, height=1200)
p
dev.off()

png("ggtree_and_gheatmap.png", width=1200, height=1200)
#  , color=NULL, offset=10, colnames_offset_y = 1, hjust=0
gheatmap(p, heatmapData2, width=0.1, 
         colnames_position="top", 
         colnames_angle=90, colnames_offset_y = 0.1, 
         hjust=0.5, font.size=4, offset = 8) + scale_fill_manual(values=heatmap.colours) + theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size = 14)) + guides(fill=guide_legend(title=""))
  #+ guides(color = guide_legend(override.aes = list(size = 10))) 
  #, breaks=0:100
  #scale_fill_manual(values=heatmap.colours2, breaks=4:15)
dev.off()  




#- The curated gene information was further loaded and plotted as a heatmap using gheatmap function with customized colors. 
#- The final figure was demonstrated in Figure 4.17.






# ###########################################################################################################
# ## finding the next strain with Phylogenetics: send both HCV231_all.png and HCV231_all.pdf to the Nicole ##
# #1, generate tree
# fasttree -gtr -nt variants_182/snippy.core.aln > snippy.core.fast.tree
# mkdir raxml-ng_182
# raxml-ng --all --model GTR+G+ASC_LEWIS --prefix raxml-ng_182/snippy.core.aln --threads 1 --msa variants_182/snippy.core.aln --bs-trees 1000 --redo
# 
# mafft --adjustdirection RSVB_0.1.fasta > RSVB_0.1.aln
# snp-sites RSVB_0.1.aln -o RSVB_0.1_.aln
# fasttree -gtr -nt RSVB_0.1_.aln > RSVB_0.1.tree
# raxml-ng --all --model GTR+G+ASC_LEWIS --prefix raxml-ng --threads 1 --msa RSVB_0.1_.aln --bs-trees 1000 --redo
# 
# 
# 
# 
# 
# #https://4va.github.io/biodatasci/r-ggtree.html
# BiocManager::install("ggtree")
# 
# library(tidyverse)
# library(ggtree)
# 
# tree <- read.tree("snippy.core.fast.tree")
# tree
# 
# #png(file="tree.png", width=1000, height=1000)
# #ggtree(tree) + 
# #  geom_tiplab() + 
# #  geom_cladelabel(node=108, label="Some random clade", 
# #                  color="red2", offset=.8, align=TRUE) + 
# #  geom_cladelabel(node=107, label="A different clade", 
# #                  color="blue", offset=.8, align=TRUE) + 
# #  theme_tree2() + 
# #  xlim(0, 70) + 
# #  theme_tree()
# #dev.off()
# 
# #https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/ggtree/inst/doc/treeManipulation.html
# png(file="tree_tiplab.png", width=5000, height=5000)
# #ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
# ggtree(tree) + geom_text2(aes(subset=isTip, label=node), hjust=-.3) + geom_tiplab()
# dev.off()
# 
# 
# #TODO: draw coverage plot also with R-package
# #https://yulab-smu.top/treedata-book/chapter4.html
# #theme_tree2() + 
# #  geom_hilight(node=88, fill="red2") + 
# #  geom_hilight(node=90, fill="blue") +
# png(file="tree.png", width=2000, height=1800)
# ggtree(tree, layout="circular") + xlim(-2, NA) +
#   geom_tiplab(align=FALSE, offset=.1, linesize=.4, size=5) + 
#   geom_cladelabel(node=108, fontsize=5, barsize=0.5, angle=0, label="Cases", 
#                   color="red2", offset=.8, align=TRUE) +
#   geom_cladelabel(node=107, fontsize=5, barsize=0.5, angle=0, label="Reference cohort", 
#                   color="blue", offset=.8, align=TRUE) +
#   geom_hilight(node=107, fill="purple")
# dev.off()
# 
# #
# png(file="tree.png", width=1000, height=900)
# #svg(file="tree.svg")
# ggtree(tree, layout="circular") + xlim(-6, NA) +
#   geom_tiplab(align=FALSE, offset=.1, linesize=.4, size=5) + 
#   geom_cladelabel(node=203, fontsize=5, barsize=0.5, label="Cases", 
#                   color="red2", offset=1.6, align=TRUE) +
#   geom_cladelabel(node=205, fontsize=5, barsize=0.5, label="Reference cohort", 
#                   color="blue", offset=1.6, align=TRUE)
# dev.off()
# 
# 
# png(file="tree.png", width=1100, height=1200)
# #svg(file="tree.svg")
# ggtree(tree) +
#   geom_tiplab(align=FALSE, offset=.0, linesize=.4, size=4,) + 
#   geom_cladelabel(node=203, fontsize=4, barsize=0.5, label="Cases", angle=58,  
#                   color="red2", offset=0.03, align=TRUE) +
#   geom_cladelabel(node=205, fontsize=4, barsize=0.5, label="Reference cohort", angle=58,  
#                   color="blue", offset=0.03, align=TRUE)
#   + theme_tree2()
# dev.off()
# 
# 
# 
# #for manuscripts_v2
# 
# library(tidyverse)
# library(ggtree)
# library(ape)
# 
# tree <- read.tree("snippy.core.aln.raxml.bestTree")
# tree
# ##https://bioconductor.statistik.tu-dortmund.de/packages/3.3/bioc/vignettes/ggtree/inst/doc/treeManipulation.html
# #ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
# 
# 
# 
# 
# ##https://stackoverflow.com/questions/44783855/ggtree-plotting-area-not-big-enough
# #png(file="tree.png", width=800, height=800)
# ##svg(file="tree.svg")
# ##pdf(file="tree.pdf", width=600, height=800)
# #ggtree(tree) +
# #  geom_tiplab(node=107, align=TRUE, offset=0.0001, linesize=.6, size=6, colour="red") + 
# #  geom_cladelabel(node=85, fontsize=6, barsize=0.8, label="Cases", angle=0,  
# #                  color="red2", offset=0.00768, align=FALSE) +
# #  geom_cladelabel(node=88, fontsize=6, barsize=0.8, label="Reference cohort", angle=0,  
# #                  color="blue", offset=0.01, align=FALSE) + theme_tree2() + xlim(0.0, 0.05)
# #dev.off()
# 
# 
# 
# #https://yulab-smu.top/treedata-book/chapter13.html
# # Define nodes for coloring later on
# tiplab <- tree$tip.label
# BXs <- tiplab[grep("^[B|X]", tiplab)]                #98
# #Xs <- tiplab[grep("^X", tiplab)]                  #20
# #others <- tiplab[grep("^[A|C|D|E|F|G|H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|Y|Z]", tiplab)]  #84
# labeltree <- groupOTU(tree, BXs)
# #svg(file="tree.svg")
# png(file="tree.png", width=1000, height=1000)
# ggtree(labeltree, aes(color=group, linetype=group), layout="circular") + xlim(-0.5, NA) +
#     scale_color_manual(values = c("#efad29", "#63bbd4")) +
#     geom_nodepoint(color="black", size=0.1) +
#     geom_tiplab(align=FALSE, offset=.1, linesize=.4, size=5, color="black") +
#     theme(
#           legend.position = "none"
#       )
# dev.off()
# 
# 
# #plot.margin = grid::unit(c(-1, -1, -1, -1), "mm")
# 
# 
# #library("ape")
# #data(chiroptera)
# #groupInfo <- split(chiroptera$tip.label, gsub("_\\w+", "", chiroptera$tip.label))
# #chiroptera <- groupOTU(chiroptera, groupInfo)
# #ggtree(chiroptera, aes(color=group), layout='circular') + geom_tiplab(size=1, aes(angle=angle))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # - the position of the case and cohort.
# # - 
# Hello Anna,
# I forgot to describe the alignment step in the figure legend. Attached is the updated version.
# Best,
# Jiabin
# 
# 
# 
# # There were a total of 2288372 positions in the final dataset.
# 
# #A bootstrap of 1000 replicates was performed, and values above 70 are shown. 
# 
# Figure 3. Molecular phylogenetic analysis of RSV subgroup B. The evolutionary history was inferred using the Maximum Likelihood method based on the General Time Reversible model. The tree with the highest log-likelihood is shown. The initial tree for the heuristic search was randomly generated. A discrete Gamma distribution was used to model evolutionary rate differences among sites with four categories. The tree is drawn to scale, with branch lengths measured in the number of substitutions per site. All positions containing gaps and missing data were eliminated.
# 
# #The analysis involved 46 nucleotide sequences and based on a 'Mafft' alignment [3]. Evolutionary analyses were conducted with RaxML-NG [4].
#  
# 1. Nei M. and Kumar S. (2000). Molecular Evolution and Phylogenetics. Oxford University Press, New York.
# 2. Felsenstein J. (1985). Confidence limits on phylogenies: An approach using the bootstrap. Evolution 39:783-791.
# 3. Mafft citation
# DELETE 3. Angiuoli,S.V. and Salzberg,S.L. (2011) Mugsy: fast multiplealignment of closely related whole genomes.Bioinformatics,27,334–342.
# 4. Alexey M. Kozlov, Diego Darriba, Tomáš Flouri, Benoit Morel, and Alexandros Stamatakis (2019) RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, btz305 doi:10.1093/bioinformatics/btz305
# 
# 
# 
# 
# 
# 
# 
# write.nexus(tree)
# 
#   geom_cladelabel(node=17, label="Some random clade", 
#                   color="red2", offset=.8, align=TRUE) + 
#   geom_cladelabel(node=21, label="A different clade", 
#                   color="blue", offset=.8, align=TRUE)
# 
# tree2 <- groupClade(tree, c(17, 21))
# p <- ggtree(tree2, aes(color=group)) + theme(legend.position='none') +
#   scale_color_manual(values=c("black", "firebrick", "steelblue"))
# scaleClade(p, node=17, scale=.1) 
# 
# 
# ggtree(tree)
# ggtree(tree, layout="roundrect")
# ggtree(tree, layout="slanted")
# ggtree(tree, layout="ellipse")
# 
# 
# ggtree(tree, layout="fan", open.angle=120)
# ggtree(tree, layout="equal_angle")
# ggtree(tree, layout="daylight")
# ggtree(tree, branch.length='none')
# ggtree(tree, layout="ellipse", branch.length="none")
# png(file="tree.png", width=600, height=1000)
# #ggtree(tree, layout='circular')
# ggtree(tree, layout="daylight")
# dev.off()
# ggtree(tree, layout="daylight", branch.length = 'none')
# 
# 
# 
# tree2 <- groupClade(tree, c(17, 21))
# p <- ggtree(tree2, aes(color=group)) + theme(legend.position='none') +
#   scale_color_manual(values=c("black", "firebrick", "steelblue"))
# scaleClade(p, node=17, scale=.1) 
# 
# 
# 
#  Trimming readpair 1:	raw_data/2029-AW_S5_L001_R1_001.fastq.gz and raw_data/2029-AW_S5_L001_R2_001.fastq.gz
#   Host reads:		0.15%
#   Fragment size:	372 (sd:0)
#   Subtracting host:	PhiX (PhiX)
#   Alignment rate:	0.14%
#   Subtracting host:	human3 (Homo_sapiens_UCSC_hg38 (dna))
#   Alignment rate:	84.88%
# 
#   
# p938-16972-nra.fasta:CTACAAACTTGCACACTCGNAAAAAAATGGGGCAAATAAGAATTTGATAAGTGCTATTTA
# p938-16972-nra.fasta:AACTTTTCAATAATTTAGCATATTGATTCCAAAATTATCATTTTAGNNNNNNGGGATTAA
# p938-16972-nra.fasta:ATAAAAGTCNAAAAC
# p953-84660-tsek.fasta:ANAAA
# p942-88507-nra.fasta:TAAAAANNNNNNNNNNNNNNNNNNNNNNAAAAATAAGGGTGAAACCAGTAACATAAATTG
# p943-98523-nra.fasta:CTACAAACTTGCACACTCNGAAAAAAATGGGGCAAATAAGAATTTGATAAGTGCTATTTA
# p948-112830-nra.fasta:CTACAAACTTGCACACTCGNAAAAAAATGGGGCAAATAAGAATTTGATAAGTGCTATTTA
# 
# ACGCGAAAAAATGCGTACTACAAACTTGCACACTCGAAAAAAAATGGGGGCAATAAGAATTTGATAAGTG
#             gcgtactacaaacttgcacactcgaaaaaaaatggggcaaataagaatttgataagtgct
# cat
# ctacaaacttgcacactcgaaaaaaaat
# ctacaaacttgcacactcgaaaaaaaatggggc
# 
# ./p944-103323-nra.fasta
# ./p947-105565-nra.fasta
# 
# cat p944-103323-nra.fasta p938-16972-nra.fasta p953-84660-tsek.fasta > c1.fasta
# mafft --adjustdirection --clustalout c1.fasta > c1.aln
# 
# cat p947-105565-nra.fasta p948-112830-nra.fasta p942-88507-nra.fasta p943-98523-nra.fasta > c2.fasta
# mafft --adjustdirection --clustalout c2.fasta > c2.aln
# 
# grep "p944-103323-nra" c1.aln > p944.fasta
# cut -f2-2 -d' ' p944.fasta > p944.fas
# sed -i -e 's/-//g' p944.fas
# grep "p938-16972-nra-" c1.aln > p938.fasta
# cut -f2-2 -d' ' p938.fasta > p938.fas
# sed -i -e 's/-//g' p938.fas
# grep "p953-84660-tsek" c1.aln > p953.fasta
# cut -f2-2 -d' ' p953.fasta > p953.fas
# sed -i -e 's/-//g' p953.fas
# 
# grep "p947-105565-nra" c2.aln > p947.fasta
# cut -f2-2 -d' ' p947.fasta > p947.fas
# sed -i -e 's/-//g' p947.fas
# grep "p948-112830-nra" c2.aln > p948.fasta
# cut -f2-2 -d' ' p948.fasta > p948.fas
# sed -i -e 's/-//g' p948.fas
# grep "p942-88507-nra-" c2.aln > p942.fasta
# cut -f2-2 -d' ' p942.fasta > p942.fas
# sed -i -e 's/-//g' p942.fas
# grep "p943-98523-nra-" c2.aln > p943.fasta
# cut -f2-2 -d' ' p943.fasta > p943.fas
# sed -i -e 's/-//g' p943.fas
# 
# seqkit seq -w 60 -u p944.fas -o p944.fa
# 
# 
# # Annotating the fasta using VAPiD
# makeblastdb -in *.fasta -dbtype nucl
# python ~/Tools/VAPiD/vapid3.py --db ~/REFs/all_virus/all_virus.fasta p938_ref.fa ~/REFs/template.sbt

library(ggtree)
library(ggplot2)

setwd("~/DATA/plotTreeHeatmap/")

info <- read.csv("typing_186.csv")
info$name <- info$Isolate
tree <- read.tree("raxml_567_186.tree")
cols <- c(infection='purple2', commensal='skyblue2')     

library(dplyr)
heatmapData2 <- info %>% select(Isolate, ST, Clonal.Complex, Phylotype)
rn <- heatmapData2$Isolate
heatmapData2$Isolate <- NULL
heatmapData2 <- as.data.frame(sapply(heatmapData2, as.character))
rownames(heatmapData2) <- rn

#https://bookdown.org/hneth/ds4psy/D-3-apx-colors-basics.html
heatmap.colours <- c("cornflowerblue","darkgreen","seagreen3","tan","red",  "navyblue", "gold",     "green","orange","pink","purple","magenta","brown", "darksalmon","chocolate4","darkkhaki", "lightcyan3", "maroon",     "blue","cyan", "skyblue2", "azure3","blueviolet","darkgoldenrod",  "tomato","mediumpurple4", 
                      "cornflowerblue","darkgreen","seagreen3","tan","red","green","orange","pink","brown","magenta",     "cornflowerblue","darkgreen","red","tan","brown",      "darkgrey")
names(heatmap.colours) <- c("1","2","3","4","5", "6","7",   "20","21","22", "28","30","33","42","52","53", "66","68",    "100","105","124","133","134","135",     "159","161",    "CC1","CC2","CC3","CC4","CC5","CC6","CC30","CC72","CC77","Singleton",    "IA1","IA2","IB","II","III",    "NA")

#circular
p <- ggtree(tree, layout='circular', branch.length='none') %<+% info + 
  geom_tippoint(aes(color=Type)) + 
  scale_color_manual(values=cols) + geom_tiplab2(aes(label=name), offset=1)
#, geom='text', align=TRUE,  linetype=NA, hjust=1.8,check.overlap=TRUE, size=3.3
#difference between geom_tiplab and geom_tiplab2?
#+ theme(axis.text.x = element_text(angle = 30, vjust = 0.5)) + theme(axis.text = element_text(size = 20))  + scale_size(range = c(1, 20))
#font.size=10, 
png("ggtree.png", width=1260, height=1260)
p
dev.off()

png("ggtree_and_gheatmap.png", width=1290, height=1000)
gheatmap(p, heatmapData2, width=0.1,colnames_position="top", colnames_angle=90, colnames_offset_y = 0.1, hjust=0.5, font.size=4, offset = 8) + scale_fill_manual(values=heatmap.colours) +  theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size = 14)) + guides(fill=guide_legend(title=""), color = guide_legend(override.aes = list(size = 5)))  
dev.off()

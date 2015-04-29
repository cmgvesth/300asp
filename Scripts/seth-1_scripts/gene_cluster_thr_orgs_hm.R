require("ggplot2")
require("reshape")
require("gplots")
require("RColorBrewer")
require("argparse")
library(ape)



args <- commandArgs(TRUE)
gc_id <- args[1]
cutoff <- args[2]
ord_gc <- args[c(-1,-2)]


dat <- read.csv(paste0("gc_m_distr_",gc_id,".csv"))
dat['X'] <- NULL
new <- cast(dat, h_clust_id ~ q_id,value = 'pident')

print(c("h_clust_id",ord_gc))



synteny <- new[c("h_clust_id",ord_gc)]


mat <- data.matrix(synteny)
rownames(mat) <- synteny[,1]
mat <- mat[,-1]
#mat[mat%%1==0] <- 1 # only use this if you have e.g. sm_protein_ids (using h_id) and you want a 1,0 heatmap
mat[is.na(mat)] <- 0


#######################
# Phylogenic ordering #
#######################

taxo <- read.tree(file="Custom.K8.FullTaxonomy.nwk")
phylo_labels <- gsub("_"," ",taxo$tip.label)
phylo_labels <- gsub("A1163","Af293",phylo_labels)
phylo_labels <- gsub("Aspka1","kawachii",phylo_labels)

phylo_order <- phylo_labels[phylo_labels %in% dat$h_clust_id]

ord.mat <- mat[phylo_order,]
# Done :)

colfunc <- brewer.pal(9,"Blues")

if(length(ord_gc)>10){
        fontsize <- 0.4
}else{
        fontsize <- 0.8
}

pdf(paste0('/home/seth/300asp/Scripts/genecluster_plots/gc_m_distr_',gc_id,'.pdf'))
par(cex.main=0.8)
heatmap.2(ord.mat,
        main = paste0("Heatmap for ", " ",gc_id,"\nwith gc above ", cutoff," members"),
        #dendrogram = 'row',
        margins=c(6,12),
        na.color = 'white',
        col = colfunc,
        key = TRUE,
        key.xlab='BLAST % identity',
        cexRow = fontsize, # Adjusts row fontsize
        trace='none',
        Colv = FALSE,
        Rowv = taxo
        )

#pheatmap(mat,
#        height = 5,
#        legend = FALSE,
#        treeheight_row = 100,
        #cellwidth = 10, scale = 'none',
        #cellheight = 10,
#        color=colfunc,
#        cluster_rows= TRUE,
#        cluster_cols = FALSE,
#        fontsize_row = 7,
#        fontsize_col = 10,
#        )

dev.off()

rm(list = ls())


# NOTES!
# phylogen <- as.factor(taxo$tip.label)


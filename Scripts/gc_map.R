require("ggplot2")
require("reshape")
require("gplots")
require("RColorBrewer")
require("argparse")
library("pheatmap")

args <- commandArgs(TRUE)

print(args)

dat <- read.csv("thmap.csv")
dat['X'] <- NULL
new <- cast(dat, h_clust_id ~ q_id,value = 'pident')

print(c("h_clust_id",args))

synteny <- new[c("h_clust_id",args)]


mat <- data.matrix(synteny)
rownames(mat) <- synteny[,1]
mat <- mat[,-1]
#mat[mat%%1==0] <- 1 # only use this if you have e.g. sm_protein_ids (using h_id) and you want a 1,0 heatmap
mat[is.na(mat)] <- 0

colfunc <- brewer.pal(9,"Blues")


pdf('trying.pdf', title = "Heatmap for 11_382281_6")
heatmap.2(mat,
        main = "HM for 11_382281_6",
        dendrogram = 'row',
        margins=c(6,10)
        na.color = 'white',
        col = colfunc,
        key = TRUE,
        key.xlab='BLAST % identity',
        cexRow = 0.5, # Adjusts row fontsize
        trace='none',
        Colv = FALSE
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
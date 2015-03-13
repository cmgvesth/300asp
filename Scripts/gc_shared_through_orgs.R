require("ggplot2")
require("reshape")
require("gplots")
require("RColorBrewer")
require("argparse")

args <- commandArgs(TRUE)
infile <- args[1]
print(infile)
description <- args[2:4] # Contains section, cutoff1 and cutoff 2
print(description)
tmpdat=read.csv(infile, header = FALSE)
#subdat <- subset(tmpdat, tmpdat$V8=="Nigri" & tmpdat$V9=="Nigri", select=c(3,4,10))
subdat <- tmpdat #[order(tmpdat$V1, tmpdat$V2), ]

#aggregate(subdat$V3, list(org=subdat$V1), max) # Find max values for each cluster to normalize

#x <- barplot(as.matrix(tmpdat$V10), xaxt='n' , ylab= "Number of best hits inside one gene cluster", main='Distribution of gc-members throughout different genomes', sub='Gene Cluster (clust_id: 22_1079950_32)')
#text(cex=.6, x=x-1, y=-2, names(tmpdat) , xpd=TRUE, srt=45) # for rotated x-labels, put xaxt='n' into plot when in use!


#qplot(tmpdat$V10, ..density.., data=tmpdat, geom="density", fill=tmpdat$V3, position="stack")
#pdf(paste0('/home/seth/Dropbox/seth-1/gc_dist_',description[1],'.pdf'))
pdf(paste0('/home/seth/Dropbox/seth-1/heatmap/density_orgs_', description[2], '_',description[3],'.pdf'))
m<- ggplot(subdat, aes(x=subdat$V3))
p<- m + geom_density(aes(fill=factor(subdat$V1)), size=2, alpha=0.3) #+ coord_cartesian(xlim=c(0, 20))
print(p)
dev.off()

dat <- cast(tmpdat, V1~V2)

mat <- data.matrix(dat)
rownames(mat) <- dat[,1]
mat <- mat[,-1]

mat[is.na(mat)] <- 0


#my_palette <- colorRampPalette(c("yellow", "red"))(n = 1000) # for costum colors
b <- c(0,0.2,0.4,0.6,0.8,1)*100

colfunc <- c("white"  ,"#74A9CF", "#3690C0" ,"#0570B0","#045A8D")
#colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
par(cex.main=0.5)
pdf(paste0('/home/seth/Dropbox/seth-1/heatmap/map_orgs_',description[2], '_',description[3],'.pdf'))
heatmap.2(mat,
        main = paste0("Shared GC among"," ",description[1], "\n", "with more than ", description[2],'% of members \n and at least ', description[3], " members"),
        dendrogram="none",
        trace = "none",
        na.color ="white",
        col=colfunc,
        Rowv=NA,
        Colv=NA,
        #srtRow = 45,
        #srtCol = 45,
        cexRow = 0.8,
        cexCol = 0.8,
        margins = c(8,8)
        )
dev.off()

pdf(paste0('/home/seth/Dropbox/seth-1/heatmap/map_orgs_clustered_', description[2], '_',description[3],'.pdf'))
heatmap.2(mat,
        main = paste0("Shared GC among", " ",description[1], "\n", "with more than ", description[2],"% of members \n and at least ", description[3], " members"),
        dendrogram="both",
        trace = "none",
        na.color ="white",
        col=colfunc,
        #srtRow = 45,
        #srtCol = 45,
        cexRow = 0.8,
        cexCol = 0.8,
        margins = c(8,8)
        )
dev.off()
require("ggplot2")
require("reshape")
require("gplots")
require("RColorBrewer")
tmpdat=read.csv("map_orgs.csv", header = FALSE)
subdat <- subset(tmpdat, tmpdat$V8=="Nigri" & tmpdat$V9=="Nigri", select=c(3,4,10))
subdat <- subdat[order(subdat$V3, subdat$V4), ]

#x <- barplot(as.matrix(tmpdat$V10), xaxt='n' , ylab= "Number of best hits inside one gene cluster", main='Distribution of gc-members throughout different genomes', sub='Gene Cluster (clust_id: 22_1079950_32)')
#text(cex=.6, x=x-1, y=-2, names(tmpdat) , xpd=TRUE, srt=45) # for rotated x-labels, put xaxt='n' into plot when in use!


#qplot(tmpdat$V10, ..density.., data=tmpdat, geom="density", fill=tmpdat$V3, position="stack")

pdf('density_orgs.pdf')
m<- ggplot(subdat, aes(x=subdat$V10))
p<- m + geom_density(aes(fill=factor(subdat$V3)), size=2, alpha=0.3) #+ coord_cartesian(xlim=c(0, 20))
print(p)
dev.off()

dat <- cast(tmpdat, V3~V4)

mat <- data.matrix(dat)
rownames(mat) <- dat[,1]
mat <- mat[,-1]

#my_palette <- colorRampPalette(c("yellow", "red"))(n = 1000) # for costum colors

pdf('map_orgs.pdf')
heatmap.2(mat,
        dendrogram="none",
        trace = "none",
        na.color ="white",
        col=brewer.pal(7, "YlOrRd"),
        Rowv=NA,
        Colv=NA)
dev.off()
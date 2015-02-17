data=read.csv("gc_asframe.csv")

data['X'] <- NULL

pdf('gc_dist.pdf')
x <- barplot(as.matrix(data), xaxt='n' , ylab= "Number of best hits inside one gene cluster", main='Distribution of gc-members throughout different genomes', sub='Gene Cluster (clust_id: 22_1079950_32)')
text(cex=.6, x=x-1, y=-2, names(data) , xpd=TRUE, srt=45) # for rotated x-labels, put xaxt='n' into plot when in use!
dev.off()
data=read.csv("new")

data['X'] <- NULL

x <- barplot(as.matrix(data), xact='n', main='Distribution of best hits of Gene Cluster members throughout different genomes')
text(cex=.6, x=x-.25, y=-1.25, names(data) , xpd=TRUE, srt=45)
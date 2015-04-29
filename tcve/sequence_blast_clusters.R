#########################################################
### A) Installing and loading required packages
#########################################################
if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }

if (!require("reshape2")) {
   install.packages("reshape2", dependencies = TRUE)
   library(reshape)
   }
if (!require("pheatmap")) {
   install.packages("pheatmap", dependencies = TRUE)
   library(pheatmap)
   }
if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2)
   }
if (!require("plyr")) {
   install.packages("plyr", dependencies = TRUE)
   library(plyr)
   }   
if (!require("cluster")) {
   install.packages("cluster", dependencies = TRUE)
   library(cluster)
   }   

# Read data from file
infile = replace_infile
tmpdat <- read.csv(infile, sep=";", header=TRUE)

afolder = replace_folderpath 
fileroot = paste(afolder, analysisNr, sep="/") 


data1 <- tmpdat
###########################################
# Analysis 1
###########################################
dat1 <- dcast(data1, q_seqkey+h_seqkey~q_org, value.var="pident", fill=0, max) # WOOOOOOOOOO 
mat <- t(data.matrix(dat1[,3:ncol(dat1)])) # NOTE transpose
colnames(mat) <- dat1$q_clustid

cvs <- unique(dat1$q_org)
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat1$q_org, function(x) cvs_color[x==cvs])

plotwidth <- 12
if (length(cvs) <= 10) {
	plotwidth <- length(cvs)+12
}
plotheight <- 10

b <- c(0,0.2,0.4,0.6,0.8,1)

#pdf(paste(fileroot, "_orgs2cluster-all-cutoff-clustered-tree.pdf", sep=""), width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, mostly white, clustered", dendrogram="both", trace = "none", 
   Rowv=TRUE, Colv=TRUE ,col=colfunc, density.info="none", margins =c(plotwidth,plotheight))

legend('bottomleft', legend=c(unique(dat1$q_org)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()


dat2 <- dcast(data1, q_seqkey+h_seqkey~q_org, value.var="pident", fill=0, max) # WOOOOOOOOOO 

# DO I need a square matrix to do dostances?
# How do I make it a square??? 
# Dont remove reciprocal hits - nope!

datdist <- as.dist(as.matrix(dat2))

wss <- (nrow(dat2)-1)*sum(apply(dat2,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(dat2, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

fit <- hclust(datdist, method="ward") 
# Method
# the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of 
# "ward", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

plot(fit) # display dendogram
#groups <- cutree(fit, k=4) # cut tree into 5 clusters
rect.hclust(fit, k=4, border="red")

clusplot(as.matrix(datdist), fit$labels, color=TRUE, shade=TRUE, labels=0, lines=0)

fit <- kmeans(datdist, 4)
clusplot(as.matrix(datdist), fit$cluster, color=TRUE, shade=TRUE, labels=4, lines=0)

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

tmpdat <- read.csv("hitsPerOrgPerClusterGene.csv", sep=";", header=TRUE)

# Number of hits, absolute values
# subdat <- subset(tmpdat, tmpdat$name=="Aspacu1" & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(2,4,7))
# subdat <- subset(tmpdat, tmpdat$name=="Aspacu1" & tmpdat$clust_size>5, select=c(2,4,7))
# subdat <- subset(tmpdat, tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(2,4,7))
# subdat <- subset(tmpdat, tmpdat$clust_size>8 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(2,4,7))
# subdat <- subset(tmpdat, tmpdat$clust_size>8 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(1,2,4,7))
# subdat <- subset(tmpdat, tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(2,4,5))
subdat <- subset(tmpdat, tmpdat$clust_size>10 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri" & tmpdat$org_id=="10", select=c(1,2,4,7))
subdat <- subset(tmpdat, tmpdat$clust_size>10 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(1,2,4,7))
subdat <- subset(tmpdat, tmpdat$clust_size>10 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri" & tmpdat$org_id=="10" |  tmpdat$org_id=="2", select=c(1,2,4,7))
# mdat <- melt(subdat, id=c("clust_id","h_org","org_id"))
# dat <- cast(mdat, h_org+clust_id~clustCov, value="clustCov", fill=0)
subdat <- subset(tmpdat, tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(1,2,4,7))

dat <- dcast(subdat, org_id+clust_id~h_org, value.var="clustCov", fill=0) # WOOOOOOOOOO 
mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
colnames(mat) <- dat$clust_id

cvs <- unique(dat$org_id)
#cvs_color <- sample(rainbow(length(cvs)), length(cvs_color))
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$org_id, function(x) cvs_color[x==cvs])

b <- c(0,0.2,0.4,0.6,0.8,1)

pdf('orgs2cluster-all-cutoff-clustered.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-all-range-clustered.pdf')
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-all-cutoff-unclustered.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-all-cutoff-orgclustered.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, organism clustered", dendrogram="row", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()










subdat <- subset(tmpdat, tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri" & tmpdat$org_id=="22", select=c(1,2,4,7))

dat <- dcast(subdat, org_id+clust_id~h_org, value.var="clustCov", fill=0) # WOOOOOOOOOO 
mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
colnames(mat) <- dat$clust_id

cvs <- unique(dat$org_id)
#cvs_color <- sample(rainbow(length(cvs)), length(cvs_color))
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$org_id, function(x) cvs_color[x==cvs])

b <- c(0,0.2,0.4,0.6,0.8,1)

pdf('orgs2cluster-Aspni7-cutoff-clustered.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Aspni7 nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-Aspni7-range-clustered.pdf')
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Aspni7 nigri clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-Aspni7-cutoff-unclustered.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Aspni7 nigri clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-Aspni7-cutoff-orgclustered.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Aspni7 nigri clusters, organism clustered", dendrogram="row", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()












subdat <- subset(tmpdat, tmpdat$clust_size>10 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(1,2,4,7))

dat <- dcast(subdat, org_id+clust_id~h_org, value.var="clustCov", fill=0) # WOOOOOOOOOO 
mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
colnames(mat) <- dat$clust_id

cvs <- unique(dat$org_id)
#cvs_color <- sample(rainbow(length(cvs)), length(cvs_color))
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$org_id, function(x) cvs_color[x==cvs])

b <- c(0,0.2,0.4,0.6,0.8,1)

pdf('orgs2cluster-10-cutoff-clustered.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-10-range-clustered.pdf')
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 nigri clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-10-cutoff-unclustered.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 nigri clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-10-cutoff-orgclustered.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 nigri clusters, organism clustered", dendrogram="row", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-color-pie.pdf')
pie(rep(1,length(cvs_color)), col=cvs_color)
dev.off()

"""
dat <- cast(subdat, h_org~clust_id+org_id, value="clustCov", fill=0)
mat <- data.matrix(dat)
rownames(mat) <- dat[,1]
mat <- mat[,-1]

heatmap.2(mat, dendrogram="row", trace = "none", col=colfunc)
heatmap.2(mat, dendrogram="both", trace = "none", col=colfunc)

# Cluster vs cluster
tmpdat <- read.csv("query_main.csv", sep=";", header=TRUE)
subdat <- subset(tmpdat, select=c(1,10,19))
dat <- cast(subdat, ta_clust_id~tb_clust_id, value="over", fill=0)
mat <- data.matrix(dat)
rownames(mat) <- dat[,1]
mat <- mat[,-1:-2]
heatmap.2(mat, dendrogram="both", trace = "none", col=colfunc)
"""
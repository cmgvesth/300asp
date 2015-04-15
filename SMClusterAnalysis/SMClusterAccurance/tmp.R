#########################################################
### A) Installing and loading required packages
#########################################################
library(ggplot2)
args <- commandArgs(TRUE)
infile <- args[1]
print (infile)

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
 library(plyr)  

#if (!require("RSVGDevice")) {
#   install.packages("RSVGDevice", dependencies = TRUE)
#   library(pheatmap)
#   }

# tmpdat$q_orgname=='Aspnov1'  |  tmpdat$q_orgname=='Aspfu1'  <- tmpdat$q_orgname=="Aspnov1" | tmpdat$q_orgname=="Aspfu1"

tmpdat <- read.csv(infile, sep=";", header=TRUE)

subdat <- subset(tmpdat, ( tmpdat$q_orgname=='Aspnov1'  |  tmpdat$q_orgname=='Aspfu1' ) & tmpdat$q_orgname!="Aspac1" & tmpdat$h_orgname!="Aspac1", select=c(12,3,6,8)) 

dat <- dcast(subdat, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
colnames(mat) <- dat$q_clustid

cvs <- unique(dat$q_orgid)
plotwidth <- length(cvs)+12
plotheight <- 10
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])

b <- c(0,0.2,0.4,0.6,0.8,1)


pdf('orgs2cluster-all-cutoff-clustered-bestcand-tree.pdf', width=plotwidth, height=plotheight)
#colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, 
	Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b,  density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
#pheatmap(mat, main="All  clusters, mostly white, clustered", show_colnames=FALSE, cluster_cols=TRUE,cluster_rows=TRUE, treeheight_col=100,color=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
dev.off()

pdf('orgs2cluster-all-cutoff-clustered-bestcand-notree.pdf', width=plotwidth, height=plotheight)
#colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, 
	Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b,  density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
#pheatmap(mat, main="All  clusters, mostly white, no tree, clustered", show_colnames=FALSE, cluster_cols=TRUE,cluster_rows=TRUE, treeheight_col=0,color=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
dev.off()

pdf('orgs2cluster-all-color-clustered-bestcand.pdf', width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, 
	Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
#pheatmap(mat, main="All  clusters, more color, clustered", show_colnames=FALSE, cluster_cols=TRUE,cluster_rows=TRUE, treeheight_col=100,color=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
dev.off()

pdf('orgs2cluster-all-color-clusteredrow-bestcand.pdf', width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, clustered",dendrogram="row", trace = "none", ColSideColor=ccol_stages, 
	Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
#pheatmap(mat, main="All  clusters, more color, clustered", show_colnames=FALSE, cluster_cols=TRUE,cluster_rows=TRUE, treeheight_col=100,color=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
dev.off()


pdf('orgs2cluster-all-cutoff-unclustered-bestcand.pdf', width=plotwidth, height=plotheight)
#colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, 
	Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
#pheatmap(mat, main="All  clusters, mostly white, unclustered", show_colnames=FALSE, cluster_cols=FALSE,cluster_rows=FALSE, color=colfunc, breaks=b, margins =c(plotwidth,plotheight))
dev.off()

# allneo <- subset(tmpdat, tmpdat$q_sec=="" & tmpdat$h_sec=="" & tmpdat$clustCov>=0.6 & tmpdat$q_orgname=="Aspneo1" & tmpdat$h_orgname!="Aspcos1") 
# neocos <- subset(tmpdat, tmpdat$q_sec=="" & tmpdat$h_sec=="" & tmpdat$clustCov>=0.6 & tmpdat$q_orgname=="Aspneo1" & tmpdat$h_orgname=="Aspcos1") 
# neocos <- subset(tmpdat, tmpdat$q_sec=="" & tmpdat$h_sec=="" & tmpdat$clustCov>=0.6 & tmpdat$q_orgname=="Aspneo1", select=c(12,3,6,8)) 

orgnames <- rownames(mat) 
fish <- rep(0, length(orgnames))
fish[match("neoniger",orgnames)] = 1
fish[match("costaricaensis",orgnames)] = 1


fish <- rep(0, length(orgnames))
fish[match("fijiensis",orgnames)] = 1
fish[match("brunneoviolaceus",orgnames)] = 1
fish[match("aculeatinus",orgnames)] = 1

for ( clust in colnames(mat)) {
	vec <- mat[,clust]
	logical <- as.numeric(vec > 0.6)
	if (paste(fish, sep="", collapse="") == paste(logical, sep="", collapse="")) {
		print(clust)
	}
}
# colsep=c(0:length(colnames(mat))), rowsep=c(0:length(rownames(mat))),

cvs <- unique(dat$q_orgid)
cvs_color <- rep(c("white", "grey50", "black"),length(cvs)/3)
ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])

pdf('orgs2cluster-all-cutoff-orgclustered-bestcand-blackwhite.pdf', width=plotwidth, height=plotheight )
#colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, organism clustered", dendrogram="row", trace = "none", 
	sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),ColSideColor=ccol_stages,
	Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight),
	xlab="Predicted clusters from all organisms", 
	distfun = function(x) dist(x,method = 'euclidean'),
	hclustfun = function(x) hclust(x,method = 'complete'))
	legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()

pdf('orgs2cluster-all-cutoff-bothclustered-bestcand-blackwhite.pdf', width=plotwidth, height=plotheight)
#colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, organism clustered", dendrogram="both", trace = "none", 
	sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),ColSideColor=ccol_stages,
	Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight),
	xlab="Predicted clusters from all organisms", 
	distfun = function(x) dist(x,method = 'euclidean'),
	hclustfun = function(x) hclust(x,method = 'complete'))
	legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()
# http://www.inside-r.org/r-doc/stats/hclust
# http://www.inside-r.org/r-doc/stats/dist
#pheatmap(mat, main="All  clusters, organism clustered", show_colnames=FALSE, cluster_cols=FALSE,cluster_rows=TRUE, treeheight_row=100,color=colfunc, breaks=b, margins =c(plotwidth,plotheight))

subdat <- subset(tmpdat, tmpdat$q_clust_size>10 & ( tmpdat$q_orgname=='Aspnov1'  |  tmpdat$q_orgname=='Aspfu1' ) & tmpdat$q_orgname!="Aspac1" & tmpdat$h_orgname!="Aspac1", select=c(12,3,6,8)) # q_orgid,

dat <- dcast(subdat, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
colnames(mat) <- dat$q_clustid
cvs <- unique(dat$q_orgid)
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])
b <- c(0,0.2,0.4,0.6,0.8,1)

pdf('orgs2cluster-10-cutoff-clustered-bestcand.pdf', width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10  clusters, mostly white, clustered", dendrogram="both", 
sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),
trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
#pheatmap(mat, main="Size >10  clusters, mostly white, clustered", annotation_colors=ccol_stages, show_colnames=FALSE, cluster_cols=TRUE,cluster_rows=TRUE, treeheight_row=100,color=colfunc, breaks=b, margins =c(plotwidth,plotheight))
dev.off()

pdf('orgs2cluster-10-range-clustered-bestcand.pdf', width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10  clusters, clustered",dendrogram="both", 
sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),
trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()

pdf('orgs2cluster-10-cutoff-unclustered-bestcand.pdf', width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10  clusters, unclustered", dendrogram="none", 
sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),
trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()

pdf('orgs2cluster-10-cutoff-orgclustered-bestcand.pdf', width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10  clusters, organism clustered", dendrogram="row", 
sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),
trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()



orgdat <- subset(tmpdat, ( tmpdat$q_orgname=='Aspnov1'  |  tmpdat$q_orgname=='Aspfu1' ) & tmpdat$q_orgname!="Aspac1" & tmpdat$h_orgname!="Aspac1") 

for (s in unique(orgdat$q_orgname)) {
	filename = paste("orgs2cluster-clustered-bestcand_", s, ".pdf", sep="")
	title = paste("Sec. met. clusters, mostly white, clustered, clusters from ", s,  sep="")

	subdat <- subset(orgdat, orgdat$q_orgname==s, select=c(12,3,6,8)) # q_orgid, q_clustid, h_realname, clustCov

	#subdat <- subset(tmpdat, tmpdat$q_realname=='fumigatus Af293' & tmpdat$q_sec=="" & tmpdat$h_sec=="", select=c(12,3,6,8)) # q_orgid, q_clustid, h_realname, clustCov

	dat <- dcast(subdat, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
	mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
	colnames(mat) <- dat$q_clustid
	cvs <- unique(dat$q_orgid)
	cvs_color <- rainbow(length(cvs))
	ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])
	b <- c(0,0.2,0.4,0.6,0.8,1)
	#heatmap.2(mat, main=title, dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))

	pdf(filename, width=plotwidth, height=plotheight)
	colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
	heatmap.2(mat, main=title, dendrogram="both", 
	sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),
	trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, density.info="none", margins =c(plotwidth,plotheight))
	legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
	dev.off()
}


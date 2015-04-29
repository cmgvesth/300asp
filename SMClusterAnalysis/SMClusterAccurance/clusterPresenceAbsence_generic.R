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

# Read data from file
infile = replace_infile
tmpdat <- read.csv(infile, sep=";", header=TRUE)

# rename columns
colnames(tmpdat)[colnames(tmpdat) == 'org_id'] <- 'q_orgid'
colnames(tmpdat)[colnames(tmpdat) == 'h_org'] <- 'h_realname'
colnames(tmpdat)[colnames(tmpdat) == 'name'] <- 'q_orgname'
colnames(tmpdat)[colnames(tmpdat) == 'clust_id'] <- 'q_clustid'
colnames(tmpdat)[colnames(tmpdat) == 'q_clust_size'] <- 'clust_size'

# Subset data for analysis 1, 2 and 3 
replace_data1
replace_data2
replace_data3

afolder = replace_folderpath 
analysisNr = replace_analysisNumber
fileroot = paste(afolder, analysisNr, sep="/") 


###########################################
# Analysis 1
###########################################
dat <- dcast(data1, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
colnames(mat) <- dat$q_clustid

cvs <- unique(dat$q_orgid)
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])

plotwidth <- 12
if (length(cvs) <= 10) {
	plotwidth <- length(cvs)+12
}
plotheight <- 10

b <- c(0,0.2,0.4,0.6,0.8,1)

# afolder
pdf(paste(fileroot, "_orgs2cluster-all-cutoff-clustered-tree.pdf", sep=""), width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, mostly white, clustered", dendrogram="row", trace = "none", 
	ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()

pdf(paste(fileroot, "_orgs2cluster-all-color-clustered.pdf", sep=""), width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, clustered",dendrogram="row", trace = "none", 
	ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()

pdf(paste(fileroot, "_orgs2cluster-all-cutoff-unclustered.pdf", sep=""), width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All clusters, unclustered", dendrogram="none", trace = "none", 
	ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()

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

# cvs <- unique(dat$q_orgid)
# cvs_color <- rep(c("white", "grey50", "black"),length(cvs)/3)
# ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])

# pdf(paste(fileroot, "_orgs2cluster-all-cutoff-orgclustered-blackwhite.pdf", sep=""), width=plotwidth, height=plotheight)
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="All clusters, organism clustered", dendrogram="row", trace = "none", 
# 	sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),ColSideColor=ccol_stages,
# 	Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight),
# 	xlab="Predicted clusters from all organisms", 
# 	distfun = function(x) dist(x,method = 'euclidean'),
# 	hclustfun = function(x) hclust(x,method = 'complete'))
# legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
# dev.off()

# pdf(paste(fileroot, "_orgs2cluster-all-cutoff-bothclustered-blackwhite.pdf", sep=""), width=plotwidth, height=plotheight)
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="All clusters, organism clustered", dendrogram="row", trace = "none", 
# 	sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),ColSideColor=ccol_stages,
# 	Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight),
# 	xlab="Predicted clusters from all organisms", 
# 	distfun = function(x) dist(x,method = 'euclidean'),
# 	hclustfun = function(x) hclust(x,method = 'complete'))
# legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
# dev.off()


###########################################
# Analysis 2
###########################################
dat <- dcast(data2, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
colnames(mat) <- dat$q_clustid
cvs <- unique(dat$q_orgid)
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])
b <- c(0,0.2,0.4,0.6,0.8,1)

pdf(paste(fileroot, "_orgs2cluster-10-cutoff-clustered.pdf", sep=""), width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 clusters, mostly white, clustered", dendrogram="both", trace = "none", 
	ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()

pdf(paste(fileroot, "_orgs2cluster-10-range-clustered.pdf", sep=""), width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 clusters, clustered",dendrogram="both", trace = "none", 
	ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()

pdf(paste(fileroot, "_orgs2cluster-10-cutoff-unclustered.pdf", sep=""), width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 clusters, unclustered", dendrogram="none", trace = "none", 
	ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()

pdf(paste(fileroot, "_orgs2cluster-10-cutoff-orgclustered.pdf", sep=""), width=plotwidth, height=plotheight)
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 clusters, organism clustered", dendrogram="row", trace = "none", 
	ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight))
legend('bottomleft', legend=c(unique(dat$q_orgid)),lty=1, col=c(unique(ccol_stages)), cex=0.8)
dev.off()


###########################################
# Analysis 3
###########################################
for (s in unique(data3$q_orgname)) {
	title = paste("Sec. met. clusters, mostly white, clustered, clusters from ", s,  sep="")
	data4 <- subset(data3, data3$q_orgname==s, select=c(q_orgid,q_clustid,h_realname,clustCov))

	dat <- dcast(data4, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
	mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
	colnames(mat) <- dat$q_clustid
	cvs <- unique(dat$q_orgid)
	cvs_color <- rainbow(length(cvs))
	ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])
	b <- c(0,0.2,0.4,0.6,0.8,1)

	pdf(paste(fileroot,"_orgs2cluster-clustered_", s, ".pdf", sep=""), width=plotwidth, height=plotheight)
	colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
	heatmap.2(mat, main=title, dendrogram="both", trace = "none", 
		ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(plotwidth,plotheight))
	dev.off()
}

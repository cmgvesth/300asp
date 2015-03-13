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
#tmpdat <- read.csv("/home/tcve/github/tcve/t_antismashLoopAntismash.csv", sep=";", header=TRUE)
tmpdat <- read.csv(infile, sep=";", header=TRUE)

subdat <- subset(tmpdat, tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri" & tmpdat$q_orgname!="Aspac1" & tmpdat$h_orgname!="Aspac1", select=c(12,3,6,8)) # q_orgid, q_clustid, h_realname, clustCov
dat <- dcast(subdat, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
colnames(mat) <- dat$q_clustid

cvs <- unique(dat$q_orgid)
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])

b <- c(0,0.2,0.4,0.6,0.8,1)

pdf('orgs2cluster-all-cutoff-clustered-bestcand-tree.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
#pheatmap(mat, main="All nigri clusters, mostly white, clustered", show_colnames=FALSE, cluster_cols=TRUE,cluster_rows=TRUE, treeheight_col=100,color=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-all-cutoff-clustered-bestcand-notree.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
#pheatmap(mat, main="All nigri clusters, mostly white, no tree, clustered", show_colnames=FALSE, cluster_cols=TRUE,cluster_rows=TRUE, treeheight_col=0,color=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-all-color-clustered-bestcand.pdf')
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
#pheatmap(mat, main="All nigri clusters, more color, clustered", show_colnames=FALSE, cluster_cols=TRUE,cluster_rows=TRUE, treeheight_col=100,color=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-all-cutoff-unclustered-bestcand.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
#pheatmap(mat, main="All nigri clusters, mostly white, unclustered", show_colnames=FALSE, cluster_cols=FALSE,cluster_rows=FALSE, color=colfunc, breaks=b, margins =c(4,10))
dev.off()

# allneo <- subset(tmpdat, tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri" & tmpdat$clustCov>=0.6 & tmpdat$q_orgname=="Aspneo1" & tmpdat$h_orgname!="Aspcos1") 
# neocos <- subset(tmpdat, tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri" & tmpdat$clustCov>=0.6 & tmpdat$q_orgname=="Aspneo1" & tmpdat$h_orgname=="Aspcos1") 
# neocos <- subset(tmpdat, tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri" & tmpdat$clustCov>=0.6 & tmpdat$q_orgname=="Aspneo1", select=c(12,3,6,8)) 

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

pdf('orgs2cluster-all-cutoff-orgclustered-bestcand-blackwhite.pdf', width=12, height=8)
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, organism clustered", dendrogram="row", trace = "none", 
	sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),ColSideColor=ccol_stages,
	Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10),
	xlab="Predicted clusters from all organisms", 
	distfun = function(x) dist(x,method = 'euclidean'),
	hclustfun = function(x) hclust(x,method = 'complete'))
dev.off()

pdf('orgs2cluster-all-cutoff-bothclustered-bestcand-blackwhite.pdf', width=12, height=8)
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, organism clustered", dendrogram="both", trace = "none", 
	sepcolor="grey50",  sepwidth=c(0.001, 0.001), rowsep=c(0:length(rownames(mat))),ColSideColor=ccol_stages,
	Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10),
	xlab="Predicted clusters from all organisms", 
	distfun = function(x) dist(x,method = 'euclidean'),
	hclustfun = function(x) hclust(x,method = 'complete'))
dev.off()
# http://www.inside-r.org/r-doc/stats/hclust
# http://www.inside-r.org/r-doc/stats/dist
#pheatmap(mat, main="All nigri clusters, organism clustered", show_colnames=FALSE, cluster_cols=FALSE,cluster_rows=TRUE, treeheight_row=100,color=colfunc, breaks=b, margins =c(4,10))


subdat <- subset(tmpdat, tmpdat$q_clust_size>10 & tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri" & tmpdat$q_orgname!="Aspac1" & tmpdat$h_orgname!="Aspac1", select=c(12,3,6,8)) # q_orgid, q_clustid, h_realname, clustCov
dat <- dcast(subdat, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
colnames(mat) <- dat$q_clustid
cvs <- unique(dat$q_orgid)
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])
b <- c(0,0.2,0.4,0.6,0.8,1)

pdf('orgs2cluster-10-cutoff-clustered-bestcand.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
#pheatmap(mat, main="Size >10 nigri clusters, mostly white, clustered", annotation_colors=ccol_stages, show_colnames=FALSE, cluster_cols=TRUE,cluster_rows=TRUE, treeheight_row=100,color=colfunc, breaks=b, margins =c(4,10))
dev.off()

pdf('orgs2cluster-10-range-clustered-bestcand.pdf')
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 nigri clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-10-cutoff-unclustered-bestcand.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 nigri clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()

pdf('orgs2cluster-10-cutoff-orgclustered-bestcand.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="Size >10 nigri clusters, organism clustered", dendrogram="row", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
dev.off()


orgdat <- subset(tmpdat, tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri" & tmpdat$q_orgname!="Aspac1" & tmpdat$h_orgname!="Aspac1") # q_orgid, q_clustid, h_realname, clustCov

for (s in unique(orgdat$q_orgname)) {
	filename = paste("orgs2cluster-clustered-bestcand_", s, ".pdf", sep="")
	title = paste("Sec. met. clusters, mostly white, clustered, clusters from ", s,  sep="")

	subdat <- subset(orgdat, orgdat$q_orgname==s, select=c(12,3,6,8)) # q_orgid, q_clustid, h_realname, clustCov

	#subdat <- subset(tmpdat, tmpdat$q_realname=='fumigatus Af293' & tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri", select=c(12,3,6,8)) # q_orgid, q_clustid, h_realname, clustCov

	dat <- dcast(subdat, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
	mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
	colnames(mat) <- dat$q_clustid
	cvs <- unique(dat$q_orgid)
	cvs_color <- rainbow(length(cvs))
	ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])
	b <- c(0,0.2,0.4,0.6,0.8,1)
	#heatmap.2(mat, main=title, dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))

	pdf(filename)
	colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
	heatmap.2(mat, main=title, dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
	dev.off()
}

	# pdf('orgs2cluster-10-range-clustered-bestcand.pdf')
	# colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
	# heatmap.2(mat, main="Size >10 nigri clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
	# dev.off()

	# pdf('orgs2cluster-10-cutoff-unclustered-bestcand.pdf')
	# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
	# heatmap.2(mat, main="Size >10 nigri clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
	# dev.off()

	# pdf('orgs2cluster-10-cutoff-orgclustered-bestcand.pdf')
	# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
	# heatmap.2(mat, main="Size >10 nigri clusters, organism clustered", dendrogram="row", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
	# dev.off()


	# filename = paste("aalength_histogram_section_", s, ".pdf", sep="")
	# title = paste("Histograms, binwidth 100 amino acids, plot per organism, section ", s,  sep="")
	# a <- ggplot(set, aes(x=len, fill=section)) 
	# a <- a + geom_bar(binwidth=100) + facet_grid(~name) + xlim(0,2000) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	# a <- a + labs(title = title, x = "Organisms, alphabetically", y = "Protein counts") 
	# ggsave(plot = a, file=filename, width=12+length(unique(set$org_id)), height=12, limitsize=FALSE)

	# filename = paste("aalength_boxplots_section_", s, ".pdf", sep="")
	# title = paste("Boxplot, orderd by median length, plot per organism, section ", s,  sep="")
	# b <- ggplot(set, aes(y=len, x=reorder(name, len, FUN=median))) 
	# b <- b + geom_boxplot() + coord_flip() 
	# b <- b + labs(title = title, x = "Organisms, ordered by median", y = "protein lengths") 
	# ggsave(plot = b, file=filename, height=12+length(unique(set$org_id)), width=12, limitsize=FALSE)
}


# Number of hits, absolute values
# subdat <- subset(tmpdat, tmpdat$name=="Aspacu1" & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(2,4,7))
# subdat <- subset(tmpdat, tmpdat$name=="Aspacu1" & tmpdat$clust_size>5, select=c(2,4,7))
# subdat <- subset(tmpdat, tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(2,4,7))
# subdat <- subset(tmpdat, tmpdat$clust_size>8 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(2,4,7))
# subdat <- subset(tmpdat, tmpdat$clust_size>8 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(1,2,4,7))
# subdat <- subset(tmpdat, tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(2,4,5))
#subdat <- subset(tmpdat, tmpdat$q_clust_size>10 & tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri" & tmpdat$q_orgid=="10", select=c(10,3,4,6))
#subdat <- subset(tmpdat, tmpdat$q_clust_size>10 & tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri" , select=c(12,3,4,8)) # q_orgid, q_clustid, h_orgname, clustCov
#subdat2 <- subset(tmpdat2, tmpdat$clust_size>10 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri" & tmpdat$org_id=="10", select=c(1,2,4,7))
# subdat <- subset(tmpdat, tmpdat$clust_size>10 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri" & tmpdat$org_id=="10", select=c(1,2,4,7))
# subdat <- subset(tmpdat, tmpdat$clust_size>10 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(1,2,4,7))
# subdat <- subset(tmpdat, tmpdat$clust_size>10 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri" & tmpdat$org_id=="10" |  tmpdat$org_id=="2", select=c(1,2,4,7))
# # mdat <- melt(subdat, id=c("clust_id","h_org","org_id"))
# # dat <- cast(mdat, h_org+clust_id~clustCov, value="clustCov", fill=0)
# subdat <- subset(tmpdat, tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(1,2,4,7))


# If the max value of a organism/cluster pair is found for several hit organisms, 
# there will be several rows with the same organism/cluster pair,
# this menas that the values are the same, and we simply select the maximum value

#dat <- dcast(subdat, q_orgid+q_clustid~h_orgname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
#mat <- (data.matrix(dat[,3:ncol(dat)]))
#cvs_color <- sample(rainbow(length(cvs)), length(cvs_color))
#heatmap.2(mat, dendrogram="both", trace = "none", Rowv=TRUE, Colv=TRUE ,density.info="none", margins =c(6,10))

# dat <- dcast(subdat, org_id+clust_id~h_org, value.var="clustCov", fill=0) # WOOOOOOOOOO 
# mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
# colnames(mat) <- dat$clust_id

# cvs <- unique(dat$org_id)
# #cvs_color <- sample(rainbow(length(cvs)), length(cvs_color))
# cvs_color <- rainbow(length(cvs))
# ccol_stages <- sapply(dat$org_id, function(x) cvs_color[x==cvs])

# b <- c(0,0.2,0.4,0.6,0.8,1)

# subdat <- subset(tmpdat, tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri" & tmpdat$org_id=="22", select=c(1,2,4,7))

# dat <- dcast(subdat, org_id+clust_id~h_org, value.var="clustCov", fill=0) # WOOOOOOOOOO 
# mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
# colnames(mat) <- dat$clust_id

# cvs <- unique(dat$org_id)
# #cvs_color <- sample(rainbow(length(cvs)), length(cvs_color))
# cvs_color <- rainbow(length(cvs))
# ccol_stages <- sapply(dat$org_id, function(x) cvs_color[x==cvs])

# b <- c(0,0.2,0.4,0.6,0.8,1)

# pdf('orgs2cluster-Aspni7-cutoff-clustered.pdf')
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="Aspni7 nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

# pdf('orgs2cluster-Aspni7-range-clustered.pdf')
# colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="Aspni7 nigri clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

# pdf('orgs2cluster-Aspni7-cutoff-unclustered.pdf')
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="Aspni7 nigri clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

# pdf('orgs2cluster-Aspni7-cutoff-orgclustered.pdf')
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="Aspni7 nigri clusters, organism clustered", dendrogram="row", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()












# subdat <- subset(tmpdat, tmpdat$clust_size>10 & tmpdat$sec1=="Nigri" & tmpdat$sec2=="Nigri", select=c(1,2,4,7))

# dat <- dcast(subdat, org_id+clust_id~h_org, value.var="clustCov", fill=0) # WOOOOOOOOOO 
# mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
# colnames(mat) <- dat$clust_id

# cvs <- unique(dat$org_id)
# #cvs_color <- sample(rainbow(length(cvs)), length(cvs_color))
# cvs_color <- rainbow(length(cvs))
# ccol_stages <- sapply(dat$org_id, function(x) cvs_color[x==cvs])

# b <- c(0,0.2,0.4,0.6,0.8,1)

# pdf('orgs2cluster-10-cutoff-clustered.pdf')
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="Size >10 nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

# pdf('orgs2cluster-10-range-clustered.pdf')
# colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="Size >10 nigri clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

# pdf('orgs2cluster-10-cutoff-unclustered.pdf')
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="Size >10 nigri clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

# pdf('orgs2cluster-10-cutoff-orgclustered.pdf')
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="Size >10 nigri clusters, organism clustered", dendrogram="row", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

# pdf('orgs2cluster-color-pie.pdf')
# pie(rep(1,length(cvs_color)), col=cvs_color)
# dev.off()

# """
# dat <- cast(subdat, h_org~clust_id+org_id, value="clustCov", fill=0)
# mat <- data.matrix(dat)
# rownames(mat) <- dat[,1]
# mat <- mat[,-1]

# heatmap.2(mat, dendrogram="row", trace = "none", col=colfunc)
# heatmap.2(mat, dendrogram="both", trace = "none", col=colfunc)

# # Cluster vs cluster
# tmpdat <- read.csv("query_main.csv", sep=";", header=TRUE)
# subdat <- subset(tmpdat, select=c(1,10,19))
# dat <- cast(subdat, ta_clust_id~tb_clust_id, value="over", fill=0)
# mat <- data.matrix(dat)
# rownames(mat) <- dat[,1]
# mat <- mat[,-1:-2]
# heatmap.2(mat, dendrogram="both", trace = "none", col=colfunc)
# """
#########################################################
### A) Installing and loading required packages
#########################################################
library(ggplot2)
args <- commandArgs(TRUE)
infile <- args[1]
print (infile)


if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
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
library(Matrix)
library(gridExtra)


tmpdat <- read.table("tab_dataset.csv", sep=",", header=TRUE, dec=".", na.strings = "NA")
subdat <- subset(tmpdat, tmpdat$Project=="300asp") # q_orgid, q_clustid, h_realname, clustCov
subdat <- subdat[order(-xtfrm(subdat$X.2)),]
#subdat$Organism.name <- factor(subdat$Organism.name, levels = subdat[order(subdat$Assembly.length),"Organism.name"])
#order(-xtfrm(factor(subdat$X.2, levels=unique(subdat$X.2))))

#pdf('FGC_assemblyLengths.pdf', width=4, height=6)
g1 <- ggplot(subdat, aes(x = X.2, y=Assembly.length/1000000)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
ggtitle("Assembly length") + 
geom_text(aes(label=trunc(Assembly.length/1000000)), vjust=+0.5, hjust=+1.2, colour="white") + 
xlab("Organism") + 
ylab("Mega bases")
#dev.off()

g2 <- ggplot(subdat, aes(x = X.2, y=Scaffolds)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
ggtitle("Number of scaffolds") + 
geom_text(aes(label=Scaffolds), vjust=+0.5, hjust=+1.2, colour="white") + 
xlab("Organism") + 
ylab("Scaffolds")

g3 <- ggplot(subdat, aes(x = X.2, y=Scaffold.N50)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
ggtitle("Number of scaffolds needed to cover 50% of genome") + 
geom_text(aes(label=Scaffold.N50), vjust=+0.5, hjust=+1.2, colour="white") + 
xlab("Organism") + 
ylab("N50")


g4 <- ggplot(subdat, aes(x = X.2, y=Genes.with.complete.CDS..has.both.start.and.stop.codons...)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
ggtitle("Genes.with.complete.CDS..has.both.start.and.stop.codons...:") + 
geom_text(aes(label=Genes.with.complete.CDS..has.both.start.and.stop.codons...), vjust=+0.5, hjust=+1.2, colour="white") + 
xlab("Organism") + 
ylab("%")

g5 <- ggplot(subdat, aes(x = X.2, y=Genes.with.complete.CDS..has.both.start.and.stop.codons.)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
ggtitle("Genes.with.complete.CDS..has.both.start.and.stop.codons.") + 
geom_text(aes(label=Genes.with.complete.CDS..has.both.start.and.stop.codons.), vjust=+0.5, hjust=+1.2, colour="white") + 
xlab("Organism") + 
ylab("Gene count")


g6 <- ggplot(subdat, aes(x = X.2, y=Genes.with.complete.CDS..has.both.start.and.stop.codons.)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
ggtitle("Genes.with.complete.CDS..has.both.start.and.stop.codons.") + 
geom_text(aes(label=Genes.with.complete.CDS..has.both.start.and.stop.codons.), vjust=+0.5, hjust=+1.2, colour="white") + 
xlab("Organism") + 
ylab("Gene count")

"""--------------------------
ANTISMASH
--------------------------"""
atmpdat <- read.table("antismash_cluster_type_counts.csv", sep=";", header=TRUE, dec=".", na.strings = "NA")
asubdat <- atmpdat[atmpdat$SpeciesName %in% subdat$X.2, ]
#asubdat <- asubdat[order(-xtfrm(asubdat$SpeciesName)),SpeciesName]
ggplot(asubdat, aes(x = asubdat$SpeciesName, y=nrPKS)) + geom_bar(stat="identity") + coord_flip()

#amdat$n <- with(amdat, relevel(SpeciesName, "violaceofuscus"))

amdat <- melt(asubdat)
#amdat <- amdat[order(xtfrm(amdat$SpeciesName)),]
amdat$n <- ordered(asubdat$SpeciesName, levels = rev(levels(asubdat$SpeciesName)))

a1 <- ggplot(amdat, aes(x = n, y=value, fill=variable)) + geom_bar(stat="identity")  + coord_flip()
+ #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=asubdat$X.2, labels=asubdat$X.2) + 
ggtitle("Assembly length") + 
geom_text(aes(label=trunc(Assembly.length/1000000)), vjust=+0.5, hjust=+1.2, colour="white") + 
xlab("Organism") + 
ylab("Mega bases")

grid.arrange(g6, a1, ncol = 2)



  



mat <- (data.matrix(tmpdat[,2:ncol(tmpdat)]))
rownames(mat) <-  tmpdat[,1]

colfunc <- colorRampPalette(c("white", "green", "blue"))(n = 10)
heatmap.2((mat), density.info="none", margins =c(10,10), trace = "none",dendrogram="none", na.color="red",Rowv=FALSE, Colv=FALSE, col=colfunc)


b <- c(seq(0, 100, by=10))
colfunc <- brewer.pal((length(b)-1),"RdBu")
heatmap.2((mat), density.info="none", margins =c(20,20), trace = "none",na.color="white",Rowv=FALSE, Colv=FALSE, col=colfunc, breaks=b)

b <- c(seq(0, 100, by=5))
colfunc <- colorRampPalette(c("white", "green", "blue"))(n = (length(b)-1))
heatmap.2((mat), density.info="none", margins =c(20,20), trace = "none",dendrogram="none", na.color="white",Rowv=FALSE, Colv=FALSE, col=colfunc, breaks=b)












tmpdat <- read.csv(infile, sep=";", header=TRUE)
tmpdat <- read.csv("allCounts_nr.csv", sep=";", header=TRUE)
tmpdat$q_per <- (tmpdat$q_count/tmpdat$count)*100
tmpdat$h_per <- (tmpdat$h_count/tmpdat$count)*100
mat <- (data.matrix(tmpdat))
heatmap.2(mat, density.info="none", margins =c(10,10), symbreaks = min(mat, na.rm=TRUE),na.color="blue")

orgvector <- unique(c(levels(tmpdat$q_org), levels(tmpdat$h_org)))

dat <- dcast(tmpdat, q_org~h_org, value.var="q_per", fill=0) # WOOOOOOOOOO 
mat <- (data.matrix(dat)) # NOTE transpose
colnames(mat) <- dat$q_org

cvs <- unique(dat$q_orgid)
cvs_color <- rainbow(length(cvs))
ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])

b <- c(0,0.2,0.4,0.6,0.8,1)

#pdf('orgs2cluster-all-cutoff-clustered-bestcand.pdf')
colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
heatmap.2(mat, main="All nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
#dev.off()

# pdf('orgs2cluster-all-range-clustered-bestcand.pdf')
# colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="All nigri clusters, clustered",dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

# pdf('orgs2cluster-all-cutoff-unclustered-bestcand.pdf')
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="All nigri clusters, unclustered", dendrogram="none", trace = "none", ColSideColor=ccol_stages, Rowv=FALSE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

# pdf('orgs2cluster-all-cutoff-orgclustered-bestcand.pdf')
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="All nigri clusters, organism clustered", dendrogram="row", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=FALSE , col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()


# subdat <- subset(tmpdat, tmpdat$q_clust_size>10 & tmpdat$q_sec=="Nigri" & tmpdat$h_sec=="Nigri" & tmpdat$q_clust_size>10, select=c(12,3,6,8)) # q_orgid, q_clustid, h_realname, clustCov
# dat <- dcast(subdat, q_orgid+q_clustid~h_realname, value.var="clustCov", fill=0, max) # WOOOOOOOOOO 
# mat <- t(data.matrix(dat[,3:ncol(dat)])) # NOTE transpose
# colnames(mat) <- dat$q_clustid
# cvs <- unique(dat$q_orgid)
# cvs_color <- rainbow(length(cvs))
# ccol_stages <- sapply(dat$q_orgid, function(x) cvs_color[x==cvs])
# b <- c(0,0.2,0.4,0.6,0.8,1)

# pdf('orgs2cluster-10-cutoff-clustered-bestcand.pdf')
# colfunc <- c("white"  ,"white", "white" ,"#2B8CBE","#084081")
# heatmap.2(mat, main="Size >10 nigri clusters, mostly white, clustered", dendrogram="both", trace = "none", ColSideColor=ccol_stages, Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, density.info="none", margins =c(4,10))
# dev.off()

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
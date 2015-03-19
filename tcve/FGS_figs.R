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

pdf('FGC_assemblyLengths.pdf', width=8, height=8)
g1 <- ggplot(subdat, aes(x = X.2, y=Assembly.length/1000000)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
geom_text(aes(label=trunc(Assembly.length/1000000)), vjust=+0.5, hjust=+1.2, colour="white") + 
theme(axis.text = element_text(size = 12, colour="grey20")) +
ggtitle("Assembly length") + 
xlab("Organism") + 
ylab("Mega bases")
g1
dev.off()

pdf('FGC_nrScaffolds.pdf', width=8, height=8)
g2 <- ggplot(subdat, aes(x = X.2, y=Scaffolds)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
geom_text(aes(label=Scaffolds), vjust=+0.5, hjust=+1.2, colour="white") + 
theme(axis.text = element_text(size = 12, colour="grey20")) +
ggtitle("Number of scaffolds") + 
xlab("Organism") + 
ylab("Scaffolds")
g2
dev.off()

pdf('FGC_n50Scaffolds.pdf', width=8, height=8)
g3 <- ggplot(subdat, aes(x = X.2, y=Scaffold.N50)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
geom_text(aes(label=Scaffold.N50), vjust=+0.5, hjust=+1.2, colour="white") + 
theme(axis.text = element_text(size = 12, colour="grey20")) +
ggtitle("Number of scaffolds needed to cover 50% of genome") + 
xlab("Organism") + 
ylab("N50")
g3
dev.off()

pdf('FGC_completeCDSpercent.pdf', width=8, height=8)
g4 <- ggplot(subdat, aes(x = X.2, y=Genes.with.complete.CDS..has.both.start.and.stop.codons...)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
geom_text(aes(label=Genes.with.complete.CDS..has.both.start.and.stop.codons...), vjust=+0.5, hjust=+1.2, colour="white") + 
theme(axis.text = element_text(size = 12, colour="grey20")) +
ggtitle("Genes.with.complete.CDS..has.both.start.and.stop.codons...:") + 
xlab("Organism") + 
ylab("%")
g4
dev.off()

pdf('FGC_completeCDScount.pdf', width=8, height=8)
g5 <- ggplot(subdat, aes(x = X.2, y=Genes.with.complete.CDS..has.both.start.and.stop.codons.)) + 
geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
coord_flip() + 
scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
geom_text(aes(label=Genes.with.complete.CDS..has.both.start.and.stop.codons.), vjust=+0.5, hjust=+1.2, colour="white") + 
theme(axis.text = element_text(size = 12, colour="grey20")) +
ggtitle("Genes.with.complete.CDS..has.both.start.and.stop.codons.") + 
xlab("Organism") + 
ylab("Gene count")
g5
dev.off()



#_______________________________________________________
# totals <- aggregate(msubdat$value, by=list(msubdat$n), FUN=sum)
# labels <- c(as.vector(totals$x) , rep("", times=length(totals$x)) )
# pos <- rep(as.vector(totals$x) +1 , each=2)
# msubdat <- data.frame(msubdat, pos=pos, labels=labels)

# ggplot(msubdat, aes(x = n, y = value)) +
# geom_bar(aes(fill = variable), stat="identity") +
# geom_text(aes(label = labels, pos=pos), vjust=+0.5, hjust=+0.5, colour="white", position = "stack") +
# scale_fill_manual(name = "Measure",values = c("#d7301f", "#2b8cbe"), labels=c("N50","Scaffolds"))
#_______________________________________________________
subdat$n <- ordered(subdat$X.2, levels = rev(levels(subdat$X.2)))
subdat$left <- (subdat$Scaffolds-subdat$Scaffold.N50)
tmsubdat <- melt(subdat)
msubdat <- subset(tmsubdat, tmsubdat$variable=="left" | tmsubdat$variable=="Scaffold.N50" )


totals <- aggregate(msubdat$value, by=list(msubdat$n), FUN=sum)
labels <- c(subset(msubdat, msubdat$variable=="Scaffold.N50")$value , as.vector(totals$x) )
msubdat$labels <- labels

pdf('FGC_scaffolds_N50.pdf', width=10, height=8)
ggplot(msubdat, aes(x = n, y = value)) +
geom_bar(aes(fill = variable), stat="identity") +
geom_text(aes(label = as.numeric(labels)), position = position_stack(width=1), vjust=+0.5, hjust=0, colour="grey20") + 
scale_fill_manual(name = "Measure",values = c("#d7301f", "#2b8cbe"), labels=c("N50","Scaffolds")) +
scale_y_continuous(limits = c(0, 530)) + 
theme(axis.text = element_text(size = 12, colour="grey20")) +
coord_flip() + 
ggtitle("Nr. scaffolds and N50") + 
xlab("Organism") + 
ylab("Scaffolds and N50")
dev.off()

# subdat$n <- ordered(subdat$X.2, levels = rev(levels(subdat$X.2)))
# subdat$left <- (subdat$X.2, levels = rev(levels(subdat$X.2)))
# msubdat <- melt(subdat)
# msubdat <- subset(msubdat, msubdat$variable=="Scaffolds" | msubdat$variable=="Scaffold.N50")

# g6 <- ggplot(msubdat, aes(x=n, y=value, fill=variable)) + 
# geom_bar(position="dodge",stat="identity") + #, color="#084081")+ 
# coord_flip() + 
# #scale_x_discrete( limits=msubdat$X.2, labels=msubdat$X.2) + 
# ggtitle("Nr. scaffolds and N50") + 
# geom_text(aes(label=value), vjust=+0.5, hjust=+1.2, colour="white", position=position_dodge(width=0.9), size=4) + 
# xlab("Organism") + 
# ylab("Scaffolds and N50") +
# scale_fill_manual(name = "Measure",values = c("#084081","#b30000"), labels=c("Scaffolds", "N50"))





# pdf('FGC_completeCDScount.pdf', width=8, height=8)
# g6 <- ggplot(subdat, aes(x = X.2, y=Genes.with.complete.CDS..has.both.start.and.stop.codons.)) + 
# geom_bar(stat="identity", fill="#2B8CBE") + #, color="#084081")+ 
# coord_flip() + 
# scale_x_discrete( limits=subdat$X.2, labels=subdat$X.2) + 
# ggtitle("Genes.with.complete.CDS..has.both.start.and.stop.codons.") + 
# geom_text(aes(label=Genes.with.complete.CDS..has.both.start.and.stop.codons.), vjust=+0.5, hjust=+1.2, colour="white") + 
# xlab("Organism") + 
# ylab("Gene count")
# g6
# dev.off()

"""--------------------------
ANTISMASH
--------------------------"""

# GnBu: {
# 		3: []string{"#e0f3db", "#a8ddb5", "#43a2ca"},
# 		4: []string{"#f0f9e8", "#bae4bc", "#7bccc4", "#2b8cbe"},
# 		5: []string{"#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac"},
# 		6: []string{"#f0f9e8", "#ccebc5", "#a8ddb5", "#7bccc4", "#43a2ca", "#0868ac"},
# 		7: []string{"#f0f9e8", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#08589e"},
# 		8: []string{"#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#08589e"},
# 		9: []string{"#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081"},

# OrRd: {
# 		3: []string{"#fee8c8", "#fdbb84", "#e34a33"},
# 		4: []string{"#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f"},
# 		5: []string{"#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"},
# 		6: []string{"#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#b30000"},
# 		7: []string{"#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"},
# 		8: []string{"#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"},
# 		9: []string{"#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000"},

atmpdat <- read.table("../antismash_cluster_type_counts.csv", sep=";", header=TRUE, dec=".", na.strings = "NA")
asubdat <- atmpdat[atmpdat$SpeciesName %in% subdat$X.2, ]
asubdat$total <- asubdat$nrPKS + asubdat$nrPKSlike + asubdat$nrNRPS + asubdat$nrNRPDlike + asubdat$nrTC + asubdat$nrHYBRID + asubdat$nrDMAT
asubdat$n <- ordered(asubdat$SpeciesName, levels = rev(levels(asubdat$SpeciesName)))

pdf('antismash_cluster_counts.pdf', width=8, height=8)
a1 <- ggplot(asubdat, aes(x = asubdat$n, y=total)) + 
geom_bar(stat="identity", fill="#2B8CBE") + 
coord_flip() +
geom_text(aes(label=total), vjust=+0.5, hjust=+1.2, colour="white") + 
theme(axis.text = element_text(size = 12, colour="grey20")) +
ggtitle("Nr. of predicted SM clusters") + 
xlab("Organism") + 
ylab("Nr. sec. metabolism clusters")
a1
dev.off()

asubdat$total <- NULL
amdat <- melt(asubdat)
amdat$n <- ordered(asubdat$SpeciesName, levels = rev(levels(asubdat$SpeciesName)))
myColors <- brewer.pal(8,"GnBu")
myColors <- c("#084081","#4eb3d3", "#2b8cbe", "#0868ac","#fc8d59", "#e34a33", "#b30000")
names(myColors) <-levels(amdat$variable)


pdf('antismash_cluster_type_counts.pdf', width=8.5, height=8)
a2 <- ggplot(amdat, aes(x = n, y=value, fill=variable)) + 
geom_bar(stat="identity") + 
coord_flip() +
scale_fill_manual(name = "Backbone type",values = myColors, labels=c("PKS", "PKS like", "NRPS", "NRPD", "TC", "Hybrid", "DMAT"))+
theme(axis.text = element_text(size = 12, colour="grey20")) +
ggtitle("Nr. and type of predicted SM clusters") + 
xlab("Organism") + 
ylab("Nr. sec. metabolism clusters")
a2
dev.off()


# mysql -u asp -p1234 -e "select metaboliteMap.* from metaboliteMap join organism using (name) where section='Nigri';" aspminedb > metaboliteMap.txt;
# mysql -u asp -p1234 -e "select metaboliteMap.* from metaboliteMap join organism using (name);" aspminedb > metaboliteMap.txt;
mtmpdat <- read.table("/home/tcve/Dropbox/FGS_asilomar/metaboliteMap.txt", sep="\t", header=TRUE, dec=".", na.strings = "NA")
#msubdat <- mtmpdat[mtmpdat$real_name %in% subdat$X.2, ]
msubdat <- mtmpdat
msubdat$status <-  as.numeric(as.character(msubdat$status))
b <- c(0,0.5,1)
colfunc <- c("white","#7f0000")

dat <- dcast(msubdat, real_name~metabolite2, value.var="status", fill=0, max)
mat1 <- t(data.matrix(dat[,2:ncol(dat)])) # NOTE transpose
colnames(mat1) <- dat$real_name

dat <- dcast(msubdat, real_name~metabolite1, value.var="status", fill=0, max)
mat2 <- t(data.matrix(dat[,2:ncol(dat)])) # NOTE transpose
colnames(mat2) <- dat$real_name

pdf('metaboliteMap_met2_clust.pdf', width=9.5, height=9)
heatmap.2(mat1, main="Metabolites in organisms", sepcolor="grey50", colsep=c(0:length(rownames(mat))), rowsep=c(0:length(colnames(mat))), sepwidth=c(0.01, 0.01), 
	dendrogram="both", trace = "none", Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, density.info="none", margins =c(10,10))
dev.off()

pdf('metaboliteMap_met2_noclust.pdf', width=9.5, height=9)
heatmap.2(mat1, main="Metabolites in organisms", sepcolor="grey50", colsep=c(0:length(rownames(mat))), rowsep=c(0:length(colnames(mat))), sepwidth=c(0.01, 0.01), 
	dendrogram="both", trace = "none", Rowv=FALSE, Colv=FALSE ,col=colfunc, breaks=b, density.info="none", margins =c(10,10))
dev.off()

pdf('metaboliteMap_met1_clust.pdf', width=9.5, height=9)
heatmap.2(mat2, main="Metabolites in organisms", sepcolor="grey50", colsep=c(0:length(rownames(mat))), rowsep=c(0:length(colnames(mat))), sepwidth=c(0.01, 0.01), 
	dendrogram="both", trace = "none", Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, density.info="none", margins =c(10,10))
dev.off()

pdf('metaboliteMap_met1_noclust.pdf', width=9.5, height=9)
heatmap.2(mat2, main="Metabolites in organisms", sepcolor="grey50", colsep=c(0:length(rownames(mat))), rowsep=c(0:length(colnames(mat))), sepwidth=c(0.01, 0.01), 
	dendrogram="both", trace = "none", Rowv=FALSE, Colv=FALSE ,col=colfunc, breaks=b, density.info="none", margins =c(10,10))
dev.off()
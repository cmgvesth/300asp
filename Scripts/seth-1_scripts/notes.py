notes.py



asubdat$total <- NULL
amdat <- melt(asubdat)
amdat$n <- ordered(asubdat$SpeciesName, levels = rev(levels(asubdat$SpeciesName)))
myColors <- brewer.pal(8,"GnBu")

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

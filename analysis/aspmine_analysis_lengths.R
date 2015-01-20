library(ggplot2)

# org_id, section, name, real_name, len
# colClasses=c("character",rep("numeric",5))
dat <- read.csv(file="tmp.csv", header = FALSE, sep=",", fill= TRUE, colClasses=c("numeric", rep("factor",3), "numeric"), col.names=c("org_id", "section", "name", "real_name", "len"))

hist_section <- ggplot(dat, aes(x=len, fill=section)) + geom_bar(binwidth=100) + facet_grid(~section) + xlim(0,2000) +
ggtitle("Histograms, binwidth 100 amino acids, plot per section, limit to 2000 amino acids")
ggsave(file="aalength_histSection.pdf", width=12, height=8)

dens_bar_section <- ggplot(dat, aes(x=len, fill=section)) + geom_bar(binwidth=100) + facet_grid(~section) + aes(y = ..density..) + xlim(0,2000) +
ggtitle("Density barplot, binwidth 100 amino acids, plot per section, limit to 2000 amino acids")
ggsave(file="aalength_densBarSection.pdf", width=12, height=8)

dens_graph_section <- ggplot(dat, aes(len, fill=section)) + geom_density(alpha = 0.2) + xlim(0,2000) + 
ggtitle("Density graph, color by section, limit to 2000 amino acids")
ggsave(file="aalength_densGraphSection.pdf", width=12, height=8)

boxplot_section <- ggplot(dat, aes(y=len, x=section, fill=section)) + geom_boxplot() + coord_flip() + 
ggtitle("Boxplot, color by section")
ggsave(file="aalength_boxplotSection.pdf", width=12, height=12)




# summary(subset(dat, section=="flavi"))
# for (i in unique(dat$section)) { 	print(summary(subset(dat, section==i))) }
for (i in 1:length(levels(dat$section)) ) { 
	set = dat[dat$section==dat$section[i],]
	print(i, summary(set$section))
	filename = paste("aalength_histogram", dat$section[i], ".pdf", sep="")
	title = paste("Histograms, binwidth 100 amino acids, plot per organism, section ", dat$section[i],  sep="")
	m <- ggplot(dat, aes(x=len)) + geom_bar(binwidth=100) + facet_grid(section~org_id) + ggtitle(title)
	ggsave(file=filename, width=12, height=12)
}

hist_org <- ggplot(dat, aes(x=len)) + geom_bar(binwidth=100) + facet_grid(section~org_id) +
ggtitle("Histograms, binwidth 100 amino acids, plot per organism")

dens_bar_org <- ggplot(dat, aes(x=len)) + geom_bar(binwidth=100) + facet_grid(section~org_id) + aes(y = ..density..) +
ggtitle("Density barplot, binwidth 100 amino acids, plot per section")


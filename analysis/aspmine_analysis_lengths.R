library(ggplot2)
args <- commandArgs(TRUE)
infile <- args[1]
print (infile)

dat <- read.csv(file=infile, header = FALSE, sep=",", fill= TRUE, colClasses=c("numeric", rep("factor",3), "numeric"), col.names=c("org_id", "section", "name", "real_name", "len"))

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

for (s in unique(dat$section)) {
	set <- subset(dat, section==s)

	filename = paste("aalength_histogram_section_", s, ".pdf", sep="")
	title = paste("Histograms, binwidth 500 amino acids, plot per organism, section ", s,  sep="")
	a <- ggplot(set, aes(x=len, fill=section)) 
	a <- a + geom_bar(binwidth=500) + facet_grid(~name) + xlim(0,2000) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	a <- a + labs(title = title, x = "Organisms, alphabetically", y = "Protein counts") 
	ggsave(plot = a, file=filename, width=12+length(unique(set$org_id)), height=12)

	filename = paste("aalength_boxplots_section_", s, ".pdf", sep="")
	title = paste("Boxplot, orderd by median length, plot per organism, section ", s,  sep="")
	b <- ggplot(set, aes(y=len, x=reorder(name, len, FUN=median))) 
	b <- b + geom_boxplot() + coord_flip() 
	b <- b + labs(title = title, x = "Organisms, ordered by median", y = "protein lengths") 
	ggsave(plot = b, file=filename, height=12+length(unique(set$org_id)), width=12)
}

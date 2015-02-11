
# R CMD BATCH --args query_result.csv ipr_PKS_plots.R test.out 
"""
head IPR_parts.csv
"ipr_id";"nr_proteins";"nr_orgs";"ipr_id";"ipr_desc";"ipr_domaindb";"ipr_domain_id";"ipr_domaindesc"
"IPR000092";1;1;"IPR000092";"Polyprenyl synthetase";"HMMPfam";"PF00348";"polyprenyl_synt"
"IPR000092";1;1;"IPR000092";"Polyprenyl synthetase";"ScanRegExp";"PS00444";"POLYPRENYL_SYNTHET_2"
"IPR000092";1;1;"IPR000092";"Polyprenyl synthetase";"ScanRegExp";"PS00723";"POLYPRENYL_SYNTHET_1"
"IPR000169";6;5;"IPR000169";"Peptidase, cysteine peptidase active site";"ScanRegExp";"PS00139";"THIOL_PROTEASE_CYS"
"IPR000169";6;5;"IPR000169";"Peptidase, cysteine peptidase active site";"ScanRegExp";"PS00639";"THIOL_PROTEASE_HIS"
"IPR000173";3;3;"IPR000173";"Glyceraldehyde 3-phosphate dehydrogenase";"HMMPfam";"PF00044";"Gp_dh_N"
"IPR000173";3;3;"IPR000173";"Glyceraldehyde 3-phosphate dehydrogenase";"HMMPfam";"PF02800";"Gp_dh_C"
"IPR000173";3;3;"IPR000173";"Glyceraldehyde 3-phosphate dehydrogenase";"HMMPIR";"PIRSF000149";"Glyceraldehyde-3-phosphate dehydrogenase"
"IPR000173";3;3;"IPR000173";"Glyceraldehyde 3-phosphate dehydrogenase";"FPrintScan";"PR00078";"G3PDHDRGNASE"

"""
library(ggplot2)
args <- commandArgs(TRUE)
IPR_single <- args[1]

# colClasses=c("numeric", rep("factor",3), "numeric"), col.names=c("org_id", "section", "name", "real_name", "len")
dat <- read.csv(file=IPR_single, header = TRUE, sep=";", fill= TRUE)
dat1 <- subset(dat,nr_proteins>2)
dat <- dat1[order(dat1$nr_proteins),]
print str(dat)


ggplot(dat, aes(x=ipr_id, y=nr_proteins)) + geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

hist_ipr <- ggplot(dat, aes(x=ipr_id, y=nr_proteins)) + geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
 + ggtitle("Histogram, nr_proteins per Interpro domain")
ggsave(file="aalength_histogram_section_all.pdf", width=4+length(unique(dat$org_id)), height=4+length(unique(dat$org_id)), limitsize=FALSE)

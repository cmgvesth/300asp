library(shiny)
library(ggplot2)
library(reshape)
library(gplots)
library(RColorBrewer)
library(argparse)
library(ape)
library(RMySQL)
library(gridExtra)
library(Hmisc)

con <- dbConnect(MySQL(), user="setd", password="1234", dbname="aspminedb", host="192.38.13.196", port = 3306)

###########
# Summary #
###########




print("How many hitclusters are hit by a querycluster?")
# query <- "select family_hits from (select *, count(*) as family_hits from (select *,  pident*members/q_size as custom_score from (select *, count(*) as members from t_genecluster_bidir_hits where pident > 50 group by q_clust,clust_id) tx where members >= 0.5*q_size and members >= 0.25*h_size and q_size > 1) ty group by q_clust) tz;"
query <- "select family_hits, q_size, count(*) as fa_hit_count from (select *, count(*) as family_hits from (select *,  pident*members/q_size as custom_score from (select *, count(*) as members from t_genecluster_bidir_hits where pident > 50 group by q_clust,clust_id) tx where members >= 0.5*q_size and members >= 0.25*h_size and q_size >= 1) ty group by q_clust) tz group by family_hits, q_size;"
hits <- dbGetQuery(con, query)

#print(summary(hits))

hits <- subset(hits, family_hits>=2) # Decided to take the singletons out because they only match themselves

print(head(hits))

sub_hits <- subset(hits, family_hits>=10)

print(head(sub_hits))

###########################################################
# Abundance of multiple backbones in bidirectional hits ? #
###########################################################
print("How many gene clusters have multiple backbones?")

query <- "select buf3 as backbone, count(*) as amount from (select *, count(*) as members, sum(buf2) as buf3 from  (select *, case when q_sm_short != 'none' or h_sm_short != 'none' then  buf+1 else buf+0 end as buf2 from t_genecluster_bidir_hits where pident > 50) tw group by q_clust,clust_id) tx where members >= 0.25*q_size and members >= 0.25*h_size and q_size > 2 group by buf3;"

multi_b <- dbGetQuery(con, query)

print(multi_b)


# Shows the same just if backbone was found or not ;)

query <- "select buf3, count(*) from (select *, count(*) as members, max(buf2) as buf3 from  (select *, case when q_sm_short != 'none' or h_sm_short != 'none' then  buf+1 else buf+0 end as buf2 from t_genecluster_bidir_hits where pident > 50) tw group by q_clust,clust_id) tx where members >= 0.25*q_size and members >= 0.25*h_size and q_size > 2 group by buf3;"

has_b <- dbGetQuery(con, query)

print(has_b)



##############
# Singletons #
##############

query <- "select real_name, singleton  from (select org_id, count(*) as singleton from (select *, count(*) as family_hits from (select *,  pident*members/q_size as custom_score from (select *, count(*) as members from t_genecluster_bidir_hits where pident > 50 group by q_clust,clust_id) tx where members >= 0.5*q_size and members >= 0.25*h_size and q_size > 1) ty group by q_clust) tz where family_hits <= 1 and q_clust = clust_id group by org_id) ta join organism using(org_id);"

singletons <- dbGetQuery(con, query)

print(singletons)


# under construction query to fetch select q_clust from (select *, count(*) as family_hits from (select *,  pident*members/q_size as custom_score from (select *, count(*) as members from t_genecluster_bidir_hits where pident > 50 group by q_clust,clust_id) tx where members >= 0.5*q_size and members >= 0.25*h_size and q_size > 1) ty group by q_clust) tz where family_hits <= 1 and q_clust = clust_id;

# select distinct q_gene, org_id from t_genecluster_bidir_hits where q_clust = '9_8731_9';

# select * from (select distinct q_gene, go_term_id, q_sm_short from t_genecluster_bidir_hits as tgbh join protein_has_go as phg on q_gene = protein_id and tgbh.org_id = phg.org_id  where q_clust = '9_8731_9') tx join go using (go_term_id) order by q_gene;


############
# Printing #
############

print("Printing to Genecluster_family_preprocessing.pdf")

pdf('Genecluster_family_preprocessing.pdf')


#qplot(family_hits,data=hits, geom='histogram', binwidth=1, main = "Hitclusters matching to queryclusters")

#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#sc <- scale_colour_gradient2(colours = myPalette(100), limits=c(1, 50))

#sc <- scale_fill_continuous(low = "#D0D1E6", high = "#023858",limits=c(1, 50), breaks = c(0,5,10,15,20,25,30,35,40,45,50))

qplot(family_hits, y=fa_hit_count, stat="identity", data=hits, geom="bar", fill=factor(q_size)) + guides(col = guide_legend(ncol = 3)) #theme(legend.text = element_text( size = 1))
#qplot(family_hits,data=sub_hits, geom='histogram', binwidth=1, main = "Hitclusters matching to queryclusters with more than 10 hits")
qplot(family_hits, y=fa_hit_count, stat="identity", data=sub_hits, geom="bar", fill=factor(q_size))

ggplot(singletons, aes(x = real_name, y = singleton))+
geom_bar(stat = "identity")+
coord_flip()

frame()

grid.table(multi_b)


dev.off()


#create table t_genecluster_testing as select *,  pident*members/q_size as custom_score from (select *, count(*) as members from t_genecluster_bidir_hits where pident > 50 group by q_clust,clust_id) tx where members >= 0.5*q_size and members >= 0.5*h_size and q_size > 2;

#create table t_gcc_testing as select * from (select *, count(*) as family_hits from (select *,  pident*members/q_size as custom_score from (select *, count(*) as members from t_genecluster_bidir_hits where pident > 50 group by q_clust,clust_id) tx where members >= 0.5*q_size and members >= 0.5*h_size and q_size > 1) ty group by q_clust) tz where family_hits >10;
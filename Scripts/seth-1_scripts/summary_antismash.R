library(shiny)
library(ggplot2)
library(reshape)
library(gplots)
library(RColorBrewer)
library(argparse)
library(ape)
library(RMySQL)


con <- dbConnect(MySQL(), user="setd", password="1234", dbname="aspminedb", host="192.38.13.196", port = 3306)

# see how many gene clusters have multiple sm_short descriptions, i.e. multiple backbones
query <- "select multi_bb, count(*) from (select concat(org_id,'_',clust_backbone,'_', clust_size), count(*) as multi_bb from antismash where sm_short != 'none' group by clust_backbone, sm_short) tx group by multi_bb;"
multi_bb <- dbGetQuery(con, query)

print(multi_bb)

query <- "select real_name, sm_short, count(*) as type from antismash join organism using (org_id) where sm_short != 'none' group by org_id, sm_short;"
clust_types <- dbGetQuery(con, query)


pdf('antismash_overview.pdf')
qplot(sm_short, type, data =clust_types, geom=c("boxplot","jitter"), fill=sm_short, main="Cluster backbone types among Aspergilli", ylab="Number of clusters with backbone")

qplot(real_name, y=type, stat="identity", data=clust_types, geom="bar", fill=sm_short) + coord_flip()
dev.off()


# In Progress:

# Find all clusters which have two or more backbones and see if they are in the same location
select *, count(*) from (select clust_id from (select *,concat(org_id,'_',clust_backbone,'_', clust_size) as clust_id, count(*) as multi_bb from antismash where sm_short != 'none' group by clust_backbone, sm_short) tx where multi_bb >= 2) tz join (select *,concat(org_id,'_',clust_backbone,'_', clust_size) as clust_id from antismash where sm_short != 'none' ) ty using (clust_id) group by clust_start;
# Result: I get 171 rows which is the number of all cases with multiple backbones I observed. Therefore, the multiple backbone is only Antsmash derived and not nature of our data

select * from (select *,concat(org_id,'_',clust_backbone,'_', clust_size) as clust_id from antismash where sm_short != 'none' ) tx where clust_id = '40_143764_25';


# For the sm_short thing:
select case when q_sm_short!='none' then custom_score +1 else custom_score +0 end as custom_score from t_genecluster_testing;

library(shiny)
library(ggplot2)
library(reshape)
library(gplots)
library(RColorBrewer)
library(argparse)
library(ape)
library(RMySQL)
library(dendextend)
library(gridExtra)


con <- dbConnect(MySQL(), user="setd", password="1234", dbname="aspminedb", host="192.38.13.196", port = 3306)

#query <- "select ta.* from t_genecluster_testing as ta join t_gcc_testing using(q_clust);"
query2 <- "select * from t_genecluster_testing_full;" 
#query <- "select *,  pident*members/q_size as custom_score from (select *, count(*) as members from t_genecluster_bidir_hits where pident > 50 group by q_clust,clust_id) tx where members >= 0.5*q_size and members >= 0.5*h_size and q_size > 1;"

# Newest version cutting hits which are too much out
query <- "select * from (select ta.* from t_genecluster_testing_full as ta join t_gcc_testing_full using(q_clust)) tb join (select distinct(q_clust) as clus from t_gcc_testing_full) tc on clus = clust_id"
# select count(distinct(q_clust)), count(distinct(clust_id)) from (select * from (select ta.* from t_genecluster_testing_full as ta join t_gcc_testing_full using(q_clust)) tb join (select distinct(q_clust) as clus from t_gcc_testing_full) tc on clus = clust_id) tz;

gcc <- dbGetQuery(con, query)

data_wide <- cast(gcc, q_clust ~ clust_id,value = 'custom_score')


data_wide[is.na(data_wide)] <-0


wssplot <- function(data_wide, nc=100, seed=1234){
               wss <- (nrow(data_wide)-1)*sum(apply(data_wide,2,var))
               for (i in 2:nc){
                    set.seed(seed)
                    wss[i] <- sum(kmeans(data_wide, centers=i)$withinss)}
                plot(1:nc, wss, type="b", xlab="Number of Clusters",
                     ylab="Within groups sum of squares")}

pdf("testing.pdf")
wssplot(data_wide)

mat <- data.matrix(data_wide)
rownames(mat) <- data_wide[,1]
mat <- mat[,-1]
mat[is.na(mat)] <- 0


    colfunc <- brewer.pal(9,"Blues")

 hv <-  heatmap.2(mat,
        main = "gcc",
        dendrogram = 'both',
        margins=c(6,12),
        na.color = 'white',
        col = colfunc,
        key = TRUE,
        #cexRow = 0.8, # Adjusts row fontsize later like in old script
        trace='none')
 
dev.off()

mytree <- hv$rowDendrogram
groups <- as.data.frame(dendextend::cutree(mytree,k=60)) # TODO fix hardcoding to make transfer easier
colnames(groups) <- "group"
print(head(groups))

groups$clust <- rownames(groups)
# rownames(groups)<-c(1:length(rownames(groups)))
row.names(groups)<- NULL
print(head(groups))
groups <- groups[order(groups$group),]
print(head(groups))

write.csv(groups,file = "gene_cluster_clustering.csv")

# Can be processed with gcc_analysis.R

#pdf("testing2.pdf")
#grid.table(groups,gpar.coretext = gpar(fontsize=6), gpar.coltext = gpar(fontsize=6), padding.h=unit(2, "mm"), padding.v=unit(2, "mm"))
#dev.off()
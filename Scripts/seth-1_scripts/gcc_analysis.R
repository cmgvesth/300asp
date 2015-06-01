library(ggplot2)
library(reshape)
library(gplots)
library(RColorBrewer)
library(argparse)
library(ape)
library(RMySQL)

# Needs gcc.R to be run first

groups <- read.csv("gene_cluster_clustering.csv")
groups$X <- NULL


con <- dbConnect(MySQL(), user="setd", password="1234", dbname="aspminedb", host="192.38.13.196", port = 3306)
pdf("single_groups.pdf")
for (i in seq(60)){
        collection <- groups[groups$group == i,]$clust

        selector <- paste(collection, collapse="','")



        query <- sprintf("select * from t_genecluster_testing where q_clust in ('%s');", selector)

        print(query)

        gcc <- dbGetQuery(con, query)


        data_wide <- cast(gcc, q_clust ~ clust_id,value = 'custom_score')


        data_wide[is.na(data_wide)] <-0




        mat <- data.matrix(data_wide)
        rownames(mat) <- data_wide[,1]
        mat <- mat[,-1]
        mat[is.na(mat)] <- 0


            colfunc <- brewer.pal(9,"Blues")

         hv <-  heatmap.2(mat,
                main = sprintf("Clustering for group %s", i),
                dendrogram = 'both',
                margins=c(6,8),
                na.color = 'white',
                col = colfunc,
                key = TRUE,
                cexRow = 0.8, # Adjusts row fontsize later like in old script
                cexCol = 0.8,
                trace='none')
}


dev.off()
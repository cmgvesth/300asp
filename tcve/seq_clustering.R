
if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
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
if (!require("pheatmap")) {
   install.packages("pheatmap", dependencies = TRUE)
   library(pheatmap)
   }
if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2)
   }
if (!require("plyr")) {
   install.packages("plyr", dependencies = TRUE)
   library(plyr)
   }   
if (!require("cluster")) {
   install.packages("cluster", dependencies = TRUE)
   library(cluster)
   }   

dat <- read.table("tmp2.tab", sep=";")
datcast <- acast(dat, V1~V2, value.var="V3")
datdist <- as.dist(as.matrix(datcast))

wss <- (nrow(datcast)-1)*sum(apply(datcast,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(datcast, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

b <- c(0,20,40,60,80,100)
colfunc <- c("white"  ,"#CCEBC5", "#7BCCC4" ,"#2B8CBE","#084081")
heatmap.2(datcast, main="Sequence clustering", dendrogram="none", trace = "none", 
	Rowv=TRUE, Colv=TRUE ,col=colfunc, breaks=b, labCol=FALSE, labRow=FALSE, density.info="none", margins =c(4,10))

fit <- hclust(datdist, method="ward") 
# Method
# the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of 
# "ward", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

plot(fit) # display dendogram
#groups <- cutree(fit, k=4) # cut tree into 5 clusters
rect.hclust(fit, k=4, border="red")

clusplot(as.matrix(datdist), fit$labels, color=TRUE, shade=TRUE, labels=0, lines=0)

fit <- kmeans(datdist, 4)
clusplot(as.matrix(datdist), fit$cluster, color=TRUE, shade=TRUE, labels=4, lines=0)

########################
# Setting clusplot lables
########################
# labels= 0,no labels are placed in the plot;
# labels= 1,points and ellipses can be identified in the plot (see identify);
# labels= 2,all points and ellipses are labelled in the plot;
# labels= 3,only the points are labelled in the plot;
# labels= 4,only the ellipses are labelled in the plot.
# labels= 5,the ellipses are labelled in the plot, and points can be identified.

########################
# Setting number of clusters
########################

#One. Look for a bend or elbow in the sum of squared error (SSE) scree plot. 
#See http://www.statmethods.net/advstats/cluster.html & http://www.mattpeeples.net/kmeans.html for more. 
#The location of the elbow in the resulting plot suggests a suitable number of clusters for the kmeans:

# wss <- (nrow(datcast)-1)*sum(apply(datcast,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(datcast, centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
# We might conclude that 4 clusters would be indicated by this method: 







# Calculate the within groups sum of squared error (SSE) for the number of cluster solutions selected by the user
n.lev <- as.integer(6)
wss <- rnorm(10)
while (prod(wss==sort(wss,decreasing=T))==0) {
wss <- (nrow(datcast)-1)*sum(apply(datcast,2,var))
for (i in 2:n.lev) wss[i] <- sum(kmeans(datcast, centers=i)$withinss)}

# Calculate the within groups SSE for 50 randomized data sets (based on the original input data)
k.rand <- function(x){
km.rand <- matrix(sample(x),dim(x)[1],dim(x)[2])
rand.wss <- as.matrix(dim(x)[1]-1)*sum(apply(km.rand,2,var))
for (i in 2:n.lev) rand.wss[i] <- sum(kmeans(km.rand, centers=i)$withinss)
rand.wss <- as.matrix(rand.wss)
return(rand.wss)}
rand.mat <- matrix(0,n.lev,50)
k.1 <- function(x) { 
for (i in 1:50) {
r.mat <- as.matrix(suppressWarnings(k.rand(datcast)))
rand.mat[,i] <- r.mat}
return(rand.mat)}

# Same function as above for data with < 3 column variables
k.2.rand <- function(x){
rand.mat <- matrix(0,n.lev,50)
km.rand <- matrix(sample(x),dim(x)[1],dim(x)[2])
rand.wss <- as.matrix(dim(x)[1]-1)*sum(apply(km.rand,2,var))
for (i in 2:n.lev) rand.wss[i] <- sum(kmeans(km.rand, centers=i)$withinss)
rand.wss <- as.matrix(rand.wss)
return(rand.wss)}
k.2 <- function(x){
for (i in 1:50) {
r.1 <- k.2.rand(datcast)
rand.mat[,i] <- r.1}
return(rand.mat)}

# Determine if the data data table has > or < 3 variables and call appropriate function above
if (dim(datcast)[2] == 2) { rand.mat <- k.2(datcast) } else { rand.mat <- k.1(datcast) }

# Plot within groups SSE against all tested cluster solutions for actual and randomized data - 1st: Log scale, 2nd: Normal scale
par(ask=TRUE)
xrange <- range(1:n.lev)
yrange <- range(log(rand.mat),log(wss))
plot(xrange,yrange, type='n', xlab='Cluster Solution', ylab='Log of Within Group SSE', main='Cluster Solutions against Log of SSE')
for (i in 1:50) lines(log(rand.mat[,i]),type='l',col='red')
lines(log(wss), type="b", col='blue')
legend('topright',c('Actual Data', '50 Random Runs'), col=c('blue', 'red'), lty=1)
par(ask=TRUE)
yrange <- range(rand.mat,wss)
plot(xrange,yrange, type='n', xlab="Cluster Solution", ylab="Within Groups SSE", main="Cluster Solutions against SSE")
for (i in 1:50) lines(rand.mat[,i],type='l',col='red')
lines(1:n.lev, wss, type="b", col='blue')
legend('topright',c('Actual Data', '50 Random Runs'), col=c('blue', 'red'), lty=1)

# Calculate the mean and standard deviation of difference between SSE of actual data and SSE of 50 randomized datasets
r.sse <- matrix(0,dim(rand.mat)[1],dim(rand.mat)[2])
wss.1 <- as.matrix(wss)
for (i in 1:dim(r.sse)[2]) {
r.temp <- abs(rand.mat[,i]-wss.1[,1])
r.sse[,i] <- r.temp}
r.sse.m <- apply(r.sse,1,mean)
r.sse.sd <- apply(r.sse,1,sd)
r.sse.plus <- r.sse.m + r.sse.sd
r.sse.min <- r.sse.m - r.sse.sd

# Plot differeince between actual SSE mean SSE from 50 randomized datasets - 1st: Log scale, 2nd: Normal scale 
par(ask=TRUE)
xrange <- range(1:n.lev)
yrange <- range(log(r.sse.plus),log(r.sse.min))
plot(xrange,yrange, type='n',xlab='Cluster Solution', ylab='Log of SSE - Random SSE', main='Cluster Solustions against (Log of SSE - Random SSE)')
lines(log(r.sse.m), type="b", col='blue')
lines(log(r.sse.plus), type='l', col='red')
lines(log(r.sse.min), type='l', col='red')
legend('topright',c('SSE - random SSE', 'SD of SSE-random SSE'), col=c('blue', 'red'), lty=1)
par(ask=TRUE)
xrange <- range(1:n.lev)
yrange <- range(r.sse.plus,r.sse.min)
plot(xrange,yrange, type='n',xlab='Cluster Solution', ylab='SSE - Random SSE', main='Cluster Solutions against (SSE - Random SSE)')
lines(r.sse.m, type="b", col='blue')
lines(r.sse.plus, type='l', col='red')
lines(r.sse.min, type='l', col='red')
legend('topright',c('SSE - random SSE', 'SD of SSE-random SSE'), col=c('blue', 'red'), lty=1)


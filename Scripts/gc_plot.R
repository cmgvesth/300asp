require("argparse")
require("ggplot2")
require("reshape")

args <- commandArgs(TRUE)
infile <- args[1]
print(infile)

description <- args[2:5] # Contains section, cutoff1 and cutoff 2
print(description)

threshold = as.numeric(description[4])*0.7
tmp = read.csv(infile)
data <- tmp[order(tmp$X),]

pdf(paste0('/home/seth/Dropbox/seth-1/special_gc/gc_dist_',description[1],'.pdf'))

DF1 <- melt(data, id.var="X")


DF1$n <- ordered(data$X, levels = rev(levels(data$X)))

seth_cols <- c("#08519C","#2171B5","#4292C6","#6BAED6","#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF")
seth_cols2 <- c("#08519C","#4292C6","#9ECAE1", "#F7FCFD", "#CCECE6", "#66C2A4")
tcveCols <- c("#084081","#4eb3d3", "#2b8cbe", "#0868ac","#fc8d59", "#e34a33", "#b30000")
ggplot(DF1, aes(x = n, y = value, fill = variable))+
geom_bar(stat = "identity")+
coord_flip()+
geom_hline(aes(yintercept=threshold), colour="#990000", linetype="dashed")+
scale_fill_manual(values=tcveCols, guide=FALSE) + xlab("Aspergillus")+
ylab("Members in gene cluster(s)")+
ggtitle(paste0("Distribution of gene cluster ",description[1]," from Aspergillus ", description[3], "\nPredicted in ", description[2]," organisms"))+
theme(axis.text = element_text(size = 8, colour="grey20"), plot.title = element_text(size=rel(0.8)))
dev.off()

#x <- barplot(as.matrix(data), xaxt='n' , ylab= "Number of best hits inside one gene cluster", main=paste0('Distribution of ',description[1],' throughout section Nigri') , sub=paste0('Found in ',description[2],' organisms')) 
#text(cex=.7, x=x-1, y=-2, names(data) , xpd=TRUE, srt=45) # for rotated x-labels, put xaxt='n' into plot when in use!

#p2<-ggplot(data,aes(x=names(data)), color=factor(vs)) + stat_summary(fun.y=mean,position="stack",geom="bar")

#ggplot(DF1, aes(x = X, y = value, fill = variable)) + geom_bar(stat = "identity")

#a1 <- ggplot(data) + geom_bar(stat="identity") + coord_flip()
#+ #, color="#084081")+
#coord_flip() +
#scale_x_discrete( limits=asubdat$X.2, labels=asubdat$X.2) +
#ggtitle("Assembly length") +
#geom_text(aes(label=trunc(Assembly.length/1000000)), vjust=+0.5, hjust=+1.2, colour="white") +
#xlab("Organism") +
#ylab("Mega bases")
#grid.arrange(g6, a1, ncol = 2)

# More colours!
# Blues"#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"

# BuGn: "#F7FCFD" "#E5F5F9" "#CCECE6" "#99D8C9" "#66C2A4" "#41AE76" "#238B45" "#006D2C" "#00441B"



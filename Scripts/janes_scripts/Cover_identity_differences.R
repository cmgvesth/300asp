#### IMPORT ####
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("lattice", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("grid", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("gridExtra", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")



####### DATA #########
# A. Carbonaris
> Aspca3_COV80_ID40_family <- read.delim("~/Aspca3_COV80_ID40_family_2015-02-25.txt", header=FALSE)
Aspca3_COV80_ID40_family$species <- "Aspca3"
Aspca3_COV80_ID40_family$coverage <- 80
Aspca3_COV80_ID40_family$identity <- 40
names(Aspca3_COV80_ID40_family)[1]<-paste("family")
names(Aspca3_COV80_ID40_family)[2]<-paste("group.size")

> Aspca3_COV80_ID50_family <- read.delim("~/Aspca3_COV80_ID50_family_2015-02-25.txt", header=FALSE)
Aspca3_COV80_ID50_family$species <- "Aspca3"
Aspca3_COV80_ID50_family$coverage <- 80
Aspca3_COV80_ID50_family$identity <- 50
names(Aspca3_COV80_ID50_family)[1]<-paste("family")
names(Aspca3_COV80_ID50_family)[2]<-paste("group.size")

> Aspca3_COV80_ID60_family <- read.delim("~/Aspca3_COV80_ID60_family_2015-02-25.txt", header=FALSE)
Aspca3_COV80_ID60_family$species <- "Aspca3"
Aspca3_COV80_ID60_family$coverage <- 80
Aspca3_COV80_ID60_family$identity <- 60
names(Aspca3_COV80_ID60_family)[1]<-paste("family")
names(Aspca3_COV80_ID60_family)[2]<-paste("group.size")

> Aspca3_COV90_ID40_family <- read.delim("~/Aspca3_COV90_ID40_family_2015-02-25.txt", header=FALSE)
Aspca3_COV90_ID40_family$species <- "Aspca3"
Aspca3_COV90_ID40_family$coverage <- 90
Aspca3_COV90_ID40_family$identity <- 40
names(Aspca3_COV90_ID40_family)[1]<-paste("family")
names(Aspca3_COV90_ID40_family)[2]<-paste("group.size")

> Aspca3_COV90_ID50_family <- read.delim("~/Aspca3_COV90_ID50_family_2015-02-25.txt", header=FALSE)
Aspca3_COV90_ID50_family$species <- "Aspca3"
Aspca3_COV90_ID50_family$coverage <- 90
Aspca3_COV90_ID50_family$identity <- 50
names(Aspca3_COV90_ID50_family)[1]<-paste("family")
names(Aspca3_COV90_ID50_family)[2]<-paste("group.size")

> Aspca3_COV90_ID60_family <- read.delim("~/Aspca3_COV90_ID60_family_2015-02-25.txt", header=FALSE)
Aspca3_COV90_ID60_family$species <- "Aspca3"
Aspca3_COV90_ID60_family$coverage <- 90
Aspca3_COV90_ID60_family$identity <- 60
names(Aspca3_COV90_ID60_family)[1]<-paste("family")
names(Aspca3_COV90_ID60_family)[2]<-paste("group.size")

# Aspergillus Nidulans
> Aspnid1_COV80_ID40_family <- read.delim("~/Aspnid1_COV80_ID40_family_2015-02-25.txt", header=FALSE)
Aspnid1_COV80_ID40_family$species <- "Aspnid1"
Aspnid1_COV80_ID40_family$coverage <- 80
Aspnid1_COV80_ID40_family$identity <- 40
names(Aspnid1_COV80_ID40_family)[1]<-paste("family")
names(Aspnid1_COV80_ID40_family)[2]<-paste("group.size")

> Aspnid1_COV80_ID50 <- read.delim("~/Aspnid1_COV80_ID50_family_2015-02-25.txt", header=FALSE)
Aspnid1_COV80_ID50$species <- "Aspnid1"
Aspnid1_COV80_ID50$coverage <- 80
Aspnid1_COV80_ID50$identity <- 50
names(Aspnid1_COV80_ID50)[1]<-paste("family")
names(Aspnid1_COV80_ID50)[2]<-paste("group.size")

> Aspnid1_COV80_ID60_family <- read.delim("~/Aspnid1_COV80_ID60_family_2015-02-25.txt", header=FALSE)
Aspnid1_COV80_ID60_family$species <- "Aspnid1"
Aspnid1_COV80_ID60_family$coverage <- 80
Aspnid1_COV80_ID60_family$identity <- 60
names(Aspnid1_COV80_ID60_family)[1]<-paste("family")
names(Aspnid1_COV80_ID60_family)[2]<-paste("group.size")

> Aspnid1_COV90_ID40_family <- read.delim("~/Aspnid1_COV90_ID40_family_2015-02-25.txt", header=FALSE)
Aspnid1_COV90_ID40_family$species <- "Aspnid1"
Aspnid1_COV90_ID40_family$coverage <- 90
Aspnid1_COV90_ID40_family$identity <- 40
names(Aspnid1_COV90_ID40_family)[1]<-paste("family")
names(Aspnid1_COV90_ID40_family)[2]<-paste("group.size")

> Aspnid1_COV90_ID50_family <- read.delim("~/Aspnid1_COV90_ID50_family_2015-02-25.txt", header=FALSE)
Aspnid1_COV90_ID50_family$species <- "Aspnid1"
Aspnid1_COV90_ID50_family$coverage <- 90
Aspnid1_COV90_ID50_family$identity <- 50
names(Aspnid1_COV90_ID50_family)[1]<-paste("family")
names(Aspnid1_COV90_ID50_family)[2]<-paste("group.size")

> Aspnid1_COV90_ID60_family <- read.delim("~/Aspnid1_COV90_ID60_family_2015-02-25.txt", header=FALSE)
Aspnid1_COV90_ID60_family$species <- "Aspnid1"
Aspnid1_COV90_ID60_family$coverage <- 90
Aspnid1_COV90_ID60_family$identity <- 60
names(Aspnid1_COV90_ID60_family)[1]<-paste("family")
names(Aspnid1_COV90_ID60_family)[2]<-paste("group.size")

# Aspergillus niger ATTC 1015
> Aspni7_COV80_ID40_family <- read.delim("~/Aspni7_COV80_ID40_family_2015-02-25.txt", header=FALSE)
Aspni7_COV80_ID40_family$species <- "Aspni7"
Aspni7_COV80_ID40_family$coverage <- 80
Aspni7_COV80_ID40_family$identity <- 40
names(Aspni7_COV80_ID40_family)[1]<-paste("family")
names(Aspni7_COV80_ID40_family)[2]<-paste("group.size")

> Aspni7_COV80_ID50_family <- read.delim("~/Aspni7_COV80_ID50_family_2015-02-25.txt", header=FALSE)
Aspni7_COV80_ID50_family$species <- "Aspni7"
Aspni7_COV80_ID50_family$coverage <- 80
Aspni7_COV80_ID50_family$identity <- 50
names(Aspni7_COV80_ID50_family)[1]<-paste("family")
names(Aspni7_COV80_ID50_family)[2]<-paste("group.size")

> Aspni7_COV80_ID60_family <- read.delim("~/Aspni7_COV80_ID60_family_2015-02-25.txt", header=FALSE)
Aspni7_COV80_ID60_family$species <- "Aspni7"
Aspni7_COV80_ID60_family$coverage <- 80
Aspni7_COV80_ID60_family$identity <- 60
names(Aspni7_COV80_ID60_family)[1]<-paste("family")
names(Aspni7_COV80_ID60_family)[2]<-paste("group.size")

> Aspni7_COV90_ID40_family <- read.delim("~/Aspni7_COV90_ID40_family_2015-02-25.txt", header=FALSE)
Aspni7_COV90_ID40_family$species <- "Aspni7"
Aspni7_COV90_ID40_family$coverage <- 90
Aspni7_COV90_ID40_family$identity <- 40
names(Aspni7_COV90_ID40_family)[1]<-paste("family")
names(Aspni7_COV90_ID40_family)[2]<-paste("group.size")

> Aspni7_COV90_ID50_family <- read.delim("~/Aspni7_COV90_ID50_family_2015-02-25.txt", header=FALSE)
Aspni7_COV90_ID50_family$species <- "Aspni7"
Aspni7_COV90_ID50_family$coverage <- 90
Aspni7_COV90_ID50_family$identity <- 50
names(Aspni7_COV90_ID50_family)[1]<-paste("family")
names(Aspni7_COV90_ID50_family)[2]<-paste("group.size")

> Aspni7_COV90_ID60_family <- read.delim("~/Aspni7_COV90_ID60_family_2015-02-25.txt", header=FALSE)
Aspni7_COV90_ID60_family$species <- "Aspni7"
Aspni7_COV90_ID60_family$coverage <- 90
Aspni7_COV90_ID60_family$identity <- 60
names(Aspni7_COV90_ID60_family)[1]<-paste("family")
names(Aspni7_COV90_ID60_family)[2]<-paste("group.size")



#### COMBIND DATA INTO ONE DATA FRAME ####
total <- rbind(Aspca3_COV80_ID40_family, Aspca3_COV80_ID50_family, Aspca3_COV80_ID60_family, 
               Aspca3_COV90_ID40_family,Aspca3_COV90_ID50_family, Aspca3_COV90_ID60_family, 
               Aspnid1_COV80_ID40_family, Aspnid1_COV80_ID50, Aspnid1_COV80_ID60_family,
               Aspnid1_COV90_ID40_family, Aspnid1_COV90_ID50_family, Aspnid1_COV90_ID60_family,
               Aspni7_COV80_ID40_family,Aspni7_COV80_ID50_family,Aspni7_COV80_ID60_family,
               Aspni7_COV90_ID40_family,Aspni7_COV90_ID50_family,Aspni7_COV90_ID60_family)
View(total)


### PLOT BARS ###
ggplot(total, aes(x=group.size, ..count.., fill = species)) + 
    geom_bar(stat="bin", position = "dodge") + 
    geom_bar(colour="black", show_guide=FALSE,stat="bin", position = "dodge") +
    facet_grid(identity ~ coverage) +
    geom_text(stat="bin", aes(y = ..count.., label=..count..),
              position = position_dodge(width=1),vjust=-0.3) +
    scale_x_discrete(limits=c(1:max(total$group.size, na.rm = TRUE)))+  ylab("# of families") + 
    xlab("Family size") 
    #+ ggtitle("Aspergillus genus") +
    #theme(plot.title=element_text(face="italic", size=20))

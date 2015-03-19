# Modules
library("plyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")

# Import data
t_complete_clustering <- read.csv("~/Desktop/t_complete_clustering.csv", sep=";")
organism <- read.csv("~/Desktop/organism.csv", sep=";")

# Count the group sizes for each organism
count_table <- count(t_complete_clustering,c("org_key","family_key"))
colnames(count_table)[3] <- "family_members"
count_table <- count(count_table, c("org_key","family_members"))
View(count_table)

# Change org_key to real names
from <- as.vector(organism$name)
to   <- as.vector(organism$real_name)

find.in.list <- function(x, y) match(TRUE, sapply(y, `%in%`, x = x))
idx.in.list  <- sapply(levels(count_table$org_key), find.in.list, from)
levels(count_table$org_key)  <- ifelse(is.na(idx.in.list), levels(count_table$org_key), to[idx.in.list])

# Sorting the x-axis (the names of the organism)
myList <- sort(levels(count_table$org_key))

# Subsetting and plotting the gene family members and frequency
subdat <- subset(count_table, count_table$family_members==1)
ggplot(subdat, aes(x = org_key, y=freq)) + 
    geom_bar(stat="identity", fill="#2B8CBE" ) + 
    coord_flip() +
    scale_x_discrete(limits=rev(myList)) +
    theme(axis.text = element_text(size = 12, colour = "grey20")) +
    ggtitle("Singleton counts") + 
    theme(plot.title = element_text(size = rel(2))) +
    geom_text(aes(label=trunc(freq)), vjust=+0.5, hjust=+1.2, colour="white") + 
    xlab("Organism") + 
    ylab("Count")


subdat <- subset(count_table, count_table$family_members==2 | count_table$family_members==3 )
ggplot(subdat, aes(x = org_key, y=freq, fill = as.factor(family_members))) + 
    guides(fill=guide_legend(title="Members")) +
    geom_bar(stat="identity",  position = "dodge") + 
    coord_flip() +
    scale_x_discrete(limits=rev(myList)) +
    scale_y_continuous(limits = c(0, 1000)) +
    theme(axis.text = element_text(size = 12, colour = "grey20")) +
    ggtitle("Gene family counts [2-3 members]") + 
    theme(plot.title = element_text(size = rel(2))) +
    geom_text(aes(y = freq, label=trunc(freq)), colour="black", position = position_dodge(width=1), hjust = -0.09) + 
    xlab("Organism") + 
    ylab("Count")

subdat <- subset(count_table, count_table$family_members==4 | count_table$family_members==5 )
ggplot(subdat, aes(x = org_key, y=freq, fill = as.factor(family_members))) + 
    guides(fill=guide_legend(title="Members")) +
    geom_bar(stat="identity",  position = "dodge") + 
    coord_flip() +
    scale_x_discrete(limits=rev(myList)) +
    scale_y_continuous(limits = c(0, 30)) +
    theme(axis.text = element_text(size = 12, colour = "grey20")) +
    ggtitle("Gene family counts [4-5 members]") + 
    theme(plot.title = element_text(size = rel(2))) +
    geom_text(aes(y = freq, label=trunc(freq)), colour="black", position = position_dodge(width=1), hjust = -1) + 
    xlab("Organism") + 
    ylab("Count")

subdat <- subset(count_table, count_table$family_members==6 | count_table$family_members==7 | count_table$family_members==8)
ggplot(subdat, aes(x = org_key, y=freq, fill = as.factor(family_members))) + 
    guides(fill=guide_legend(title="Members")) +
    geom_bar(stat="identity",  position = "dodge") + 
    coord_flip() +
    scale_x_discrete(limits=rev(myList)) +
    scale_y_continuous(limits = c(0, 8)) +
    theme(axis.text = element_text(size = 12, colour = "grey20")) +
    ggtitle("Gene family counts [6-8 members]") + 
    theme(plot.title = element_text(size = rel(2))) +
    geom_text(aes(y = freq, label=trunc(freq)), colour="black", position = position_dodge(width=1), hjust = -1) + 
    xlab("Organism") + 
    ylab("Count")


subdat <- subset(count_table, count_table$family_members==9 | count_table$family_members==10 | count_table$family_members==11 
                 | count_table$family_members==12 | count_table$family_members==13 | count_table$family_members==14
                 | count_table$family_members==15 | count_table$family_members==16 | count_table$family_members==17
                 | count_table$family_members==18)
ggplot(subdat, aes(x = org_key, y=freq, fill = as.factor(family_members))) + 
    guides(fill=guide_legend(title="Members")) +
    geom_bar(stat="identity",  position = "dodge") + 
    coord_flip() +
    scale_x_discrete(limits=rev(myList)) +
    scale_y_continuous(limits = c(0, 5)) +
    theme(axis.text = element_text(size = 12, colour = "grey20")) +
    ggtitle("Gene family counts [9-18 members]") + 
    theme(plot.title = element_text(size = rel(2))) +
    geom_text(aes(y = freq, label=trunc(freq)), colour="black", position = position_dodge(width=1), hjust = -1) + 
    xlab("Organism") + 
    ylab("Count")

subdat <- subset(count_table, count_table$family_members==20 | count_table$family_members==24 | count_table$family_members==29 
                 | count_table$family_members==32)
ggplot(subdat, aes(x = org_key, y=freq, fill = as.factor(family_members))) + 
    guides(fill=guide_legend(title="Members")) +
    geom_bar(stat="identity",  position = "dodge") + 
    coord_flip() +
    scale_x_discrete(limits=rev(myList)) +
    scale_y_continuous(limits = c(0, 3)) +
    theme(axis.text = element_text(size = 12, colour = "grey20")) +
    ggtitle("Gene family counts [20-32 members]") + 
    theme(plot.title = element_text(size = rel(2))) +
    geom_text(aes(y = freq, label=trunc(freq)), colour="black", position = position_dodge(width=1), hjust = -1) + 
    xlab("Organism") + 
    ylab("Count")


# Extract rows containing a specific value
look.for <- c(18,20,24,29,32)
look.for.18.to.32 <- t_complete_clustering[t_complete_clustering$group_size %in% look.for, ]

# Order by org_id and then family_members
look.for.18.to.32 <- look.for.18.to.32[with(look.for.18.to.32, order(org_id, group_size)), ]
View(look.for.18.to.32)

look.for <- c(24)
look.for.24 <- t_complete_clustering[t_complete_clustering$group_size %in% look.for, ]
View(look.for.24)

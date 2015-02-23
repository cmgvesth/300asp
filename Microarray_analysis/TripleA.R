###########################
## R-script file
## Automated Array Analysis (TripleA) for MicroArray data
## Developed for R 2.5.1 or above
###########################

## Initial loadings
library(affy); # Used for loading and normalization of the data
library(limma); # Used for statistical analysis
library(makecdfenv)
NigerDSMcdf <- make.cdf.env("dsmM_ANIGERa_anColl.CDF")


Exp1Data <- ReadAffy(cdfname="NigerDSMcdf")


## Show summary of data for verification
show(Exp1Data)


## 2: Make boxplots of the variations prior to normalization
X11(); par(mfrow = c(1,1)); boxplot(Exp1Data)

## 3: Plots of intensities prior to normalization
X11(); hist(Exp1Data)

## Running the preprocessing
ExpSet <- expresso(Exp1Data, bgcorrect.method="rma", normalize.method="quantiles", pmcorrect.method="pmonly",summary.method="medianpolish")


# QUALITY-check: Show values before and after:
pdf(file="./NormalizationQC.pdf", width=40, height=15)
par(mfrow = c(1,2));
boxplot(Exp1Data);
boxplot(as.data.frame(exprs(ExpSet)));
dev.off()

## Dumping a comma-separated file to Excel for those so inclined
write.exprs(ExpSet,file="NormDataV2.csv", sep=",")

###########################
## Getting a list of genes for further analysis:
###########################

## Extracting list of specific genes
probelist <- read.table("UniqueGeneList.txt")
index <- pmatch(as.character(probelist[,1]), rownames(exprs(ExpSet)))
SubsetOnly <- exprs(ExpSet)[index,]


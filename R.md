
#R SCRIPT 
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("edgeR", "Homo.sapiens"))

# Set working directory & load packages 
#working directory (all post QC files should be here) 
setwd("U:/Research/Projects/ihbi/grc/grc_general/Stem Cell Group/Sofia/Transcriptomics")
#create folders for R QC and output per below 
dir.create("3-QC")
dir.create("4-Output")
options(digits=3)
getwd()

#load libraries 
library(tximport)
library(readr)

# tximport of pre-constructed tx2gene table - created from bioMart Ensembl site. Ensure this is in folder, otherwise will not be run 
tx2gene <- read.csv("2-Input/Homo_sapiens.GRCh38.91_tx2gene.csv")
head(tx2gene, 5)

# create links to salmon quant folders to import salmon files. "<-" defines the folders/directories that we want to look at --- file format is *.csv
folder <- c("2-Input/cell_quants")
salmon.dir <- as.matrix(read.csv(file="2-Input/quant_filenames.csv", sep=",", header=F))

salmon.files <- file.path(folder, salmon.dir, "quant.sf")
names(salmon.files) <- as.matrix(read.csv(file="2-Input/names.csv", sep=",", header=F))
all(file.exists(salmon.files))

# tximport
# countsFromAbundance default is "no" just returns counts not TPM (transcripts-per-million)
# scaledTPM is TPM scaled up to library sizewhile lengthScaledTPM
# lengthScaledTPM first multiplies TPM by feature length and then scales up to library size
# Both are then quantities that are on the same scale as original counts,
# except no longer correlated with feature length across samples.

#gene level summary and salmon index building .. quantification file built using salmon as above, this generates an offset matrix for downstream gene-level differential analysis of count matrices ---- a matrix is a collection of elements of the same data type (numeric, character, or logical) arranged into a fixed number of rows and columns
# lengthScaledTPM first multiplies TPM by feature length and then scales up to library size
# tximport is able to import counts produced by different software, and different workflows are described for each in the tximport vignette.
txi <- tximport(salmon.files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = TRUE, countsFromAbundance="lengthScaledTPM") 
names(txi)
head(txi$counts, 3)
#write.csv = save summaries of partitioned breeding values to CSV files on disk for further analyses of processing with other software or just for saving (backing up) results
write.csv(txi, file = "4-Output/geneSMART_tximport_matrix.csv")

# prepare data
#start by loading libraries. #Use edgeR package to import, organise, filter and normalise the data
library(edgeR)
library(Homo.sapiens)
library(limma)

# edgeR used to read in count data (DGEList -- Creates a DGEList object from a table of counts (rows=features, columns=samples), group indicator for each column, library size (optional) and a table of feature annotation (optional); dim -- Retrieve or set the dimension of an object.)
y <- DGEList(txi$counts)
dim(y)

# Read in sample information table
csvfile <- file.path("2-Input/Cell_sample_table.csv")
sampleTable <- read.csv(csvfile, row.names=1)
y$samples$names <- sampleTable$name
y$samples$filename <- sampleTable$filename
y$samples$Sample_ID <- sampleTable$Sample_ID
y$samples$cell_line <- sampleTable$cell_line
y$samples$Passage <- sampleTable$Passage
y$samples$Day <- sampleTable$Day
y$samples$Treatment <- sampleTable$Treatment
y$samples$run_date <- sampleTable$run_date
y$samples$day <- sampleTable$Day
y$samples$id <- sampleTable$Sample_ID
y$samples$run_date <- sampleTable$run_date
y$samples$cell_line <- sampleTable$cell_line
y$samples

#gives you all the column names for all the samples
colnames1 <- colnames(y$samples)
colnames1

# Add gene annotation
geneid <- rownames(y)
columns(Homo.sapiens)
genes <- select(Homo.sapiens, key=geneid, keytype="ENSEMBL",
                columns=c("SYMBOL","GENEID","GENENAME","TXCHROM"),
                multiVals="first")

# Add gene annotation
library("biomaRt")
geneid <- rownames(y)
ensembl91 <- useMart(host="dec2017.archive.ensembl.org", 
                     biomart="ENSEMBL_MART_ENSEMBL", 
                     dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl91)
attributes[1:20,]

genes <- select(ensembl91, keys=geneid, keytype="ensembl_gene_id",
                columns=c("ensembl_gene_id", "external_gene_name",
                          "description","entrezgene","chromosome_name","gene_biotype"))

(colnames(genes) <- c("ENSEMBL","SYMBOL","GENENAME","GENEID","TXCHROM","BIOTYPE"))

#remove duplicated genes
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes
head(genes,5)
dim(genes)

#This is where the actual analysis begins, above is QC/cleanup 
#change point
# Barplot of library sizes -- this will create a barplot of library sizes (i.e. how big the reads are)
png("3-QC/Barplot of library sizes.png", width = 90, height = 30, units = 'cm', res = 300)
col <- as.numeric(y$sample$timepoint)
dt <- colSums((y$counts)*1e-6)
barplot(dt, names=colnames(dt), col=c("lightsalmon","lightcoral","steelblue1","cornsilk")[col], las=2, cex.names=0.8)
abline(h=5,col="black",lty=3)
abline(h=10,col="black",lty=3)
abline(h=15,col="black",lty=3)
title(main="Barplot of library sizes",ylab="Library size (millions)")
dev.off()

# Summary of total mapped read counts
summary(dt)

# Add all reads
CountMeans <- rowMeans(y$counts)
TotalCounts <- sum(CountMeans)

# Average count for each gene in all samples, sorted from highest to lowest
PercentReads <- (CountMeans/TotalCounts)*100
AvgCounts <- rowMeans(y$counts)
Symbol <- y$genes$SYMBOL
GeneName <- y$genes$GENENAME
AvgCountsTable <- (data.frame(Symbol,GeneName,AvgCounts,PercentReads))
AvgCountsTable <- AvgCountsTable[order(-AvgCountsTable$AvgCounts),]
head(AvgCountsTable,20)
write.csv(AvgCountsTable, file ="4-Output/Average_Counts.csv")

# Top 5, 10, 25 transcripts
top5 <- ((sum(head(AvgCountsTable, 5)$AvgCounts))/TotalCounts)*100
top10 <- ((sum(head(AvgCountsTable, 10)$AvgCounts))/TotalCounts)*100
top25 <- ((sum(head(AvgCountsTable, 25)$AvgCounts))/TotalCounts)*100

top5; top10; top25

# Summary of the lengthscaledTPM count data
dim(y)
summary(rowMeans(y$counts))

# Histogram of count distribution - this will create a plot of read count distributions BEFORE FILTERING
# pseudocount of 0.25 added to values to prevent logging zero values
# therefore O CPM = lcpm -6.9 --- counts per million reads mapped (CPM) = the count of sequenced fragments mapping to the feature 
lcpm.AvgCounts <- cpm(AvgCounts, log = TRUE)
png("3-QC/Read Count Distribution-B4.png", width = 45, height = 25, units = 'cm', res = 600)
par(mar=c(5,6,4,1)+.1)
hist(lcpm.AvgCounts, col="salmon", border="salmon",
     cex.lab=2.5, cex.main=3, cex.axis=2,
     xlab="Median log2-CPM", ylab="No. of Transcripts",
     breaks=100, xlim=c(-10,20), main ="Read Count Distribution (before)")
dev.off()

# number of genes with zero counts across all 211 samples -- this will spit out the info but not save it 
table(AvgCounts>=1)
table(AvgCounts>=10)
table(AvgCounts>=100)
table(AvgCounts>=1000)
table(AvgCounts>=10000)
table(AvgCounts>=100000)

# Filter unexpressed and very low expressed genes
# Count genes with zero counts across all 211 samples
table(rowSums(y$counts==0)==211)

# FILTERING
# Based on median log2 CPM
# Expect approximately 9K-12K genes to remain after filtering
# https://www.biostars.org/p/211954/
median_cpm <- apply(cpm(y), 1, median)
expr_cutoff <- 0.5 # in cpm
sum(median_cpm > expr_cutoff)
y.Filt <- y[median_cpm > expr_cutoff, ]
lcpm.Filt <- cpm(y.Filt$counts, log = TRUE)

log.cutoff <- log2(expr_cutoff)

# Summary of filtered count data
summary(rowMeans(y.Filt$counts))

# Density of count values (code modified from 'RNAseq 1-2-3') - this will create graphs of raw and filtered data -- density vs log of CPM
png("3-QC/Density of count values.png", width = 10, height = 20, units = 'cm', res = 600)
nsamples <- ncol(y)
col <- rainbow(nsamples)
par(mfrow=c(2,1))
lcpm.Raw <- cpm(y$counts, log=TRUE)
plot(density(lcpm.Raw[,1]), col=col[1], xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="")
for (i in 2:nsamples){
  den <-density(lcpm.Raw[,i])
  lines(den$x, den$y, col=col[i])
}
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="")
title("Raw data",xlab="log2-CPM")

lcpm.Filt <- cpm(y.Filt$counts, log=TRUE)
plot(density(lcpm.Filt[,1]), col=col[1], xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="")
for (i in 2:nsamples){
  den <-density(lcpm.Filt[,i])
  lines(den$x, den$y, col=col[i])
}
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="")
title("Filtered data (median CPM > 0.5)",xlab="log2-CPM")
dev.off()


# Histogram of count distribution  - this will create a plot of read count distributions AFTER FILTERING
# pseudocount of 0.25 added to values to prevent logging zero values
# therefore O CPM = lcpm -6.9
AvgCounts <- rowMeans(y.Filt$counts)
lcpm.AvgCounts <- cpm(AvgCounts, log = TRUE)
png("3-QC/Read Count Distribution-AFTER.png", width = 45, height = 25, units = 'cm', res = 600)
par(mar=c(5,6,4,1)+.1)
hist(lcpm.AvgCounts, col="salmon", border="salmon",
     cex.lab=2.5, cex.main=3, cex.axis=2,
     xlab="Median log2-CPM", ylab="No. of Transcripts",
     breaks=100, xlim=c(-10,20), main ="Read Count Distribution (after)")
dev.off()


# Heatmap of samples -- this is useful for visualising the expression of genes across the MSC samples
png("3-QC/Sample heatmap.png", width = 60, height = 60, units = 'cm', res = 300)
par(mfrow=c(1,2))
lcpm.Raw <- cpm(y$counts, log = TRUE)
heatmap(cor(lcpm.Raw))
title("Raw data")
lcpm.Filt <- cpm(y.Filt$counts, log = TRUE)
heatmap(cor(lcpm.Filt))
title("Filtered data")
dev.off() 


#################
##  limma --- provides an integrated solution for analysing data from gene expression experiments.
#load limma package
#################
library(limma)

# Normalisation (to eliminate systematic experimental bias and technical variation while preserving biological variation)
y.Norm <- calcNormFactors(y.Filt, method="TMM")
y.Norm$samples$norm.factors
summary(y.Norm$samples$norm.factors)

#boxplot BEFORE normalisation
png("3-QC/BoxplotBefore.png", width = 90, height = 45, units = 'cm', res = 300)
par(mfrow=c(2,1))
lcpm.Filt <- cpm(y.Filt, log=TRUE)
col <- as.numeric(y$sample$tDay)
boxplot(lcpm.Filt,las=2,
        col=c("lightsalmon","lightcoral","steelblue1","cornsilk")[col],main="")
abline(h=median(lcpm.Filt),col="black",lty=3)
title(main="Unnormalised data",ylab="log-counts")

#boxplot AFTER normalisation (limma gains statistical power by modelling the variances AKA normalising data)
png("3-QC/BoxplotAfter.png", width = 90, height = 45, units = 'cm', res = 300)
par(mfrow=c(2,1))
lcpm.Norm <- cpm(y.Norm, log=TRUE)
boxplot(lcpm.Norm,las=2,
        col=c("lightsalmon","lightcoral","steelblue1","cornsilk")[col],main="")
abline(h=median(lcpm.Norm),col="black",lty=3)
title(main="TMM-Normalised data",ylab="log-counts")
dev.off()



# UNSUPERVISED CLUSTERING OF SAMPLES
# MDS plots - MDS arranges the points on the plot so that the distances among each pair of points correlates as best as possible to the dissimilarity between those two samples. The values on the two axes tell you nothing about the variables for a given sample - the plot is just a two dimensional space to arrange the points.

dir.create("5-Glimma")
library(Glimma)
# Glimma MDS Plot --- Draws an interactive MD plot from a DGEList object with distances calculated from most variable genes.
glMDSPlot(y.Norm, top=500, labels=y.Norm$samples$subject,
          groups=y.Norm$samples[,c(1:9)], launch=TRUE,
          path="5-Glimma", folder="glMDSPlot", html="MDS_plot")


##### MODIFY CODE ######

# MDS plot of sample names
png("3-QC/MDS Plot_names.png", width = 15, height = 30, units = 'cm', res = 600)
par(mfrow=c(2,1))
plotMDS(y.Norm, top=500, cex=0.8, labels=y.Norm$samples$cell_line, col=as.numeric(y.Norm$samples$Treatment), main="MDS Plot (SampleID)")
legend("topright", legend=c("hMSC-20176", "hMSC-21558"), cex=0.8, col=1:16, pch=16)

plotMDS(y.Norm, top=500, cex=1, pch=21, col="white", bg=as.numeric(y.Norm$samples$Treatment), main="MDS Plot (Symbol)")
pch=x; 0= open square; 1= open circle; 15= solid square, 16=solid circle, 21=filled circle/border
legend("topright", legend=c("hMSC-20176", "hMSC-21558"), cex=0.8, col=1:16, pch=16)
dev.off()




# Mean Difference Plots - XY scatter plot that compares the disagreement, or differences, between two quantitative measurements
# Library size-adjusted log-fold change between two libraries (the difference)
# against the average log-expression across those libraries (the mean).
# The following command produces an MD plot that compares

# MD plot before and after normalisation
png("3-QC/MDPlots-1LI00.png", width=11, height=20, units='cm', res=200)
par(mfrow=c(2,1))
plotMD(y.Filt, column=1, main="First sample (raw)"); abline(h=0, col="red", lty=2, lwd=2)
plotMD(y.Norm, column=1, main="First sample (TMM-normalised)"); abline(h=0, col="red", lty=2, lwd=2)
dev.off()

# Limma design matrix
data <- y$samples
# create a design matrix of no intercept, only treated vs untreated
design2 <- model.matrix(~0+data$Treatment)
# set the names of columns in the design matrix accordingly
colnames(design2) <- c("Treated", "Untreated")
# create a contrast matrix comparing the treated vs untreated groups in this design matrix
contrast.matrix <- makeContrasts("Treated-Untreated", levels=design2)
# do some magic to grab the data from the RNA sequences corresponding to the regression we want to run so that it's ready for the linear model
v2 <- voom(y.Norm, design2)
# fit a linear model using the design matrix only
fit2 <- lmFit(v2, design2)
# fit the contrasts to get info between groups
fit2c <- contrasts.fit(fit2, contrast.matrix)
# compute statistics of the linear model using empirical bayes
fit2c <- eBayes(fit2c)
# print table of top predictors
topTable(fit2c)



fit2 <- lmFit(data, design2)

group <- y$samples$Day
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

# Paired design using DuplicateCorrelation
#https://support.bioconductor.org/p/52920/  https://support.bioconductor.org/p/59700/
Group <- y.Norm$samples$Sample_ID
design <- model.matrix(~0 + Group)

v <- voom(y.Norm, design)
colnames(design) <- c("hMSC20176P13D3untreated", "hMSC20176P13D3treated", "hMSC20176P13D5untreated", "hMSC20176P13D5treated", "hMSC20176P5D3untreated", "hMSC20176P5D3treated", "hMSC20176P5D5untreated", "hMSC20176P5D5treated", "hMSC20176P7D3untreated", "hMSC20176P7D3treated", "hMSC20176P7D5untreated", "hMSC20176P7D5treated", "hMSC21558P13D3untreated", "hMSC21558P13D3treated", "hMSC21558P13D5untreated", "hMSC21558P13D5treated", "hMSC21558P5D3untreated", "hMSC21558P5D3treated", "hMSC21558P5D5untreated", "hMSC21558P5D5treated", "hMSC21558P7D3untreated", "hMSC21558P7D3treated", "hMSC21558P7D5untreated", "hMSC21558P7D5treated")
corfit <- duplicateCorrelation(v, design, block=y.Norm$samples$cell_line)
v <- voom(y.Norm, design, block=y.Norm$samples$Sample_ID, correlation=corfit$consensus)

save(v, file="4-Output/v.rda")


# Use the SVD functon in ChAMP to assess which covariates correlate with the top components
# Covariate === an independent variable that can influence the outcome of a given statistical trial, but which is not of direct interest.
library(ChAMP)

#THE ISSUE IS HERE!!!
pd <- v$targets[,c(colnames1)]
colnames(pd)[1] <- y$samples$filename
champ.SVD(beta=v$E, pd=pd, PDFplot=TRUE, Rplot=FALSE, resultsDir="./3-QC/")

glMDSPlot(v, top=500, labels=v$targets$cell_line,
          groups=v$targets[,c(1:9)], launch=TRUE,
          path="5-Glimma", folder="glMDSPlot", html="MDS_plot_voom")

# MD plot before and after normalisation
png("3-QC/MDPlots-sample54.png", width=20, height=40, units='cm', res=300)
par(mfrow=c(3,1))
plotMD(y.Filt, column=54, main="sample 54 (raw)"); abline(h=0, col="red", lty=2, lwd=2)
plotMD(y.Norm, column=54, main="sample 54 (TMM-normalised)"); abline(h=0, col="red", lty=2, lwd=2)
plotMD(v, column=54, main="sample 54 (voom)"); abline(h=0, col="red", lty=2, lwd=2)
dev.off()


# Linear model fitting (eBayes or TREAT)
#lfc <- log2(1.1) # use for TREAT only

#include group before??????
# Design matrix=	Used to define the form of a statistical model and to store observed values of the explanatory variable(s). Used in the computation process to estimate model parameters. 
#make a model matrix instead of whatever is happening here ????? 
#Contrast matrix=	Used in conjunction with a design matrix to calculate specific values of interest between estimated parameters.
#Levels	Unique = values within a factor, e.g. wildtype or mutant.
fit <- lmFit(v, design, block=y.Norm$samples$cell_line, correlation=corfit$consensus)
cont.matrix = makeContrasts("20176_P5_D3_treated-20176_P13_D3_treated" = hMSC20176P5D3treated - hMSC20176P13D3treated, "20176_P5_D3_untreated-20176_P13_D3_untreated" - hMSC20176P5D3untreated - hMSC20176P13D3untreated, "hMSC_P5_D3_treated-hMSC_P13_D3_treated" = (hMSC20176P5D3treated - hMSC21558P5D3treated) - (hMSC20176P13D3treated - hMSC2155813D3treated), "20176_P5_D3_untreated-treated" = hMSC20176P5D3treated - hMSC20176P5D3untreated, levels=design)

cont.matrix = makeContrasts("20176_P13_D3_treated-Untreated" = hMSC20176P13D3treated - hMSC20176P13D3untreated, levels=design)

("20176_P13_D3_treated-Untreated" = hMSC20176P13D3treated - hMSC20176P13D3untreated,
  "20176_P13_D5_untreated-treated" = hMSC20176P13D5treated - hMSC20176P13D5untreated,
  "20176_P13_D5-D3treated" = (hMSC20176P13D5treated - hMSC20176P13D5untreated) - (hMSC20176P13D3treated - hMSC20176P13D3untreated),
  "20176_P5_D3_untreated-treated" = hMSC20176P5D3treated - hMSC20176P5D3untreated,
  "20176_P5_D5_untreated-treated" = hMSC20176P5D5treated - hMSC20176P5D5untreated,
  "20176_P5_D5-D3treated" = (hMSC20176P5D5treated - hMSC20176P5D5untreated) - (hMSC20176P5D3treated - hMSC20176P5D3untreated),
  "20176_P7_D3_untreated-treated" = hMSC20176P7D3treated - hMSC20176P7D3untreated,
  "20176_P7_D5_untreated-treated" = hMSC20176P7D5treated - hMSC20176P7D5untreated,
  "20176_P7_D5-D3treated" = (hMSC20176P7D5treated - hMSC20176P7D5untreated) - (hMSC20176P7D3treated - hMSC20176P7D3untreated),
  "20176_P13-P5treated" = ((hMSC20176P13D5treated - hMSC20176P13D5untreated) - (hMSC20176P13D3treated - hMSC20176P13D3untreated)) - ((hMSC20176P5D5treated - hMSC20176P5D5untreated) - (hMSC20176P5D3treated - hMSC20176P5D3untreated)),
  "20176_P7-P5treated" = ((hMSC20176P7D5treated - hMSC20176P7D5untreated) - (hMSC20176P7D3treated - hMSC20176P7D3untreated)) - ((hMSC20176P5D5treated - hMSC20176P5D5untreated) - (hMSC20176P5D3treated - hMSC20176P5D3untreated)),
  "21558_P13_D3_treated-Untreated" = hMSC21558P13D3treated - hMSC21558P13D3untreated,
  "21558_P13_D5_untreated-treated" = hMSC21558P13D5treated - hMSC21558P13D5untreated,
  "21558_P13_D5-D3treated" = (hMSC21558P13D5treated - hMSC21558P13D5untreated) - (hMSC21558P13D3treated - hMSC21558P13D3untreated),
  "21558_P5_D3_untreated-treated" = hMSC21558P5D3treated - hMSC21558P5D3untreated,
  "21558_P5_D5_untreated-treated" = hMSC21558P5D5treated - hMSC21558P5D5untreated,
  "21558_P5_D5-D3treated" = (hMSC21558P5D5treated - hMSC21558P5D5untreated) - (hMSC21558P5D3treated - hMSC21558P5D3untreated),
  "21558_P7_D3_untreated-treated" = hMSC21558P7D3treated - hMSC21558P7D3untreated,
  "21558_P7_D5_untreated-treated" = hMSC21558P7D5treated - hMSC21558P7D5untreated,
  "21558_P7_D5-D3treated" = (hMSC21558P7D5treated - hMSC21558P7D5untreated) - (hMSC21558P7D3treated - hMSC21558P7D3untreated),
  "21558_P13-P5treated" = ((hMSC21558P13D5treated - hMSC21558P13D5untreated) - (hMSC21558P13D3treated - hMSC21558P13D3untreated)) - ((hMSC21558P5D5treated - hMSC21558P5D5untreated) - (hMSC21558P5D3treated - hMSC21558P5D3untreated)),
  "21558_P7-P5treated" = ((hMSC21558P7D5treated - hMSC21558P7D5untreated) - (hMSC21558P7D3treated - hMSC21558P7D3untreated)) - ((hMSC21558P5D5treated - hMSC21558P5D5untreated) - (hMSC21558P5D3treated - hMSC21558P5D3untreated)),
  "20176_P13-21558_P13treated" = ((hMSC20176P13D5treated - hMSC20176P13D5untreated) - (hMSC20176P13D3treated - hMSC20176P13D3untreated)) - ((hMSC21558P13D5treated - hMSC21558P13D5untreated) - (hMSC21558P13D3treated - hMSC21558P13D3untreated)),
  "20176_P7-21558_P7treated" = ((hMSC20176P7D5treated - hMSC20176P7D5untreated) - (hMSC20176P7D3treated - hMSC20176P7D3untreated)) - ((hMSC21558P7D5treated - hMSC21558P7D5untreated) - (hMSC21558P7D3treated - hMSC21558P7D3untreated)),
  "20176_P5-21558_P5treated" = ((hMSC20176P5D5treated - hMSC20176P5D5untreated) - (hMSC20176P5D3treated - hMSC20176P5D3untreated)) - ((hMSC21558P5D5treated - hMSC21558P5D5untreated) - (hMSC21558P5D3treated - hMSC21558P5D3untreated)),
  levels=design)
# check the matrix 
cont.matrix

#contrast.fit == Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
fit2 = contrasts.fit(fit, cont.matrix)
#non-Bayesian analysis??
#try different matrices, see if that makes a differenve ---- Try a different analysis style too. Lol 
fit2 = eBayes(fit2)
#fit2 <- treat(fit2, lfc=lfc)

# SA plots
png("3-QC/SA Plot.png", width = 30, height = 15, units = 'cm', res = 600)
par(mfrow=c(1,2))
v <- voom(y.Norm, design, plot=TRUE)
plotSA(fit2, main="Final model: Mean???variance trend")
dev.off()


# Differential Expression 
# p value and method in here --- p<0.05 , log fold change min 0 
results <- decideTests(fit2)
summary(results)
vennCounts(results)
# Total number of significant genes
ngenes <-  which(results[,1] != 0 | results[,2] != 0 | results[,3] != 0)
ngenes <- length(ngenes)

# Venn Diagram
png("4-Output/Venn Diagram.png", width = 30, height = 15, units = 'cm', res = 600)
par(mfrow=c(1,2))
vennDiagram(results,include=c("both"), circle.col=c("blue","yellow","green","orange"), counts.col=c("blue3"), cex=c(1,0.8,0.8))
vennDiagram(results,include=c("up","down"), circle.col=c("blue","yellow","green","orange"), counts.col=c("red","green3"), cex=c(1,0.7,0.7))
mtext("Cell hMSC (eBayes)", side = 3, line = -2, outer = TRUE, cex=1.8)
mtext(paste0("(", ngenes," genes)"), side = 3, line = -3.5, outer = TRUE, cex=1.5)
dev.off()

# MD Plots
png("4-Output/MDplot-Average.png", width = 20, height = 30, units = 'cm', res = 600)
par(mfrow=c(3,1))
plotMD(fit2, coef=1, status=results[,"P0-PRE"], values = c(-1, 1), main=colnames(fit2)[1])
plotMD(fit2, coef=2, status=results[,"P3-PRE"], values = c(-1, 1), main=colnames(fit2)[2])
plotMD(fit2, coef=3, status=results[,"P4-PRE"], values = c(-1, 1), main=colnames(fit2)[3])
dev.off()

glMDPlot(fit2, coef=1, status=results[,"P0-PRE"], main=colnames(fit2)[1],
         side.main="SYMBOL", counts=y.Norm$counts, groups=v$targets$chip,
         launch=TRUE, path="5-Glimma", folder="glMDPlots", html="P0-PRE")

glMDPlot(fit2, coef=2, status=results[,"P3-PRE"], main=colnames(fit2)[2],
         side.main="SYMBOL", counts=y.Norm$counts, groups=v$targets$timepoint,
         launch=TRUE, path="5-Glimma", folder="glMDPlots", html="P3-PRE")

glMDPlot(fit2, coef=3, status=results[,"P4-PRE"], main=colnames(fit2)[3],
         side.main="SYMBOL", counts=y.Norm$counts, groups=v$targets$timepoint,
         launch=TRUE, path="5-Glimma", folder="glMDPlots", html="P4-PRE")


#Volcano Plots
png("4-Output/Volcano_plots.png", width = 20, height = 30, units = 'cm', res = 600)
par(mfrow=c(3,1))
cutoff = -log10(0.05)
volcanoplot(fit2,coef=1,highlight=25,names=v$genes$SYMBOL,main=colnames(fit2)[1])
volcanoplot(fit2,coef=2,highlight=25,names=v$genes$SYMBOL,main=colnames(fit2)[2])
volcanoplot(fit2,coef=3,highlight=25,names=v$genes$SYMBOL,main=colnames(fit2)[3])
dev.off()


#https://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/Glimma.pdf

glXYPlot(x=fit2$coef[,"P0-PRE"], y=fit2$lod[,"P0-PRE"], xlab="logFC", ylab="logodds",
         status=results[,"P0-PRE"], main=colnames(fit2)[1], side.main="SYMBOL",
         anno=fit2$genes, counts=y.Norm$counts, groups=v$targets$timepoint,
         path="5-Glimma", folder="glVolcano", html="P0-PRE")

glXYPlot(x=fit2$coef[,"P3-PRE"], y=fit2$lod[,"P3-PRE"], xlab="logFC", ylab="logodds",
         status=results[,"P3-PRE"], main=colnames(fit2)[2], side.main="SYMBOL",
         anno=fit2$genes, counts=y.Norm$counts, groups=v$targets$timepoint,
         path="5-Glimma", folder="glVolcano", html="P3-PRE")

glXYPlot(x=fit2$coef[,"P4-PRE"], y=fit2$lod[,"P4-PRE"], xlab="logFC", ylab="logodds",
         status=results[,"P4-PRE"], main=colnames(fit2)[3], side.main="SYMBOL",
         anno=fit2$genes, counts=y.Norm$counts, groups=v$targets$timepoint,
         path="5-Glimma", folder="glVolcano", html="P4-PRE")


###############################################################

# which genes respond to heparin treatment -- P5 D3 vs D5?
# which genes respond to heparin treatment -- P5 D3 vs P13 D3?
# which genes respond to heparin treatment -- 20176s vs 21558s?
r1 <- topTable(fit2, adjust="BH", coef=1, n=Inf)
write.csv(r1, file="4-Output/topTable_P0-PRE.csv")
sig <- r1$adj.P.Val <0.05
cat("No.Sig.Genes.P0-PRE:", length(which(sig==1)))

# which genes respond to exercise - P3 relative to PRE?
r2 <- topTable(fit2, adjust="BH", coef=2, n=Inf)
write.csv(r2, file="4-Output/topTable_P3-PRE.csv")
sig <- r2$adj.P.Val <0.05
cat("No.Sig.Genes.P3-PRE:", length(which(sig==1)))

# which genes respond to exercise - P4 relative to PRE?
r3 <- topTable(fit2, adjust="BH", coef=3, n=Inf)
write.csv(r3, file="4-Output/topTable_P4-PRE.csv")
sig <- r3$adj.P.Val <0.05
cat("No.Sig.Genes.P4-PRE:", length(which(sig==1)))

sessionInfo()


######### END #######################


##MARTINAS CONTRIBUTION
TRANSCRIPTOME
------------------------------------------------LOAD PACKAGES----------------------------------------------------
library(tximport) 
library(readr)
library(edgeR)
library(Homo.sapiens)
library(limma)
library(ChAMP)
library(biomaRt)
library(goseq)
library(qusage)
library(maditr)
------------------------------------------------------------------------------------------------------------------
setwd("U:/Research/Projects/ihbi/grc/grc_general/Stem Cell Group/Sofia/Transcriptomics")
tx2gene <- read.csv("2-Input/Homo_sapiens.GRCh38.91_tx2gene.csv") 
head(tx2gene, 5)
folder <- c("2-Input/cell_quants") 
salmon.dir <- as.matrix(read.csv(file="2-Input/quant_filenames.csv", sep=",", header=F))
salmon.files <- file.path(folder, salmon.dir, "quant.sf") 
names(salmon.files) <- as.matrix(read.csv(file="2-Input/names.csv", sep=",", header=F)) 
all(file.exists(salmon.files))
txi <- tximport(salmon.files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = TRUE, countsFromAbundance="lengthScaledTPM")  ### NOTE : 4545 transcripts missing ####
names(txi)
head(txi$counts, 3) 
write.csv(txi, file = "4-Output/geneSMART_tximport_matrix.csv")

--------------------------------------------------------------------------------------------------------------------

y <- DGEList(txi$counts)
dim(y)
csvfile <- file.path("2-Input/Cell_sample_table.csv") 
sampleTable <- read.csv(csvfile, row.names=1) 
y$samples$names <- sampleTable$name 
y$samples$filename <- sampleTable$filename 
y$samples$Sample_ID <- sampleTable$Sample_ID 
y$samples$cell_line <- sampleTable$cell_line 
y$samples$Passage <- sampleTable$Passage 
y$samples$Day <- sampleTable$Day 
y$samples$Treatment <- sampleTable$Treatment 
y$samples$day <- sampleTable$Day 
y$samples$id <- sampleTable$Sample_ID 
y$samples$run_date <- sampleTable$run_date 
y$samples
geneid <- rownames(y) 
columns(Homo.sapiens) 
genes <- select(Homo.sapiens, key=geneid, keytype="ENSEMBL", columns=c("SYMBOL","GENEID","GENENAME","TXCHROM"), multiVals="first")
geneid <- rownames(y) 
ensembl91 <- useMart(host="dec2017.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl") 
attributes <- listAttributes(ensembl91) 
attributes[1:20,]
genes <- select(ensembl91, keys=geneid, keytype="ensembl_gene_id", columns=c("ensembl_gene_id", "external_gene_name", "description","entrezgene","chromosome_name","gene_biotype"))
(colnames(genes) <- c("ENSEMBL","SYMBOL","GENENAME","GENEID","TXCHROM","BIOTYPE"))
genes <- genes[!duplicated(genes$ENSEMBL),] 
y$genes <- genes 
head(genes,5) 
dim(genes)
------------------------------------------PLOT FOR SIZES-------------------------------------------------------------------------

png("3-QC/Barplot of library sizes.png", width = 90, height = 30, units = 'cm', res = 300) 
col <- as.numeric(y$sample$timepoint)
dt <- colSums((y$counts)*1e-6) 
barplot(dt, names=colnames(dt), 
col=c("lightsalmon","lightcoral","steelblue1","cornsilk")[col], las=2, cex.names=0.8)
abline(h=5,col="black",lty=3) 
abline(h=10,col="black",lty=3) 
abline(h=15,col="black",lty=3) 
title(main="Barplot of library sizes",ylab="Library size (millions)") 
dev.off()
-----------------------------------------------------------------------------------------------------------------------

summary(dt)
CountMeans <- rowMeans(y$counts) 
TotalCounts <- sum(CountMeans)
PercentReads <- (CountMeans/TotalCounts)*100
AvgCounts <- rowMeans(y$counts)
Symbol <- y$genes$SYMBOL 
GeneName <- y$genes$GENENAME
AvgCountsTable <- (data.frame(Symbol,GeneName,AvgCounts,PercentReads))
AvgCountsTable <- AvgCountsTable[order(-AvgCountsTable$AvgCounts),] 
head(AvgCountsTable,20) 
write.csv(AvgCountsTable, file ="4-Output/Average_Counts.csv")
top5 <- ((sum(head(AvgCountsTable, 5)$AvgCounts))/TotalCounts)*100
top10 <- ((sum(head(AvgCountsTable, 10)$AvgCounts))/TotalCounts)*100
top25 <- ((sum(head(AvgCountsTable, 25)$AvgCounts))/TotalCounts)*100
dim(y) 
summary(rowMeans(y$counts))
lcpm.AvgCounts <- cpm(AvgCounts, log = TRUE) 
----------------------------------------------------PLOT-------------------------------------------------------------------------------
png("3-QC/Read Count Distribution-B4.png", width = 45, height = 25, units = 'cm', res = 600) 
par(mar=c(5,6,4,1)+.1) 
hist(lcpm.AvgCounts, col="salmon", border="salmon", cex.lab=2.5, cex.main=3, cex.axis=2, xlab="Median log2-CPM", ylab="No. of Transcripts", breaks=100, xlim=c(-10,20), 
main ="Read Count Distribution (before)") 
------------------------------------------------------------------------------------------------------------------------------------
dev.off()

table(AvgCounts>=1) 
table(AvgCounts>=10) 
table(AvgCounts>=100) 
table(AvgCounts>=1000) 
table(AvgCounts>=10000) 
table(AvgCounts>=100000)
table(rowSums(y$counts==0)==211)

median_cpm <- apply(cpm(y), 1, median)
expr_cutoff <- 0.5 # in cpm sum(median_cpm > expr_cutoff)
y.Filt <- y[median_cpm > expr_cutoff, ]
lcpm.Filt <- cpm(y.Filt$counts, log = TRUE)
log.cutoff <- log2(expr_cutoff)
summary(rowMeans(y.Filt$counts))
-------------------------------------------------------------------------------------------------------
png("3-QC/Density of count values.png", width = 10, height = 20, units = 'cm', res = 600) 
nsamples <- ncol(y) 
col <- rainbow(nsamples) 
par(mfrow=c(2,1)) 
lcpm.Raw <- cpm(y$counts, log=TRUE) 
plot(density(lcpm.Raw[,1]), col=col[1], 
xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="") 
for (i in 2:nsamples){ den <-density(lcpm.Raw[,i]) 
lines(den$x, den$y, col=col[i]) } 
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="") 
title("Raw data",xlab="log2-CPM")

lcpm.Filt <- cpm(y.Filt$counts, log=TRUE) 
plot(density(lcpm.Filt[,1]), col=col[1], xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="") 
for (i in 2:nsamples){ den <-density(lcpm.Filt[,i]) 
lines(den$x, den$y, col=col[i]) } 
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="") 
title("Filtered data (median CPM > 0.5)",xlab="log2-CPM") 
dev.off()
---------------------------------------------------------------------------------------------------------
AvgCounts <- rowMeans(y.Filt$counts) 
lcpm.AvgCounts <- cpm(AvgCounts, log = TRUE) 

png("3-QC/Read Count Distribution-AFTER.png", width = 45, height = 25, units = 'cm', res = 600) 
par(mar=c(5,6,4,1)+.1) 
hist(lcpm.AvgCounts, col="salmon", border="salmon", cex.lab=2.5, cex.main=3, cex.axis=2, xlab="Median log2-CPM", ylab="No. of Transcripts", breaks=100, xlim=c(-10,20), 
main ="Read Count Distribution (after)") 
dev.off()

png("3-QC/Sample heatmap.png", width = 60, height = 60, units = 'cm', res = 300) 
par(mfrow=c(1,2)) 
lcpm.Raw <- cpm(y$counts, log = TRUE) 
heatmap(cor(lcpm.Raw)) 
title("Raw data") 
lcpm.Filt <- cpm(y.Filt$counts, log = TRUE) 
heatmap(cor(lcpm.Filt)) 
title("Filtered data") 
dev.off()
---------------------------------------------------------------------------------------------------------------
y.Norm <- calcNormFactors(y.Filt, method="TMM") 
y.Norm$samples$norm.factors 
summary(y.Norm$samples$norm.factors)

png("3-QC/BoxplotBefore.png", width = 90, height = 45, units = 'cm', res = 300) 
par(mfrow=c(2,1)) 
lcpm.Filt <- cpm(y.Filt, log=TRUE) 
col <- as.numeric(y$sample$tDay) 
boxplot(lcpm.Filt,las=2, col=c("lightsalmon","lightcoral","steelblue1","cornsilk")[col],main="") 
abline(h=median(lcpm.Filt),col="black",lty=3) 
title(main="Unnormalised data",ylab="log-counts")

png("3-QC/BoxplotAfter.png", width = 90, height = 45, units = 'cm', res = 300) 
par(mfrow=c(2,1)) 
lcpm.Norm <- cpm(y.Norm, log=TRUE) 
boxplot(lcpm.Norm,las=2, col=c("lightsalmon","lightcoral","steelblue1","cornsilk")[col],main="") 
abline(h=median(lcpm.Norm),col="black",lty=3) 
title(main="TMM-Normalised data",ylab="log-counts") 
dev.off()
-----------------------------------------------------THIS DOESNT WORK BUT I DONT THINK IT MATTERS-------------------------------------------------------------

dir.create("5-Glimma") 
glMDSPlot(y.Norm, top=500, labels=y.Norm$samples$subject, groups=y.Norm$samples[,c(1:9)], launch=TRUE, path="5-Glimma", folder="glMDSPlot", html="MDS_plot")

png("3-QC/MDS Plot_names.png", width = 15, height = 30, units = 'cm', res = 600) 
par(mfrow=c(2,1)) 
plotMDS(y.Norm, top=500, cex=0.8, labels=y.Norm$samples$cell_line, col=as.numeric(y.Norm$samples$Treatment), 
main="MDS Plot (SampleID)") 
legend("topright", legend=c("hMSC-20176", "hMSC-21558"), cex=0.8, col=1:16, pch=16)

plotMDS(y.Norm, top=500, cex=1, pch=21, col="white", bg=as.numeric(y.Norm$samples$Treatment), main="MDS Plot (Symbol)") pch=x; 0= open square; 1= open circle; 15= solid square, 
16=solid circle, 21=filled circle/border 
legend("topright", legend=c("hMSC-20176", "hMSC-21558"), cex=0.8, col=1:16, pch=16) dev.off()
png("3-QC/MDPlots-1LI00.png", width=11, height=20, units='cm', res=200) par(mfrow=c(2,1)) plotMD(y.Filt, column=1, main="First sample (raw)"); abline(h=0, col="red", lty=2, lwd=2) plotMD(y.Norm, column=1, main="First sample (TMM-normalised)"); abline(h=0, col="red", lty=2, lwd=2) dev.off()

-------------------------------------------------------------------------------------------------------------------------------------------------------------

data <- y$samples
group <- paste(data$Treatment,data$Day,sep=".")
group2 <- paste(data$cell_line,data$Treatment,sep=".")
group3 <- paste(data$cell_line,data$Treatment,data$Day,sep=".")
group4 <- paste(data$Treatment,data$Passage,sep=".")
group5 <- paste(data$cell_line,data$Treatment,data$Day,data$Passage,sep=".")

-------------------------------- TREATED VS UNTREATED * POOLED PASSAGE * POOLED DAYS * POOLED POPULATION------------------------------------------------------
design <- model.matrix(~0+data$Treatment)
colnames(design) <- c("Treated", "Untreated")

---------------------------------TREATED VS UNTREATED * POOLED POPULATIONS * POOLED PASSAGE * DIFFERENT DAYS---------------------------------------------------
design <- model.matrix(~0+group)
colnames(design2) <- c("Treated_D3", "Treated_D5", "Untreated_D3", "Untreated_D5")
contrast.matrix <- makeContrasts("Treated_D3-Untreated_D3", levels=design2)

---------------------------------TREATED VS UNTREATED * SEPARATE POPULATIONS * POOLED PASSAGES * POOLED DAYS-------------------------------------------------
design <- model.matrix(~0+group2)
colnames(design3) <- c("hMSC 220176 Treated", "hMSC 220176 Untreated", "hMSC 221558 Treated", "hMSC 221558 Untreated")

-------------------------------_TREATED VS UNTREATED * SEP POPULATIONS * POOLED PASSAGES * SEP DAYS------------------------------------------------------
design <- model.matrix(~0+group3)
colnames(design4) <- c("hMSC 220176 Treated Day 3", "hMSC 220176 Treated Day 5", "hMSC 220176 Untreated Day 3", 
"hMSC 220176 Untreated Day 5","hMSC 21558 Treated Day 3", "hMSC 21558 Treated Day 5", "hMSC 21558 Untreated Day 3", "hMSC 21558 Untreated Day 5")

-------------------------------TREATED VS UNTREATED * POOLED POP * SEP PASSAGES * POOLED DAYS------------------------------------------------------------

design <- model.matrix(~0+group4)
colnames(design) <- c("TreatedP13", "TreatedP5", "TreatedP7", 
"UntreatedP13", "UntreatedP5", "UntreatedP7")

---------------------------------EACH SAMPLE SEPARATE (DAYS, POP, PASSAGE, HEP NO HEP)-------------------------------------------------------------------

design6 <- model.matrix(~0+group5)
#### Leaving names as it is... Theres so much....

------------------------------------------------------------------------------------------------------------------------------------------------------


-----------------------------------------------FOLLOWING NICK-----------------------------------------------------------------------
png("3-QC/SA Plot.png", width = 30, height = 15, units = 'cm', res = 600) 
par(mfrow=c(1,2)) 
v <- voom(y.Norm, design, plot=TRUE) 
plotSA(fit2, main="Final model: Mean???variance trend") 
dev.off()

v2 <- voom(y.Norm, design)
corfit <- duplicateCorrelation(v2, design, block=y.Norm$samples$cell_line)
v <- voom(y.Norm, design, block = y.Norm$samples$cell_line, correlation =
corfit$consensus)
fit2 <- lmFit(v, design, block=y.Norm$samples$cell_line, correlation=corfit$consensus)
contrast.matrix <- makeContrasts("TreatedP5-TreatedP13", levels=design)
fit2c <- contrasts.fit(fit2, contrast.matrix)
names(fit2c)

METHOD1
fit2c <- eBayes(fit2c)
Then use TopTable

METHOD2
lfc <- log2(1.1)
fit2c <- treat(fit2c, lfc=lfc)
save(fit2c, file="4-Output/v.rda")
Then use topTreat

results <- decideTests(fit2c) 
summary(results) 
vennCounts(results)
ngenes <- which(results[,1] != 0)
ngenes <- length(ngenes)

png("4-Output/Venn Diagram.png", width = 30, height = 15, units = 'cm', res = 600) 
par(mfrow=c(1,2)) 
vennDiagram(results,include=c("both"), circle.col=c("blue","yellow","green","orange"), counts.col=c("blue3"), cex=c(1,0.8,0.8)) 
vennDiagram(results,include=c("up","down"), circle.col=c("blue","yellow","green","orange"), counts.col=c("red","green3"), 
cex=c(1,0.7,0.7)) 
mtext("Cell hMSC", side = 3, line = -2, outer = TRUE, cex=1.8) 
mtext(paste0("(", ngenes," genes)"), side = 3, line = -3.5, outer = TRUE, cex=1.5) 
dev.off()
-----------------------------------------------FOLLOWING ANOTHER PROTOCOL IF NEEDED ----------------------------------------------------------
design <- model.matrix(~group4)
dgeObj <- estimateCommonDisp(y.Norm, verbose=TRUE)
Disp = 0.30772 , BCV = 0.5547 
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)
plotBCV(dgeObj, xlab="Average log CPM", ylab="Biological coefficient of variation", pch=16, cex=0.2, col.common="red", col.trend="blue", col.tagwise="black")
fitd <- glmFit(dgeObj, design)
names(fitd)
head(coef(fitd))
lrt.BvsL <- glmLRT(fitd) 
topTags(lrt.BvsL)

PvsV <- makeContrasts(group4Treated.P5-group4Untreated.P5, levels=design)
lrt.pVsV <- glmLRT(fitd, contrast=PvsV)
topTags(lrt.pVsV)
results <- as.data.frame(topTags(lrt.pVsV,n = Inf))
dim(results)
summary(de <- decideTestsDGE(lrt.pVsV))
detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.pVsV, de.tags=detags)
------------------------------------------------------------------------------------------------------------------------------------------------------------

png("4-Output/MDplot-Average.png", width = 20, height = 30, units = 'cm', res = 600)  
plotMD(fit2c, coef=1, status=results[,1], values = c(-1, 1), main=colnames(fit2c)[1]) 
dev.off()
glMDPlot(fit2c, coef=1, status=results[,1], main=colnames(fit2c)[1], side.main="SYMBOL", counts=y.Norm$counts, groups=v$targets$chip, launch=TRUE, path="5-Glimma", folder="glMDPlots", html="Treated P5 - Treated P13")

###This shows an error but seems to be working####

--------------------------------------------------CAN ALSO PLOT MORE THAN ONE COMPARISON BUT ADD PAR--------------------------------------------------------
par(mfrow=c(3,1)) <--- change accordingly
plotMD(fit2c, coef=2, status=results[,2], values = c(-1, 1), main=colnames(fit2)[2]) 
plotMD(fit2, coef=3, status=results[,3], values = c(-1, 1), main=colnames(fit2)[3]) 
glMDPlot(fit2c, coef=2, status=results[,2], main=colnames(fit2c)[1], side.main="SYMBOL", counts=y.Norm$counts, groups=v$targets$chip, launch=TRUE, path="5-Glimma", folder="glMDPlots", html="Treated P5 - Treated P13")
glMDPlot(fit2c, coef=3, status=results[,3], main=colnames(fit2c)[1], side.main="SYMBOL", counts=y.Norm$counts, groups=v$targets$chip, launch=TRUE, path="5-Glimma", folder="glMDPlots", html="Treated P5 - Treated P13")

png("4-Output/Volcano_plots.png", width = 20, height = 30, units = 'cm', res = 600) 
cutoff = -log10(0.05) 
volcanoplot(fit2c,coef=1,highlight=25,names=v$genes$SYMBOL,main=colnames(fit2c)[1]) 
dev.off()

glXYPlot(x=fit2c$coef[,1], y=fit2c$p.value[,1], xlab="logFC", ylab="p-value", status=results[,1], main=colnames(fit2c)[1], 
side.main="SYMBOL", anno=fit2c$genes, counts=y.Norm$counts, groups=v$targets$timepoint, path="5-Glimma", folder="glVolcano", html="Treated P5 - Treated P13")



---------------------------------------------------ADD MORE-----------------------------------------------------------------------------------------------
par(mfrow=c(3,1))
volcanoplot(fit2,coef=2,highlight=25,names=v$genes$SYMBOL,main=colnames(fit2)[2]) 
volcanoplot(fit2,coef=3,highlight=25,names=v$genes$SYMBOL,main=colnames(fit2)[3]) 
----------------------------------------------------------------------------------------------------------------------------------------------------------
r1 <- topTreat(fit2c, adjust="BH", coef=1, n=Inf) 
write.csv(r1, file="4-Output/topTable_TvsUT.csv") 
sig <- r1$adj.P.Val <0.05 
cat("No.Sig.Genes.Treated P5 - Treated P13:", length(which(sig==1)))

---------------------------------------------------GENE SET TESTING-------------------------------------------------------------------------------------------------------

results <- as.data.frame(r1)
genes <- as.list(r1$adj.P.Val<0.05)
names(genes) <- results$SYMBOL
result <- as.data.frame(genes[1:265])
result <- names(result)
gene.df <- bitr(result, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID" ),
                OrgDb = org.Hs.eg.db)
ggo <- groupGO(gene     = result,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)
ggo_df <- data.frame(ggo)
write.csv(ggo_df, file="4-Output/GO_Analysis.csv") 


results.ord <- r1[ order(-r1[,"logFC"]), ]
head(results.ord)
ranks <- results.ord$logFC
names(ranks) <- results.ord$SYMBOL
head(ranks)

pathways.hallmark <- gmtPathways("2-Input/h.all.v7.0.symbols.gmt")
head(pathways.hallmark)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) 
gene.in.pathway <- pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest(cols = c(SYMBOL)) %>% 
  inner_join(r1, by="SYMBOL")


dev.new()
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = factor(padj))) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Hallmark pathways Enrichment Score from GSEA")
dev.off()


plotEnrichment(pathway = pathways.hallmark[["HALLMARK_HEDGEHOG_SIGNALING"]], ranks)
plotGseaTable(pathways.hallmark[fgseaRes$pathway[fgseaRes$padj < 0.05]], ranks, fgseaRes, 
                gseaParam=0.5)

sig.path <- fgseaResTidy$pathway[fgseaResTidy$adjPvalue == "significant"]
sig.gen <- unique(na.omit(gene.in.pathway$SYMBOL[gene.in.pathway$pathway %in% sig.path]))
h.dat <- dcast(gene.in.pathway[, c(1,2)], value.var="SYMBOL", SYMBOL~pathway, fun.aggregate=NULL)
rownames(h.dat) <- h.dat$SYMBOL
h.dat <- h.dat[, -1]


h.dat <- h.dat[rownames(h.dat) %in% sig.gen, ]
h.dat <- h.dat[, colnames(h.dat) %in% sig.path]

table(data.frame(rowSums(as.data.frame(h.dat))))
h.dat <- h.dat[data.frame(rowSums(h.dat)) >= 3, ]
topTable <- results.ord[results.ord$SYMBOL %in% rownames(h.dat), ]
rownames(topTable) <- topTable$SYMBOL

topTableAligned <- topTable[which(rownames(topTable) %in% rownames(h.dat)),]
head(topTableAligned) <- topTableAligned[match(rownames(h.dat), rownames(topTableAligned)),]
all(rownames(topTableAligned) == rownames(h.dat))







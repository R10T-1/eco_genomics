# WGCNA 
#script for analzing and visualizing gene correlation networks

library(DESeq2) 
library(ggplot2)
library(WGCNA); options(stringAsFactors = FALSE);
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)
setwd("~/projects/eco_genomics/transcriptomics/")

##import counts data from 01
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header=TRUE, row.names = 1)
tail(countsTable)
dim(countsTable)

#DESeq2 doesnt like decimals, but salmon uses them, so we need to round 
countsTableRound <- round(countsTable) 

tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header=TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

tail(countsTable)
dim(countsTableRound)

countsTableRound <- round(countsTable)
tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", 
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

traitData = read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header = T, row.names = 1)

###1) subset data to just baseline samples 

filtered_count_matrix_BASEonly <- countsTable[, conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE",]
rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)

###2) detecting outliers
#find outlier genes 
gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg)

table(gsg$goodGenes)
#FALSE  TRUE 
#37235 82203 

#filter out bad genes 
data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes==TRUE,]
dim(data_WGCNA)
#[1] 82203     7

#use clustering with a tree dendrogram to identify outlier samples 
htree <- hclust(dist(t(data_WGCNA)), method="average")
plot(htree)

#PCA - outlier detection method
pca <- prcomp(t(data_WGCNA))
pca_data <- pca$x
# make a data frame 
pca_data <- as.data.frame(pca_data)

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

ggplot(pca_data, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2:', pca.var.percent[2]))

###3) normalization
colData <- row.names(filtered_sample_metadata_BASEonly)

#run deseq2 without model defined, this is going to let wgcna cluster without any knowledge of groups
dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA,
                                    colData = filtered_sample_metadata_BASEonly,
                                    design=~1) #there are no specified groups

#filter again getting out low counts 
dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >=6,]
nrow(dds_WGCNA_75) #filtered down to 299559 transcripts

#use variance normalization 
dds_norm <- vst(dds_WGCNA_75) #perform variance stabilization 

#extract normalized counts 
norm.counts <- assay(dds_norm) %>%
  t()

###4 network construction 
#choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

#call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
#"signed" means we focus on transcripts that are positively correlated

sft.data <- sft$fitIndices

#plot to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = "Power", y = "Scale free topology model, signed R^2") +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

#WGCNA wants us to pick a degree of relatedness 
#discription of network, some clusters are closely related and some are not
#power thresholds are the level of the correlation
#so basically it says that, when we pick different relateness levels how does that affect what we predict using the natural model of gene expression 
#as you increase power(strength of relationship) you lower the correlation 
#you need to chose what power is going to maximize the biological relevance while maintaining connectivity
#we will chose a threshold around 24-26
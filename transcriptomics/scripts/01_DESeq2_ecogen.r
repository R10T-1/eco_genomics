### code for analyzing RNASeq data using DESeq2 

# load libraries 
library(DESeq2) 
library(ggplot2)


setwd("~/projects/eco_genomics/transcriptomics/")

#import counts matrix 
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header=TRUE, row.names = 1)
tail(countsTable)
dim(countsTable)

#DESeq2 doesnt like decimals, but salmon uses them, so we need to round 
countsTableRound <- round(countsTable) 

tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header=TRUE, stringsAsFactors = TRUE, row.names = 1)
conds

############################################
# explore counts matrix
############################################

#lets see how many reads we got for each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound),
        cex.names = 0.5, las = 2, ylim = c(0,30000000))
abline(h=mean(colSums(countsTableRound)), col = "blue4", lwd=2)

#the average number of counts per gene 
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) #3244.739
median(rowSums(countsTableRound)) #64
#the stark difference means that there is a lot of overdispersion

apply(countsTableRound,2,mean) #gives a sense of variation in sequencing effort across samples 

#########################################
#start DESeq2 analysis
########################################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ DevTemp + FinalTemp)
dim(dds)

#start filtering 
dds <- dds[rowSums(counts(dds) >= 10) >= 15,] #we are filtering the transcripts by taking the sums of the rows (transcripts)
#if when you sum, theres not more than 10 reads across 15 samples, we filter it out 
nrow(dds) #35,527 = number of transcripts with more than 10  reads and in greater than or equal to 15 samples 

# run the DESeq2 model to test for global differential gene expressiom 
dds <- DESeq(dds)

# list the results
resultsNames(dds)
# [1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28" 
# [4] "FinalTemp_BASE_vs_A28"

#visualize our global gene expression patterns using pca
#first we need to tranform the data for plotting using variance stabilization 

vsd <- vst(dds, blind=FALSE)

#build plot 
pcaData <- plotPCA(vsd, intgroup=c("DevTemp","FinalTemp"), returnData = TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))

final_temp_colors <- c("BASE" = "grey", "A28" = "hotpink", "A33" = 'red')
shapes_choose <- c("D18" = 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, shape = DevTemp)) + 
  geom_point(size = 5) +
  scale_shape_manual(values = shapes_choose) +
  scale_color_manual(values = final_temp_colors) + 
  labs(x = paste0("PC1: ", percentVar[1], ' %'),
       y = paste0("PC2: ", percentVar[2], " %")) +
  theme_bw(base_size = 16)
p
  
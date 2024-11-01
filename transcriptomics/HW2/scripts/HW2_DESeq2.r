###getting set up for DESeq2 

setwd("~/projects/eco_genomics/transcriptomics/")
library(DESeq2)
library(ggplot2)

#import the counts matrix
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, row.names = 1)
#round the counts for DESeq2
countsTableRound <- round(countsTable)

tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, row.names = 1)
conds

##before I do any plotting I have to set my margins- 
##I dont know why I always have an issue with margins in ggplot
par(mar=c(1,1,1,1))
####################### explore the counts matrix #########################

########reads per sample######## 

colSums(countsTableRound) 
# N1C3     N2C2    sN1L3     N1C4    sN2H2     N1L2     N1H2     N2H1     N1H3     N1H4    sN1L1     N1C2     N2L4 
# 16599225 18041715 30040661 16944319 16495945 16937756 22077351 18706970 17993067 16797749 18255512 17360936 18659738 
# N2C3    sN2C5     N2L2     N1L4     N2L3     N2H3     N2H4     N2C1 
# 24349215 12265167 19248965 20929974 16166346 17165334 16782595 15726568 

mean(colSums(countsTableRound)) #18454529

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound),
        cex.names = 0.5, las = 2, ylim = c(0,30000000))
abline(h=mean(colSums(countsTableRound)), col = "mediumblue", lwd = 2)

########mean counts per gene######## 

rowSums(countsTableRound)
mean(rowSums(countsTableRound)) #3244.739
median(rowSums(countsTableRound)) #64

apply(countsTableRound, 2, mean)
# N1C3     N2C2    sN1L3     N1C4    sN2H2     N1L2     N1H2     N2H1     N1H3     N1H4    sN1L1     N1C2     N2L4 
# 138.9778 151.0551 251.5168 141.8671 138.1130 141.8121 184.8436 156.6249 150.6478 140.6399 152.8451 145.3552 156.2295 
# N2C3    sN2C5     N2L2     N1L4     N2L3     N2H3     N2H4     N2C1 
# 203.8649 102.6907 161.1628 175.2371 135.3535 143.7175 140.5130 131.6714 


###################### Start DESeq 2 for real #################################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds,
                              design = ~ DevTemp + FinalTemp)
dim(dds) #[1] 119438     21

#start filtering#
dds <- DESeq(dds)

#list the results 
resultsNames(dds) 
#[1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

#variance stabilization so we can start visualizing global data 
vsd <- vst(dds, blind = FALSE)

#build plot 
pcaData <- plotPCA(vsd, intgroup=c("DevTemp", "FinalTemp"), returnData = TRUE)
percentVar <- round(100*attr(pcaData, "percentVAR"))

final_temp_colors <- c("BASE" = "azure4", "A28" = "chartreuse2", "A33" = "olivedrab")
shapes_choose <- c("D18" = 16, "D22" = 18) #where do we get these values?

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, chape = DevTemp)) +
    geom_point(size = 5) +
    scale_shape_manual(values = shapes_choose) +
    scale_color_manual(values = final_temp_colors) +
    labs(x = paste0("PC1: ", percentVar[1], '%'),
         y = paste0("PC2: ", percentVar[2], '%')) +
    theme_bw(base_size = 16)
p
#ggsave("~/projects/eco_genomics/transcriptomics/HW2/figures/HW2_PCA_BASE_28_33.pdf")

dev.off()

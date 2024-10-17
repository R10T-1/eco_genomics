##relies on running both previous scripts 
library(eulerr)

#set up groups with DESeq object 
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp)) #setup factors
design(dds) <- ~ group #set up the groups using the factors
dds <- DESeq(dds)
dim(dds)
resultsNames(dds) #these are the results we will refer back to
#[1] "Intercept"               "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28"
#[4] "group_D22A28_vs_D18A28"  "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

#make separate files for each contrast so that we can compare what is differentially expressed within the groups 
#we are comparing within groups adn between groups (compare devtemp vs final temp or compare the levels within the groups)

#actually do the contrasts with the groups

  #1. compare baseline gene expression between developmental treatment groups
res_D18_BASE_D22_BASE <- results(dds, contrast=c("group", "D18BASE", "D22BASE"), alpha = 0.05)
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_[!is.na(res_D18_BASE_D22_$padj),] #getting rid of empty data points and adjusting pvalues
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),] #order by adjusted p value
summary(res_D18_BASE_D22_BASE)
      # out of 32768 with nonzero total read count
      # adjusted p-value < 0.05
      # LFC > 0 (up)       : 556, 1.7%
      # LFC < 0 (down)     : 1379, 4.2%
      # outliers [1]       : 0, 0%
      # low counts [2]     : 0, 0%
      # (mean count < 18)
      # [1] see 'cooksCutoff' argument of ?results
      # [2] see 'independentFiltering' argument of ?results 

      #make a list of whcih genes in our comparison of interest are diffirentially expressed (list of DEGs) (remember in deseq each row is a gene)
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj<0.05,])

plotMA(res_D18_BASE_D22_BASE, ylim=c(-4,4))

  # 2. compare gene expression between dev temp treatment groups at a28

res_D18_A28_D22_A28 <- results(dds, contrast=c("group", "D18A28", "D22A28"), alpha = 0.05)
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),] #getting rid of empty data points and adjusting pvalues
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),] #order by adjusted p value
summary(res_D18_A28_D22_A28)

degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj<0.05,])

plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))

  # 3. compare gene expression between dev temp treatment groups at A33

res_D18_A33_D22_A33 <- results(dds, contrast=c("group", "D18A33", "D22A33"), alpha = 0.05)
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),] #getting rid of empty data points and adjusting pvalues
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),] #order by adjusted p value
summary(res_D18_A33_D22_A33)

degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj<0.05,])

plotMA(res_D18_A33_D22_A33, ylim=c(-4,4)) 

#we have our lists of DEGs for each contrast so now we build euler plot 
length(degs_D18_BASE_D22_BASE) #1935 
length(degs_D18_A28_D22_A28) #296
length(degs_D18_A33_D22_A33) #78 

#look at the overlaps in which genes are DEGs in multiple contrasts
length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A28_D22_A28)) #107
length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A33_D22_A33)) #44
length(intersect(degs_D18_A33_D22_A33,degs_D18_A28_D22_A28)) #29

length(intersect(degs_D18_BASE_D22_BASE,(intersect(degs_D18_A28_D22_A28, degs_D18_A33_D22_A33)))) #23 
  #you could also name a variable nested_intersection <- [intersection of two of them] and then use that as an argument in an intersect with the last group


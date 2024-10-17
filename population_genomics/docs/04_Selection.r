#this script will test fro selection outliers using the PCAdapt method 
#this is an extension of the Fst outlier method, developed by Duforest-Fribourg et. al 

library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)

#first we use PCAdapt to read in our (uncompressed) VCF data from class 
vcf <- read.pcadapt("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf",
                    type="vcf")

#we also need the compressed version to subset the metadata file with 
vcfR <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf.gz")

#NOTE they ran into issue in class and used the 2 files above stored in class dir
#use the paths for whatever file you are using (this one or ur own) 
#uncompressed files go in the read.pcadapt and the compressed go in read.vcfr

#now, read in the metadata and subset it for the indvs in vcfr objects 
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(vcfR@gt[,1]),]

#run pcadapt, options annotated 
pcadapt.pca <- pcadapt(vcf,
                       K=2,
                       method="componentwise",
                       min.maf=0.01,
                       LD.clumping = list(size=500, thr=0.2))

#K: number of PC axes you want to test for 
#method: selecting "componentwise" returns a separate test for selection on each PC axis; otherwise it reutrns a global selection test across al PS axes together
#min.maf: the minor allele frequency value to use when testing -- SNPs below this value wont get tested 
#LD.clumping: removes loci in strong LD when initially fitting the PCs, then tests for selection using all loci; size is distance bwteeen adjacemt SNPs in bp; thr = threshold value of the correlation coefficient (r^2, a measure of LD)between adjacent loci

summary(pcadapt.pca)

#quick manhattan plot 
plot(pcadapt.pca, K=2)

#to make the manhattan plot where we can actually see the outliers as a function
#of chromosomes and position, we need ti get the SNP info out of the "fix" -
#region of the vcfR object

#we can get a quick view of that region by combining view and head 
View(head(vcfR@fix))

#this shows us that the chromosome and position info are in the first few columns
#lets grab that info
vcfR.fix <- as.data.frame(vcfR@fix[,1:2])

#we need to create a new variable that relabels each chromosome to a numeric value
#here we'll jusr keep the major 8 chromosomes and not label small scaffolds 
chr.main <- unique(vcfR.fix$CHROM)[1:8]


#make a data fram by binding together the unique chroms w info
chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))

#now extract pvalues
# note: rows number exaclty the same # of loci to bind properly 
Pval <- pcadapt.pca$pvalues 
pcadapt.MHplot <- cbind(vcfR.fix, Pval) 
names (pcadapt.MHplot) [4:5] <- c("pPC1", "pPC2") 

#now use tidyr:left_join to assign the  hrom numbers to loci 
pcadapt.MHplot <- left_join(chrnum, pcadapt.MHplot, join_by(chr.main==CHROM))



#create SNP variables for qqman
pcadapt.MHplot <- pcadapt.MHplot %>%
  mutate(SNP=paste0(chr.main, "_", POS))

#set variables as numeric that contain numbers
pcadapt.MHplot$V2 = as.numeric(pcadapt.MHplot$V2)
pcadapt.MHplot$POS = as.numeric(pcadapt.MHplot$POS)
pcadapt.MHplot$pPC1 = as.numeric(pcadapt.MHplot[,4])
pcadapt.MHplot$pPC2 = as.numeric(pcadapt.MHplot[,5])

#drop loci with NA, these got filitered by the min.maf 
pcadapt.MHplot <- pcadapt.MHplot %>% drop_na(pPC1)

#finially, finish the manhattan plot 
#lets make seperate once for selection outliers on PC1 and PC2 
manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log10 p-value",
          genomewideline = F,
          main="PCAdapt genome scan for selection (PC1)")
manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC2",
          col=c("blue3","orange2"),
          logP=T,
          ylab="-log10 p-value",
          genomewideline = F,
          main="PCAdapt genome scan for selection (PC2)")

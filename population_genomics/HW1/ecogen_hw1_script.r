#sof's ecogen HW1 

##Start by loading all the libraries and set wd to avoid issues later on 
library(vcfR)
library(tidyverse)
library(qqman)
library(SNPfiltR)
library(LEA)
library(ggplot2)
library(pcadapt)
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
install.packages("ape")

##note: plots kept getting error that margins were too large so i used par(mar)
par("mar")
par(mar=c(1,1,1,1))

#this fixed my margins and allowed me to plot 

list.files("variants/")
list.files("reference/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")
#dna2 = dna$CM058040.1 Centaurea solstitialis isolate CAN-66 chromosome 1, whole genome shotgun sequence

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")
#gff2 = grep("CM058040.1",gff)

vcf
head(vcf)

#chromR object 
chromR1 <- create.chromR(name="Centaurea Chromosome 1", vcf=vcf, seq=dna, ann=gff)

plot(chromR1)

#make a pdf with the plot
pdf(file="~/projects/eco_genomics/population_genomics/HW1/chromoPlot_chromR1.pdf")
chromoqc(chromR1, xlim=c(1e1, 1.1e8))
dev.off()

#mask poor quality variants
chromR1_masked <- masker(chromR1, 
                         min_QUAL = 50,
                         min_DP = 1000,
                         max_DP = 10000,
                         min_MQ = 30)
plot(chromR1_masked)

chromoqc(chromR1_masked, xlim=c(1e1, 1.1e8))

#process the object with proc
chromR1_proc <- proc.chromR(chromR1_masked, win.size = 1e5)
plot(chromR1_proc)
chromoqc(chromR1_proc, xlim=c(1e1, 1.1e8))
head(chromR1_proc)

DP <- extract.gt(vcf, element = "DP", as.numeric = T)

quantile(DP)

DP[DP==0] <- NA #set 0 depth genos to 'NA'

quantile(DP, na.rm=T) #much better 
dim(DP) #ensure loci are in rows 

#average DP per individual 
avgDP = colMeans(DP, na.rm=T)
summary(avgDP)
hist(avgDP, breaks=50)

#missingness heatmap 

pdf(file="~/projects/eco_genomics/population_genomics/HW1/chromR1_missingness.pdf")
heatmap.bp(DP[1:5000,], rlabels=F, clabels=F)
dev.off()

#set individual genotypes 
vcf@gt[,-1][is.na(DP)==TRUE] <- NA
vcf

#time to filter 

meta <- read.csv("metadata/meta4vcf.csv", header=T)

meta2=meta[,c(1,4)]
names(meta2) = c("id","pop")
meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

#mean depth per individual 
hard_filter(vcf)
vcf.filt <- hard_filter(vcf, 
                        depth=3)
#look at allele balance (note these are autotetraploid ... 0.25, etc.)
vcf.filt <- filter_allele_balance(vcf.filt,
                                  min.ratio = 0.15,
                                  max.ratio = 0.85)
vcf.file <- max_depth(vcf.filt,
                     maxdepth = 60)
###this is at 25
vcf.filt.indMiss.25 <- missing_by_sample(vcf.filt,
                                      popmap=meta2, 
                                      cutoff=0.25)

#subset popmap to only retained individuals
meta <- meta[meta$id %in% colnames(vcf.filt.indMiss.25@gt),]

#get rid of monomorphic or multi-allelic sites
vcf.filt.indMiss.25 <- filter_biallelic(vcf.filt.indMiss.25)
vcf.filt.indMiss.25 <- min_mac(vcf.filt.indMiss.25, min.mac = 1)
vcf.filt.indSNPMiss.25 <- missing_by_snp(vcf.filt.indMiss.25, cutoff=0.5)
DP2.25 <- extract.gt(vcf.filt.indSNPMiss.25,
                  element="DP",
                  as.numeric=T)
#use PCA to verify missingness is not driving clustering 
install.packages("adgenet")

missPCA.25 <- assess_missing_data_pca(vcfR=vcf.filt.indMiss.25,
                                   popmap = meta, 
                                   thresholds = 0.1,
                                   clustering = FALSE)
write.vcf(vcf.filt.indSNPMiss.25, "~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.25.vcf.gz")

###this is at 50

vcf.filt.indMiss.50 <- missing_by_sample(vcf.filt,
                                      popmap=meta2, 
                                      cutoff=0.50)

#subset popmap to only retained individuals
meta <- read.csv("metadata/meta4vcf.csv", header=T)

meta2=meta[,c(1,4)]
names(meta2) = c("id","pop")
meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

meta <- meta[meta$id %in% colnames(vcf.filt.indMiss.50@gt),]

#get rid of monomorphic or multi-allelic sites
vcf.filt.indMiss.50 <- filter_biallelic(vcf.filt.indMiss.50)
vcf.filt.indMiss.50 <- min_mac(vcf.filt.indMiss.50, min.mac = 1)
vcf.filt.indSNPMiss.50 <- missing_by_snp(vcf.filt.indMiss.50, cutoff=0.5)
DP2.50 <- extract.gt(vcf.filt.indSNPMiss.50,
                  element="DP",
                  as.numeric=T)
#use PCA to verify missingness is not driving clustering 
install.packages("adgenet")

missPCA <- assess_missing_data_pca(vcfR=vcf.filt.indMiss.50,
                                   popmap = meta, 
                                   thresholds = 0.1,
                                   clustering = FALSE)
write.vcf(vcf.filt.indSNPMiss.50, "~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.50.vcf.gz")

###75

vcf.filt.indMiss.75 <- missing_by_sample(vcf.filt,
                                      popmap=meta2, 
                                      cutoff=0.75)

#subset popmap to only retained individuals

meta <- read.csv("metadata/meta4vcf.csv", header=T)

meta2=meta[,c(1,4)]
names(meta2) = c("id","pop")
meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

meta <- meta[meta$id %in% colnames(vcf.filt.indMiss.75@gt),]

#get rid of monomorphic or multi-allelic sites
vcf.filt.indMiss.75 <- filter_biallelic(vcf.filt.indMiss.75)
vcf.filt.indMiss.75 <- min_mac(vcf.filt.indMiss.75, min.mac = 1)
vcf.filt.indSNPMiss.75 <- missing_by_snp(vcf.filt.indMiss.75, cutoff=0.5)
DP2 <- extract.gt(vcf.filt.indSNPMiss.75,
                  element="DP",
                  as.numeric=T)
#use PCA to verify missingness is not driving clustering 
install.packages("adgenet")

missPCA <- assess_missing_data_pca(vcfR=vcf.filt.indMiss.75,
                                   popmap = meta, 
                                   thresholds = 0.1,
                                   clustering = FALSE)
write.vcf(vcf.filt.indSNPMiss.75, "~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.75.vcf.gz")

#end of filtering, time to do diversity and differentiation

setwd("~/projects/eco_genomics/population_genomics/")
### at 25 
X11.options(type="cairo")
options(bitmapType = "cairo")

vcf25 <- read.vcfR("~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.25.vcf.gz")

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta25 <- meta[meta$id %in% colnames(vcf25@gt[,-1]),] 

vcf25.div <- genetic_diff(vcf25,
                        pops=as.factor(meta25$region),
                        method="nei")
#visualize the data 
str(vcf25.div)
unique(vcf25.div$CHROM)

chr.main.25 <- unique(vcf25.div$CHROM)[1:8] 
chrnum.25 <- as.data.frame(cbind(chr.main.25, seq(1,8,1)))

#merge vcf25.div with chrnum25
vcf25.div.MHplot <- left_join(chrnum.25, vcf25.div, join_by(chr.main.25==CHROM))

vcf25.div.MHplot <- vcf25.div.MHplot %>% 
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main.25,"_",POS))

vcf25.div.MHplot$V2 = as.numeric(vcf25.div.MHplot$V2)
vcf25.div.MHplot$POS = as.numeric(vcf25.div.MHplot$POS)

#make manhattan plot 
manhattan(vcf25.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf25.div.MHplot$Gst,0.999))

write.csv(vcf25.div.MHplot, "~/projects/eco_genomics/population_genomics/HW1/Genetic_diff_by_Region_HW1_.25.csv",
          quote=F,
          row.names = F)
#looking at Hs values 
names(vcf25.div.MHplot)

vcf25.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title="Genome-wide expected heterozygosity (Hs)",fill="Regions",
       x="Gene diversity within Regions", y="Counts of SNPs")

ggsave("Histogram_GenomDiversity_byRegion_HW1_.25.pdf",
       path="~/projects/eco_genomics/population_genomics/HW1/")

vcf25.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>% 
  filter(value!=0) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())

### at 50 
X11.options(type="cairo")
options(bitmapType = "cairo")

vcf50 <- read.vcfR("~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.50.vcf.gz")

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta50 <- meta[meta$id %in% colnames(vcf50@gt[,-1]),] 

vcf50.div <- genetic_diff(vcf,
                         pops=as.factor(meta50$region),
                         method="nei")
#visualize the data 
str(vcf50.div)
unique(vcf50.div$CHROM)

chr.main.50 <- unique(vcf50.div$CHROM)[1:8] 
chrnum.50 <- as.data.frame(cbind(chr.main.50, seq(1,8,1)))

#merge vcf2.div with chrnum
vcf50.div.MHplot <- left_join(chrnum.50, vcf50.div, join_by(chr.main.50==CHROM))

vcf50.div.MHplot <- vcf50.div.MHplot %>% 
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main.50,"_",POS))

vcf50.div.MHplot$V2 = as.numeric(vcf50.div.MHplot$V2)
vcf50.div.MHplot$POS = as.numeric(vcf50.div.MHplot$POS)

#make manhattan plot 
manhattan(vcf50.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf50.div.MHplot$Gst,0.999))

write.csv(vcf50.div.MHplot, "~/projects/eco_genomics/population_genomics/HW1/Genetic_diff_by_Region_HW1_.50.csv",
          quote=F,
          row.names = F)
#looking at Hs values 
names(vcf50.div.MHplot)

vcf50.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title="Genome-wide expected heterozygosity (Hs)",fill="Regions",
       x="Gene diversity within Regions", y="Counts of SNPs")

ggsave("Histogram_GenomDiversity_byRegion_HW1_.50.pdf",
       path="~/projects/eco_genomics/population_genomics/HW1/")

vcf50.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>% 
  filter(value!=0) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())

##at 75 
X11.options(type="cairo")
options(bitmapType = "cairo")

vcf75 <- read.vcfR("~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.75.vcf.gz")

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta75 <- meta[meta$id %in% colnames(vcf75@gt[,-1]),] 

vcf75.div <- genetic_diff(vcf,
                          pops=as.factor(meta75$region),
                          method="nei")
#visualize the data 
str(vcf75.div)
unique(vcf75.div$CHROM)

chr.main.75 <- unique(vcf75.div$CHROM)[1:8] 
chrnum.75 <- as.data.frame(cbind(chr.main.75, seq(1,8,1)))

#merge vcf2.div with chrnum
vcf75.div.MHplot <- left_join(chrnum.75, vcf75.div, join_by(chr.main.75==CHROM))

vcf75.div.MHplot <- vcf75.div.MHplot %>% 
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main.75,"_",POS))

vcf75.div.MHplot$V2 = as.numeric(vcf75.div.MHplot$V2)
vcf75.div.MHplot$POS = as.numeric(vcf75.div.MHplot$POS)

#make manhattan plot 
manhattan(vcf75.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf75.div.MHplot$Gst,0.999))

write.csv(vcf75.div.MHplot, "~/projects/eco_genomics/population_genomics/HW1/Genetic_diff_by_Region_HW1_.75.csv",
          quote=F,
          row.names = F)
#looking at Hs values 
names(vcf75.div.MHplot)

vcf75.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title="Genome-wide expected heterozygosity (Hs)",fill="Regions",
       x="Gene diversity within Regions", y="Counts of SNPs")

ggsave("Histogram_GenomDiversity_byRegion_HW1_.75.pdf",
       path="~/projects/eco_genomics/population_genomics/HW1/")

vcf75.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>% 
  filter(value!=0) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n()) 

#vcf50.div.MHplot %>%
 # as_tibble()%>%
  #pivot_longer(cols = c(4:9), names_to = "region", values_to = "Hs") %>%
  #group_by(region)%>%
  #summarise(Hs_zero = sum(Hs == 0, na.rm = TRUE),
   #         Hs_greater_than_zero = sum (Hs > 0, na.rm = TRUE))
# used this to understand the Hs=0 and Hs>0 part. 

#PCA and Admixture 

## starting with 25 
options(bitmapType="cairo")
#read in vcf
vcf.25 <- read.vcfR("~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.25.vcf.gz") 
#thin the vcf
vcf.thin.25 <- distance_thin(vcf.25, min.distance = 500)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

meta.25 <- meta[meta$id %in% colnames(vcf.thin.25@gt[,-1]),]

dim(meta.25)

write.vcf(vcf.thin.25, "~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_thinned_hw1_.25.vcf.gz")

#uncompress the file 
system("gunzip -c ~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_thinned_hw1_.25.vcf.gz > ~/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.25.vcf")

geno25 <- vcf2geno(input.file= "/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.25.vcf",
                 output.file= "/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.25.geno")
CentPCA <- LEA::pca("~/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.25.geno", scale=TRUE)
plot(CentPCA$projections,
     col=as.factor(meta.25$region))
legend("bottomright", legend=as.factor(unique(meta.25$region)),
       fill=as.factor(unique(meta.25$region)))

CentPCA<- load.pcaProject("vcf_final_filtered_thinned_.25.pcaProject")

show(CentPCA)

plot(CentPCA)

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta.25$region, shape=meta.25$continent))+
  geom_point(alpha=0.5) +
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent") +
  xlim(-10,10) + 
  ylim(-10,10)
ggsave("~/projects/eco_genomics/population_genomics/HW1/Cent_PCA_PC1vPC2_HW1_25.pdf", width=6, height=6, units="in")

#admixture 
CentAdmix <- snmf("HW1/vcf_final_filtered_thinned_.25.geno",
                  K=1:10,
                  entropy=T,
                  repetitions =3,
                  project="new")

#compare evidence for different levels of K (or PCs) using the cross-entropy
par("mar")
par(mar=c(1,1,1,1))

par(mfrow=c(2,1)) #sets up multi-panel plot
plot(CentAdmix, col="blue4", main="SNMF")
plot(CentPCA$eigenvalues[1:10], ylab="Eigenvalues", xlab="Number of PCs", col="blue4",main="PCA")
dev.off()

#set k value 
myK=4

#calculate the cross entropy (=model fit; lower values are better)
CE = cross.entropy(CentAdmix, K=myK)
best = which.min(CE)

#extract the ancestry coefficients (the q scores) 
myKQ= Q(CentAdmix, K=myK, run=best)

#cbind to metadata
myKQmeta = cbind(myKQ, meta.25)

#set up color panel 
my.colors = c("blue4", "gold", "tomato", "lightblue", "olivedrab")

#sort the entire dataset by features of interest 
#first group by continent. then sort by region and pop within 
myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region,pop, .by_group = TRUE)

#lastly, make an ancestry plot and save it 
pdf("HW1/Admixture_K4_HW1_.25.pdf", width=10, height=5)
barplot(as.matrix(t(myKQmeta[,1:myK])),
        border=NA,
        space=0,
        col=my.colors[1:myK],
        xlab = "Geographic Regions",
        ylab = "Ancestry Proportions",
        main=paste0("Ancestry matrix K=",myK))
axis(1,
     at=1: length(myKQmeta$region),
     labels = myKQmeta$region,
     cex.axis=0.5,
     las=3)
dev.off()

## 50 
options(bitmapType="cairo")
#read in vcf
vcf.50 <- read.vcfR("~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.50.vcf.gz") 
#thin the vcf
vcf.thin.50 <- distance_thin(vcf.50, min.distance = 500)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

meta.50 <- meta[meta$id %in% colnames(vcf.thin.50@gt[,-1]),]

dim(meta.50)

write.vcf(vcf.thin.50, "~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_thinned_hw1_.50.vcf.gz")

#uncompress the file 
system("gunzip -c ~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_thinned_hw1_.50.vcf.gz > ~/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.50.vcf")

geno50 <- vcf2geno(input.file= "/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.50.vcf",
                 output.file= "/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.50.geno")
CentPCA <- LEA::pca("/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.50.geno", scale=TRUE)
plot(CentPCA$projections,
     col=as.factor(meta.50$region))
legend("bottomright", legend=as.factor(unique(meta.50$region)),
       fill=as.factor(unique(meta.50$region)))

CentPCA<- load.pcaProject("vcf_final_filtered_thinned_.50.pcaProject")

show(CentPCA)

plot(CentPCA)

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta.50$region, shape=meta.50$continent))+
  geom_point(alpha=0.5) +
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent") +
  xlim(-10,10) + 
  ylim(-10,10)
ggsave("~/projects/eco_genomics/population_genomics/HW1/Cent_PCA_PC1vPC2_HW1_50.pdf", width=6, height=6, units="in")

#admixture 
CentAdmix <- snmf("HW1/vcf_final_filtered_thinned_.50.geno",
                  K=1:10,
                  entropy=T,
                  repetitions =3,
                  project="new")

#compare evidence for different levels of K (or PCs) using the cross-entropy
par("mar")
par(mar=c(1,1,1,1))

par(mfrow=c(2,1)) #sets up multi-panel plot
plot(CentAdmix, col="blue4", main="SNMF")
plot(CentPCA$eigenvalues[1:10], ylab="Eigenvalues", xlab="Number of PCs", col="blue4",main="PCA")
dev.off()

#set k value 
myK=4

#calculate the cross entropy (=model fit; lower values are better)
CE = cross.entropy(CentAdmix, K=myK)
best = which.min(CE)

#extract the ancestry coefficients (the q scores) 
myKQ= Q(CentAdmix, K=myK, run=best)

#cbind to metadata
myKQmeta = cbind(myKQ, meta.50)

#set up color panel 
my.colors = c("blue4", "gold", "tomato", "lightblue", "olivedrab")

#sort the entire dataset by features of interest 
#first group by continent. then sort by region and pop within 
myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region,pop, .by_group = TRUE)

#lastly, make an ancestry plot and save it 
pdf("HW1/Admixture_K4_HW1_.50.pdf", width=10, height=5)
barplot(as.matrix(t(myKQmeta[,1:myK])),
        border=NA,
        space=0,
        col=my.colors[1:myK],
        xlab = "Geographic Regions",
        ylab = "Ancestry Proportions",
        main=paste0("Ancestry matrix K=",myK))
axis(1,
     at=1: length(myKQmeta$region),
     labels = myKQmeta$region,
     cex.axis=0.5,
     las=3)
dev.off()

## 75 
options(bitmapType="cairo")
#read in vcf
vcf.75 <- read.vcfR("~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.75.vcf.gz") 
#thin the vcf
vcf.thin.75 <- distance_thin(vcf.75, min.distance = 500)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

meta.75 <- meta[meta$id %in% colnames(vcf.thin.75@gt[,-1]),]

dim(meta.75)

write.vcf(vcf.thin.75, "~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_thinned_hw1_.75.vcf.gz")

#uncompress the file 
system("gunzip -c ~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_thinned_hw1_.75.vcf.gz > ~/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.75.vcf")

geno75 <- vcf2geno(input.file= "/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.75.vcf",
                   output.file= "/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.75.geno")
CentPCA <- LEA::pca("/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_final_filtered_thinned_.75.geno", scale=TRUE)
plot(CentPCA$projections,
     col=as.factor(meta.75$region))
legend("bottomright", legend=as.factor(unique(meta.75$region)),
       fill=as.factor(unique(meta.75$region)))

CentPCA<- load.pcaProject("vcf_final_filtered_thinned_.75.pcaProject")

show(CentPCA)

plot(CentPCA)

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta.75$region, shape=meta.75$continent))+
  geom_point(alpha=0.5) +
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent") +
  xlim(-10,10) + 
  ylim(-10,10)
ggsave("~/projects/eco_genomics/population_genomics/HW1/Cent_PCA_PC1vPC2_HW1_75.pdf", width=6, height=6, units="in")

#admixture 
CentAdmix <- snmf("HW1/vcf_final_filtered_thinned_.75.geno",
                  K=1:10,
                  entropy=T,
                  repetitions =3,
                  project="new")

#compare evidence for different levels of K (or PCs) using the cross-entropy
par("mar")
par(mar=c(1,1,1,1))

par(mfrow=c(2,1)) #sets up multi-panel plot
plot(CentAdmix, col="blue4", main="SNMF")
plot(CentPCA$eigenvalues[1:10], ylab="Eigenvalues", xlab="Number of PCs", col="blue4",main="PCA")
dev.off()

#set k value 
myK=4

#calculate the cross entropy (=model fit; lower values are better)
CE = cross.entropy(CentAdmix, K=myK)
best = which.min(CE)

#extract the ancestry coefficients (the q scores) 
myKQ= Q(CentAdmix, K=myK, run=best)

#cbind to metadata
myKQmeta = cbind(myKQ, meta.75)

#set up color panel 
my.colors = c("blue4", "gold", "tomato", "lightblue", "olivedrab")

#sort the entire dataset by features of interest 
#first group by continent. then sort by region and pop within 
myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region,pop, .by_group = TRUE)

#lastly, make an ancestry plot and save it 
pdf("HW1/Admixture_K4_HW1_.75.pdf", width=10, height=5)
barplot(as.matrix(t(myKQmeta[,1:myK])),
        border=NA,
        space=0,
        col=my.colors[1:myK],
        xlab = "Geographic Regions",
        ylab = "Ancestry Proportions",
        main=paste0("Ancestry matrix K=",myK))
axis(1,
     at=1: length(myKQmeta$region),
     labels = myKQmeta$region,
     cex.axis=0.5,
     las=3)
dev.off()

#selection 

## starting with 25 
system("gunzip -c ~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.25.vcf.gz > ~/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.25.vcf")

#read in uncompressed file 
vcf_25 <- read.pcadapt("/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.25.vcf",
                    type="vcf")
#read in compressed one too 
vcfR_25 <- read.vcfR("/gpfs1/home/s/a/santelo/projects/eco_genomics/population_genomics/HW1/vcf_filtered_hw1_.25.vcf.gz")

#read in metadata
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta_25 <- meta[meta$id %in% colnames(vcfR_25@gt[,1]),]

#run pcadapt 
pcadapt.pca_25 <- pcadapt(vcf_25,
                       K=2,
                       method="componentwise",
                       min.maf=0.01,
                       LD.clumping = list(size=500, thr=0.2))
summary(pcadapt.pca_25)

#quick manhattan plot 
plot(pcadapt.pca_25, K=2)

View(head(vcfR_25@fix))
vcfR_25.fix <- as.data.frame(vcfR_25@fix[,1:2])

chr.main_25 <- unique(vcfR_25.fix$CHROM)[1:8]

chrnum_25 <- as.data.frame(cbind(chr.main_25, seq(1,8,1)))

#extract p values 
Pval25 <- pcadapt.pca_25$pvalues 
pcadapt_25.MHplot <- cbind(vcfR_25.fix, Pval25) 
names (pcadapt_25.MHplot) [4:5] <- c("pPC1", "pPC2") #attribute must be the same length as vector

#create SNP variables for qqman
pcadapt_25.MHplot <- pcadapt_25.MHplot %>%
  mutate(SNP=paste0(chr.main_25, "_", POS))

#set variables as numeric that contain numbers
pcadapt_25.MHplot$V2 = as.numeric(pcadapt_25.MHplot$V2) #(*tmp*, V2, value=numeric(0))
pcadapt_25.MHplot$POS = as.numeric(pcadapt_25.MHplot$POS)
pcadapt_25.MHplot$pPC1 = as.numeric(pcadapt_25.MHplot[,4])
pcadapt_25.MHplot$pPC2 = as.numeric(pcadapt_25.MHplot[,5])


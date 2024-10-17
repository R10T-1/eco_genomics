library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics")

list.files()
list.files("variants/")

install.packages("ape")

library(vcfR) #program used to read the vcf files


vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz") 
vcf
head(vcf) #this is creating a list-summary of the data 
dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta") #reading in DNA
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="") #annotation data
chr1 <- create.chromR(name="Chromosome1", vcf=vcf, seq=dna, ann=gff) #we want to look at the data for chromosome 1,use the vcf file from before,using dna,the annotation data is in the gff

plot(chr1) #fun little graphing feature 

pdf(file="~/Projects/eco_genomics_2/Population Genomics/Figures/ChromoPlot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()
##this is where we left off last week## 

#using extract allows us to specify an element and have the program extract all the data from all the observations for that element

DP <- extract.gt(vcf, element="DP", as.numeric = T)

dim(DP)

DP[1:5,1:10]

quantile(DP)

DP[DP==0] <- NA
quantile(DP, na.rm=T)

#visualize the matrix of depth and missingness in our vcf file#

heatmap.bp(DP)
heatmap.bp(DP[1:1000,], rlabels = F, clabels = F)

#start filtering on depth
library(SNPfiltR)

vcf.filt <- hard_filter(vcf, depth=3) #setting the depth to 3 means that if you dont have a minimum of three reads at a snp it will be filtered as NA

vcf.filt <- max_depth(vcf.filt, maxdepth = 60) #filter out genotypes with >60 reads/snps

meta <- read.csv("metadata/meta4vcf.csv", header=T)

#for now, the snp filter wants a single column for ids and a single column for grouping filters. What we want to do is extract two columns, one of which will be a regional grouping, because theres less regions it will be easier to read#

meta2 <- meta[,c(1,4)]
head(meta2)

names(meta2) <- c("id","pop")

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta$pop)

#time to look at missingness, start w new object that looks at individual missingness

vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                     popmap=meta2, 
                                     cutoff=0.75)

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff = .5)

DP2 <- extract.gt(vcf.filt.indSNPMiss,
                  element="DP",
                  as.numeric=T)

heatmap.bp(DP2[5001:15000,],
           rlabels=F, clables=F)
write.vcf(vcf.filt.indSNPMiss,
          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

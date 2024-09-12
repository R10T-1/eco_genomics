library(vcfR)
install.packages("ape")
library(vcfR) #program used to read the vcf files

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz") 
head(vcf)
dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta") #reading in DNA
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="") #annotation data
chr1 <- create.chromR(name="Chromosome1", vcf=vcf, seq=dna, ann=gff) #we want to look at the data for chromosome 1,use the vcf file from before,using dna,the annotation data is in the gff

plot(chr1) #fun little graphing feature 

pdf(file="~/projects/eco_genomics/population_genomics/figures/ChromoPlot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()
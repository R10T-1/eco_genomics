library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType="cairo")

setwd("~/projects/eco_genomics/population_genomics/")

vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz") 

#we need to thin the SNPs for LD (linkage disequilibrium) before we run PCA and 
#Admixture analyses to satusfy the assumptions of indepedance among loci 

vcf.thin <- distance_thin(vcf, min.distance = 500)

# we went from 15454 to 3646, this tells us that alot of our SNPs are clustered
#in fragments, this is why when we cut down how many SNPs per base pair distance
#we see such a stark decrease in reads

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[,-1]),]

dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

#the next thing we wsnt to do is to uncompress the file, but store it outside of 
#github repo because the file is too big and will make github break

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

geno <- vcf2geno(input.file= "/gpfs1/home/s/a/santelo/vcf_final.filtered.thinned.vcf",
                 output.file= "outputs/vcf_final.filtered.thinned.geno")
CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)
plot(CentPCA$projections,
     col=as.factor(meta2$region))
legend("bottomright", legend=as.factor(unique(meta2$region)),
       fill=as.factor(unique(meta2$region)))

CentPCA<- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)

plot(CentPCA)
#screeplot, magnitude of PC values in decending
#shpwing the eigenvalue -> axis put through the dots 
#pca dots get smaller with each consecutive plot

ggplot(as.data.frame(CentPCA$projections),
      aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent))+
      geom_point(alpha=0.5) +
      labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent") +
      xlim(-10,10) + 
      ylim(-10,10)
ggsave("~/projects/eco_genomics/population_genomics/figures/Cent_PCA_PC1vPC2.pdf", width=6, height=6, units="in") #didnt work 

#run admixture analyses and create plots
#admixture, were going to use the LEA R package 
#function inside LEA s called snmf 
CentAdmix <- snmf("outputs/vcf_final.filtered.thinned.geno",
                  K=1:10,
                  entropy=T,
                  repetitions =3,
                  project="new")
# we can compare evidence fro different levels of K (or PCs) using the cross-entropy

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
myKQmeta = cbind(myKQ, meta2)

#set up color panel 
my.colors = c("blue4", "gold", "tomato", "lightblue", "olivedrab")

#sort the entire dataset by features of interest 
#first group by continent. then sort by region and pop within 
myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region,pop, .by_group = TRUE)

#lastly, make an ancestry plot and save it 
pdf("figures/Admixture_K4.pdf", width=10, height=5)
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

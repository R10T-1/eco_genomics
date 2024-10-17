#estimatng diversity and genomic differentiation in the filtered centaurea data 

library(vcfR)
library(tidyverse)
library(qqman)

#helps solve plotting issues 
X11.options(type="cairo")
options(bitmapType = "cairo")

#read in vcf file from oputputs directory 
vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

#read in metadata (info on population of origin, what region they come from, what continent, etc.)
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) #vcf files have 593
dim(meta) #meta has 629 inds

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),] 

dim(meta2)

#calculate diversity stats using the genetic_diff fxn in vcfR
vcf.div <- genetic_diff(vcf,
                        pops=as.factor(meta2$region),
                        method="nei")

#visualize the data 
str(vcf.div) #tells you the structure of the data frame as opposed to head showing you the first few lines

#data wrangling- if you peak into the first slot of vcf.div, there is an entry for every single snp, how do we get the unique values only?
unique(vcf.div$CHROM)

chr.main <- unique(vcf.div$CHROM)[1:8] #this lets us look at our "main" chromosomes aka the first 8 unique ones

chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))

#merge vcf.div with chrnum to be able to assign the numbers to all the rows in vcf.div - it says that we will use chrnum as the left frame an vcf.div as the right, and it is joining using the chr.main vairable from the left which is CHROM in the right frame
vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))

vcf.div.MHplot <- vcf.div.MHplot %>% 
                    filter(Gst>0) %>%
                    mutate(SNP=paste0(chr.main,"_",POS))

vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2) #this is to make sure it knows its a number
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)

#make manhattan plot 
manhattan(vcf.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst,0.999))

#manhattan plot is like THE ecological genomics plot, each dot is a SNP, Fst is the frequency of the SNPs along the chromosome among the populations
#the dip in the middle of the chromosomes on the plot represents the centromere 
#to plot the median use the 50th quantile
manhattan(vcf.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst,0.5))

write.csv(vcf.div.MHplot, "~/projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion.csv",
          quote=F,
          row.names = F)
#look at the Hs (diversity) values
names(vcf.div.MHplot)
#lets look at it as more of a histogram or density plot
vcf.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title="Genome-wide expected heterozygosity (Hs)",fill="Regions",
       x="Gene diversity within Regions", y="Counts of SNPs")

ggsave("Histogram_GenomDiversity_byRegion.pdf",
       path="~/projects/eco_genomics/population_genomics/firgures/")

vcf.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>% 
  filter(value!=0) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())
#sample size (N) = number of SNPs
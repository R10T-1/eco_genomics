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
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),] #getting rid of empty data points and adjusting pvalues
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

#calculate the number of unique genes in each portion of the euler plot
#we can do that based on the values above 
1935-107-44+23 #number of genes in A28 not overlapping=1807 genes diffirentially expressed uniquely at baseline v btwn 18 vs 22
296-107-29+23 #183= uniquely expressed when exposed to 28
78-44-29+23 #28 uniquely ex;ressed when exosed to 33

107-23 #84 unique to BASE & A28 
44-23 #21 unque to BASE & A33
29-23 #6 unique to A28 & A23

myEuler <- euler(c("BASE"=1807, "A28"=183, "A33"=28, "BASE&A28"=84, "BASE&A33"=21, "A28&A33"=6, "BASE&A28&A33"=23))

plot(myEuler, lty=1:3, quatities=TRUE)

###########################
#
#make a scatter plot of responses to A28/33 when copepods develop at 18 vs 22 

#contrast  D18_A28vsBASE
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast =c("group","D18BASE","D18A28"), alpha = 0.05))

#contrast D22_A28vsBASE
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast =c("group","D22BASE","D22A28"), alpha = 0.05))

#merge data frames 
res_df28 <- merge(res_D18_BASEvsA28, res_D22_BASEvsA28, by = "row.names", suffixes = c(".18", ".22"))
rownames(res_df28) <- res_df28$Row.names
res_df28 <- res_df28[, -1]

library(dplyr)
library(tidyr)

#define color mapping logic with mutate function 

res_df28 <- res_df28 %>%
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise2", 
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta1",
    padj.22 < 0.05 & stat.22 < 0 ~ "blue2",
    padj.22 < 0.05 & stat.22 > 0 ~ "red"
  ))

# Count the number of points per fill color 
color_counts <- res_df28 %>%
  group_by(fill) %>%
  summarise(count = n())

#create a data frame for the labels
label_positions <- data.frame(
  fill = c("blue2","magenta1", "red", "turquoise2"),
  x_pos = c(1,5,0,-7.5),
  y_pos = c(-5,0,9,3)
)

label_data28 <- merge(color_counts, label_positions, by = "fill")

#plot 
plot28 <- ggplot(res_df28, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  geom_text(data = label_data28, aes( x = x_pos, y = y_pos, label = count, color = fill),
            size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") +
  xlim(-10,10) + ylim(-10,10) +
  labs(x = "Log2FoldChange 28 vs BASE at 18",
       y = "Log2FoldChange 28 vs BASE at 22", 
       title = "How does response to 28 C vary by DevTemp?") +
  theme_minimal()
plot28

####repeat for A33
#make a scatter plot of responses to A28/33 when copepods develop at 18 vs 22 

#contrast  D18_A33vsBASE
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast =c("group","D18BASE","D18A33"), alpha = 0.05))

#contrast D22_A28vsBASE
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast =c("group","D22BASE","D22A33"), alpha = 0.05))

#merge data frames 
res_df33 <- merge(res_D18_BASEvsA33, res_D22_BASEvsA33, by = "row.names", suffixes = c(".18", ".22"))
rownames(res_df33) <- res_df33$Row.names
res_df33 <- res_df33[, -1]

library(dplyr)
library(tidyr)

#define color mapping logic with mutate function 

res_df33 <- res_df33 %>%
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "cadetblue", 
    padj.18 < 0.05 & stat.18 > 0 ~ "deeppink",
    padj.22 < 0.05 & stat.22 < 0 ~ "cornflowerblue",
    padj.22 < 0.05 & stat.22 > 0 ~ "firebrick2"
  ))

# Count the number of points per fill color 
color_counts <- res_df33 %>%
  group_by(fill) %>%
  summarise(count = n())

#create a data frame for the labels
label_positions <- data.frame(
  fill = c("deepskyblue","deeppink", "firebrick2", "cornflowerblue"),
  x_pos = c(1,5,0,-7.5),
  y_pos = c(-5,0,9,3)
)

label_data33 <- merge(color_counts, label_positions, by = "fill")

#plot 
plot33 <- ggplot(res_df33, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  geom_text(data = label_data33, aes( x = x_pos, y = y_pos, label = count, color = fill),
            size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") +
  labs(x = "Log2FoldChange 33 vs BASE at 18",
       y = "Log2FoldChange 33 vs BASE at 22", 
       title = "How does response to 33 C vary by DevTemp?") +
  theme_minimal()
plot33

#put the plots togethr into a two pnel plot 

library(gridExtra)
combined_plot <- grid.arrange(plot28, plot33, ncol = 2)

#for some reason, the labels get messed up sometimes and will have the same numbers on each figure, i found that if i run the whole thing as a chunk it works

#save plot 
ggsave("~/projects/eco_genomics/transcriptomics/figures/combined_scatter_plot.png", combined_plot, width = 12, height = 6)

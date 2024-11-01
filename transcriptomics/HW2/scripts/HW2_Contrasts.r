## requires the outputs from steps one and two##
library(eulerr)

#setting up my groups
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)
dim (dds) #[1] 119438     21
resultsNames(dds) #[1] "Intercept"               "group_D18A33_vs_D18A28" 
#[3] "group_D18BASE_vs_D18A28" "group_D22A28_vs_D18A28" 
#[5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

###run contrasts###

##18_BASEvs18_A28
HW2_res_D18_BASE_D18_A28 <- results(dds, contrast=c("group", "D18BASE", "D18A28"), alpha = 0.05)
HW2_res_D18_BASE_D18_A28 <- HW2_res_D18_BASE_D18_A28[!is.na(HW2_res_D18_BASE_D18_A28$padj),]
HW2_res_D18_BASE_D18_A28 <- HW2_res_D18_BASE_D18_A28[order(HW2_res_D18_BASE_D18_A28$padj),]
#i have no removed empty data points, adjusted the pvalues, and ordered them by their adjusted p values
summary(HW2_res_D18_BASE_D18_A28)
# out of 60738 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 33, 0.054%
# LFC < 0 (down)     : 62, 0.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

HW2_degs_D18_BASE_D18_A28 <- row.names(HW2_res_D18_BASE_D18_A28[HW2_res_D18_BASE_D18_A28$padj<0.05,]) #list of DEGs

plotMA(HW2_res_D18_BASE_D18_A28, ylim=c(-4,4))

##18_BASE vs 18_A33
HW2_res_D18_BASE_D18_A33 <- results(dds, contrast=c("group", "D18BASE", "D18A33"), alpha = 0.05)
HW2_res_D18_BASE_D18_A33 <- HW2_res_D18_BASE_D18_A33[!is.na(HW2_res_D18_BASE_D18_A33$padj),]
HW2_res_D18_BASE_D18_A33 <- HW2_res_D18_BASE_D18_A33[order(HW2_res_D18_BASE_D18_A33$padj),]

summary(HW2_res_D18_BASE_D18_A33)
# out of 54320 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 132, 0.24%
# LFC < 0 (down)     : 350, 0.64%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 5)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

HW2_degs_D18_BASE_D18_A33 <- row.names(HW2_res_D18_BASE_D18_A33[HW2_res_D18_BASE_D18_A33$padj<0.05,])

plotMA(HW2_res_D18_BASE_D18_A33, ylim=c(-4,4))

##D18 Euler Plot 
length(HW2_degs_D18_BASE_D18_A28) #95
length(HW2_degs_D18_BASE_D18_A33) #482

length(intersect(HW2_degs_D18_BASE_D18_A28, HW2_degs_D18_BASE_D18_A33)) #61
#because i am ding 2 comparisons, there is only one overlap(61), so subtract that or unique ones
#95-61 = 34
#482-61 = 421

HW2_euler_D18 <- euler(c("A28"=34, "A33"=421, "A28&A33"=61))
plot(HW2_euler_D18, lty=1:3, quantities=TRUE, fill = c("paleturquoise", "palevioletred1", "springgreen"))
##22_BASEvs18_A28
HW2_res_D22_BASE_D22_A28 <- results(dds, contrast=c("group", "D22BASE", "D22A28"), alpha = 0.05)
HW2_res_D22_BASE_D22_A28 <- HW2_res_D22_BASE_D22_A28[!is.na(HW2_res_D22_BASE_D22_A28$padj),]
HW2_res_D22_BASE_D22_A28 <- HW2_res_D22_BASE_D22_A28[order(HW2_res_D22_BASE_D22_A28$padj),]
#i have no removed empty data points, adjusted the pvalues, and ordered them by their adjusted p values
summary(HW2_res_D22_BASE_D22_A28)
# out of 49523 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 689, 1.4%
# LFC < 0 (down)     : 49, 0.099%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 6)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
HW2_degs_D22_BASE_D22_A28 <- row.names(HW2_res_D22_BASE_D22_A28[HW2_res_D22_BASE_D22_A28$padj<0.05,]) #list of DEGs

plotMA(HW2_res_D22_BASE_D22_A28, ylim=c(-4,4))


##22_BASE vs 18_A33
HW2_res_D22_BASE_D22_A33 <- results(dds, contrast=c("group", "D22BASE", "D22A33"), alpha = 0.05)
HW2_res_D22_BASE_D22_A33 <- HW2_res_D22_BASE_D22_A33[!is.na(HW2_res_D22_BASE_D22_A33$padj),]
HW2_res_D22_BASE_D22_A33 <- HW2_res_D22_BASE_D22_A33[order(HW2_res_D22_BASE_D22_A33$padj),]

summary(HW2_res_D22_BASE_D22_A33)
# out of 46319 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1094, 2.4%
# LFC < 0 (down)     : 481, 1%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 8)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

HW2_degs_D22_BASE_D22_A33 <- row.names(HW2_res_D22_BASE_D22_A33[HW2_res_D22_BASE_D22_A33$padj<0.05,])

plotMA(HW2_res_D22_BASE_D22_A33, ylim=c(-4,4))

##D22 Euler Plot 
length(HW2_degs_D22_BASE_D22_A28) #738
length(HW2_degs_D22_BASE_D22_A33) #1575

length(intersect(HW2_degs_D22_BASE_D22_A28, HW2_degs_D22_BASE_D22_A33)) #196
#because i am ding 2 comparisons, there is only one overlap(61), so subtract that or unique ones
# 738-196 = 542
# 1575-196 = 1379
HW2_euler_D22 <- euler(c("A28"=542, "A33"=1379, "A28&A33"=196))
plot(HW2_euler_D22, lty=1:3, quantities=TRUE, fill = c("cyan", "maroon1", "palegreen"))

###make scatterplot###

##contrast D18 base v A28 & A33
HW2_res_D18_BASE_D18_A28 <- as.data.frame(results(dds, contrast =c("group", "D18BASE", "D18A28"), alpha = 0.05))
HW2_res_D18_BASE_D18_A33 <- as.data.frame(results(dds, contrast =c("group", "D18BASE", "D18A33"), alpha = 0.05))
#merge
HW2_res_d18 <- merge(HW2_res_D18_BASE_D18_A28, HW2_res_D18_BASE_D18_A33, by = "row.names", suffixes = c(".28", ".33"))
rownames(HW2_res_d18) <- HW2_res_d18$Row.names
HW2_res_d18 <- HW2_res_d18[, -1]

library(dplyr)
library(tidyr)

#define color mapping logic 
HW2_res_d18 <- HW2_res_d18 %>%
  mutate(fill=case_when(
  padj.28 < 0.05 & stat.28 < 0 ~ "orchid1",
  padj.28 < 0.05 & stat.28 > 0 ~ "cornflowerblue",
  padj.33 < 0.05 & stat.33 < 0 ~ "palegreen2",
  padj.33 < 0.05 & stat.33 > 0 ~ "aquamarine"
  ))

#count number of points per fill color
color_counts <- HW2_res_d18 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions <- data.frame(
  fill = c("deepskyblue", "forestgreen", "red2", "orange"), 
  x_pos = c(1,5,0,-7.5),
  y_pos = c(-5,0,9,3))

label_data18 <- merge(color_counts, label_positions)

#plot 
plot18 <- ggplot(HW2_res_d18, aes(x = log2FoldChange.28, y = log2FoldChange.33, collor = fill)) +
  geom_point(alpha = 0.08) +
  scale_color_identity() +
  geom_text(data = label_data18, aes(x = x_pos, y = y_pos, label = count, color = fill), size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") +
  xlim(-10,10) + ylim(-10,10) +
  labs(x = "Log2FoldChange D18 BASE vs D18 A28",
       y = "log2FoldChange D18 BASE vs D18 A33",
       title = "DevTemp 18 C respond to 28 C and 33 C?") +
  theme_minimal()

plot18
ggsave("~/projects/eco_genomics/transcriptomics/HW2/figures/HW2_D18_scatterplot.pdf")

#repeat for D22

###make scatterplot###

##contrast D22 base v A28 & A33
HW2_res_D22_BASE_D22_A28 <- as.data.frame(results(dds, contrast =c("group", "D22BASE", "D22A28"), alpha = 0.05))
HW2_res_D22_BASE_D22_A33 <- as.data.frame(results(dds, contrast =c("group", "D22BASE", "D22A33"), alpha = 0.05))
#merge
HW2_res_d22 <- merge(HW2_res_D22_BASE_D22_A28, HW2_res_D22_BASE_D22_A33, by = "row.names", suffixes = c(".28", ".33"))
rownames(HW2_res_d22) <- HW2_res_d22$Row.names
HW2_res_d22 <- HW2_res_d22[, -1]

library(dplyr)
library(tidyr)

#define color mapping logic 
HW2_res_d22 <- HW2_res_d22 %>%
  mutate(fill=case_when(
    padj.28 < 0.05 & stat.28 < 0 ~ "orchid1",
    padj.28 < 0.05 & stat.28 > 0 ~ "cornflowerblue",
    padj.33 < 0.05 & stat.33 < 0 ~ "palegreen2",
    padj.33 < 0.05 & stat.33 > 0 ~ "aquamarine"
  ))

#count number of points per fill color
color_counts <- HW2_res_d22 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions <- data.frame(
  fill = c("deepskyblue", "forestgreen", "red2", "orange"), 
  x_pos = c(1,5,0,-7.5),
  y_pos = c(-5,0,9,3))

label_data22 <- merge(color_counts, label_positions)

#plot 
plot22 <- ggplot(HW2_res_d22, aes(x = log2FoldChange.28, y = log2FoldChange.33, collor = fill)) +
  geom_point(alpha = 0.08) +
  scale_color_identity() +
  geom_text(data = label_data22, aes(x = x_pos, y = y_pos, label = count, color = fill), size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") +
  xlim(-10,10) + ylim(-10,10) +
  labs(x = "Log2FoldChange D22 BASE vs D22 A28",
       y = "log2FoldChange D22 BASE vs D22 A33",
       title = "DevTemp 22 C response to 28 C and 33 C?") +
  theme_minimal()

plot22
ggsave("~/projects/eco_genomics/transcriptomics/HW2/figures/HW2_D22_scatterplot.pdf")

#combine plots 

library(gridExtra)
combined_plot_hw2 <- grid.arrange(plot18, plot22, ncol = 2)
#couldnt use ggsave so i did export 

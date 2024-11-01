###picking up from DESEq2, have WD and packages from DESeq2 loaded###
library(pheatmap)

options(bitmapType = "cairo")

#bring in the rsult names from deseq object
resultsNames(dds) 
#[1] "Intercept" "DevTemp_D22_vs_D18"  "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

#mnake a results file

results_D22vsD18_HW2 <- results(dds, name="DevTemp_D22_vs_D18", alpha=0.05)

#order by significance 
results_D22vsD18_HW2 <- results_D22vsD18_HW2[order(results_D22vsD18_HW2$padj),]
head(results_D22vsD18_HW2) #first few genes in results, ordered by significance of adjusted pvalue -> significant difference of expression

#look at counts of a specific top gene that we're interested in to validate that the model is working 
d <- plotCounts(dds, gene="TRINITY_DN140854_c0_g5_i2", int=(c("DevTemp", "FinalTemp")), returnData = TRUE)
d 
#make the plot 
p <- ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp)) + 
  theme_minimal() + theme(text=element_text(size=20), panel.grid.major=element_line(colour = "grey"))

p  <- p + geom_point(position=position_jitter(w=0.2,h=0),size=3)
p
#the direction was -1.1 for this gene, and we compared from 22 to 18, which tells us that expression decreased from 22 to 18

#MA plot, similar to manhattan plot. m=logfoldchange A=average. so its the expression change vs the average 
MAplot <- plotMA(results_D22vsD18_HW2, ylim=c(-4,4))
#blue dots are the more significant 
#we saw a lot of upregulation in development at 22 compared to 18

#volcano plot 
#convert our deseq results object into a data frame 
res_df_HW <- as.data.frame(results_D22vsD18_HW2)

#add a column to denote if a gene is being differentially expressed or not
res_df_HW$Significant <- ifelse(res_df_HW$padj < 0.05 & abs(res_df_HW$log2FoldChange) > 1, "Significant", "Not Significant")

#plot 
ggplot(res_df_HW, aes(x = log2FoldChange, y = -log10(padj), color = Significant))+
  geom_point(alpha = 0.8)+
  scale_color_manual(values = c('darkseagreen',"lightcyan3"))+
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", title = "Volcano Plot")+
  theme_minimal()+
  theme(legend.position = "top")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "plum4")+ 
  geom_vline(xintercept = c(-1,1),linetype = "dashed", color = "plum4")
# ggsave("~/projects/eco_genomics/transcriptomics/HW2/figures/volcanoplot_HW2.pdf")
 
#another way to look at the expression data is through a heatmap 
vsd <- vst(dds, blind=FALSE)

topgenes <- head(rownames(results_D22vsD18_HW2), 20) #selecting the top 20 genes
mat <- assay(vsd)[topgenes,] #making matrix of all the gene expression data across the dataset 
df <- as.data.frame(colData(dds)[,c("DevTemp","FinalTemp")])

library(RColorBrewer) #i want to make my heatmap pretty 

pallete <- colorRampPalette(brewer.pal(11,"PRGn")) (100)

pheatmap(mat, annotation_col=df, show_rownames=FALSE, cluster_cols=T, cluster_rows=T, color = pallete, filename = "~/projects/eco_genomics/transcriptomics/HW2/figures/heatmap_HW.pdf")
#columns are samples, rows are genes, color correlates with the magnitude of expression 

dev.off()

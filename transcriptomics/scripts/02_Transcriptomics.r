#load in transcriptomics day one data#

#make sure that DESeq2 and ggplot are loaded, load if not
library(pheatmap)

#this is to help us make a heatmap

options(bitmapType = "cairo")

#bring in the rsult names from deseq object
resultsNames(dds) 
#[1] "Intercept" "DevTemp_D22_vs_D18"  "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

#take the results from the one we are interested in and put them in a results file

res_D22vsD18 <- results(dds, name="DevTemp_D22_vs_D18", alpha=0.05) #create the object so we can look into it 

#order by significance 
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
head(res_D22vsD18) #first few genes in results, ordered by significance of adjusted pvalue -> significant difference of expression

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
plotMA(res_D22vsD18, ylim=c(-4,4))
#blue dots are the more significant 
#we saw a lot of upregulation in development at 22 compared to 18

#volcano plot 
#convert our deseq results object into a data frame 
res_df <- as.data.frame(res_D22vsD18)

#add a column to denote if a gene is being differentially expressed or not
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

#plot 
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant))+
  geom_point(alpha = 0.8)+
  scale_color_manual(values = c('slateblue',"tomato"))+
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", title = "Volcano Plot")+
  theme_minimal()+
  theme(legend.position = "top")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange")+ 
  geom_vline(xintercept = c(-1,1),linetype = "dashed", color = "orange")

#another way to look at the expression data is through a heatmap 
vsd <- vst(dds, blind=FALSE)

topgenes <- head(rownames(res_D22vsD18), 20) #selecting the top 20 genes
mat <- assay(vsd)[topgenes,] #making matrix of all the gene expression data across the dataset 
df <- as.data.frame(colData(dds)[,c("DevTemp","FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=FALSE, cluster_cols=T, cluster_rows=T)
#columns are samples, rows are genes, color correlates with the magnitude of expression 


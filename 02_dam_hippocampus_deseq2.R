#Gestational bisphenol project
#Hannah Lapp PhD
#Champagne Lab, Fall 2021

#Dam hippocampus deseq2 analysis and figures

library(tidyverse)
library("DESeq2")
library(pheatmap)
library(glue)
source("./analysis/DeSeq2_functions.R") #function: get_LFC_results 

#-----------------Load data -
brain_region <- "Hippocampus"
brain_region_short <- "hipp"
group <- "dam"
cts <- read.csv("./data/dam_hipp_raw.csv", row.names=1, check.names=FALSE)
coldata <- read.csv("./data/dam_condition_lg_data.csv", row.names=1)

#-----------------DEseq2 object--
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition + lg_centered)

keep <- rowSums(counts(dds)) > 10 #Prefiltering: remove rows that have less than 10 reads total 
dds <- dds[keep,]
dds2 <- estimateSizeFactors(dds)
filt_norm_counts <- counts(dds2, normalized = TRUE) # df with filtered normalized counts, use for pheatmap

ddsMF <-  dds

ddsMF$condition <- factor(ddsMF$condition, levels = c("vehicle", "mixed")) #set condition levels as factors 
design(ddsMF) <-  formula(~ lg_centered + condition)
dds <- DESeq(ddsMF)  

#----------------Create gene lists from deseq2 results object (table 1)
#mix vs vehicle************************
res_con <- results(dds, contrast = c("condition", "mixed", "vehicle"))
res_con_p <- subset(res_con, pvalue <0.05)
res_con_padj <-  subset(res_con, padj <0.05) 
cat(paste("DEGs Mixed vs. Vehicle in", brain_region,
          "\nunadj p ", nrow(filter(as.data.frame(res_con), pvalue <.05 & log2FoldChange > 0)), "upregualted genes", 
          "\nunadj p ", nrow(filter(as.data.frame(res_con), pvalue <.05 & log2FoldChange < 0)), "downregulated genes", 
          "\nadj p ", nrow(filter(as.data.frame(res_con), padj <.05 & log2FoldChange > 0)), "upregualted genes",
          "\nadj p ", nrow(filter(as.data.frame(res_con), padj <.05 & log2FoldChange < 0)), "downregulated genes"))

#lg attendance***********************
res_lg <-  results(dds, name = "lg_centered") #pulling out results for DEG for lg attendance
res_lg_p <- subset(res_lg, pvalue <.05)
res_lg_padj <-  subset(res_lg, padj <0.05) 
cat(paste("DEGs associated with L/G in", brain_region,
          "\nunadj p ", nrow(filter(as.data.frame(res_lg), pvalue <.05 & log2FoldChange > 0)), "upregualted genes", 
          "\nunadj p ", nrow(filter(as.data.frame(res_lg), pvalue <.05 & log2FoldChange < 0)), "downregulated genes", 
          "\nadj p ", nrow(filter(as.data.frame(res_lg), padj <.05 & log2FoldChange > 0)), "upregualted genes",
          "\nadj p ", nrow(filter(as.data.frame(res_lg), padj <.05 & log2FoldChange < 0)), "downregulated genes"))

#------------------Save gene lists --
file_folder <- "./results/"

#all gene results for BP vs. Vehicle
all_genes_bp <- as.data.frame(res_con)
all_genes_file1 <- paste(file_folder, group,"-", brain_region_short,"-", "all-results-BP",".csv", sep = "")
write.csv(all_genes_bp, all_genes_file1)

#all gene results for LG
all_genes_lg <- as.data.frame(res_lg)
all_genes_file2 <- paste(file_folder, group,"-", brain_region_short,"-", "all-results-LG",".csv", sep = "")
write.csv(all_genes_lg, all_genes_file2)

# Create pheatmap with top 100 most variable genes- (fig 2)
sd_counts <- transform(filt_norm_counts, stdev = apply(filt_norm_counts, 1, sd)) #find standard dev by row (gene)
colnames(sd_counts) <- str_replace((colnames(sd_counts)), "X", "")
counts <-  sd_counts[order(-sd_counts$stdev),] #put in order by decesnding stardard dev
topvar <- counts[1:100,]  #take top 100
topvar <- subset( topvar, select = -stdev) # remove column with standard dev.

exposures <- coldata %>%  dplyr::select(condition, lg_centered)
colnames(exposures) <- c("Condition",  "Licking/grooming")

cal_z_score <- function(x){ #to scale heatmap
  (x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(topvar, 1, cal_z_score))

ann_colors = list(
  Condition = c(mixed = "coral2", vehicle = "darkslategray2"),
  
  `Licking/grooming` = c("darkorchid4", "plum1"))

h <- pheatmap(data_subset_norm, 
              main = brain_region,
              annotation_col = exposures,
              show_rownames = FALSE,
              fontsize_row = 6, 
              treeheight_row = 5, # dendrogram is hidden
              show_colnames = FALSE, 
              treeheight_col = 1,
              border_color = NA,
              annotation_colors = ann_colors,
              legend = FALSE,
              annotation_legend = FALSE,
              fontsize = 9,
              col = colorRampPalette(c("navy", "white", "darkorange2"))(50)) 
h
ggsave(h, filename = paste(brain_region, "_pheatmap _dam.png", sep = ""), path = "./figures_tables/", dpi= 600, height = 5, width = 2 )

#-------------PCA plot -
vsd <- vst(dds, blind=FALSE)
z <- plotPCA(vsd, intgroup=c( "condition"))
pca_plot <- z + 
  ggtitle(paste(brain_region)) +
  theme_bw(base_size = 9)+
  coord_fixed(ratio = 1.7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        legend.text = element_text(size=3),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid"),
  )
pca_plot
ggsave(pca_plot, filename = paste(brain_region, "_pca_dam.png", sep = ""), path = "./figures_tables/", dpi= 600, height = 2, width = 2 )

#--------------QQ plot--
con <- as.data.frame(res_con$pvalue)
rank_con <- con %>% arrange(res_con$pvalue)
lg <- as.data.frame(res_lg$pvalue)
rank_lg <- lg %>% arrange(res_lg$pvalue)
nrow(rank_con) == nrow(rank_lg) #check is True
num_row <- as.numeric(nrow(rank_con))
exp = -log10((1:num_row)/(num_row+1))
df <-  cbind( expd = exp,
              LG = -log10(rank_lg$`res_lg$pvalue`),
              BP = -log10(rank_con$`res_con$pvalue`)) 
df <- as.data.frame(df)
df1 <- df %>% gather(key, value, 2:3) %>% 
  mutate(key= factor(key, levels = c("LG", "BP")))

d <- df1 %>% 
  ggplot(aes(expd, value, color = key))+
  geom_point(shape = 21, size = 1, alpha = .7)+
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype ="dashed")+
  theme_bw(base_size = 9)+
  labs(x = expression(-log[10]~p-values~(expected)),
       y = expression(-log[10]~p-values~(observed)),
       title = brain_region,
       color = "",
       fill = "")+
  theme(legend.position = c(0.85,0.15),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size=5),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("darkorchid4","goldenrod2" ))
d
ggsave(d, filename = paste(brain_region, "_QQ_dam.png", sep = ""), path = "./figures_tables/", dpi= 600, height = 2, width = 2 )










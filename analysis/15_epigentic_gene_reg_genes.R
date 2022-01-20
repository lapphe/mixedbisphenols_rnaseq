#Gestational bisphenol project
#Hannah Lapp PhD
#Champagne Lab, Fall 2021

#hypothesis-driven analysis for epigentic regualtion of gene expression gene list

library(tidyverse)
library(msigdbr)
library(fgsea)
library(purrr)

C5 <- msigdbr(species = "Rattus norvegicus", category = "C5", subcategory = "BP") 
C5_E <- C5 %>% filter(gs_name == "GOBP_REGULATION_OF_GENE_EXPRESSION_EPIGENETIC") %>% 
  dplyr::select(gene_symbol,ensembl_gene)

format_res <- function(region, comp, age){
  path <- paste("./results/pup-", region, "-all-results-", comp, ".csv", sep = "")
  
  df<- read.csv(path)
  df <- dplyr::rename(df, "ensembl_gene" = "X")
  df1 <- left_join(C5_E, df) %>% 
    mutate(region = region) %>% 
    mutate(comparison = comp)
  df2 <- df1[!is.na(df1$padj),]
  df3 <- df2 %>% filter(pvalue < .05) #filter out to only unadj p<.05
  
  df4 <- df3 %>% 
    dplyr::select(gene_symbol, log2FoldChange)
  colnames(df4) <- c("gene_symbol", paste(age, region, comp, sep = "_"))
  
  return(df4)
}

p_amy_bp <- format_res("amy", "BP" , "pup")
p_hipp_bp <- format_res("hipp", "BP", "pup")
p_hypo_bp <- format_res("hypo", "BP", "pup")
p_nac_bp <- format_res("nac", "BP", "pup")
p_pl_bp <- format_res("pl", "BP", "pup")

p_amy_lg <- format_res("amy", "LG", "pup")
p_hipp_lg <- format_res("hipp", "LG", "pup")
p_hypo_lg <- format_res("hypo", "LG", "pup")
p_nac_lg <- format_res("nac", "LG", "pup")
p_pl_lg <- format_res("pl", "LG", "pup")

format_res_dam <- function(region, comp, age){
  path <- paste("./results/dam-", region, "-all-results-", comp, ".csv", sep = "")
  
  df<- read.csv(path)
  df <- dplyr::rename(df, "ensembl_gene" ="X")
  df1 <- left_join(C5_E, df) %>% 
    mutate(region = region) %>% 
    mutate(comparison = comp)
  df2 <- df1[!is.na(df1$padj),]
  df3 <- df2 %>% filter(pvalue < .05)
  
  df4 <- df3 %>% 
    dplyr::select(gene_symbol, log2FoldChange)
  colnames(df4) <- c("gene_symbol", paste(age, region, comp, sep = "_"))
  
  return(df4)
}

d_amy_bp <- format_res_dam("amy", "BP" , "dam")
d_hipp_bp <- format_res_dam("hipp", "BP", "dam")
d_hypo_bp <- format_res_dam("hypo", "BP", "dam")
d_nac_bp <- format_res_dam("nac", "BP", "dam")
d_pl_bp <- format_res_dam("pl", "BP", "dam")

d_amy_lg <- format_res_dam("amy", "LG", "dam")
d_hipp_lg <- format_res_dam("hipp", "LG", "dam")
d_hypo_lg <- format_res_dam("hypo", "LG", "dam")
d_nac_lg <- format_res_dam("nac", "LG", "dam")
d_pl_lg <- format_res_dam("pl", "LG", "dam")

all_data <- list(p_amy_bp,
                 p_hipp_bp,
                 p_hypo_bp, 
                 p_nac_bp, 
                 p_pl_bp,
                 p_amy_lg, 
                 p_hipp_lg,
                 p_hypo_lg,
                 p_nac_lg, 
                 p_pl_lg,
                 d_amy_bp,
                 d_hipp_bp,
                 d_hypo_bp, 
                 d_nac_bp, 
                 d_pl_bp,
                 d_amy_lg, 
                 d_hipp_lg,
                 d_hypo_lg,
                 d_nac_lg, 
                 d_pl_lg) %>% 
  purrr::reduce(full_join, by = "gene_symbol")

all_data1 <- as.data.frame(all_data[,-1])
rownames(all_data1) <- all_data$gene_symbol

all_data_long <- all_data %>% pivot_longer(cols = !gene_symbol, names_to = "dataset", values_to = "log2FC")
all_long_filt <- all_data_long %>% filter(!is.na(log2FC)) %>% 
  mutate(Age = str_extract(dataset, "pup")) %>% 
  mutate(Age = ifelse(is.na(Age), "dam", Age)) %>% 
  mutate(Comparison = str_extract(dataset, "BP")) %>% 
  mutate(Comparison = ifelse( Comparison == "BP", "Bisphenol", Comparison)) %>% 
  mutate(Comparison = ifelse(is.na(Comparison), "Licking/grooming", Comparison)) %>% 
  mutate(Region = str_extract(dataset, "amy")) %>% 
  mutate(Region = ifelse(is.na(Region), str_extract(dataset, "hipp"), Region)) %>% 
  mutate(Region = ifelse(is.na(Region), str_extract(dataset, "hypo"), Region)) %>% 
  mutate(Region = ifelse(is.na(Region), str_extract(dataset, "nac"), Region)) %>% 
  mutate(Region = ifelse(is.na(Region), str_extract(dataset, "pl"), Region)) %>% 
  mutate(Region = str_to_sentence(Region)) %>% 
  mutate(Age = str_to_sentence(Age)) %>% 
  select(dataset, Age, Comparison, Region, gene_symbol, log2FC)

all_long_filt %>% 
  ggplot(aes(x = gene_symbol, y = log2FC, label = log2FC))+
  geom_segment(aes(y=0, 
                   yend = `log2FC`,
                   x= `gene_symbol`,
                   xend= `gene_symbol`, 
                   color = `Comparison`), 
               size = 1, 
               alpha = .7)+
  geom_point(stat= "identity", aes(col = Comparison, shape = `Region`), size = 5, alpha = .7, stroke = 1.5)+
  scale_shape_manual(values=c(0,1,2,6,5))+
  scale_color_manual(values= c("darkorange", "blueviolet"))+
  geom_hline(yintercept= 0)+
  coord_flip()+
  ylab("Log2 fold change")+
  xlab("")+
  ggtitle("Epigentic gene regulation genes")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size=8),
        legend.background = element_rect(color = "black", linetype = "solid"),
        plot.title = element_text(size=13, hjust = 0.5),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        axis.text = element_text(size = 11))+
  facet_grid(~ `Age`)

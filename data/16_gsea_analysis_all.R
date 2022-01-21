#Gestational bisphenol project
#Hannah Lapp PhD
#Champagne Lab, Fall 2021

#gsea all analysis

library(msigdbr)
library(fgsea)
library(tidyverse)

C5 <- msigdbr(species = "Rattus norvegicus", category = "C5", subcategory = "BP") 

C5.ensembl.ls <-  C5 %>% #Biological processes gene list
  dplyr::select(gs_name, ensembl_gene) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(ensembl_gene))) %>% 
  deframe()

#load and format deseq2 results for each dataset----
pup_amy_bp.res <- read.csv("./results/pup-amy-all-results-BP.csv")
pup_amy_bp.arrange <- pup_amy_bp.res %>% arrange(desc(log2FoldChange))
pup_amy_bp.v <- pup_amy_bp.arrange$log2FoldChange
names(pup_amy_bp.v) <- pup_amy_bp.arrange$X
head(pup_amy_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_amy_bp.v) #two tailed
rm(pup_amy_bp.res, pup_amy_bp.arrange)

pup_amy_lg.res <- read.csv("./results/pup-amy-all-results-LG.csv")
pup_amy_lg.arrange <- pup_amy_lg.res %>% arrange(desc(log2FoldChange))
pup_amy_lg.v <- pup_amy_lg.arrange$log2FoldChange
names(pup_amy_lg.v) <- pup_amy_lg.arrange$X
head(pup_amy_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_amy_lg.v) #two tailed
rm(pup_amy_lg.arrange, pup_amy_lg.res)

pup_hipp_bp.res <- read.csv("./results/pup-hipp-all-results-BP.csv")
pup_hipp_bp.arrange <- pup_hipp_bp.res %>% arrange(desc(log2FoldChange))
pup_hipp_bp.v <- pup_hipp_bp.arrange$log2FoldChange
names(pup_hipp_bp.v) <- pup_hipp_bp.arrange$X
head(pup_hipp_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_hipp_bp.v) #two tailed
rm(pup_hipp_bp.arrange, pup_hipp_bp.res)

pup_hipp_lg.res <- read.csv("./results/pup-hipp-all-results-LG.csv")
pup_hipp_lg.arrange <- pup_hipp_lg.res %>% arrange(desc(log2FoldChange))
pup_hipp_lg.v <- pup_hipp_lg.arrange$log2FoldChange
names(pup_hipp_lg.v) <- pup_hipp_lg.arrange$X
head(pup_hipp_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_hipp_lg.v) #two tailed
rm(pup_hipp_lg.arrange, pup_hipp_lg.res)

pup_hypo_bp.res <- read.csv("./results/pup-hypo-all-results-BP.csv")
pup_hypo_bp.arrange <- pup_hypo_bp.res %>% arrange(desc(log2FoldChange))
pup_hypo_bp.v <- pup_hypo_bp.arrange$log2FoldChange
names(pup_hypo_bp.v) <- pup_hypo_bp.arrange$X
head(pup_hypo_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_hypo_bp.v) #two tailed
rm(pup_hypo_bp.arrange, pup_hypo_bp.res)

pup_hypo_lg.res <- read.csv("./results/pup-hypo-all-results-LG.csv")
pup_hypo_lg.arrange <- pup_hypo_lg.res %>% arrange(desc(log2FoldChange))
pup_hypo_lg.v <- pup_hypo_lg.arrange$log2FoldChange
names(pup_hypo_lg.v) <- pup_hypo_lg.arrange$X
head(pup_hypo_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_hypo_lg.v) #two tailed
rm(pup_hypo_lg.arrange, pup_hypo_lg.res)

pup_nac_bp.res <- read.csv("./results/pup-nac-all-results-BP.csv")
pup_nac_bp.arrange <- pup_nac_bp.res %>% arrange(desc(log2FoldChange))
pup_nac_bp.v <- pup_nac_bp.arrange$log2FoldChange
names(pup_nac_bp.v) <- pup_nac_bp.arrange$X
head(pup_nac_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_nac_bp.v) #two tailed
rm(pup_nac_bp.arrange, pup_nac_bp.res)

pup_nac_lg.res <- read.csv("./results/pup-nac-all-results-LG.csv")
pup_nac_lg.arrange <- pup_nac_lg.res %>% arrange(desc(log2FoldChange))
pup_nac_lg.v <- pup_nac_lg.arrange$log2FoldChange
names(pup_nac_lg.v) <- pup_nac_lg.arrange$X
head(pup_nac_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_nac_lg.v) #two tailed
rm(pup_nac_lg.res, pup_nac_lg.arrange)

pup_pl_bp.res <- read.csv("./results/pup-pl-all-results-BP.csv")
pup_pl_bp.arrange <- pup_pl_bp.res %>% arrange(desc(log2FoldChange))
pup_pl_bp.v <- pup_pl_bp.arrange$log2FoldChange
names(pup_pl_bp.v) <- pup_pl_bp.arrange$X
head(pup_pl_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_pl_bp.v) #two tailed
rm(pup_pl_bp.arrange, pup_pl_bp.res)

pup_pl_lg.res <- read.csv("./results/pup-pl-all-results-LG.csv")
pup_pl_lg.arrange <- pup_pl_lg.res %>% arrange(desc(log2FoldChange))
pup_pl_lg.v <- pup_pl_lg.arrange$log2FoldChange
names(pup_pl_lg.v) <- pup_pl_lg.arrange$X
head(pup_pl_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(pup_pl_lg.v) #two tailed
rm(pup_pl_lg.arrange, pup_pl_lg.res)

#dam datasets
dam_amy_bp.res <- read.csv("./results/dam-amy-all-results-BP.csv")
dam_amy_bp.arrange <- dam_amy_bp.res %>% arrange(desc(log2FoldChange))
dam_amy_bp.v <- dam_amy_bp.arrange$log2FoldChange
names(dam_amy_bp.v) <- dam_amy_bp.arrange$X
head(dam_amy_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_amy_bp.v) #two tailed
rm(dam_amy_bp.arrange, dam_amy_bp.res)

dam_amy_lg.res <- read.csv("./results/dam-amy-all-results-LG.csv")
dam_amy_lg.arrange <- dam_amy_lg.res %>% arrange(desc(log2FoldChange))
dam_amy_lg.v <- dam_amy_lg.arrange$log2FoldChange
names(dam_amy_lg.v) <- dam_amy_lg.arrange$X
head(dam_amy_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_amy_lg.v) #two tailed
rm(dam_amy_lg.arrange, dam_amy_lg.res)

dam_hipp_bp.res <- read.csv("./results/dam-hipp-all-results-BP.csv")
dam_hipp_bp.arrange <- dam_hipp_bp.res %>% arrange(desc(log2FoldChange))
dam_hipp_bp.v <- dam_hipp_bp.arrange$log2FoldChange
names(dam_hipp_bp.v) <- dam_hipp_bp.arrange$X
head(dam_hipp_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_hipp_bp.v) #two tailed
rm(dam_hipp_bp.arrange, dam_hipp_bp.res)

dam_hipp_lg.res <- read.csv("./results/dam-hipp-all-results-LG.csv")
dam_hipp_lg.arrange <- dam_hipp_lg.res %>% arrange(desc(log2FoldChange))
dam_hipp_lg.v <- dam_hipp_lg.arrange$log2FoldChange
names(dam_hipp_lg.v) <- dam_hipp_lg.arrange$X
head(dam_hipp_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_hipp_lg.v) #two tailed
rm(dam_hipp_lg.arrange, dam_hipp_lg.res)

dam_hypo_bp.res <- read.csv("./results/dam-hypo-all-results-BP.csv")
dam_hypo_bp.arrange <- dam_hypo_bp.res %>% arrange(desc(log2FoldChange))
dam_hypo_bp.v <- dam_hypo_bp.arrange$log2FoldChange
names(dam_hypo_bp.v) <- dam_hypo_bp.arrange$X
head(dam_hypo_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_hypo_bp.v) #two tailed
rm(dam_hypo_bp.arrange, dam_hypo_bp.res)

dam_hypo_lg.res <- read.csv("./results/dam-hypo-all-results-LG.csv")
dam_hypo_lg.arrange <- dam_hypo_lg.res %>% arrange(desc(log2FoldChange))
dam_hypo_lg.v <- dam_hypo_lg.arrange$log2FoldChange
names(dam_hypo_lg.v) <- dam_hypo_lg.arrange$X
head(dam_hypo_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_hypo_lg.v) #two tailed
rm(dam_hypo_lg.arrange, dam_hypo_lg.res)

dam_nac_bp.res <- read.csv("./results/dam-nac-all-results-BP.csv")
dam_nac_bp.arrange <- dam_nac_bp.res %>% arrange(desc(log2FoldChange))
dam_nac_bp.v <- dam_nac_bp.arrange$log2FoldChange
names(dam_nac_bp.v) <- dam_nac_bp.arrange$X
head(dam_nac_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_nac_bp.v) #two tailed
rm(dam_nac_bp.arrange, dam_nac_bp.res)

dam_nac_lg.res <- read.csv("./results/dam-nac-all-results-LG.csv")
dam_nac_lg.arrange <- dam_nac_lg.res %>% arrange(desc(log2FoldChange))
dam_nac_lg.v <- dam_nac_lg.arrange$log2FoldChange
names(dam_nac_lg.v) <- dam_nac_lg.arrange$X
head(dam_nac_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_nac_lg.v) #two tailed
rm(dam_nac_lg.arrange, dam_nac_lg.res)

dam_pl_bp.res <- read.csv("./results/dam-pl-all-results-BP.csv")
dam_pl_bp.arrange <- dam_pl_bp.res %>% arrange(desc(log2FoldChange))
dam_pl_bp.v <- dam_pl_bp.arrange$log2FoldChange
names(dam_pl_bp.v) <- dam_pl_bp.arrange$X
head(dam_pl_bp.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_pl_bp.v) #two tailed
rm(dam_pl_bp.arrange, dam_pl_bp.res)

dam_pl_lg.res <- read.csv("./results/dam-pl-all-results-LG.csv")
dam_pl_lg.arrange <- dam_pl_lg.res %>% arrange(desc(log2FoldChange))
dam_pl_lg.v <- dam_pl_lg.arrange$log2FoldChange
names(dam_pl_lg.v) <- dam_pl_lg.arrange$X
head(dam_pl_lg.v) #check that they are ordered and that highest is positive and lowest is negative
tail(dam_pl_lg.v) #two tailed
rm(dam_pl_lg.arrange, dam_pl_lg.res)

gc()
#run fgsea analysis-------
scoreType <- "std" #all are two tailed tests

#function to run fgsea and collapse results. returns df of mainpathways
get.C5gsea <- function(df, analysis_name, region, age, comparison){
  set.seed(42)
  fgseaRes <- fgseaMultilevel(pathways = C5.ensembl.ls,
                              stats = df,
                              scoreType = scoreType, 
                              maxSize = 250, 
                              minSize = 15, 
                              sampleSize = 101,
                              eps = 0, 
                              nproc = 1)
  set.seed(42)
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][pval < 0.01], 
                                        C5.ensembl.ls, 
                                        df)
  
  mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
    order(-padj), pathway]
  
  keep_rows <- fgseaRes[, which((pathway %in% mainPathways) ==TRUE)]
  Main_pathways <- fgseaRes[keep_rows,] 
  Main_pathways <- Main_pathways %>% 
    arrange(padj) %>% 
    mutate(pathway = gsub("GOBP_", "", pathway)) %>% 
    mutate(pathway = gsub("_", " ", pathway)) %>% 
    mutate(pval = round(pval, digits = 4)) %>% 
    mutate(NES = round(pval, digits = 4)) %>% 
    mutate(padj = round(padj, digits =  4)) %>% 
    mutate( Age = age) %>% 
    mutate( Region = region) %>% 
    mutate(Comparison = comparison) %>% 
    mutate(FDR = padj) %>% 
    dplyr::select(Age, Region, Comparison, pathway, pval, FDR, NES, leadingEdge)
  
  library(org.Rn.eg.db)
  
  Main_pathways <- Main_pathways[, leadingEdge := mapIdsList(
    x=org.Rn.eg.db, 
    keys=leadingEdge,
    keytype="ENSEMBL", 
    column="SYMBOL")]
  
  library(data.table)
  fwrite(Main_pathways, file= paste("./results/", analysis_name, "_gsea.csv", sep = "") , sep="\t", sep2=c("", " ", ""))
  
  return(Main_pathways)
}


#run analysis for
#run analysis for pups
p_amy_bp <-  get.C5gsea(pup_amy_bp.v,   "pup_amy_bp",  region = "Amygdala",  age = "Pup", comparison = "Bisphenol")
p_hipp_bp <-  get.C5gsea(pup_hipp_bp.v, "pup_hipp_bp", region = "Hippocampus", age = "Pup", comparison = "Bisphenol")
p_hypo_bp <-  get.C5gsea(pup_hypo_bp.v, "pup_hypo_bp", region = "Hypothalamus", age = "Pup", comparison = "Bisphenol")
p_nac_bp <-  get.C5gsea(pup_nac_bp.v,   "pup_nac_bp",  region = "Nucleus accumbens",  age = "Pup", comparison = "Bisphenol")
p_pl_bp <-  get.C5gsea(pup_pl_bp.v,     "pup_pl_bp",   region = "Prelimbic cortex",   age = "Pup", comparison = "Bisphenol")

p_amy_lg <-  get.C5gsea(pup_amy_lg.v,   "pup_amy_lg",  region = "Amygdala",  age = "Pup", comparison = "Licking/grooming")
p_hipp_lg <-  get.C5gsea(pup_hipp_lg.v, "pup_hipp_lg", region = "Hippocampus", age = "Pup", comparison = "Licking/grooming")
p_hypo_lg <-  get.C5gsea(pup_hypo_lg.v, "pup_hypo_lg", region = "Hypothalamus", age = "Pup", comparison = "Licking/grooming")
p_nac_lg <-  get.C5gsea(pup_nac_lg.v,   "pup_nac_lg",  region = "Nucleus accumbens",  age = "Pup", comparison = "Licking/grooming")
p_pl_lg <-  get.C5gsea(pup_pl_lg.v,     "pup_pl_lg",   region = "Prelimbic cortex",   age = "Pup", comparison = "Licking/grooming")

#run analysis for dams
d_amy_bp <-  get.C5gsea(dam_amy_bp.v,   "dam_amy_bp",  region = "Amygdala",  age = "dam", comparison = "Bisphenol")
d_hipp_bp <-  get.C5gsea(dam_hipp_bp.v, "dam_hipp_bp", region = "Hippocampus", age = "dam", comparison = "Bisphenol")
d_hypo_bp <-  get.C5gsea(dam_hypo_bp.v, "dam_hypo_bp", region = "Hypothalamus", age = "dam", comparison = "Bisphenol")
d_nac_bp <-  get.C5gsea(dam_nac_bp.v,   "dam_nac_bp",  region = "Nucleus accumbens",  age = "dam", comparison = "Bisphenol")
d_pl_bp <-  get.C5gsea(dam_pl_bp.v,     "dam_pl_bp",   region = "Prelimbic cortex",   age = "dam", comparison = "Bisphenol")

d_amy_lg <-  get.C5gsea(dam_amy_lg.v,   "dam_amy_lg",  region = "Amygdala",  age = "dam", comparison = "Licking/grooming")
d_hipp_lg <-  get.C5gsea(dam_hipp_lg.v, "dam_hipp_lg", region = "Hippocampus", age = "dam", comparison = "Licking/grooming")
d_hypo_lg <-  get.C5gsea(dam_hypo_lg.v, "dam_hypo_lg", region = "Hypothalamus", age = "dam", comparison = "Licking/grooming")
d_nac_lg <-  get.C5gsea(dam_nac_lg.v,   "dam_nac_lg",  region = "Nucleus accumbens",  age = "dam", comparison = "Licking/grooming")
d_pl_lg <-  get.C5gsea(dam_pl_lg.v,     "dam_pl_lg",   region = "Prelimbic cortex",   age = "dam", comparison = "Licking/grooming")

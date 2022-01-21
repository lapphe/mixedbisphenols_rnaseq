# bisphenol manuscript supplemental tables
library(tidyverse)
setwd("./Bisphenols pilot study 2019-2020/R/RNA_seq/clean_data/")

#Suppl table 1: all dam deseq2 results

amy_BP <- read.csv("./deseq2-results-dams/dam-amy-all-results-BP.csv") %>% 
  mutate(Region = "Amygdala") %>% 
  mutate(comparison = "Bisphenol")

hipp_BP <- read.csv("./deseq2-results-dams/dam-hipp-all-results-BP.csv") %>% 
  mutate(Region = "Hippocampus") %>% 
  mutate(comparison = "Bisphenol")

hypo_BP <- read.csv("./deseq2-results-dams/dam-hypo-all-results-BP.csv") %>% 
  mutate(Region = "Hypothalamus") %>% 
  mutate(comparison = "Bisphenol")

nac_BP <- read.csv("./deseq2-results-dams/dam-nac-all-results-BP.csv") %>% 
  mutate(Region = "Nucleus Accumbens") %>% 
  mutate(comparison = "Bisphenol")

pl_BP <- read.csv("./deseq2-results-dams/dam-pl-all-results-BP.csv") %>% 
  mutate(Region = "Prelimbic Cortex") %>% 
  mutate(comparison = "Bisphenol")

amy_LG <- read.csv("./deseq2-results-dams/dam-amy-all-results-LG.csv") %>% 
  mutate(Region = "Amygdala") %>% 
  mutate(comparison = "Licking/grooming")

hipp_LG <- read.csv("./deseq2-results-dams/dam-hipp-all-results-LG.csv") %>% 
  mutate(Region = "Hippocampus") %>% 
  mutate(comparison = "Licking/grooming")

hypo_LG <- read.csv("./deseq2-results-dams/dam-hypo-all-results-LG.csv") %>% 
  mutate(Region = "Hypothalamus") %>% 
  mutate(comparison = "Licking/grooming")

nac_LG <- read.csv("./deseq2-results-dams/dam-nac-all-results-LG.csv") %>% 
  mutate(Region = "Nucleus Accumbens") %>% 
  mutate(comparison = "Licking/grooming")

pl_LG <- read.csv("./deseq2-results-dams/dam-pl-all-results-Lg.csv") %>% 
  mutate(Region = "Prelimbic Cortex") %>% 
  mutate(comparison = "Licking/grooming")

dam_dfs <- list(amy_BP, hipp_BP, hypo_BP, nac_BP, pl_BP,
                amy_LG, hipp_LG, hypo_LG, nac_LG, pl_LG)

df <- do.call("rbind", dam_dfs)

write.csv(df, "supp_table_1_dam_deseq2_res.csv")

#----table 2 pup deseq2 results------

amy_BP <- read.csv("./deseq2-results/pup-amy-all-results-BP.csv") %>% 
  mutate(Region = "Amygdala") %>% 
  mutate(comparison = "Bisphenol")

hipp_BP <- read.csv("./deseq2-results/pup-hipp-all-results-BP.csv") %>% 
  mutate(Region = "Hippocampus") %>% 
  mutate(comparison = "Bisphenol")

hypo_BP <- read.csv("./deseq2-results/pup-hypo-all-results-BP.csv") %>% 
  mutate(Region = "Hypothalamus") %>% 
  mutate(comparison = "Bisphenol")

nac_BP <- read.csv("./deseq2-results/pup-nac-all-results-BP.csv") %>% 
  mutate(Region = "Nucleus Accumbens") %>% 
  mutate(comparison = "Bisphenol")

pl_BP <- read.csv("./deseq2-results/pup-pl-all-results-BP.csv") %>% 
  mutate(Region = "Prelimbic Cortex") %>% 
  mutate(comparison = "Bisphenol")

amy_LG <- read.csv("./deseq2-results/pup-amy-all-results-LG.csv") %>% 
  mutate(Region = "Amygdala") %>% 
  mutate(comparison = "Licking/grooming")

hipp_LG <- read.csv("./deseq2-results/pup-hipp-all-results-LG.csv") %>% 
  mutate(Region = "Hippocampus") %>% 
  mutate(comparison = "Licking/grooming")

hypo_LG <- read.csv("./deseq2-results/pup-hypo-all-results-LG.csv") %>% 
  mutate(Region = "Hypothalamus") %>% 
  mutate(comparison = "Licking/grooming")

nac_LG <- read.csv("./deseq2-results/pup-nac-all-results-LG.csv") %>% 
  mutate(Region = "Nucleus Accumbens") %>% 
  mutate(comparison = "Licking/grooming")

pl_LG <- read.csv("./deseq2-results/pup-pl-all-results-Lg.csv") %>% 
  mutate(Region = "Prelimbic Cortex") %>% 
  mutate(comparison = "Licking/grooming")

pup_dfs <- list(amy_BP, hipp_BP, hypo_BP, nac_BP, pl_BP,
                amy_LG, hipp_LG, hypo_LG, nac_LG, pl_LG)

df <- do.call("rbind", pup_dfs)

write.csv(df, "supp_table_2_pup_deseq2_res.csv")

# suppl table 3 dam GSEA results------
setwd("~/Bisphenols pilot study 2019-2020/R/RNA_seq/gsea results")

amy_BP <- read.csv("dam_amy_bp.csv") %>% 
  mutate(Region = "Amygdala") %>% 
  mutate(comparison = "Bisphenol")

hipp_BP <- read.csv("dam_hipp_bp.csv") %>% 
  mutate(Region = "Hippocampus") %>% 
  mutate(comparison = "Bisphenol")

hypo_BP <- read.csv("dam_hypo_bp.csv") %>% 
  mutate(Region = "Hypothalamus") %>% 
  mutate(comparison = "Bisphenol")

nac_BP <- read.csv("dam_nac_bp.csv") %>% 
  mutate(Region = "Nucleus Accumbens") %>% 
  mutate(comparison = "Bisphenol")

pl_BP <- read.csv("dam_pl_bp.csv") %>% 
  mutate(Region = "Prelimbic Cortex") %>% 
  mutate(comparison = "Bisphenol")

amy_LG <- read.csv("dam_amy_lg.csv") %>% 
  mutate(Region = "Amygdala") %>% 
  mutate(comparison = "Licking/grooming")

hipp_LG <- read.csv("dam_hipp_lg.csv") %>% 
  mutate(Region = "Hippocampus") %>% 
  mutate(comparison = "Licking/grooming")

hypo_LG <- read.csv("dam_hypo_lg.csv") %>% 
  mutate(Region = "Hypothalamus") %>% 
  mutate(comparison = "Licking/grooming")

nac_LG <- read.csv("dam_nac_lg.csv") %>% 
  mutate(Region = "Nucleus Accumbens") %>% 
  mutate(comparison = "Licking/grooming")

pl_LG <- read.csv("pup_pl_lg.csv") %>% 
  mutate(Region = "Prelimbic Cortex") %>% 
  mutate(comparison = "Licking/grooming")

dam_dfs <- list(amy_BP, hipp_BP, hypo_BP, nac_BP, pl_BP,
                amy_LG, hipp_LG, hypo_LG, nac_LG, pl_LG)

df <- do.call("rbind", dam_dfs)

write.csv(df, "supp_table_3_dam_gsea_res.csv")


# suppl table 4 pup GSEA results------
setwd("~/Bisphenols pilot study 2019-2020/R/RNA_seq/gsea results")

amy_BP <- read.csv("pup_amy_bp.csv") %>% 
  mutate(Region = "Amygdala") %>% 
  mutate(comparison = "Bisphenol")

hipp_BP <- read.csv("pup_hipp_bp.csv") %>% 
  mutate(Region = "Hippocampus") %>% 
  mutate(comparison = "Bisphenol")

hypo_BP <- read.csv("pup_hypo_bp.csv") %>% 
  mutate(Region = "Hypothalamus") %>% 
  mutate(comparison = "Bisphenol")

nac_BP <- read.csv("pup_nac_bp.csv") %>% 
  mutate(Region = "Nucleus Accumbens") %>% 
  mutate(comparison = "Bisphenol")

pl_BP <- read.csv("pup_pl_bp.csv") %>% 
  mutate(Region = "Prelimbic Cortex") %>% 
  mutate(comparison = "Bisphenol")

amy_LG <- read.csv("pup_amy_lg.csv") %>% 
  mutate(Region = "Amygdala") %>% 
  mutate(comparison = "Licking/grooming")

hipp_LG <- read.csv("pup_hipp_lg.csv") %>% 
  mutate(Region = "Hippocampus") %>% 
  mutate(comparison = "Licking/grooming")

hypo_LG <- read.csv("pup_hypo_lg.csv") %>% 
  mutate(Region = "Hypothalamus") %>% 
  mutate(comparison = "Licking/grooming")

nac_LG <- read.csv("pup_nac_lg.csv") %>% 
  mutate(Region = "Nucleus Accumbens") %>% 
  mutate(comparison = "Licking/grooming")

pl_LG <- read.csv("pup_pl_lg.csv") %>% 
  mutate(Region = "Prelimbic Cortex") %>% 
  mutate(comparison = "Licking/grooming")

pup_dfs <- list(amy_BP, hipp_BP, hypo_BP, nac_BP, pl_BP,
                amy_LG, hipp_LG, hypo_LG, nac_LG, pl_LG)

df <- do.call("rbind", pup_dfs)

write.csv(df, "supp_table_4_pup_gsea_res.csv")

#BP project RNAseq functions

#------------------------------------------------------------------------------
#function to get results from dds objects, all three comparisons
get_LFC_results <- function(dds, LFC_threshold = 0.2, pvalue_threshold = 0.05) 
{
  temp1 <- results(dds, contrast=c("condition","bpa","vehicle")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "bpa - vehicle") 
  
  temp2 <- results(dds, contrast=c("condition","mixed","vehicle")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "mixed - vehicle") 
  
  temp3 <- results(dds, contrast=c("condition","bpa", "mixed")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "mixed - bpa") 
  
  rbind(temp1, temp2, temp3) %>% 
    arrange(desc(abs(log2FoldChange))) -> result
  
  return(result)
}


#-------------------------------------------------------------------------------
#for use with DeSeq2 analysis to generate normalized gene counts
number_the_genes <- function(df){
  
  result <- df %>% 
    left_join(DEG_results_final %>% 
                dplyr::rename(real_LFC = log2FoldChange)) %>% 
    dplyr::rename(LFC_with_NAs = log2FoldChange) %>% 
    dplyr::rename(log2FoldChange = real_LFC) %>% 
    dplyr::select(ensgene, log2FoldChange,contrast) %>%  
    spread(contrast, log2FoldChange, fill = 0) %>% 
    arrange(`bpa - vehicle`,`mixed - vehicle`) %>%
    mutate(x_var = row_number()) %>% 
    gather(contrast,log2FoldChange,2:4) %>% 
    mutate(contrast = factor(contrast, 
                             levels = rev(c('bpa - vehicle','mixed - vehicle')))) 
  
  return(result)
  
}
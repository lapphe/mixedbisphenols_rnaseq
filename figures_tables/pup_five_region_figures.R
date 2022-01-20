#Gestational bisphenol project
#Hannah Lapp PhD
#Champagne Lab, Fall 2021

#BP pilot pups RNAseq DEGs for 5 regions figure for unadjusted p val <.05

library("readr")
library(apeglm)
library(glue)
library(tidyverse)

amy1 <- read.csv("./results/pup-amy-all-results-BP.csv") %>% 
  filter(pvalue < .05)
hipp1 <- read.csv("./results/pup-hipp-all-results-BP.csv") %>% 
  filter(pvalue < .05)
pl1 <- read.csv("./results/pup-pl-all-results-BP.csv") %>% 
  filter(pvalue < .05)
nac1 <- read.csv("./results/pup-nac-all-results-BP.csv") %>% 
  filter(pvalue < .05)
hypo1 <- read.csv("./results/pup-hypo-all-results-BP.csv") %>% 
  filter(pvalue < .05)

amy1$region <- "Amygdala"
hipp1$region <- "Hippocampus"
pl1$region <- "Prelimbic Cortex"
nac1$region <- "Nucleus Accumbens"
hypo1$region <- "Hypothalamus"

DEG_results_final <- rbind(amy1, hipp1, pl1, nac1, hypo1)

names(DEG_results_final)[1] <- "ensgene"

#split copied df into groups by region. create df with only amygdala rows
DEG_results_final %>% 
  split(.$region) %>% 
  .$Amy -> df

df <- amy1 %>% 
  filter(is.na(log2FoldChange))

number_the_genes <- function(df){
  
  result <- df %>% 
    dplyr::left_join(DEG_results_final %>% 
                       dplyr::rename(real_LFC = log2FoldChange)) %>% #add new column with lo2fc called real_flc
    dplyr::rename(LFC_with_NAs = log2FoldChange) %>% 
    dplyr::rename(log2FoldChange = real_LFC) %>% #change name back( why???)
    dplyr::select(ensgene, log2FoldChange,region) %>%  
    spread(region, log2FoldChange, fill = 0) %>% 
    arrange(Amygdala,
            `Prelimbic Cortex`,
            `Nucleus Accumbens`,
            Hypothalamus,
            Hippocampus
    ) %>%
    mutate(x_var = row_number()) %>% 
    gather(region,log2FoldChange,2:6) %>% 
    mutate(region = factor(region, 
                           levels = rev(c('Amygdala', 
                                          'Prelimbic Cortex',
                                          'Nucleus Accumbens',
                                          'Hypothalamus',
                                          'Hippocampus')))) 
  
  return(result)
  
}
DEG_result_for_plot <- number_the_genes(DEG_results_final)


#find number of genes per region
DEG_results_final_temp<-DEG_results_final%>%
  filter(region == "Amygdala")

#make table for number of DEGs by brain region
table(DEG_results_final$region)%>%
  as.data.frame()%>%
  dplyr::rename(region = Var1) %>% 
  mutate(facet_it = glue('{region} ({Freq})')) %>% 
  dplyr::select(region, facet_it) -> freq1



DEG_result_for_plot %>% 
  filter(is.na(log2FoldChange))

DEG_result_for_plot_replace <- DEG_result_for_plot%>%
  mutate(log2FoldChange = (ifelse(.$log2FoldChange > 2, 2, .$log2FoldChange)))%>%
  mutate(log2FoldChange = (ifelse(.$log2FoldChange < -2, -2, .$log2FoldChange)))





DEG_result_for_plot_replace$region <- factor(DEG_result_for_plot_replace$region, 
                                             levels = c("Amygdala", "Prelimbic Cortex",
                                                        "Nucleus Accumbens",
                                                        "Hypothalamus",
                                                        "Hippocampus"))

df <- DEG_result_for_plot_replace %>% 
  left_join(freq1)

desired_order <- c("Hippocampus",
                   "Hypothalamus",
                   "Nucleus Accumbens",
                   "Prelimbic Cortex",
                   "Amygdala")
df$region <- factor(as.character(df$region), levels = desired_order)
freq1$region <- factor(freq1$region, levels = desired_order)

a <- grep("^Amy", freq1$facet_it, value = TRUE) 
b <- grep("^Pre", freq1$facet_it, value = TRUE)
c <- grep("^Nuc", freq1$facet_it, value = TRUE)
d <- grep("^Hypo", freq1$facet_it, value = TRUE)
e <- grep("^Hipp", freq1$facet_it, value = TRUE)
y_axis <- c(e,d,c,b,a)

df%>%
  ggplot(aes(x =x_var, y = region, fill = log2FoldChange))+
  geom_tile(color = NA, size = 20)+
  #facet_grid(facet_it ~., switch = "y")+
  labs(x = "",
       y = "",       ,
       fill = "Log2 Fold Change")+
  scale_fill_gradient2(low="seagreen", 
                       mid="white", 
                       high="darkgoldenrod1", #colors in the scale
                       midpoint=0,    #same midpoint for plots 
                       limits=c(-2, 2),#same limits for plots
                       breaks=seq(-2,2,1),   #breaks in the scale bar
                       na.value = "white")+
  
  scale_y_discrete(labels = y_axis)+
  ggtitle("Differentially expressed genes for bisphenol exposure in pups")+
  theme_bw(base_size = 9)+
  theme(legend.position = c(-0.18,-0.04),
        
        legend.direction = "horizontal",
        axis.text.x= element_blank(),
        #axis.text.y= element_blank(),
        text = element_text(size = 10),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        plot.title = element_text(size = 10, face = "bold"),
        plot.background =  element_rect(color = "black"),
        panel.background = element_rect(fill = "transparent",colour = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 12)) -> ticktick
ticktick

ggsave(ticktick, filename = "pup_bp_unionpheatmap.png", path = "./figures_tables/", dpi= 600, height = 2, width = 6 )


#----------------LG figure 5 region figure--
#filter all licking/grooming deseq2 results by pvalue <.05
amy2 <- read.csv("./results/pup-amy-all-results-LG.csv") %>% 
  filter(pvalue < .05)
hipp2 <- read.csv("./results/pup-hipp-all-results-LG.csv") %>% 
  filter(pvalue < .05)
pl2<- read.csv("./results/pup-pl-all-results-LG.csv") %>% 
  filter(pvalue < .05)
nac2<- read.csv("./results/pup-nac-all-results-LG.csv") %>% 
  filter(pvalue < .05)
hypo2 <- read.csv("./results/pup-hypo-all-results-LG.csv") %>% 
  filter(pvalue < .05)

amy2$region <- "Amygdala"
hipp2$region <- "Hippocampus"
pl2$region <- "Prelimbic Cortex"
nac2$region <- "Nucleus Accumbens"
hypo2$region <- "Hypothalamus"

DEG_results_final <- rbind(amy2, hipp2, pl2, nac2, hypo2)

names(DEG_results_final)[1] <- "ensgene"

#split copied df into groups by region. create df with only amygdala rows
DEG_results_final %>% 
  split(.$region) %>% 
  .$Amy -> df

df <- amy2 %>% 
  filter(is.na(log2FoldChange))

DEG_result_for_plot <- number_the_genes(DEG_results_final)

#find number of genes per region
DEG_results_final_temp<-DEG_results_final%>%
  filter(region == "Amygdala")

#make table for number of DEGs by brain region
table(DEG_results_final$region)%>%
  as.data.frame()%>%
  dplyr::rename(region = Var1) %>% 
  mutate(facet_it = glue('{region} ({Freq})')) %>% 
  dplyr::select(region, facet_it) -> freq1

DEG_result_for_plot %>% 
  filter(is.na(log2FoldChange))

DEG_result_for_plot_replace <- DEG_result_for_plot%>%
  mutate(log2FoldChange = (ifelse(.$log2FoldChange > 2, 2, .$log2FoldChange)))%>%
  mutate(log2FoldChange = (ifelse(.$log2FoldChange < -2, -2, .$log2FoldChange)))

DEG_result_for_plot_replace$region <- factor(DEG_result_for_plot_replace$region, 
                                             levels = c("Amygdala", "Prelimbic Cortex",
                                                        "Nucleus Accumbens",
                                                        "Hypothalamus",
                                                        "Hippocampus"))

df <- DEG_result_for_plot_replace %>% 
  left_join(freq1)

desired_order <- c("Hippocampus",
                   "Hypothalamus",
                   "Nucleus Accumbens",
                   "Prelimbic Cortex",
                   "Amygdala")
df$region <- factor(as.character(df$region), levels = desired_order)
freq1$region <- factor(freq1$region, levels = desired_order)

a <- grep("^Amy", freq1$facet_it, value = TRUE) 
b <- grep("^Pre", freq1$facet_it, value = TRUE)
c <- grep("^Nuc", freq1$facet_it, value = TRUE)
d <- grep("^Hypo", freq1$facet_it, value = TRUE)
e <- grep("^Hipp", freq1$facet_it, value = TRUE)
y_axis <- c(e,d,c,b,a)

df%>%
  ggplot(aes(x =x_var, y = region, fill = log2FoldChange))+
  geom_tile(color = NA, size = 20)+
  #facet_grid(facet_it ~., switch = "y")+
  labs(x = "",
       y = "",       ,
       fill = "Log2 Fold Change")+
  scale_fill_gradient2(low="seagreen", 
                       mid="white", 
                       high="darkgoldenrod1", #colors in the scale
                       midpoint=0,    #same midpoint for plots 
                       limits=c(-2, 2),#same limits for plots
                       breaks=seq(-2,2,1),   #breaks in the scale bar
                       na.value = "white")+
  
  scale_y_discrete(labels = y_axis)+
  ggtitle("Differentially expressed genes for licking/grooming in pups")+
  theme_bw(base_size = 9)+
  theme(legend.position = c(-0.18,-0.04),
        
        legend.direction = "horizontal",
        axis.text.x= element_blank(),
        #axis.text.y= element_blank(),
        text = element_text(size = 10),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        plot.title = element_text(size = 10, face = "bold"),
        plot.background =  element_rect(color = "black"),
        panel.background = element_rect(fill = "transparent",colour = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 12)) -> ticktick2
ticktick2

ggsave(ticktick2, filename = "pup_lg_unionpheatmap.png", path = "./figures_tables/", dpi= 600, height = 2, width = 6 )

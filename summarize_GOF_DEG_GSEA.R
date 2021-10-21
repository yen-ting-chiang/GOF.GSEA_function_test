
summarise.GOF.DEG.GSEA.result <- function()
{
  
}

getwd()
setwd("C:/Users/dannyj/Documents/Rproject/GOF.GSEA_function_test/C5_positive")
files <- list.files(path=
                      "C:/Users/dannyj/Documents/Rproject/GOF.GSEA_function_test/C5_positive")
results_bind <- sapply(files, function(x) mget(read.csv(x)), simplify = TRUE)

#summarize DEG data-------------------------------------------------------

DEG_GBM_TP53_IMPACT_MODERATE_vs_WT_positive <- 
  read.csv("DEG_GBM_TP53_IMPACT_MODERATE_vs_WT_positive.csv")
DEG_PAAD_TP53_IMPACT_MODERATE_vs_WT_positive <- 
  read.csv("DEG_PAAD_TP53_IMPACT_MODERATE_vs_WT_positive.csv")
DEG_ESCA_TP53_IMPACT_MODERATE_vs_WT_positive <- 
  read.csv("DEG_ESCA_TP53_IMPACT_MODERATE_vs_WT_positive.csv")

GBM_M = DEG_GBM_TP53_IMPACT_MODERATE_vs_WT_positive %>% 
  mutate(cancer_type = "GBM") %>% 
  distinct(ENSEMBL, .keep_all = TRUE)
PAAD_M = DEG_PAAD_TP53_IMPACT_MODERATE_vs_WT_positive %>% 
  mutate(cancer_type = "PAAD")%>% 
  distinct(ENSEMBL, .keep_all = TRUE)
ESCA_M = DEG_ESCA_TP53_IMPACT_MODERATE_vs_WT_positive %>% 
  mutate(cancer_type = "ESCA")%>% 
  distinct(ENSEMBL, .keep_all = TRUE)

combine1 = bind_rows(GBM_M, PAAD_M, ESCA_M)

DEG_number_ENTREZID = combine1 %>% 
  group_by(ENTREZID) %>%
  summarise(n()) %>% 
  filter(is.na(ENTREZID) == FALSE) %>% 
  arrange(desc(n()))

DEG_number_SYMBOL = combine1 %>% 
  group_by(SYMBOL) %>%
  summarise(n()) %>% 
  filter(is.na(SYMBOL) == FALSE) %>% 
  arrange(desc(n()))


#summarize GSEA data-------------------------------------------------------

GBM_TP53_IMPACT_MODERATE_vs_WT_C5_positive <- 
  read.csv("GBM_TP53_IMPACT_MODERATE_vs_WT_C5_positive.csv")
PAAD_TP53_IMPACT_MODERATE_vs_WT_C5_positive <- 
  read.csv("PAAD_TP53_IMPACT_MODERATE_vs_WT_C5_positive.csv")
ESCA_TP53_IMPACT_MODERATE_vs_WT_C5_positive <- 
  read.csv("ESCA_TP53_IMPACT_MODERATE_vs_WT_C5_positive.csv")

GBM_M = GBM_TP53_IMPACT_MODERATE_vs_WT_C5_positive %>% 
  mutate(cancer_type = "GBM")
PAAD_M = PAAD_TP53_IMPACT_MODERATE_vs_WT_C5_positive %>% 
  mutate(cancer_type = "PAAD")
ESCA_M = ESCA_TP53_IMPACT_MODERATE_vs_WT_C5_positive %>% 
  mutate(cancer_type = "ESCA")

combine1 = bind_rows(GBM_M, PAAD_M, ESCA_M)

GSEA_number = combine1 %>% 
  group_by(ID) %>%
  summarise(n())%>% 
  arrange(desc(n()))

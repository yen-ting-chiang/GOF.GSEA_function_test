# #scripts
# filenames  <- list.files(pattern = "*_TP53_IMPACT_MODERATE_vs_WT_C5_positive.csv" )
# combo_data <- purrr::map_df(filenames, 
#                             ~read.csv(.x, stringsAsFactors = FALSE) %>% 
#                               mutate(filename = .x))
# save(combo_data, 
#      file = "combo_TP53_IMPACT_MODERATE_vs_WT_C5_positive.rdata")
# wirte.csv(combo_data, 
#      file = "combo_TP53_IMPACT_MODERATE_vs_WT_C5_positive.csv")
# # combo_data2 <- do.call(rbind, lapply(filenames, function(x) 
# #   cbind(read.csv(x, stringsAsFactors = FALSE), filename = x)))
# 
# GSEA_number = combo_data %>% 
#   group_by(ID) %>%
#   summarise(n())%>% 
#   arrange(desc(n()))
# 
# save(GSEA_number, 
#      file = "summary_TP53_IMPACT_MODERATE_vs_WT_C5_positive.rdata")
# wirte.csv(GSEA_number, 
#           file = "summary_TP53_IMPACT_MODERATE_vs_WT_C5_positive.csv")


library(dplyr)
setwd("C:/Users/dannyj/Documents/Rproject/GOF.GSEA_function_test/GOF_result")
getwd()

summarise.GOF.result <- function(hugo_gene_name,
                                 grouping_type, 
                                 conpairing_type, 
                                 reference_type,
                                 msigdb_type,
                                 pos_or_neg)
{
  filenames  <- list.files(pattern = sprintf("*_%s_%s_%s_vs_%s_%s_%s.csv",
                                             hugo_gene_name,
                                             grouping_type, 
                                             conpairing_type, 
                                             reference_type,
                                             msigdb_type,
                                             pos_or_neg))
  combo_data <- purrr::map_df(filenames, 
                              ~read.csv(.x, stringsAsFactors = FALSE) %>% 
                                mutate(filename = .x))
  save(combo_data, 
       file = sprintf("combo_%s_%s_%s_vs_%s_%s_%s.rdata",
                      hugo_gene_name,
                      grouping_type, 
                      conpairing_type, 
                      reference_type,
                      msigdb_type,
                      pos_or_neg))
  write.csv(combo_data, 
            file = sprintf("combo_%s_%s_%s_vs_%s_%s_%s.csv",
                           hugo_gene_name,
                           grouping_type, 
                           conpairing_type, 
                           reference_type,
                           msigdb_type,
                           pos_or_neg))
  
  GSEA_number = combo_data %>% 
    group_by(ID) %>%
    summarise(n())%>% 
    arrange(desc(n()))
  
  save(GSEA_number, 
       file = sprintf("summary_%s_%s_%s_vs_%s_%s_%s.rdata",
                      hugo_gene_name,
                      grouping_type, 
                      conpairing_type, 
                      reference_type,
                      msigdb_type,
                      pos_or_neg))
  write.csv(GSEA_number, 
            file = sprintf("summary_%s_%s_%s_vs_%s_%s_%s.csv",
                           hugo_gene_name,
                           grouping_type, 
                           conpairing_type, 
                           reference_type,
                           msigdb_type,
                           pos_or_neg))
}

summarise.GOF.result(hugo_gene_name = "TP53",
                     grouping_type = "IMPACT", 
                     conpairing_type = "MODERATE", 
                     reference_type = "WT",
                     msigdb_type = "C3",
                     pos_or_neg = "positive")

msigdb_list <- c("H", "C1", "C2", "C3", "C4", "C5", "C6","C7","C8")

lapply(msigdb_list,
       summarise.GOF.result,
       hugo_gene_name = "TP53",
       grouping_type = "IMPACT", 
       conpairing_type = "MODERATE", 
       reference_type = "WT",
       pos_or_neg = "positive")

lapply(msigdb_list,
       summarise.GOF.result,
       hugo_gene_name = "TP53",
       grouping_type = "IMPACT", 
       conpairing_type = "MODERATE", 
       reference_type = "WT",
       pos_or_neg = "negative")








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

library(dplyr)
setwd("C:/Users/dannyj/Documents/Rproject/GOF.GSEA_function_test/GOF_result2")
getwd()

#summarize GSEA data-------------------------------------------------------
summarise.GSEA.result <- function(hugo_gene_name,
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
                              ~read.csv(.x, stringsAsFactors = FALSE,
                                        colClasses = "character") %>%
                                mutate(filename = .x))

  
  save(combo_data, 
       file = sprintf("%s_%s_%s_vs_%s_%s_%s_combo.rdata",
                      hugo_gene_name,
                      grouping_type, 
                      conpairing_type, 
                      reference_type,
                      msigdb_type,
                      pos_or_neg))
  write.csv(combo_data, 
            file = sprintf("%s_%s_%s_vs_%s_%s_%s_combo.csv",
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
       file = sprintf("%s_%s_%s_vs_%s_%s_%s_summary.rdata",
                      hugo_gene_name,
                      grouping_type, 
                      conpairing_type, 
                      reference_type,
                      msigdb_type,
                      pos_or_neg))
  write.csv(GSEA_number, 
            file = sprintf("%s_%s_%s_vs_%s_%s_%s_summary.csv",
                           hugo_gene_name,
                           grouping_type, 
                           conpairing_type, 
                           reference_type,
                           msigdb_type,
                           pos_or_neg))
}


msigdb_list <- c("H", "C1", "C2", "C3", "C4", "C5", "C6","C7","C8")

lapply(msigdb_list,
       summarise.GSEA.result,
       hugo_gene_name = "TP53",
       grouping_type = "IMPACT", 
       conpairing_type = "MODERATE", 
       reference_type = "WT",
       pos_or_neg = "positive")

lapply(msigdb_list,
       summarise.GSEA.result,
       hugo_gene_name = "TP53",
       grouping_type = "IMPACT", 
       conpairing_type = "MODERATE", 
       reference_type = "WT",
       pos_or_neg = "negative")

#summarize DEG data-------------------------------------------------------

summarise.DEG.result <- function(hugo_gene_name,
                                 grouping_type, 
                                 conpairing_type, 
                                 reference_type,
                                 pos_or_neg)
{
  filenames  <- list.files(pattern = sprintf("*_%s_%s_%s_vs_%s_%s.csv",
                                             hugo_gene_name,
                                             grouping_type, 
                                             conpairing_type, 
                                             reference_type,
                                             pos_or_neg))
  combo_data <- purrr::map_df(filenames,
                              ~read.csv(.x, stringsAsFactors = FALSE,
                                        colClasses = "character") %>%
                                mutate(filename = .x) %>% 
                                distinct(ENSEMBL, .keep_all = TRUE))
  
  save(combo_data, 
       file = sprintf("%s_%s_%s_vs_%s_%s_DEG_combo.rdata",
                      hugo_gene_name,
                      grouping_type, 
                      conpairing_type, 
                      reference_type,
                      pos_or_neg))
  
  write.csv(combo_data, 
            file = sprintf("%s_%s_%s_vs_%s_%s_DEG_combo.csv",
                           hugo_gene_name,
                           grouping_type, 
                           conpairing_type, 
                           reference_type,
                           pos_or_neg))
  
  DEG_ENTREZID_number = combo_data %>% 
    group_by(ENTREZID) %>%
    summarise(n()) %>% 
    filter(is.na(ENTREZID) == FALSE) %>% 
    arrange(desc(n()))
  
  save(DEG_ENTREZID_number, 
       file = sprintf("DEG_summary_%s_%s_%s_vs_%s_ENTREZID_%s.rdata",
                      hugo_gene_name,
                      grouping_type, 
                      conpairing_type, 
                      reference_type,
                      pos_or_neg))
  
  write.csv(DEG_ENTREZID_number, 
            file = sprintf("DEG_summary_%s_%s_%s_vs_%s_ENTREZID_%s.csv",
                           hugo_gene_name,
                           grouping_type, 
                           conpairing_type, 
                           reference_type,
                           pos_or_neg))
  
  
  DEG_SYMBOL_number = combo_data %>% 
    group_by(SYMBOL) %>%
    summarise(n()) %>% 
    filter(is.na(SYMBOL) == FALSE) %>% 
    arrange(desc(n()))
  
  save(DEG_SYMBOL_number, 
       file = sprintf("DEG_summary_%s_%s_%s_vs_%s_SYMBOL_%s.rdata",
                      hugo_gene_name,
                      grouping_type, 
                      conpairing_type, 
                      reference_type,
                      pos_or_neg))
  write.csv(DEG_SYMBOL_number, 
            file = sprintf("DEG_summary_%s_%s_%s_vs_%s_SYMBOL_%s.csv",
                           hugo_gene_name,
                           grouping_type, 
                           conpairing_type, 
                           reference_type,
                           pos_or_neg))
  
}


pos_neg_list <- c("positive", "negative")

lapply(pos_neg_list,
       summarise.DEG.result,
       hugo_gene_name = "TP53",
       grouping_type = "IMPACT", 
       conpairing_type = "MODERATE", 
       reference_type = "WT")

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

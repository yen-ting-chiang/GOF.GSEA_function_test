library(tidyr)
library(dplyr)

setwd("C:/Users/dannyj/Documents/Rproject/GOF.GSEA_function_test")

summarise.TP53.mutation.condition <- function(project_name, 
                                               hugo_gene_name)
{
  load(sprintf("%s.%s.ColData.rdata",
               project_name,
               hugo_gene_name))
  ColData_summary = ColData %>% 
    group_by(IMPACT) %>% 
    summarise(n())
  write.csv(ColData_summary,
            file = sprintf("%s.%s.mutation.condition.summary.csv",
                           project_name,
                           hugo_gene_name))
}

TCGA_project_list <- c("ACC", "BLCA", "LGG", "BRCA",
                       "CESC", "CHOL", "COAD", "LAML",
                       "ESCA", "GBM", "HNSC", "KICH",
                       "KIRC", "KIRP", "LIHC", "LUAD", "LUSC",
                       "DLBC", "MESO", "OV", "PAAD",
                       "PCPG", "PRAD", "READ", "SARC", "SKCM",
                       "STAD", "TGCT", "THYM", "THCA", "UCS",
                       "UCEC", "UVM")

a <- lapply(TCGA_project_list,
       summarise.TP53.mutation.condition,
       hugo_gene_name = "TP53")

save(a, file = "summarise.TP53.mutation.condition.rdata")
write.table(a, file = "summarise.TP53.mutation.condition.csv")

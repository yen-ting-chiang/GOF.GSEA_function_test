library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)


# TCGA_project_list <- c("ACC", "BLCA", "LGG", "BRCA", 
#                        "CESC", "CHOL", "COAD", "LAML",
#                        "ESCA", "GBM", "HNSC", "KICH", 
#                        "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", 
#                        "DLBC", "MESO", "OV", "PAAD", 
#                        "PCPG", "PRAD", "READ", "SARC", "SKCM", 
#                        "STAD", "TGCT", "THYM", "THCA", "UCS", 
#                        "UCEC", "UVM")

TCGA_project_list <- c("DLBC", "GBM")

GOF.TCGA.MakeColRow <- function(project_name, 
                                hugo_gene_name)
{
  
  #TCGAbiolinks(TCGA RNAseq data download)------------------------------
  
  query <- GDCquery(project = sprintf("TCGA-%s", project_name),
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",
                    experimental.strategy = "RNA-Seq",
                    legacy = FALSE)
  
  GDCdownload(query,
              method = "api",
              files.per.chunk = 10)
  
  RNAseq_data <- GDCprepare(query,
                            save = TRUE,
                            save.filename = sprintf("%s_RNAseq_counts.rdata",
                                                    project_name))
  
  #deseq2 RowData prepare----------------------------------------------------
  
  RNAseq_data_matrix <- as.data.frame(assay(RNAseq_data))
  
  RNAseq_data_matrix_sample <- RNAseq_data_matrix
  colnames(RNAseq_data_matrix_sample)=
    gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
         "",
         colnames(RNAseq_data_matrix_sample))
  
  #TCGAbiolinks(TCGA mutation data download)-------------------------
  
  muse.maf <- GDCquery_Maf(sprintf("%s", project_name), 
                           pipelines = "muse",
                           save.csv=TRUE)
  
  
  muse.maf_target_gene <- muse.maf %>% 
    filter(Hugo_Symbol==sprintf("%s", hugo_gene_name)) %>% 
    mutate(Tumor_Case_Barcode=gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
                                   "",
                                   Tumor_Sample_Barcode))%>%
    dplyr::select(Tumor_Case_Barcode,
                  Protein_position,
                  HGVSp_Short,
                  Exon_Number,
                  Variant_Classification,
                  Consequence,
                  IMPACT,
                  SIFT,
                  PolyPhen
    )
  
  #create list of all case IDs with mutation data-------------------------------
  muse.maf_case_list=muse.maf%>%
    mutate(Tumor_Case_Barcode=gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
                                   "",
                                   Tumor_Sample_Barcode))%>%
    dplyr::select(Tumor_Case_Barcode)
  
  muse.maf_case_list_unique <- muse.maf_case_list
  muse.maf_case_list_unique=unique(muse.maf_case_list_unique)
  # write.csv(PAAD.muse.maf_case_list_unique, 
  # file = "PAAD.muse.maf_case_list_unique.csv")
  
  #create the RowData that matches the ColData--------------------------------------------------------
  
  #delete duplicated RNAseq data
  RNAseq_data_matrix_sample_delete_duplicated <- 
    RNAseq_data_matrix_sample[, !duplicated(colnames(RNAseq_data_matrix_sample))]
  
  mutation_list <- muse.maf_case_list_unique[,1, drop = TRUE]
  
  #select samples that have mutation data
  RowData <- 
    RNAseq_data_matrix_sample_delete_duplicated[, colnames(RNAseq_data_matrix_sample_delete_duplicated)%in%
                                                  (mutation_list)]
  
  save(RowData,
       file = sprintf("%s.%s.RowData.rdata", project_name, hugo_gene_name))
  save(muse.maf_target_gene,
       file = sprintf("%s.muse.maf_%s.rdata", project_name, hugo_gene_name))
  save(muse.maf_case_list_unique,
       file = sprintf("%s.muse.maf_case_list_unique.rdata", project_name))
  # write.csv(RowData, file = "RowData.csv")
  
  # load(sprintf("%s.%s.RowData.rdata", project_name, hugo_gene_name))
  # load(sprintf("%s.muse.maf_%s.rdata", project_name, hugo_gene_name))
  # load(sprintf("%s.muse.maf_case_list_unique.rdata", project_name))
  
  #create the WT sample Coldata -----------
  
  list_of_samples_containing_MT_target_gene <- muse.maf_target_gene[,1, drop = TRUE]
  list_of_samples_containing_MT_target_gene_unique <- 
    unique(list_of_samples_containing_MT_target_gene)
  
  
  DF1 <- data.frame(mutation_list)
  DF1.1 <- DF1 %>% filter(mutation_list %in% 
                            list_of_samples_containing_MT_target_gene_unique == FALSE) %>% 
    mutate(Protein_position = "WT", HGVSp_Short = "WT", IMPACT = "WT")
  colnames(DF1.1)[1] <- "Tumor_Case_Barcode"
  
  #create the multiple-mutation sample Coldata -----------
  
  DF2 <- data.frame(mutation_list)
  DF2.1 <- muse.maf_target_gene %>% 
    filter(duplicated(Tumor_Case_Barcode) == TRUE)%>% 
    dplyr::select(Tumor_Case_Barcode, Protein_position, HGVSp_Short, IMPACT)%>% 
    mutate(Protein_position = "multiple", 
           HGVSp_Short = "multiple", 
           IMPACT = "multiple")
  
  #create the single-mutation sample Coldata -----------
  DF3 <-muse.maf_target_gene%>% 
    dplyr::select(Tumor_Case_Barcode, Protein_position, HGVSp_Short, IMPACT) %>% 
    filter(Tumor_Case_Barcode %in% DF2.1[,1, drop = TRUE] == FALSE)
  
  
  #combine ColDatas----------------------------------------------------
  DF4 <- full_join(DF1.1, DF2.1, 
                   by = c("Tumor_Case_Barcode", 
                          "Protein_position", 
                          "HGVSp_Short", 
                          "IMPACT"))
  
  DF5 <- full_join(DF4, DF3, 
                   by = c("Tumor_Case_Barcode", 
                          "Protein_position", 
                          "HGVSp_Short", 
                          "IMPACT"))
  
  ColData <- DF5 %>% 
    filter(Tumor_Case_Barcode %in% colnames(RowData) == TRUE)
  
  save(ColData, file = sprintf("%s.%s.ColData.rdata", 
                               project_name, 
                               hugo_gene_name))
  write.csv(ColData, file = sprintf("%s.%s.ColData.csv", 
                                    project_name,
                                    hugo_gene_name))
  
  
  
}


lapply(TCGA_project_list,
       GOF.TCGA.MakeColRow,
       hugo_gene_name = "TP53")



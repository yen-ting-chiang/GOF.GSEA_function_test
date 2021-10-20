
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
browseVignettes("TCGAbiolinks")

BiocManager::install("SummarizedExperiment")
install.packages("DT")
install.packages("dplyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")


BiocManager::install("AnnotationHub")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("msigdbr")

library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(DESeq2)

library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(msigdbr)

# Check memory limit
memory.limit()
# Change memory limit
memory.limit(size = 28000)

GOF.TCGA <- function(project_name, 
                     hugo_gene_name, 
                     grouping_type, 
                     conpairing_type, 
                     reference_type, 
                     design_formula)
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
  write.csv(RowData, file = "RowData.csv")
  
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
  DF2.1 <- unique(DF2.1)
  
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
  
  
  #launch DESeq2 analysis---------------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(RowData),
                                colData = ColData,
                                design = design_formula)
  
  dds$IMPACT <- relevel(dds$IMPACT, ref = reference_type)
  dds2 <- DESeq(dds)
  
  res <- results(dds2, 
                 name = sprintf("%s_%s_vs_%s", 
                                grouping_type, 
                                conpairing_type, 
                                reference_type))
  
  save(res, file = sprintf("deseq2_%s_%s_%s_%s_vs_%s.rdata", 
                           project_name,
                           hugo_gene_name,
                           grouping_type, 
                           conpairing_type, 
                           reference_type))
  
  write.csv(as.data.frame(res),
            file = sprintf("deseq2_%s_%s_%s_%s_vs_%s.csv", 
                           project_name,
                           hugo_gene_name,
                           grouping_type, 
                           conpairing_type, 
                           reference_type))
  
  
  #Launch GSEA analysis----------------------------------
  
  D <- as.data.frame(res)
  
  D2=D%>%
    dplyr::filter(grepl('ENSG', rownames(D)))%>%
    dplyr::filter(is.na(log2FoldChange)==FALSE)
  
  D2_ENTREZID <- bitr(rownames(D2), fromType = "ENSEMBL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db,
                      drop = TRUE)
  
  D2 <- tibble::rownames_to_column(D2, "ENSEMBL")
  
  D3=right_join(D2,D2_ENTREZID, by=c("ENSEMBL"="ENSEMBL"))
  
  ## feature 1: numeric vector
  geneList_ENTREZ <- D3[,3]
  
  ## feature 2: named vector
  names(geneList_ENTREZ) <- D3[,8]
  
  ## feature 3: decreasing order
  geneList_ENTREZ <- sort(geneList_ENTREZ, decreasing = TRUE)
  
  collections <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
  for(i in collections)
  {
    m_t2g <- msigdbr(species = "Homo sapiens", 
                     category = i) %>% 
      dplyr::select(gs_name, entrez_gene)
    
    run.GSEA=GSEA(geneList_ENTREZ, 
                  exponent = 1,
                  minGSSize = 25,
                  maxGSSize = 500,
                  eps = 1e-10, 
                  pvalueCutoff = 0.25,
                  pAdjustMethod = "BH",
                  TERM2GENE= m_t2g,
                  TERM2NAME = NA,
                  verbose = TRUE,
                  seed = FALSE,
                  by = "fgsea")
    
    save(run.GSEA,
         file = sprintf("%s_%s_%s_%s_vs_%s_%s.rdata",
                        project_name,
                        hugo_gene_name,
                        grouping_type, 
                        conpairing_type, 
                        reference_type,
                        i))
    
    pos_NES <- as.data.frame(run.GSEA@result) %>% 
      filter(NES >= 0) %>% 
      arrange(desc(NES))
    
    write.csv(pos_NES, 
              file = sprintf("%s_%s_%s_%s_vs_%s_%s_positive.csv",
                             project_name,
                             hugo_gene_name,
                             grouping_type, 
                             conpairing_type, 
                             reference_type,
                             i))
    
    neg_NES <- as.data.frame(run.GSEA@result) %>% 
      filter(NES < 0) %>% 
      arrange(NES)
    
    write.csv(neg_NES, 
              file = sprintf("%s_%s_%s_%s_vs_%s_%s_negative.csv",
                             project_name,
                             hugo_gene_name,
                             grouping_type, 
                             conpairing_type, 
                             reference_type,
                             i))
    
  }
  
  # launch DEG analysis-----------------------------------------------------
  
  deg <- as.data.frame(res)
  
  #positive
  deg2=deg%>%
    dplyr::filter(grepl('ENSG', rownames(deg)))%>%
    dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
    dplyr::filter(is.na(padj)==FALSE)%>%
    dplyr::filter(log2FoldChange>=0.58)%>%
    dplyr::filter(padj<0.05)
  deg2_ENTREZID <- bitr(rownames(deg2), fromType = "ENSEMBL",
                        toType = c("ENTREZID"),
                        OrgDb = org.Hs.eg.db,
                        drop = TRUE)
  deg2_SYMBOL <- bitr(rownames(deg2), fromType = "ENSEMBL",
                      toType = c("SYMBOL"),
                      OrgDb = org.Hs.eg.db,
                      drop = TRUE)
  deg2 <- tibble::rownames_to_column(deg2, "ENSEMBL")
  deg3=left_join(deg2,deg2_ENTREZID, by=c("ENSEMBL"="ENSEMBL"))
  deg4=left_join(deg3,deg2_SYMBOL, by=c("ENSEMBL"="ENSEMBL"))
  deg4<- deg4[ , c(1,8, 9, 2:7)]
  deg5=deg4%>%
    arrange(desc(log2FoldChange))
  
  save(deg5, file = sprintf("DEG_%s_%s_%s_%s_vs_%s_positive.rdata",
                            project_name,
                            hugo_gene_name,
                            grouping_type, 
                            conpairing_type, 
                            reference_type))
  write.csv(deg5,
            file=sprintf("DEG_%s_%s_%s_%s_vs_%s_positive.csv",
                         project_name,
                         hugo_gene_name,
                         grouping_type, 
                         conpairing_type, 
                         reference_type))
  
  #negative
  deg2n=deg%>%
    dplyr::filter(grepl('ENSG', rownames(deg)))%>%
    dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
    dplyr::filter(is.na(padj)==FALSE)%>%
    dplyr::filter(log2FoldChange <= -0.58)%>%
    dplyr::filter(padj<0.05)
  deg2n_ENTREZID <- bitr(rownames(deg2n), fromType = "ENSEMBL",
                         toType = c("ENTREZID"),
                         OrgDb = org.Hs.eg.db,
                         drop = TRUE)
  deg2n_SYMBOL <- bitr(rownames(deg2n), fromType = "ENSEMBL",
                       toType = c("SYMBOL"),
                       OrgDb = org.Hs.eg.db,
                       drop = TRUE)
  deg2n <- tibble::rownames_to_column(deg2n, "ENSEMBL")
  deg3n=left_join(deg2n,deg2n_ENTREZID, by=c("ENSEMBL"="ENSEMBL"))
  deg4n=left_join(deg3n,deg2n_SYMBOL, by=c("ENSEMBL"="ENSEMBL"))
  deg4n<- deg4n[ , c(1,8, 9, 2:7)]
  deg5n=deg4n%>%
    arrange(log2FoldChange)
  
  save(deg5n, file = sprintf("DEG_%s_%s_%s_%s_vs_%s_negative.rdata",
                             project_name,
                             hugo_gene_name,
                             grouping_type, 
                             conpairing_type, 
                             reference_type))
  write.csv(deg5n,
            file=sprintf("DEG_%s_%s_%s_%s_vs_%s_negative.csv",
                         project_name,
                         hugo_gene_name,
                         grouping_type, 
                         conpairing_type, 
                         reference_type))
  
}

GOF.TCGA(project_name = "GBM", 
         hugo_gene_name = "TP53",
         grouping_type = "IMPACT", 
         conpairing_type = "MODERATE", 
         reference_type = "WT",
         design_formula = ~ IMPACT)



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

DEG_number_SYMBOL = combine1 %>% 
  group_by(SYMBOL) %>%
  summarise(n()) %>% 
  filter(is.na(SYMBOL) == FALSE) %>% 
  arrange(desc(n()))






load("ColData.rdata")
load(sprintf("%s.%s.RowData.rdata", project_name, hugo_gene_name))
load("%s.muse.maf_%s.rdata")
load("%s.muse.maf_case_list_unique.rdata")



#deseq2 RowData prepare----------------------------------------------------

# load("PAAD_RNAseq_counts.rda")

RNAseq_data_matrix= as.data.frame(assay(data))
save(RNAseq_data_matrix,
     file = "PAAD_RNAseq_data_matrix.rdata")
write.csv(RNAseq_data_matrix, 
          file="PAAD_RNAseq_data_matrix.csv")

# load("PAAD_RNAseq_data_matrix.rdata")

RNAseq_data_matrix_sample <- RNAseq_data_matrix
colnames(RNAseq_data_matrix_sample)=
  gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
       "",
       colnames(RNAseq_data_matrix_sample))


write.csv(RNAseq_data_matrix_sample, 
          file="RNAseq_data_matrix_sample.csv",
          quote=FALSE)



#TCGAbiolinks(TCGA mutation data download)------

PAAD.muse.maf <- GDCquery_Maf("PAAD", 
                              pipelines = "muse",
                              save.csv=TRUE)
save(PAAD.muse.maf,
     file = "PAAD.muse.maf.rdata")


load("PAAD.muse.maf.rdata")
#PAAD.muse.maf=read.csv("TCGA.PAAD.muse.59a84472-27d4-497c-8f37-8bc447ff9374.DR-10.0.somatic.maf.csv",  header=FALSE, stringsAsFactors = FALSE)


PAAD.muse.maf_TP53_clean=PAAD.muse.maf %>% 
  filter(Hugo_Symbol=="TP53") %>% 
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
         PolyPhen,
  )

write.csv(PAAD.muse.maf_TP53_clean, 
          file="PAAD.muse.maf_TP53_clean.csv")

#create list of all case IDs with mutation data-------------------------------
PAAD.muse.maf_case_list=PAAD.muse.maf%>%
  mutate(Tumor_Case_Barcode=gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
                                 "",
                                 Tumor_Sample_Barcode))%>%
  dplyr::select(Tumor_Case_Barcode)

PAAD.muse.maf_case_list_unique <- PAAD.muse.maf_case_list
PAAD.muse.maf_case_list_unique=unique(PAAD.muse.maf_case_list_unique)
write.csv(PAAD.muse.maf_case_list_unique, file = "PAAD.muse.maf_case_list_unique.csv")



#create the RowData that matches the ColData--------------------------------------------------------

#delete duplicated RNAseq data
RNAseq_data_matrix_sample_delete_duplicated <- 
  RNAseq_data_matrix_sample[, !duplicated(colnames(RNAseq_data_matrix_sample))]

mutation_list <- PAAD.muse.maf_case_list_unique[,1, drop = TRUE]

#select samples that have mutation data
RowData <- 
  RNAseq_data_matrix_sample_delete_duplicated[, colnames(RNAseq_data_matrix_sample_delete_duplicated)%in%
                                                (mutation_list)]

# save(RowData,
#      file = "RowData.rdata")
# save(PAAD.muse.maf_TP53_clean,
#      file = "PAAD.muse.maf_TP53_clean.rdata")
# save(PAAD.muse.maf_case_list_unique,
#      file = "PAAD.muse.maf_case_list_unique.rdata")
# write.csv(RowData, file = "RowData.csv")

load("RowData.rdata")
load("PAAD.muse.maf_TP53_clean.rdata")
laod("PAAD.muse.maf_case_list_unique.rdata")


#create the TP53-WT sample Coldata -----------

samples_with_TP53_MT_list <- PAAD.muse.maf_TP53_clean[,1, drop = TRUE]
samples_with_TP53_MT_list_unique <- unique(samples_with_TP53_MT_list)


DF1 <- data.frame(mutation_list)
DF1.1 <- DF1 %>% filter(mutation_list %in% 
                          samples_with_TP53_MT_list_unique == FALSE) %>% 
  mutate(Protein_position = "WT", HGVSp_Short = "WT", IMPACT = "WT")
colnames(DF1.1)[1] <- "Tumor_Case_Barcode"

#create the TP53-double-mutation sample Coldata -----------

DF2 <- data.frame(mutation_list)
DF2.1 <- PAAD.muse.maf_TP53_clean %>% 
  filter(duplicated(Tumor_Case_Barcode) == TRUE)%>% 
  dplyr::select(Tumor_Case_Barcode, Protein_position, HGVSp_Short, IMPACT)%>% 
  mutate(Protein_position = "DB", HGVSp_Short = "DB", IMPACT = "DB")

#create the TP53-single-mutation sample Coldata -----------
DF3 <-  PAAD.muse.maf_TP53_clean%>% 
  select(Tumor_Case_Barcode, Protein_position, HGVSp_Short, IMPACT) %>% 
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

# save(ColData, file = "ColData.rdata")
# write.csv(ColData, file = "ColData.csv")
load("ColData.rdata")


#launch DESeq2 analysis----------------------------------------------

load("RowData.rdata")
load("ColData.rdata")

library("DESeq2")
#cts <- RowData_match




GOF.deseq2 <- function(grouping_type, 
                       conpairing_type, 
                       reference_type,
                       design_formula,
                       dds_relevel)
{
  
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(RowData),
                                colData = ColData,
                                design = design_formula)
  
  dds_relevel <- relevel(dds_relevel, ref = reference_type)
  dds2 <- DESeq(dds)
  
  
  
  res <- results(dds2, 
                 name = sprintf("%s_%s_vs_%s", 
                                grouping_type, 
                                conpairing_type, 
                                reference_type))
  
  save(res, file = sprintf("%s_%s_vs_%s.rdata", 
                           grouping_type, 
                           conpairing_type, 
                           reference_type))
  
  # write.csv(as.data.frame(res),                                                                 
  #           file="IMPACT_HIGH_vs_WT.csv")
  
}


GOF.deseq2(grouping_type = "IMPACT", 
           conpairing_type = "MODERATE", 
           reference_type = "HIGH",
           design_formula = ~ IMPACT,
           dds_relevel = dds$IMPACT)


dds <- DESeqDataSetFromMatrix(countData = as.matrix(RowData),
                              colData = ColData,
                              design= ~ IMPACT)

dds$IMPACT <- relevel(dds$IMPACT, ref = "WT")
dds2 <- DESeq(dds)


resultsNames(dds2) # lists the coefficients
res <- results(dds2, name = "IMPACT_HIGH_vs_WT")

save(res, file = "IMPACT_HIGH_vs_WT.rdata")

write.csv(as.data.frame(res1),                                                                 
          file="IMPACT_HIGH_vs_WT.csv")






res2 <- results(dds2, name="IMPACT_MODERATE_vs_WT")

save(res2, file = "IMPACT_MODERATE_vs_WT.rdata")
write.csv(as.data.frame(res2), 
          file="IMPACT_MODERATE_vs_WT.csv")

res3 <- results(dds2, name="IMPACT_MODERATE_vs_WT")

write.csv(as.data.frame(res3),
          file="IMPACT_MODERATE_vs_WT.csv")

#launch GSEA analysis----------------------------------------------------


library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(msigdbr)

load("IMPACT_MODERATE_vs_HIGH.rdata")


GOF.GSEA <- function(GSEA_collections)
{
  D <- as.data.frame(res)
  
  D2=D%>%
    dplyr::filter(grepl('ENSG', rownames(D)))%>%
    dplyr::filter(is.na(log2FoldChange)==FALSE)
  
  D2_ENTREZID <- bitr(rownames(D2), fromType = "ENSEMBL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db,
                      drop = TRUE)
  
  D2 <- tibble::rownames_to_column(D2, "ENSEMBL")
  
  D3=right_join(D2,D2_ENTREZID, by=c("ENSEMBL"="ENSEMBL"))
  
  ## feature 1: numeric vector
  geneList_ENTREZ <- D3[,3]
  
  ## feature 2: named vector
  names(geneList_ENTREZ) <- D3[,8]
  
  ## feature 3: decreasing order
  geneList_ENTREZ <- sort(geneList_ENTREZ, decreasing = TRUE)
  
  m_t2g <- msigdbr(species = "Homo sapiens", 
                   category = GSEA_collections) %>% 
    dplyr::select(gs_name, entrez_gene)
  
  run.GSEA=GSEA(geneList_ENTREZ, 
                   exponent = 1,
                   minGSSize = 25,
                   maxGSSize = 500,
                   eps = 1e-10, 
                   pvalueCutoff = 0.25,
                   pAdjustMethod = "BH",
                   TERM2GENE= m_t2g,
                   TERM2NAME = NA,
                   verbose = TRUE,
                   seed = FALSE,
                   by = "fgsea")
  
  save(run.GSEA,
       file = sprintf("GBM_IMPACT_MODERATE_vs_HIGH_%s.rdata",
                      GSEA_collections))
  
  pos_NES <- as.data.frame(run.GSEA@result) %>% 
    filter(NES >= 0) %>% 
    arrange(desc(NES))
  
  write.csv(pos_NES, 
            file = sprintf("GBM_IMPACT_MODERATE_vs_HIGH_%s_positive.csv", 
                           GSEA_collections))
  
  neg_NES <- as.data.frame(run.GSEA@result) %>% 
    filter(NES < 0) %>% 
    arrange(NES)
  
  write.csv(neg_NES, 
            file = sprintf("GBM_IMPACT_MODERATE_vs_HIGH_%s_negative.csv", 
                           GSEA_collections))
}

GOF.GSEA(GSEA_collections = "C6")


D <- as.data.frame(res)

D2=D%>%
  dplyr::filter(grepl('ENSG', rownames(D)))%>%
  dplyr::filter(is.na(log2FoldChange)==FALSE)

D2_ENTREZID <- bitr(rownames(D2), fromType = "ENSEMBL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db,
                   drop = TRUE)

D2 <- tibble::rownames_to_column(D2, "ENSEMBL")

D3=right_join(D2,D2_ENTREZID, by=c("ENSEMBL"="ENSEMBL"))


## feature 1: numeric vector
geneList_ENTREZ <- D3[,3]
## feature 2: named vector
names(geneList_ENTREZ) <- D3[,8]
## feature 3: decreasing order
geneList_ENTREZ <- sort(geneList_ENTREZ, decreasing = TRUE)
head(geneList_ENTREZ)
str(geneList_ENTREZ)

keytypes(org.Hs.eg.db)

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

H=GSEA(
  geneList_ENTREZ,
  exponent = 1,
  minGSSize = 25,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.25,
  pAdjustMethod = "BH",
  TERM2GENE= m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)

write.csv(as.data.frame(H@result), 
            file="GBM_IMPACT_MODERATE_vs_HIGH_H.csv")


# launch DEG analysis-----------------------------------------------------
library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

deg <- read.csv("IMPACT_HIGH_vs_WT.csv")
deg <- as.data.frame(res)

names(deg)[1] <- "ENSEMBL_ID"

#positive
deg2=deg%>%
  dplyr::filter(grepl('ENSG', rownames(deg)))%>%
  dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
  dplyr::filter(is.na(padj)==FALSE)%>%
  dplyr::filter(log2FoldChange>=0.58)%>%
  dplyr::filter(padj<0.05)
deg2_ENTREZID <- bitr(rownames(deg2), fromType = "ENSEMBL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db,
                   drop = TRUE)
deg2_SYMBOL <- bitr(rownames(deg2), fromType = "ENSEMBL",
                 toType = c("SYMBOL"),
                 OrgDb = org.Hs.eg.db,
                 drop = TRUE)
deg2 <- tibble::rownames_to_column(deg2, "ENSEMBL")
deg3=left_join(deg2,deg2_ENTREZID, by=c("ENSEMBL"="ENSEMBL"))
deg4=left_join(deg3,deg2_SYMBOL, by=c("ENSEMBL"="ENSEMBL"))
deg4<- deg4[ , c(1,8, 9, 2:7)]
deg5=deg4%>%
  arrange(desc(log2FoldChange))

save(deg5, file = )
write.csv(deg5,
          file="PAAD_IMPACT_HIGH_vs_WT_DEG_pos_unfilter.csv")

#negative
deg2n=deg%>%
  dplyr::filter(grepl('ENSG', rownames(deg)))%>%
  dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
  dplyr::filter(is.na(padj)==FALSE)%>%
  dplyr::filter(log2FoldChange <= -0.58)%>%
  dplyr::filter(padj<0.05)
deg2n_ENTREZID <- bitr(rownames(deg2n), fromType = "ENSEMBL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db,
                      drop = TRUE)
deg2n_SYMBOL <- bitr(rownames(deg2n), fromType = "ENSEMBL",
                    toType = c("SYMBOL"),
                    OrgDb = org.Hs.eg.db,
                    drop = TRUE)
deg2n <- tibble::rownames_to_column(deg2n, "ENSEMBL")
deg3n=left_join(deg2n,deg2n_ENTREZID, by=c("ENSEMBL"="ENSEMBL"))
deg4n=left_join(deg3n,deg2n_SYMBOL, by=c("ENSEMBL"="ENSEMBL"))
deg4n<- deg4n[ , c(1,8, 9, 2:7)]
deg5n=deg4n%>%
  arrange(log2FoldChange)
write.csv(deg5n,
          file="PAAD_IMPACT_HIGH_vs_WT_DEG_pos_unfilter.csv")

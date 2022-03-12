library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(DESeq2)

library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(DO.db)
library(msigdbr)

LUAD.varscan2.maf <- GDCquery_Maf("LUAD", 
                                  pipelines = "varscan2",
                                  save.csv=TRUE)
save(LUAD.mutect2.maf,
     file = "LUAD.mutect2.maf.rdata")




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
  
  varscan2.maf <- GDCquery_Maf(sprintf("%s", project_name),
                               pipelines = "varscan2",
                               save.csv=TRUE)
  
  
  varscan2.maf_target_gene <- varscan2.maf %>%
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
  varscan2.maf_case_list=varscan2.maf%>%
    mutate(Tumor_Case_Barcode=gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$",
                                   "",
                                   Tumor_Sample_Barcode))%>%
    dplyr::select(Tumor_Case_Barcode)
  
  varscan2.maf_case_list_unique <- varscan2.maf_case_list
  varscan2.maf_case_list_unique=unique(varscan2.maf_case_list_unique)
  # write.csv(PAAD.varscan2.maf_case_list_unique,
  # file = "PAAD.varscan2.maf_case_list_unique.csv")
  
  #create the RowData that matches the ColData--------------------------------------------------------
  
  #delete duplicated RNAseq data
  RNAseq_data_matrix_sample_delete_duplicated <-
    RNAseq_data_matrix_sample[, !duplicated(colnames(RNAseq_data_matrix_sample))]
  
  mutation_list <- varscan2.maf_case_list_unique[,1, drop = TRUE]
  
  #select samples that have mutation data
  RowData <-
    RNAseq_data_matrix_sample_delete_duplicated[, colnames(RNAseq_data_matrix_sample_delete_duplicated)%in%
                                                  (mutation_list)]
  
  save(RowData,
       file = sprintf("%s.%s.RowData.rdata", project_name, hugo_gene_name))
  save(varscan2.maf_target_gene,
       file = sprintf("%s.varscan2.maf_%s.rdata", project_name, hugo_gene_name))
  save(varscan2.maf_case_list_unique,
       file = sprintf("%s.varscan2.maf_case_list_unique.rdata", project_name))
  write.csv(RowData, file = "RowData.csv")
  
  # load(sprintf("%s.%s.RowData.rdata", project_name, hugo_gene_name))
  # load(sprintf("%s.varscan2.maf_%s.rdata", project_name, hugo_gene_name))
  # load(sprintf("%s.varscan2.maf_case_list_unique.rdata", project_name))
  
  #create the WT sample Coldata -----------
  
  list_of_samples_containing_MT_target_gene <- varscan2.maf_target_gene[,1, drop = TRUE]
  list_of_samples_containing_MT_target_gene_unique <-
    unique(list_of_samples_containing_MT_target_gene)
  
  
  DF1 <- data.frame(mutation_list)
  DF1.1 <- DF1 %>% filter(mutation_list %in%
                            list_of_samples_containing_MT_target_gene_unique == FALSE) %>%
    mutate(Protein_position = "WT", HGVSp_Short = "WT", IMPACT = "WT")
  colnames(DF1.1)[1] <- "Tumor_Case_Barcode"
  
  #create the multiple-mutation sample Coldata -----------
  
  DF2 <- data.frame(mutation_list)
  DF2.1 <- varscan2.maf_target_gene %>%
    filter(duplicated(Tumor_Case_Barcode) == TRUE)%>%
    dplyr::select(Tumor_Case_Barcode, Protein_position, HGVSp_Short, IMPACT)%>%
    mutate(Protein_position = "multiple",
           HGVSp_Short = "multiple",
           IMPACT = "multiple")
  DF2.1 <- unique(DF2.1)
  
  #create the single-mutation sample Coldata -----------
  DF3 <-varscan2.maf_target_gene%>%
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
  ColData[is.na(ColData)] <- "NA"
  
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
  
  dds$HGVSp_Short <- relevel(dds$HGVSp_Short, ref = reference_type)
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

GOF.TCGA(project_name = "LUAD", 
         hugo_gene_name = "EGFR",
         grouping_type = "HGVSp_Short", 
         conpairing_type = "p.E746_A750del", 
         reference_type = "WT",
         design_formula = ~ HGVSp_Short)

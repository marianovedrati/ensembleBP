
setwd("/Users/Maria/Desktop/giorgio/ensembleBP/Codes/")

#' @description DEG:
# This function conducts Differential Expression Analysis (DEA) using DESeq2 for RNA-Seq data.
# It accepts count data and phenotype data as input, performs DEA based on specified parameters,
# and saves the resulting list of the name of differentially expressed genes (DEGs) as an RDS object.
# If batch correction is needed, a batch column can be specified in the phenotype data.
# The function provides flexibility to select DEGs based on different criteria such as all DEGs,
# upregulated DEGs, or downregulated DEGs.

# Parameters:
#' @param pathdf Path to the count data file. genes on rows, samples on columns
#' @param pathclin relative path to dfPheno: 2 columns should be named ID and class
#' @param which_genes relative to which genes to keep  to select ("all", "up", or "down") (default: "all").
#' @param batchcol Specifies the name of the batch column in the phenotype data for batch correction (default: "").
#' @returns list of dge that will be save in rds 

# Libraries Required:
# DESeq2: For conducting Differential Expression Analysis.

DEG <- function(pathdf = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Count.csv",
                pathclin = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Pheno.csv",
                which_genes ="all",batchcol="") {
  
  library(DESeq2) # Load DESeq2 library
  
  # Import the count dataframe
  df <- read.csv2(pathdf, row.names = 1)
  
  # Import the descriptive table
  df_pheno <- read.csv2(pathclin, row.names = 1)
  
  
  #use the transponse of df to have genes in rows and sample in the column
  df <- as.data.frame(t(df)) 
  
  
  # Match the order of df and df_pheno
  m<-match(colnames(df), rownames(df_pheno))
  df_pheno<-df_pheno[m,]
  
  
  # Convert class column to factor
  df_pheno$class <- as.factor(df_pheno$class)
  
  # If batch column is specified, include it in the df_pheno
  if (batchcol != "") {
    
    df_pheno$batch <- as.factor(df_pheno[, batchcol])
    
    # Create the DESeq dataset with batch adjustment
    dds <- DESeqDataSetFromMatrix(countData = df,
                                  colData = df_pheno,
                                  design = ~ batch + class)  
    
  } else {
    # Create the DESeq dataset without batch adjustment
    dds <- DESeqDataSetFromMatrix(countData = df,
                                  colData = df_pheno,
                                  design = ~ class)
  }
  
  # Filter out low count rows
  keep <- rowSums(counts(dds)>10) > 10
  dds <- dds[keep,]
  
  # Run DESeq analysis
  dds <- DESeq(dds)
  
  # lists the coefficient
  resultsNames(dds)
  
  # Extract DESeq results making the comparison between class 1 and 0
  res <- results(dds,alpha = 0.05, contrast =c("class","1","0"))
  
  
  # Based on which_genes argument, select genes accordingly
  if (which_genes == "all") {
    list_dge <- list(rownames(res[res$padj < 0.05 , ]))
  } else if (which_genes == "up") {
    list_dge <- list(rownames(res[res$padj < 0.05 & res$log2FoldChange >0, ]))
  } else if (which_genes == "down") {
    list_dge <- list(rownames(res[res$padj < 0.05 & res$log2FoldChange < 0, ]))
  }
  
  # Save list of gene names as an RDS object
  saveRDS(list_dge, paste0("../Results/list_dge.rds"))

  
}


#run example

DEG(pathdf = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Count.csv",
    pathclin = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Pheno.csv",
    which_genes ="up",batchcol = "patient.gender")

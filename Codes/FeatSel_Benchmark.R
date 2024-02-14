#setwd("/Users/giorgiomontesi/Desktop/Universita_di_Siena/A_PhD_Project/Biomarker_Prediction/Review_Master/Codes")

## Questo codice esegue un Benchmark tra embedded methods di Biomarker discovery
## l'obiettivo è quello di comparare diversi metodi di classificazione ML-based
## in termini di Accuracy e altre metriche su diversi Dataset pubblici
## in modo da compararne le metriche al variare delle dimensioni dei df in input




# Import first df
df <- read.csv2("../Data/ACC_Adrenocortical_Carcinoma/ACC_Count.csv", row.names = 1)
df_pheno <- read.csv2("../Data/ACC_Adrenocortical_Carcinoma/ACC_Pheno.csv", row.names = 1)

# Select from df_pheno the only col we are interested in:
df_pheno <- df_pheno[,c(1,9)]
# transform alive status into factor
# 0: alive
# 1: dead
df_pheno$patient.vital_status <- ifelse(df_pheno$patient.vital_status == "alive", 0, 1)




#' Normalize raw count dataframe using DeSeq2 and order clinical dataframe
#'
#' DeSeq_Norm returns the count df normalized and its matched clinical dataframe
#'
#' @param df_count df of raw counts with genes on cols and samples on rows
#' @param df_pheno df of clinical observations for each sample of df_count. It
#' should have a column named 'patient.vital_status' containing dependent variables info
#' @returns df_count with genes on cols and samples on rows normalized with DeSeq2 
#' and the ordered matched df_pheno
DeSeq_Norm <- function(df_count, df_pheno){
  
  # Rotate df_count
  df_count <- as.data.frame(t(df_count))
  # match df_count and df_pheno
  m <- match(colnames(df_count), rownames(df_pheno))
  df_pheno_matched <- df_pheno[m, ]
  
  print(rownames(df_pheno_matched) == colnames(df_count))
  
  library(DESeq2)
  
  # Build DDS matrix for DeSeq
  dds <- DESeqDataSetFromMatrix(countData = df_count, 
                                colData = df_pheno_matched, 
                                design =~patient.vital_status)
  # Apply DeSeq to DDS matrix
  dds <- DESeq(dds)
  # Compute Normalized Count Data
  norm_df <- t(as.data.frame(assay(vst(dds))))
  
  return(list(norm_df, df_pheno_matched))
  
}






#' Main function used to simulate the benchmark
#'
#' 
#'
#' @param df_count df of raw counts with genes on cols and samples on rows
#' @param df_pheno df of clinical observations for each sample of df_count. It
#' should have a column named 'patient.vital_status' containing dependent variables info
#' @returns df_count with genes on cols and samples on rows normalized with DeSeq2 
#' and the ordered matched df_pheno
featSel_Bench <- function(df, df_pheno, split = 0.7){
  
  DeSeq_Norm_list <- DeSeq_Norm(df, df_pheno)
  df_norm <- DeSeq_Norm_list[[1]]
  df_desc <- DeSeq_Norm_list[[2]]
  
}





























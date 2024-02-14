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




#' @description
#' Function that given a df_count of input split it into Train and Test
#' 
#' @param df_count df of raw counts with genes on cols and samples on rows
#' @param df_pheno df of clinical observations for each sample of df_count. It
#' should have a column named 'patient.vital_status' containing dependent variables info
#' @returns df_count_train with genes on cols and samples on rows and the ordered matched df_pheno_train
#' @returns df_count_test with genes on cols and samples on rows and the ordered matched df_pheno_test
Split_train_test <- function(df_count, df_pheno, split = 0.8, seed = 123){
  
  library(caret)
  library(dplyr)
  
  # match df_count and df_pheno
  m <- match(rownames(df_count), rownames(df_pheno))
  df_pheno <- df_pheno[m, ]
  
  # Create data partition on df_pheno outcome to create a balanced split
  set.seed(seed)
  train_samples <- df_pheno$patient.vital_status %>%
    createDataPartition(p = split, list = FALSE)
  #print(train_samples)

  df_count_train <- df_count[train_samples, ] 
  df_pheno_train <- df_pheno[train_samples, ]
  
  df_count_test <- df_count[-train_samples, ]
  df_pheno_test <- df_pheno[-train_samples, ]
  
  return(list(df_count_train, df_pheno_train, 
              df_count_test, df_pheno_test))
  
}





#' @description
#' PPLasso pipeline
#' @param df_count_train df of raw counts with genes on cols and samples on rows
#' @param df_pheno_train df of clinical observations for each sample of df_count. It
#' should have a column named 'patient.vital_status' containing dependent variables info
#' @param df_count_test df of raw counts with genes on cols and samples on rows
#' @param df_pheno_test df of clinical observations for each sample of df_count. It
#' should have a column named 'patient.vital_status' containing dependent variables info
#' @returns PPLasso pipeline results on test set
pplasso <- function(df_norm_train, df_pheno_matched_train, df_norm_test, df_pheno_matched_test){
  
  library(PPLasso)
  
}










#' @description
#' Main function used to simulate the benchmark
#' @param df df of raw counts with genes on cols and samples on rows
#' @param df_pheno df of clinical observations for each sample of df_count. It
#' should have a column named 'patient.vital_status' containing dependent variables info
#' @returns pipeline results
featSel_Bench <- function(df, df_pheno){
  
  print("############################################")
  print("Started splitting df into train and test ...")
  # Split df into train and test
  Splitted_df <- Split_train_test(df, df_pheno)
  df_count_train <- Splitted_df[[1]]
  df_pheno_train <- Splitted_df[[2]]
  df_count_test <- Splitted_df[[3]]
  df_pheno_test <- Splitted_df[[4]]
  
  print("Starting Normalization with DeSeq2")
  
  ## Apply Normalization on Train and Test independently
  # Train Normalization
  DeSeq_Norm_list_tr <- DeSeq_Norm(df_count_train, df_pheno_train)
  df_norm_train <- as.data.frame(DeSeq_Norm_list_tr[[1]])
  df_pheno_matched_train <- DeSeq_Norm_list_tr[[2]]
  # Test Normalization
  DeSeq_Norm_list_tst <- DeSeq_Norm(df_count_test, df_pheno_test)
  df_norm_test <- as.data.frame(DeSeq_Norm_list_tst[[1]])
  df_pheno_matched_test <- DeSeq_Norm_list_tst[[2]]
  
  print("Finished successfully")
  print("############################################")
  
  return(list(df_norm_train, df_pheno_matched_train, df_norm_test, df_pheno_matched_test))
  
}




a <- featSel_Bench(df, df_pheno)
























# Import first df
df <- read.csv2("../Data/ACC_Adrenocortical_Carcinoma/ACC_Count.csv", row.names = 1)
df_pheno <- read.csv2("../Data/ACC_Adrenocortical_Carcinoma/ACC_Pheno.csv", row.names = 1)

# Select from df_pheno the only col we are interested in:
df_pheno <- df_pheno[,c(1,9)]
# transform alive status into factor
# 0: alive
# 1: dead
df_pheno$patient.vital_status <- as.factor(ifelse(df_pheno$patient.vital_status == "alive", 0, 1))
colnames(df_pheno) <- c("ID", "condition")





#' @description
#' Function that given a df_count of input split it into Train and Test
#' 
#' @param df_count df of raw counts with genes on cols and samples on rows
#' @param df_pheno df of clinical observations for each sample of df_count. It
#' should have a column named 'Class' containing dependent variables info
#' @returns df_count_train with genes on cols and samples on rows and the ordered matched df_pheno_train
#' @returns df_count_test with genes on cols and samples on rows and the ordered matched df_pheno_test
Split_train_test <- function(df_count, df_pheno, split = 0.7, seed = 123){
  
  library(caret)
  library(dplyr)
  
  # match df_count and df_pheno
  m <- match(rownames(df_count), rownames(df_pheno))
  df_pheno <- df_pheno[m, ]
  
  # Create data partition on df_pheno outcome to create a balanced split
  set.seed(seed)
  train_samples <- df_pheno$condition %>%
    createDataPartition(p = split, list = FALSE)
  #print(train_samples)
  
  df_count_train <- df_count[train_samples, ] 
  df_pheno_train <- df_pheno[train_samples, ]
  
  df_count_test <- df_count[-train_samples, ]
  df_pheno_test <- df_pheno[-train_samples, ]
  
  return(list(df_count_train, df_pheno_train, 
              df_count_test, df_pheno_test))
  
}


#' @description Build S4 objects one for Train one for Test to pass MLSeq
#'
#' @param df_count_train df of raw counts with genes on cols and samples on rows
#' @param df_pheno_train df of clinical observations for each sample of df_count_train. It
#' should have a column named 'Class' containing dependent variables info
#' @param df_count_test df of raw counts with genes on cols and samples on rows
#' @param df_pheno_test df of clinical observations for each sample of df_count_test. It
#' should have a column named 'Class' containing dependent variables info
#' @returns data.trainS4 and data.testS4 to feed MLSeq
returnS4 <- function(df_count_train, df_pheno_train, df_count_test, df_pheno_test, k=10, l=10){
  
  ## Train
  # Rotate df_count
  df_count_train <- as.data.frame(t(df_count_train))
  # match df_count and df_pheno
  m <- match(colnames(df_count_train), rownames(df_pheno_train))
  df_pheno_matched_train <- df_pheno_train[m, ]
  df_pheno_matched_train <- DataFrame(condition = df_pheno_matched_train)
  
  ## Test
  # Rotate df_count
  df_count_test <- as.data.frame(t(df_count_test))
  # match df_count and df_pheno
  n <- match(colnames(df_count_test), rownames(df_pheno_test))
  df_pheno_matched_test <- df_pheno_test[n, ] 
  df_pheno_matched_test <- DataFrame(condition = df_pheno_matched_test)
  
  library(DESeq2)
  
  ## Train
  print("... Building data.trainS4 ...")
  # Build DDS matrix for DeSeq
  data.trainS4 <- DESeqDataSetFromMatrix(countData = as.matrix(df_count_train), 
                                         colData = df_pheno_matched_train,
                                         design = formula(~condition))
  keep <- rowSums(counts(data.trainS4) > k) > l
  data.trainS4 <- data.trainS4[keep, ]
  
  ## Test
  print("... Building data.testS4 ...")
  data.testS4 <- DESeqDataSetFromMatrix(countData = as.matrix(df_count_test), 
                                        colData = df_pheno_matched_test, 
                                        design = formula(~condition))
  # Apply filter computed on Train set
  data.testS4 <- data.testS4[keep, ]
  
  return(list(data.trainS4, data.testS4))
  
}



#' @description
#' This function computes Confusion Matrix of voom-based classifiers
#' @param data.trainS4 S4 object previously created
#' @param data.testS4 S4 object previously created
#' @param cv number of Cross-Validation for hyper parameter tuning
#' @returns Confusion Matrix for each tested algorithm
voomBased <- function(data.trainS4, data.testS4, cv=1){
  
  library(MLSeq)
  library(edgeR)
  
  set.seed(1234)
  # Define control parameters for voom-based algorithms
  ctrl.voom <- voomControl(method = "repeatedcv", number = cv, repeats = 2,
                           tuneLength = 10)
  # Fit train for voomDLDA
  print("Training voomDLDA ...")
  fit.voomDLDA <- classify(data = data.trainS4, method = "voomDLDA",
                           normalize = "none" ,ref = "1", control = ctrl.voom)
  # Fit train for voomDQDA
  print("Training voomDQDA ...")
  fit.voomDQDA <- classify(data = data.trainS4, method = "voomDQDA",
                           normalize = "deseq", ref = "1", control = ctrl.voom)
  # Fit train for voomNSC
  print("Training voomNSC ...")
  fit.voomNSC <- classify(data = data.trainS4, method = "voomNSC",
                          normalize = "deseq", ref = "1", control = ctrl.voom)
  
  return(list(trained(fit.voomDLDA), trained(fit.voomDQDA), trained(fit.voomNSC)))
  
}










################################################################################
################################# MAIN #########################################
################################################################################

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
  
  print("Building S4 objects ...")
  
  S4_objects <- returnS4(df_count_train, df_pheno_train, 
                         df_count_test, df_pheno_test, k=10, l=10)
  data.trainS4 <- S4_objects[[1]]
  data.testS4 <- S4_objects[[2]]
  #prova <- S4_objects[[3]]
  
  # a <- estimateSizeFactors(data.trainS4)
  # a <- estimateDispersions(a)
  # avst <- getVarianceStabilizedData(a)
  
  print(data.trainS4)
  
  print("... Train and Test correctly created!")
  print("Daje de voom")
  voom <- voomBased(data.trainS4, data.testS4)
  
  
  print("Finished successfully")
  print("############################################")
  
  #return(list(data.trainS4, data.testS4, voom))
  return(S4_objects)
  
}

################################################################################
################################################################################
################################################################################


main_return <- featSel_Bench(df, df_pheno)
dtrains4 <- main_return[[1]]

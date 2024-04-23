setwd("/Users/giorgiomontesi/Desktop/Universita_di_Siena/A_PhD_Project/Biomarker_Prediction/ensembleBP/Codes")

options(expressions = 5e5)

library(MLSeq)
library(DESeq2)
library(edgeR)
library(VennDiagram)
library(pamr)
library(caret)
library(devtools)
#install_github("enriquea/feseR", force = T)
library(feseR)


#' @description Import dfCount and dfPheno given their paths
#' @param pathdf relative path to dfCount. genes on rows, samples on columns
#' @param pathclin relative path to dfPheno 2 columns should be named ID and class
#' @returns df: dfCount with samples on cols and genes on rows
#' @returns class: S4 DF dfPheno matched to df with only one col of response variable named condition 
dfs.import <- function(pathdf = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Count.csv",
                       pathclin = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Pheno.csv"){
  
  # Import first df
  df <- read.csv2(pathdf, row.names = 1)
  df_pheno <- read.csv2(pathclin, row.names = 1)
  
  df <- as.data.frame(t(df))
  # Select from df_pheno the only col we are interested in:
  df_pheno <- df_pheno[ ,c("ID", "class")]
  rownames(df_pheno) <- df_pheno$ID
  # transform alive status into factor
  # L: alive
  # D: dead
  # match df_count and df_pheno
  m <- match(colnames(df), rownames(df_pheno))
  df_pheno <- df_pheno[m, ]
  df_pheno$class <- as.factor(df_pheno$class)
  
  df_pheno$class <- as.factor(ifelse(df_pheno$class == 1, "L", "D"))
  df_pheno <- DataFrame(condition = df_pheno$class)
  class <- df_pheno
  
  return(list(df, class))
}


#' @description Split dfCount and Class into train and test according split ratio
#' @param df dfCount as preprocessed from dfs.import
#' @param class S4 dfPheno matched and preprocessed as from dfs.import
#' @param ratio split ratio of test set
#' @param mincorr correlation threshold for soft filter
#' @returns data.trainS4
#' @returns data.testS4
#' @returns classts: real test labels
trainTest.split <- function(df, class, ratio = 0.3, mincorr = 0.1, seed = 123){
  
  set.seed(seed)
  
  data <- df
  nTest <- ceiling(ncol(data) * ratio)
  ind <- sample(ncol(data), nTest, FALSE)
  
  # Minimum count is set to 1 in order to prevent 0 division problem within
  # classification models.
  data.train <- as.matrix(data[ ,-ind] + 1)
  data.test <- as.matrix(data[ ,ind] + 1)
  classtr <- DataFrame(condition = class[-ind, ])
  classts <- DataFrame(condition = class[ind, ])
  
  # Apply very basic correlation filter to train df
  classtr.num <- c(1,0)[classtr$condition]
  dtr <- filter.corr(scale(t(data.train), center = T, scale = T), 
                     classtr.num, mincorr = mincorr)
  dtr <- (t(dtr))
  data.train <- data.train[rownames(data.train) %in% rownames(dtr), ]
  # Apply results to test dataset
  data.test <- data.test[rownames(data.test) %in% rownames(data.train), ]
  
  # Define S4 objects
  data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,
                                        design = formula(~condition))
  data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,
                                       design = formula(~condition))
  
  return(list(data.trainS4, data.testS4, classts))
  
}


#' @description Computes elbow point given a list of named coefficients
#' @param coefficients a Named list of Coefficients to compute elbow point
#' @returns a list of selected features above elbow point 
elbow_comp <- function(coefficients){
  
  library(dplyr)
  library(segmented)
  
  # coefficients are saved in a named vector called coefficients e.g.
  # coefficients <- abs(rnorm(100))  # coeff example
  # names(coefficients) <- as.character(1:100)
  coefficients <- coefficients
  
  # 1. order coefficients in decreasing order
  sorted_coefficients <- sort(coefficients, decreasing = TRUE)
  
  # 2. fit a curve over ordered coefficients
  # build a df with oordered data
  data <- data.frame(index = seq_along(sorted_coefficients), coefficient = sorted_coefficients)
  
  # fit the curve with library(segmented)
  fit <- segmented(lm(coefficient ~ index, data = data), seg.Z = ~ index, psi = list(index = 5))
  
  # 3. Find elbow point over fitted curve
  elbow_point <- round(summary(fit)$psi[[2]])
  
  # Plot data of fitted curve
  #plot(data$index, data$coefficient, type = "l", main = "Elbow Method", xlab = "Index", ylab = "Coefficient")
  plot(sorted_coefficients)
  lines(data$index, predict(fit), col = "indianred")
  abline(v = elbow_point, col = "lightblue", lty = 2)
  legend("topright", legend = c("Data", "Fitted Curve", "Elbow Point"), col = c("black", "red", "blue"), lty = c(1, 1, 2))
  
  # print elbow point
  print(paste("Elbow Point:", elbow_point))
  selected_features <- names(sorted_coefficients)[1:elbow_point]
  
  return(selected_features)
  
}


#' @description Performs SVM-based svm.Radial classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
svm.Radial <- function(data.trainS4, data.testS4, classts, 
                        tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all svm.based classifiers
  svmControl <- trainControl(method = "repeatedcv", number = n,
                             repeats = r, classProbs = TRUE)
  
  print("Fitting SVM-Radial")
  set.seed(1510)
  # Support vector machines with radial basis function kernel
  fit.svmRadial <- classify(data = data.trainS4, method = "svmRadial",
                            preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                            control = svmControl)
  
  #Predicted class labels
  pred.svmRadial <- predict(fit.svmRadial, data.testS4)
  pred.svmRadial <- relevel(pred.svmRadial, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblRadial <- table(Predicted = pred.svmRadial, Actual = actual)
  svmRadial.cm <- confusionMatrix(tblRadial, positive = "D")
  # Compute important genes
  genes_svmRadial_SVMBased <- list(colnames(fit.svmRadial@modelInfo@trainedModel$trainingData[, c(fit.svmRadial@modelInfo@trainedModel$finalModel@SVindex)]))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(svmRadial.cm, genes_svmRadial_SVMBased))
  
}


#' @description Performs SVM-based svm.Poly classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
svm.Poly <- function(data.trainS4, data.testS4, classts, 
                      tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all svm.based classifiers
  svmControl <- trainControl(method = "repeatedcv", number = n,
                             repeats = r, classProbs = TRUE)
  
  print("Fitting SVM-Poly")
  set.seed(1510)
  # Support vector machines with poly basis function kernel
  fit.svmPoly <- classify(data = data.trainS4, method = "svmPoly",
                          preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                          control = svmControl)
  
  #Predicted class labels
  pred.svmPoly <- predict(fit.svmPoly, data.testS4)
  pred.svmPoly <- relevel(pred.svmPoly, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblPoly <- table(Predicted = pred.svmPoly, Actual = actual)
  svmPoly.cm <- confusionMatrix(tblPoly, positive = "D")
  
  # Compute important genes
  genes_svmPoly_SVMBased <- list(colnames(fit.svmPoly@modelInfo@trainedModel$trainingData[, c(fit.svmPoly@modelInfo@trainedModel$finalModel@SVindex)]))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(svmPoly.cm, genes_svmPoly_SVMBased))
  
}


#' @description Performs SVM-based svm.Linear classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
svm.Linear <- function(data.trainS4, data.testS4, classts, 
                        tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all svm.based classifiers
  svmControl <- trainControl(method = "repeatedcv", number = n,
                             repeats = r, classProbs = TRUE)
  
  print("Fitting SVM-Linear")
  set.seed(1510)
  # Support vector machines with linear basis function kernel
  fit.svmLinear <- classify(data = data.trainS4, method = "svmLinear",
                            preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                            control = svmControl)
  
  #Predicted class labels
  pred.svmLinear <- predict(fit.svmLinear, data.testS4)
  pred.svmLinear <- relevel(pred.svmLinear, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblLinear <- table(Predicted = pred.svmLinear, Actual = actual)
  svmLinear.cm <- confusionMatrix(tblLinear, positive = "D")
  
  # Compute important genes
  genes_svmLinear_SVMBased <- list(colnames(fit.svmLinear@modelInfo@trainedModel$trainingData[, c(fit.svmLinear@modelInfo@trainedModel$finalModel@SVindex)]))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(svmLinear.cm, genes_svmLinear_SVMBased))
  
}


#' @description Tests VOOM-based voom.DLDA classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
voom.DLDA <- function(data.trainS4, data.testS4, classts, 
                       tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all voom.based classifiers
  voomControl <- voomControl(method = "repeatedcv", number = n, repeats = r,
                             tuneLength = tL)
  
  print("Fitting voom-DLDA")
  set.seed(1510)
  
  # voomDLDA
  fit.voomDLDA <- classify(data = data.trainS4, method = "voomDLDA",
                           normalize = "deseq", ref = "D",
                           control = voomControl)
  
  #Predicted class labels
  pred.voomDLDA <- predict(fit.voomDLDA, data.testS4)
  pred.voomDLDA <- relevel(pred.voomDLDA, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblDLDA <- table(Predicted = pred.voomDLDA, Actual = actual)
  voomDLDA.cm <- confusionMatrix(tblDLDA, positive = "D")
  
  # Compute elbow genes
  coeff_voomDLDA <- fit.voomDLDA@modelInfo@trainedModel@finalModel[["model"]][["weightedStats"]][["weightedMean"]]
  names(coeff_voomDLDA) <- rownames(fit.voomDLDA@modelInfo@trainedModel@finalModel[["model"]][["weightedStats"]][["weightedMean.C"]])
  genes_voomDLDA_voomBased <- list(elbow_comp(coefficients = coeff_voomDLDA))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(voomDLDA.cm, genes_voomDLDA_voomBased))
  
}
  
  
#' @description Tests VOOM-based voom.DQDA classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
voom.DQDA <- function(data.trainS4, data.testS4, classts, 
                      tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  print("Fitting voom-DQDA")
  
  # Define control function for all voom.based classifiers
  voomControl <- voomControl(method = "repeatedcv", number = n, repeats = r,
                             tuneLength = tL)
  
  set.seed(123)
  # voomDQDA
  fit.voomDQDA <- classify(data = data.trainS4, method = "voomDQDA",
                           normalize = "deseq", ref = "D",
                           control = voomControl)
  
  #Predicted class labels
  pred.voomDQDA <- predict(fit.voomDQDA, data.testS4)
  pred.voomDQDA <- relevel(pred.voomDQDA, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblDQDA <- table(Predicted = pred.voomDQDA, Actual = actual)
  voomDQDA.cm <- confusionMatrix(tblDQDA, positive = "D")
  
  # Compute elbow genes
  coeff_voomDQDA <- fit.voomDQDA@modelInfo@trainedModel@finalModel[["model"]][["weightedStats"]][["weightedMean"]]
  names(coeff_voomDQDA) <- rownames(fit.voomDQDA@modelInfo@trainedModel@finalModel[["model"]][["weightedStats"]][["weightedMean.C"]])
  genes_voomDQDA_voomBased <- list(elbow_comp(coefficients = coeff_voomDQDA))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(voomDQDA.cm, genes_voomDQDA_voomBased))
  
}
  

#' @description Tests VOOM-based voom.NSC classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
voom.NSC <- function(data.trainS4, data.testS4, classts, 
                      tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  print("Fitting voom-NSC")
  
  # Define control function for all voom.based classifiers
  voomControl <- voomControl(method = "repeatedcv", number = n, repeats = r,
                             tuneLength = tL)
  
  set.seed(1510)
  # voomNSC
  fit.voomNSC <- classify(data = data.trainS4, method = "voomNSC",
                          normalize = "deseq", ref = "D",
                          control = voomControl)
  
  #Predicted class labels
  pred.voomNSC <- predict(fit.voomNSC, data.testS4)
  pred.voomNSC <- relevel(pred.voomNSC, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblNSC <- table(Predicted = pred.voomNSC, Actual = actual)
  voomNSC.cm <- confusionMatrix(tblNSC, positive = "D")
  
  # Compute elbow genes
  coeff_voomNSC <- fit.voomNSC@modelInfo@trainedModel@finalModel[["model"]][["weightedMean"]]
  names(coeff_voomNSC) <- rownames(fit.voomNSC@modelInfo@trainedModel@finalModel[["model"]][["weightedMean.C"]])
  genes_voomNSC_voomBased <- list(elbow_comp(coefficients = coeff_voomNSC))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(voomNSC.cm, genes_voomNSC_voomBased))
  
}
  
 
#' @description Perform linear-based PLDA classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
lin.PLDA <- function(data.trainS4, data.testS4, classts, 
                     tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all linear.based classifiers
  linearControl <- discreteControl(method = "repeatedcv", number = n, repeats = r,
                                   tuneLength = tL)
  
  print("Fitting PLDA")
  set.seed(1510)
  # PLDA
  fit.PLDA <- classify(data = data.trainS4, method = "PLDA",
                       normalize = "deseq", ref = "D",
                       control = linearControl)
  
  #Predicted class labels
  pred.PLDA <- predict(fit.PLDA, data.testS4)
  pred.PLDA <- relevel(pred.PLDA, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblPLDA <- table(Predicted = pred.PLDA, Actual = actual)
  PLDA.cm <- confusionMatrix(tblPLDA, positive = "D")
  genes_PLDA_LDABased <- list(selectedGenes(fit.PLDA))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(PLDA.cm, genes_PLDA_LDABased))
  
}
  
 
#' @description Perform linear-based PLDA2 classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
lin.PLDA2 <- function(data.trainS4, data.testS4, classts, 
                      tL = 5, n = 5, r = 5){ 
  
  start.time <- Sys.time()
  
  # Define control function for all linear.based classifiers
  linearControl <- discreteControl(method = "repeatedcv", number = n, repeats = r,
                                   tuneLength = tL)
  
  print("Fitting PLDA2")
  set.seed(1510)
  # PLDA2
  fit.PLDA2 <- classify(data = data.trainS4, method = "PLDA2",
                        normalize = "deseq", ref = "D",
                        control = linearControl)
  
  #Predicted class labels
  pred.PLDA2 <- predict(fit.PLDA2, data.testS4)
  pred.PLDA2 <- relevel(pred.PLDA2, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblPLDA2 <- table(Predicted = pred.PLDA2, Actual = actual)
  PLDA2.cm <- confusionMatrix(tblPLDA2, positive = "D")
  genes_PLDA2_LDABased <- list(selectedGenes(fit.PLDA2))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(PLDA2.cm, genes_PLDA2_LDABased))
  
}
  

#' @description Performs sparseLDA classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
sparse.LDA <- function(data.trainS4, data.testS4, classts, 
                         tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  library(sdwd)
  library(sparseLDA)
  library(spls)
  
  # Define control function for all sparse.based classifiers
  sparseControl <- trainControl(method = "repeatedcv", number = n,
                                repeats = r, classProbs = TRUE)
  
  
  print("Fitting sparseLDA")
  set.seed(1510)
  # Sparse linear discriminant analysis
  fit.sparseLDA <- classify(data = data.trainS4, method = "sparseLDA",
                            preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                            control = sparseControl)
  
  #Predicted class labels
  pred.sparseLDA <- predict(fit.sparseLDA, data.testS4)
  pred.sparseLDA <- relevel(pred.sparseLDA, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblsparseLDA <- table(Predicted = pred.sparseLDA, Actual = actual)
  sparseLDA.cm <- confusionMatrix(tblsparseLDA, positive = "D")
  
  # Compute important genes
  genes_sparseLDA_LDAbased <- list(genes_sparseLDA_LDAbased = fit.sparseLDA@modelInfo@trainedModel[["finalModel"]][["varNames"]])
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(sparseLDA.cm, genes_sparseLDA_LDAbased))
  
}


#' @description Performs tree-based rpart classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
tree.rpart <- function(data.trainS4, data.testS4, classts, 
                       tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all tree.based classifiers
  treeControl <- trainControl(method = "repeatedcv", number = n,
                              repeats = r, classProbs = TRUE, 
                              savePredictions = "all", returnData = T)
  
  print("Fitting rpart")
  set.seed(1510)
  # Sparse Distance Weighted Discrimination
  fit.rpart <- classify(data = data.trainS4, method = "rpart",
                        preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                        control = treeControl)
  
  #Predicted class labels
  pred.rpart <- predict(fit.rpart, data.testS4)
  pred.rpart <- relevel(pred.rpart, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblrpart <- table(Predicted = pred.rpart, Actual = actual)
  rpart.cm <- confusionMatrix(tblrpart, positive = "D")
  genes_rpart_treeBased <- list(fit.rpart@modelInfo@trainedModel[["finalModel"]][["variable.importance"]])
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(rpart.cm, genes_rpart_treeBased))
  
}
  
  
#' @description Performs tree-based cforest classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
tree.cforest <- function(data.trainS4, data.testS4, classts, 
                          tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all tree.based classifiers
  treeControl <- trainControl(method = "repeatedcv", number = n,
                              repeats = r, classProbs = TRUE, 
                              savePredictions = "all", returnData = T)
  
  print("Fitting cforest")
  set.seed(1510)
  # Sparse linear discriminant analysis
  fit.cforest <- classify(data = data.trainS4, method = "cforest",
                          preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                          control = treeControl)
  
  
  #Predicted class labels
  pred.cforest <- predict(fit.cforest, data.testS4)
  pred.cforest <- relevel(pred.cforest, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblcforest <- table(Predicted = pred.cforest, Actual = actual)
  cforest.cm <- confusionMatrix(tblcforest, positive = "D")
  
  # Compute elbow genes
  coeff_cforest <- varImp(fit.cforest@modelInfo@trainedModel[["finalModel"]])$Overall
  names(coeff_cforest) <- rownames(varImp(fit.cforest@modelInfo@trainedModel[["finalModel"]]))
  genes_cforest_treeBased <- list(elbow_comp(coefficients = coeff_cforest))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(cforest.cm, genes_cforest_treeBased))
  
}
  
  
#' @description Performs tree-based cforest classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
tree.ctree <- function(data.trainS4, data.testS4, classts, 
                        tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all tree.based classifiers
  treeControl <- trainControl(method = "repeatedcv", number = n,
                              repeats = r, classProbs = TRUE, 
                              savePredictions = "all", returnData = T)
  
  print("Fitting ctree")
  set.seed(1510)
  # Sparse partial least squares
  fit.ctree <- classify(data = data.trainS4, method = "ctree",
                        preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                        control = treeControl)
  
  #Predicted class labels
  pred.ctree <- predict(fit.ctree, data.testS4)
  pred.ctree <- relevel(pred.ctree, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblctree <- table(Predicted = pred.ctree, Actual = actual)
  ctree.cm <- confusionMatrix(tblctree, positive = "D")
  
  # compute elbow genes
  coeff_ctree <- varImp(fit.ctree@modelInfo@trainedModel)[["importance"]]$D
  names(coeff_ctree) <- rownames(varImp(fit.ctree@modelInfo@trainedModel)[["importance"]])
  genes_ctree_treeBased <- list(elbow_comp(coefficients = coeff_ctree))
  # fit.ctree@modelInfo@trainedModel[["finalModel"]]@tree
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(ctree.cm, genes_ctree_treeBased))
  
}
  
  
#' @description Performs tree-based Random Forest classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
tree.rf <- function(data.trainS4, data.testS4, classts, 
                    tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all tree.based classifiers
  treeControl <- trainControl(method = "repeatedcv", number = n,
                              repeats = r, classProbs = TRUE, 
                              savePredictions = "all", returnData = T) 
  
  print("Fitting rf")
  set.seed(1510)
  # Sparse partial least squares
  fit.rf <- classify(data = data.trainS4, method = "rf",
                     preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                     control = treeControl)
  
  #Predicted class labels
  pred.rf <- predict(fit.rf, data.testS4)
  pred.rf <- relevel(pred.rf, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblrf <- table(Predicted = pred.rf, Actual = actual)
  rf.cm <- confusionMatrix(tblrf, positive = "D")
  
  # Compute elbow genes
  coeff_rf <- as.data.frame(fit.rf@modelInfo@trainedModel[["finalModel"]][["importance"]])$MeanDecreaseGini
  names(coeff_rf) <- rownames(as.data.frame(fit.rf@modelInfo@trainedModel[["finalModel"]][["importance"]]))
  genes_rf_treeBased <- list(elbow_comp(coefficients = coeff_rf))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(rf.cm, genes_rf_treeBased))
  
}


#' @description Performs bagging-based AdaBag classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
bagg.AdaBag <- function(data.trainS4, data.testS4, classts, 
                        tL = 3, n = 3, r = 3){
  
  start.time <- Sys.time()
  
  library(adabag)
  library(earth)
  
  # Define control function for all bagg.based classifiers
  baggControl <- trainControl(method = "repeatedcv", number = n,
                              repeats = r, classProbs = TRUE,
                              savePredictions = "all", returnData = T,
                              allowParallel = F)
  
  print("Fitting AdaBag")
  set.seed(1510)
  # Sparse Distance Weighted Discrimination
  fit.AdaBag <- classify(data = data.trainS4, method = "AdaBag",
                         preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                         control = baggControl)
  
  #Predicted class labels
  pred.AdaBag <- predict(fit.AdaBag, data.testS4)
  pred.AdaBag <- relevel(pred.AdaBag, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblAdaBag <- table(Predicted = pred.AdaBag, Actual = actual)
  AdaBag.cm <- confusionMatrix(tblAdaBag, positive = "D")
  genes_AdaBag_baggedBased <- list(fit.AdaBag@modelInfo@trainedModel[["finalModel"]][["importance"]])
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(AdaBag.cm, genes_AdaBag_baggedBased))
  
}
  

#' @description Performs bagging-based treebag classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
bagg.treebag <- function(data.trainS4, data.testS4, classts, 
                        tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  library(adabag)
  library(earth)
  
  # Define control function for all bagg.based classifiers
  baggControl <- trainControl(method = "repeatedcv", number = n,
                              repeats = r, classProbs = TRUE,
                              savePredictions = "all", returnData = T)
  
  print("Fitting treebag")
  set.seed(1510)
  # Sparse linear discriminant analysis
  fit.treebag <- classify(data = data.trainS4, method = "treebag",
                          preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                          control = baggControl)
  
  #Predicted class labels
  pred.treebag <- predict(fit.treebag, data.testS4)
  pred.treebag <- relevel(pred.treebag, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tbltreebag <- table(Predicted = pred.treebag, Actual = actual)
  treebag.cm <- confusionMatrix(tbltreebag, positive = "D")
  
  # Compute elbow genes
  coeff_treebag <- varImp(fit.treebag@modelInfo@trainedModel[["finalModel"]])$Overall
  names(coeff_treebag) <- rownames(varImp(fit.treebag@modelInfo@trainedModel[["finalModel"]]))
  genes_treebag_baggedBased <- list(elbow_comp(coefficients = coeff_treebag))
  # fit.treebag@modelInfo@trainedModel[["finalModel"]][["mtrees"]]
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(treebag.cm, genes_treebag_baggedBased))
  
}

 
#' @description Performs bagging-based bagFDA classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
bagg.bagFDA <- function(data.trainS4, data.testS4, classts, 
                        tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  library(adabag)
  library(earth)
  
  # Define control function for all bagg.based classifiers
  baggControl <- trainControl(method = "repeatedcv", number = n,
                              repeats = r, classProbs = TRUE,
                              savePredictions = "all", returnData = T)
  
  print("Fitting bagFDA")
  set.seed(1510)
  # Sparse partial least squares
  fit.bagFDA <- classify(data = data.trainS4, method = "bagFDA",
                         preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                         control = baggControl)
  
  #Predicted class labels
  pred.bagFDA <- predict(fit.bagFDA, data.testS4)
  pred.bagFDA <- relevel(pred.bagFDA, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblbagFDA <- table(Predicted = pred.bagFDA, Actual = actual)
  bagFDA.cm <- confusionMatrix(tblbagFDA, positive = "D")
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(bagFDA.cm))
  
}


#' @description Performs boost-based gamboost classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
boost.gamboost <- function(data.trainS4, data.testS4, classts, 
                          tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all boost.based classifiers
  boostControl <- trainControl(method = "repeatedcv", number = n,
                               repeats = r, classProbs = TRUE)
  
  print("Fitting gamboost")
  set.seed(1510)
  # Sparse Distance Weighted Discrimination
  fit.gamboost <- classify(data = data.trainS4, method = "gamboost",
                           preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                           control = boostControl)
  
  #Predicted class labels
  pred.gamboost <- predict(fit.gamboost, data.testS4)
  pred.gamboost <- relevel(pred.gamboost, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblgamboost <- table(Predicted = pred.gamboost, Actual = actual)
  gamboost.cm <- confusionMatrix(tblgamboost, positive = "D")
  genes_gamboost_boostBased <- list(names(coef(fit.gamboost@modelInfo@trainedModel[["finalModel"]])))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(gamboost.cm, genes_gamboost_boostBased))
  
}
  
  
#' @description Performs boost-based bstSm classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
boost.bstSm <- function(data.trainS4, data.testS4, classts, 
                        tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all boost.based classifiers
  boostControl <- trainControl(method = "repeatedcv", number = n,
                               repeats = r, classProbs = TRUE)  
  
  print("Fitting bstSm")
  set.seed(1510)
  # Sparse linear discriminant analysis
  fit.bstSm <- classify(data = data.trainS4, method = "bstSm",
                        preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                        control = boostControl)
  
  #Predicted class labels
  pred.bstSm <- predict(fit.bstSm, data.testS4)
  pred.bstSm <- relevel(pred.bstSm, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblbstSm <- table(Predicted = pred.bstSm, Actual = actual)
  bstSm.cm <- confusionMatrix(tblbstSm, positive = "D")
  genes_bstSm_boostBased <- list(fit.bstSm@modelInfo@trainedModel[["finalModel"]][["xNames"]][c(fit.bstSm@modelInfo@trainedModel[["finalModel"]][["xselect"]])])
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(bstSm.cm, genes_bstSm_boostBased))
  
}
  
 
#' @description Performs boost-based bstSm classifier
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
boost.bstTree <- function(data.trainS4, data.testS4, classts, 
                          tL = 3, n = 3, r = 3){
  
  start.time <- Sys.time()
  
  # Define control function for all boost.based classifiers
  boostControl <- trainControl(method = "repeatedcv", number = n,
                               repeats = r, classProbs = TRUE)  
  
  print("Fitting bstTree")
  set.seed(1510)
  # Sparse partial least squares
  fit.bstTree <- classify(data = data.trainS4, method = "bstTree",
                          preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                          control = boostControl)
  
  #Predicted class labels
  pred.bstTree <- predict(fit.bstTree, data.testS4)
  pred.bstTree <- relevel(pred.bstTree, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblbstTree <- table(Predicted = pred.bstTree, Actual = actual)
  bstTree.cm <- confusionMatrix(tblbstTree, positive = "D")
  genes_bstTree_boostBased <- list(fit.bstTree@modelInfo@trainedModel[["finalModel"]][["xNames"]][c(fit.bstTree@modelInfo@trainedModel[["finalModel"]][["xselect"]])])
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(bstTree.cm, genes_bstTree_boostBased))
  
}


#' @description Perform pls-based gpls classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
pls.gpls <- function(data.trainS4, data.testS4, classts, 
                      tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  library(gpls)
  
  # Define control function for all pls.based classifiers
  plsControl <- trainControl(method = "repeatedcv", number = n,
                             repeats = r, classProbs = TRUE)
  
  print("Fitting gpls")
  set.seed(1510)
  # generalized partial least squares
  fit.gpls <- classify(data = data.trainS4, method = "gpls",
                       preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                       control = plsControl)
  
  #Predicted class labels
  pred.gpls <- predict(fit.gpls, data.testS4)
  pred.gpls <- relevel(pred.gpls, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblgpls <- table(Predicted = pred.gpls, Actual = actual)
  gpls.cm <- confusionMatrix(tblgpls, positive = "D")
  # compute elbow genes
  coeff_gpls <- abs(fit.gpls@modelInfo@trainedModel[["finalModel"]][["coefficients"]][-1])
  genes_gpls_plsBased <- list(elbow_comp(coefficients = coeff_gpls))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(gpls.cm, genes_gpls_plsBased))
  
}
  
 
#' @description Perform pls-based pls classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
pls.pls <- function(data.trainS4, data.testS4, classts, 
                    tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all pls.based classifiers
  plsControl <- trainControl(method = "repeatedcv", number = n,
                             repeats = r, classProbs = TRUE) 
  
  print("Fitting pls")
  set.seed(1510)
  # partial least squares
  fit.pls <- classify(data = data.trainS4, method = "pls",
                      preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                      control = plsControl)
  
  #Predicted class labels
  pred.pls <- predict(fit.pls, data.testS4)
  pred.pls <- relevel(pred.pls, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblpls <- table(Predicted = pred.pls, Actual = actual)
  pls.cm <- confusionMatrix(tblpls, positive = "D")
  
  #compute elbow genes
  coeff_pls <- abs(fit.pls@modelInfo@trainedModel[["finalModel"]][["coefficients"]][,1,1])
  genes_pls_plsBased <- list(elbow_comp(coefficients = coeff_pls))
  # genes_pls_plsBased <- list(c("ADA", "Gino"))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(pls.cm, genes_pls_plsBased))
  
}
  
  
#' @description Perform pls-based pls classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list containing a Confusion Matrix and a list of important genes
pls.spls <- function(data.trainS4, data.testS4, classts, 
                     tL = 5, n = 5, r = 5){
  
  start.time <- Sys.time()
  
  # Define control function for all pls.based classifiers
  plsControl <- trainControl(method = "repeatedcv", number = n,
                             repeats = r, classProbs = TRUE)   
  
  print("Fitting SPLS")
  set.seed(1510)
  # Sparse partial least squares
  fit.spls <- classify(data = data.trainS4, method = "spls",
                       preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                       control = plsControl)
  
  #Predicted class labels
  pred.spls <- predict(fit.spls, data.testS4)
  pred.spls <- relevel(pred.spls, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  
  tblspls <- table(Predicted = pred.spls, Actual = actual)
  spls.cm <- confusionMatrix(tblspls, positive = "D")
  # compute elbow genes
  coeff_spls <- abs(fit.spls@modelInfo@trainedModel[["finalModel"]][["normx"]])
  genes_spls_plsBased <- list(elbow_comp(coefficients = coeff_spls))
  
  # Computation time
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("Accomplished in ", time.taken, "secs"))
  
  return(list(spls.cm, genes_spls_plsBased))
  
}


seed=123
#' @description
#' Function to perform classification prediction non-CV
#' @param seed seed to inizialize analysis
#' @param mincorr correlation threshold to filter genes
#' @param pathdf relative path to count table
#' @param pathclin relative path to desc table
#' @param method list of methods to use for the analysis
crossVal.1layer <- function(seed, i, mincorr = 0.39, 
                            pathdf = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Count.csv",
                            pathclin = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Pheno.csv",
                            methods = c("all")){
  
  print("Importing specified datasets")
  dfsImport <- dfs.import(pathdf = pathdf, pathclin = pathclin)
  df <- as.data.frame(t(dfsImport[[1]]))
  class <- dfsImport[[2]]
  
  print("Performing mini features trimming")
  # Keep only features with at least 10 counts in 1/3 of samples
  keep <- rowSums(df > 10) > round(ncol(df)/3)
  df <- df[keep, ]
  
  print("Splitting datasets into train and test")
  tts <- trainTest.split(df, class, mincorr = mincorr, seed = seed)
  data.trainS4 <- tts[[1]]
  data.testS4 <- tts[[2]]
  classts <- tts[[3]]
  
  # mini-check to see if filtered genes are the same among train and test
  print(sum(rownames(assay(data.trainS4)) == rownames(assay(data.testS4))))
  
  print("Starting classification tasks")
  ## Start running Classifiers
  
  if ("svmRadial" %in% methods | "svmBased" %in% methods | "all" %in% methods){ 
    svmR <- svm.Radial(data.trainS4, data.testS4, classts)
    svmRadial <- svmR[[1]]
    genes_svmRadial_SVMBased <- svmR[[2]]
  }
  
  if ("svmPoly" %in% methods | "svmBased" %in% methods | "all" %in% methods){
    svmP <- svm.Poly(data.trainS4, data.testS4, classts)
    svmPoly <- svmP[[1]]
    genes_svmPoly_SVMBased <- svmP[[2]]
  }
  
  if ("svmLinear" %in% methods | "svmBased" %in% methods | "all" %in% methods){
    svmL <- svm.Linear(data.trainS4, data.testS4, classts)
    svmLinear <- svmL[[1]]
    genes_svmLinear_SVMBased <- svmL[[2]]
  }
  
  if ("voomDLDA" %in% methods | "voomBased" %in% methods | "all" %in% methods){
    vDLDA <- voom.DLDA(data.trainS4, data.testS4, classts)
    voomDLDA <- vDLDA[[1]]
    genes_voomDLDA_voomBased <- vDLDA[[2]]
  }
  
  if ("voomDQDA" %in% methods | "voomBased" %in% methods | "all" %in% methods){
    vDQDA <- voom.DQDA(data.trainS4, data.testS4, classts)
    voomDQDA <- vDQDA[[1]]
    genes_voomDQDA_voomBased <- vDQDA[[2]]
  }
  
  if ("voomNSC" %in% methods | "voomBased" %in% methods | "all" %in% methods){
    vNSC <- voom.NSC(data.trainS4, data.testS4, classts)
    voomNSC <- vNSC[[1]]
    genes_voomNSC_voomBased <- vNSC[[2]]
  }
  
  if ("PLDA" %in% methods | "LDABased" %in% methods | "all" %in% methods){
    linPLDA <- lin.PLDA(data.trainS4, data.testS4, classts)
    PLDA <- linPLDA[[1]]
    genes_PLDA_LDABased <- linPLDA[[2]]
  }
  
  if ("PLDA2" %in% methods | "LDABased" %in% methods | "all" %in% methods){
    linPLDA2 <- lin.PLDA2(data.trainS4, data.testS4, classts)
    PLDA2 <- linPLDA2[[1]]
    genes_PLDA2_LDABased <- linPLDA2[[2]]
  }

  if ("sparseLDA" %in% methods | "LDABased" %in% methods | "all" %in% methods){
    sparse <- sparse.LDA(data.trainS4, data.testS4, classts) # <-- too slow!!
    sparseLDA <- sparse[[1]]
    genes_sparseLDA_LDABased <- sparse[[2]]
  }

  if ("rpart" %in% methods | "treeBased" %in% methods | "all" %in% methods){
    treerpart <- tree.rpart(data.trainS4, data.testS4, classts)
    rpart <- treerpart[[1]]
    genes_rpart_treeBased <- treerpart[[2]]
    genes_rpart_treeBased <- list(names(genes_rpart_treeBased[[1]]))
  }
  
  if ("cforest" %in% methods | "treeBased" %in% methods | "all" %in% methods){
    treecforest <- tree.cforest(data.trainS4, data.testS4, classts)
    cforest <- treecforest[[1]]
    genes_cforest_treeBased <- treecforest[[2]]
  }
  
  if ("ctree" %in% methods | "treeBased" %in% methods | "all" %in% methods){
    treectree <- tree.ctree(data.trainS4, data.testS4, classts)
    ctree <- treectree[[1]]
    genes_ctree_treeBased <- treectree[[2]]
  }
  
  if ("rf" %in% methods | "treeBased" %in% methods | "all" %in% methods){
    treerf <- tree.rf(data.trainS4, data.testS4, classts)
    RF <- treerf[[1]]
    genes_rf_treeBased <- treerf[[2]]
  }
  
  if ("AdaBag" %in% methods | "baggedBased" %in% methods | "all" %in% methods){
    bagAda <- bagg.AdaBag(data.trainS4, data.testS4, classts)
    AdaBag <- bagAda[[1]]
    genes_AdaBag_baggedBased <- bagAda[[2]]
    genes_AdaBag_baggedBased <- list(names(genes_AdaBag_baggedBased[[1]][genes_AdaBag_baggedBased[[1]]>0]))
  }
  
  if ("treebag" %in% methods | "baggedBased" %in% methods | "all" %in% methods){
    bagtree <- bagg.treebag(data.trainS4, data.testS4, classts)
    treebag <- bagtree[[1]]
    genes_treebag_baggedBased <- bagtree[[2]]
  }
  
  if ("bagFDA" %in% methods | "baggedBased" %in% methods | "all" %in% methods){
    baggfda <- bagg.bagFDA(data.trainS4, data.testS4, classts)
    BagFDA <- baggfda[[1]]
  }

  if ("gamboost" %in% methods | "boostBased" %in% methods | "all" %in% methods){
    bstgam <- boost.gamboost(data.trainS4, data.testS4, classts)
    gamboost <- bstgam[[1]]
    genes_gamboost_boostBased <- bstgam[[2]]
    genes_gamboost_boostBased[[1]] <- sub("bbs\\(([^,]+),.*", "\\1", genes_gamboost_boostBased[[1]])
  }
  
  if ("bstSm" %in% methods | "boostBased" %in% methods | "all" %in% methods){
    bstbstsm <- boost.bstSm(data.trainS4, data.testS4, classts)
    bstSm <- bstbstsm[[1]]
    genes_bstSm_boostBased <- bstbstsm[[2]]
  }
  
  if ("bstTree" %in% methods | "boostBased" %in% methods | "all" %in% methods){
    bsttree <- boost.bstTree(data.trainS4, data.testS4, classts)
    bstTree <- bsttree[[1]]
    genes_bstTree_boostBased <- bsttree[[2]]
  }
 
  if ("gpls" %in% methods | "plsBased" %in% methods | "all" %in% methods){
    plsgpls <- pls.gpls(data.trainS4, data.testS4, classts)
    gpls <- plsgpls[[1]]
    genes_gpls_plsBased <- plsgpls[[2]]
  }
  
  if ("pls" %in% methods | "plsBased" %in% methods | "all" %in% methods){
    plspls <- pls.pls(data.trainS4, data.testS4, classts)
    pls <- plspls[[1]]
    genes_pls_plsBased <- plspls[[2]]
  }
  
  if ("spls" %in% methods | "plsBased" %in% methods | "all" %in% methods){
    plsspls <- pls.spls(data.trainS4, data.testS4, classts)
    spls <- plsspls[[1]]
    genes_spls_plsBased <- plsspls[[2]]
  }
  
  
  # Build Performance-Metrics Table
  acc.df <- data.frame(
    svmRadial = if (exists("svmRadial")) c(svmRadial$overall, svmRadial$byClass) else NA,
    svmPoly = if (exists("svmPoly")) c(svmPoly$overall, svmPoly$byClass) else NA,
    svmLinear = if (exists("svmLinear")) c(svmLinear$overall, svmLinear$byClass) else NA,
    voomDLDA = if (exists("voomDLDA")) c(voomDLDA$overall, voomDLDA$byClass) else NA,
    voomDQDA = if (exists("voomDQDA")) c(voomDQDA$overall, voomDQDA$byClass) else NA,
    voomNSC = if (exists("voomNSC")) c(voomNSC$overall, voomNSC$byClass) else NA,
    PLDA = if (exists("PLDA")) c(PLDA$overall, PLDA$byClass) else NA,
    PLDA2 = if (exists("PLDA2")) c(PLDA2$overall, PLDA2$byClass) else NA,
    sparseLDA = if (exists("sparseLDA")) c(sparseLDA$overall, sparseLDA$byClass) else NA,
    rpart = if (exists("rpart")) c(rpart$overall, rpart$byClass) else NA,
    cforest = if (exists("cforest")) c(cforest$overall, cforest$byClass) else NA,
    ctree = if (exists("ctree")) c(ctree$overall, ctree$byClass) else NA,
    rf = if (exists("RF")) c(RF$overall, RF$byClass) else NA,
    AdaBag = if (exists("AdaBag")) c(AdaBag$overall, AdaBag$byClass) else NA,
    treebag = if (exists("treebag")) c(treebag$overall, treebag$byClass) else NA,
    bagFDA = if (exists("BagFDA")) c(BagFDA$overall, BagFDA$byClass) else NA,
    gamboost = if (exists("gamboost")) c(gamboost$overall, gamboost$byClass) else NA,
    bstSm = if (exists("bstSm")) c(bstSm$overall, bstSm$byClass) else NA,
    bstTree = if (exists("bstTree")) c(bstTree$overall, bstTree$byClass) else NA,
    gpls = if (exists("gpls")) c(gpls$overall, gpls$byClass) else NA,
    pls = if (exists("pls")) c(pls$overall, pls$byClass) else NA,
    spls = if (exists("spls")) c(spls$overall, spls$byClass) else NA
  )
  
  # Save Metrics-Performance Table
  write.csv2(t(acc.df), paste0("../Results/eboplusAccuracyTable_",i,".csv"))

  
  # Build Selected genes list of list
  list_genes <- list(
    genes_svmRadial_SVMBased = if (exists("genes_svmRadial_SVMBased")) genes_svmRadial_SVMBased[[1]] else NA,
    genes_svmPoly_SVMBased = if (exists("genes_svmPoly_SVMBased")) genes_svmPoly_SVMBased[[1]] else NA,
    genes_svmLinear_SVMBased = if (exists("genes_svmLinear_SVMBased")) genes_svmLinear_SVMBased[[1]] else NA,
    genes_voomDLDA_voomBased = if (exists("genes_voomDLDA_voomBased")) genes_voomDLDA_voomBased[[1]] else NA,
    genes_voomDQDA_voomBased = if (exists("genes_voomDQDA_voomBased")) genes_voomDQDA_voomBased[[1]] else NA,
    genes_voomNSC_voomBased = if (exists("genes_voomNSC_voomBased")) genes_voomNSC_voomBased[[1]] else NA,
    genes_PLDA_LDABased = if (exists("genes_PLDA_LDABased")) genes_PLDA_LDABased[[1]] else NA,
    genes_PLDA2_LDABased = if (exists("genes_PLDA2_LDABased")) genes_PLDA2_LDABased[[1]] else NA,
    genes_sparseLDA_LDABased = if (exists("genes_sparseLDA_LDABased")) genes_sparseLDA_LDABased[[1]] else NA,
    genes_rpart_treeBased = if (exists("genes_rpart_treeBased")) genes_rpart_treeBased[[1]] else NA,
    genes_cforest_treeBased = if (exists("genes_cforest_treeBased")) genes_cforest_treeBased[[1]] else NA,
    genes_ctree_treeBased = if (exists("genes_ctree_treeBased")) genes_ctree_treeBased[[1]] else NA,
    genes_rf_treeBased = if (exists("genes_rf_treeBased")) genes_rf_treeBased[[1]] else NA,
    genes_AdaBag_baggedBased = if (exists("genes_AdaBag_baggedBased")) genes_AdaBag_baggedBased[[1]] else NA,
    genes_treebag_baggedBased = if (exists("genes_treebag_baggedBased")) genes_treebag_baggedBased[[1]] else NA,
    genes_gamboost_boostBased = if (exists("genes_gamboost_boostBased")) genes_gamboost_boostBased[[1]] else NA,
    genes_bstSm_boostBased = if (exists("genes_bstSm_boostBased")) genes_bstSm_boostBased[[1]] else NA,
    genes_bstTree_boostBased = if (exists("genes_bstTree_boostBased")) genes_bstTree_boostBased[[1]] else NA,
    genes_gpls_plsBased = if (exists("genes_gpls_plsBased")) genes_gpls_plsBased[[1]] else NA,
    genes_pls_plsBased = if (exists("genes_pls_plsBased")) genes_pls_plsBased[[1]] else NA,
    genes_spls_plsBased = if (exists("genes_spls_plsBased")) genes_spls_plsBased[[1]] else NA
  )
  
  
  # Save list of list of genes as RDS object
  saveRDS(list_genes, paste0("../Results/ebopluslist_genes_",i,".rds"))
  
  return(list(acc.df, list_genes))
}

#' @description
#' Function to perform n-fold CV classification prediction
#' @param seed seed to inizialize analysis
#' @param mincorr correlation threshold to filter genes
#' @param pathdf relative path to count table
#' @param pathclin relative path to desc table
#' @param method list of methods to use for the analysis
#' @param cv number of cross validation to perform
ensembleBP <- function(mincorr = 0.4, cv = 10, 
                       pathdf = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Count.csv",
                       pathclin = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Pheno.csv",
                       methods = c("all")){
  
  
  
  for (i in c(1:cv)) {
    
    print(paste0("Performing Cross-Validation of ",i," layer"))
    crossVal.1layer(seed = i+5, i = i, mincorr = mincorr, 
                    pathdf = pathdf,
                    pathclin = pathclin,
                    methods = methods)
    
  }
  
  
}




## EXAMPLE USAGE:

ensembleBP(mincorr = 0.4, cv = 5, 
           pathdf = "../Data/Eboplus/count.day1.csv",
           pathclin = "../Data/Eboplus/desc_ab_class.csv",
           methods = c("svmBased", "plsBased", "rpart",
                       "voomBased", "rf"))

a <- readRDS("../Results/provalist_genes_1.rds")
b <- readRDS("../Results/provalist_genes_2.rds")










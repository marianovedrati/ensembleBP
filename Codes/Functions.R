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
#' @param pathdf relative path to dfCount
#' @param pathclin relative path to dfPheno
#' @returns df: dfCount with samples on cols and genes on rows
#' @returns class: S4 DF dfPheno matched to df with only one col of response variable named condition 
dfs.import <- function(pathdf = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Count.csv",
                       pathclin = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Pheno.csv"){
  
  # Import first df
  df <- read.csv2(pathdf, row.names = 1)
  df_pheno <- read.csv2(pathclin, row.names = 1)
  
  df <- as.data.frame(t(df))
  # Select from df_pheno the only col we are interested in:
  df_pheno <- df_pheno[,c(1,9)]
  # transform alive status into factor
  # L: alive
  # D: dead
  # match df_count and df_pheno
  m <- match(colnames(df), rownames(df_pheno))
  df_pheno <- df_pheno[m, ]
  
  df_pheno$patient.vital_status <- as.factor(ifelse(df_pheno$patient.vital_status == "alive", "L", "D"))
  df_pheno <- DataFrame(condition = df_pheno$patient.vital_status)
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


#' #' @description Tests several NNet-based classifiers
#' #' @param data.trainS4
#' #' @param data.testS4
#' #' @param classts
#' #' @param tL tune Length
#' #' @param n number of CV
#' #' @param r number of repeats for CV
#' #' @returns A list of Confusion Matrices one for each NNet-based method
# nnet.based <- function(data.trainS4, data.testS4, classts,
#                          tL = 2, n = 2, r = 2){
#   
#   start.time <- Sys.time()
#   # Define control function for all NNet.based classifiers
#   nnetControl <- trainControl(method = "repeatedcv", number = n,
#                                 repeats = r, classProbs = TRUE)
# 
#   print("Fitting nnet")
#   # set.seed(1510)
#   # # Sparse Distance Weighted Discrimination
#   # fit.nnet <- classify(data = data.trainS4, method = "nnet",
#   #                      preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
#   #                      control = nnetControl)
#   #
#   # #Predicted class labels
#   # pred.nnet <- predict(fit.nnet, data.testS4)
#   # pred.nnet <- relevel(pred.nnet, ref = "D")
#   actual <- relevel(classts$condition, ref = "D")
# 
#   # tblnnet <- table(Predicted = pred.nnet, Actual = actual)
#   # nnet.cm <- confusionMatrix(tblnnet, positive = "D")
#   nnet.cm=1
# 
#   print("Fitting mlp")
#   set.seed(1510)
#   # Sparse linear discriminant analysis
#   fit.mlp <- classify(data = data.trainS4, method = "mlp",
#                             preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
#                             control = nnetControl)
# 
#   #Predicted class labels
#   pred.mlp <- predict(fit.mlp, data.testS4)
#   pred.mlp <- relevel(pred.mlp, ref = "D")
# 
#   tblmlp <- table(Predicted = pred.mlp, Actual = actual)
#   mlp.cm <- confusionMatrix(tblmlp, positive = "D")
# 
#   print("Fitting mlpML")
#   set.seed(1510)
#   # Sparse partial least squares
#   fit.mlpML <- classify(data = data.trainS4, method = "mlpML",
#                        preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
#                        control = nnetControl)
# 
#   #Predicted class labels
#   pred.mlpML <- predict(fit.mlpML, data.testS4)
#   pred.mlpML <- relevel(pred.mlpML, ref = "D")
# 
#   tblmlpML <- table(Predicted = pred.mlpML, Actual = actual)
#   mlpML.cm <- confusionMatrix(tblmlpML, positive = "D")
# 
#   # print("Fitting avNNet")
#   # set.seed(1510)
#   # # Sparse partial least squares
#   # fit.avNNet <- classify(data = data.trainS4, method = "avNNet",
#   #                       preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
#   #                       control = nnetControl)
#   #
#   # #Predicted class labels
#   # pred.avNNet <- predict(fit.avNNet, data.testS4)
#   # pred.avNNet <- relevel(pred.avNNet, ref = "D")
#   #
#   # tblavNNet <- table(Predicted = pred.avNNet, Actual = actual)
#   # avNNet.cm <- confusionMatrix(tblavNNet, positive = "D")
#   avNNet.cm = 1
# 
#   print("Successfully accomplished NNet-based methods")
#
#   # Computation time
#   end.time <- Sys.time()
#   time.taken <- round(end.time - start.time,2)
#   print(paste0("Accomplished in ", time.taken, "secs"))
#
#   return(list(nnet.cm, mlp.cm, mlpML.cm, avNNet.cm))
# 
# }


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
crossVal.1layer <- function(seed, i, mincorr = 0.4, 
                            pathdf = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Count.csv",
                            pathclin = "../Data/ACC_Adrenocortical_Carcinoma/ACC_Pheno.csv"){
  
  print("Importing specified datasets")
  dfsImport <- dfs.import(pathdf = pathdf, pathclin = pathclin)
  df <- dfsImport[[1]]
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
  svmR <- svm.Radial(data.trainS4, data.testS4, classts)
  svmRadial <- svmR[[1]]
  genes_svmRadial_SVMBased <- svmR[[2]]
  
  svmP <- svm.Poly(data.trainS4, data.testS4, classts)
  svmPoly <- svmP[[1]]
  genes_svmPoly_SVMBased <- svmP[[2]]
  
  svmL <- svm.Linear(data.trainS4, data.testS4, classts)
  svmLinear <- svmL[[1]]
  genes_svmLinear_SVMBased <- svmL[[2]]
  
  vDLDA <- voom.DLDA(data.trainS4, data.testS4, classts)
  voomDLDA <- vDLDA[[1]]
  genes_voomDLDA_voomBased <- vDLDA[[2]]
  
  vDQDA <- voom.DQDA(data.trainS4, data.testS4, classts)
  voomDQDA <- vDQDA[[1]]
  genes_voomDQDA_voomBased <- vDQDA[[2]]
  
  vNSC <- voom.NSC(data.trainS4, data.testS4, classts)
  voomNSC <- vNSC[[1]]
  genes_voomNSC_voomBased <- vNSC[[2]]
  
  linPLDA <- lin.PLDA(data.trainS4, data.testS4, classts)
  PLDA <- linPLDA[[1]]
  genes_PLDA_LDABased <- linPLDA[[2]]
  
  linPLDA2 <- lin.PLDA2(data.trainS4, data.testS4, classts)
  PLDA2 <- linPLDA2[[1]]
  genes_PLDA2_LDABased <- linPLDA2[[2]]

  sparse <- sparse.LDA(data.trainS4, data.testS4, classts) # <-- too slow!!
  sparseLDA <- sparse[[1]]
  genes_sparseLDA_LDABased <- sparse[[2]]
  
  # net <- nnet.based(data.trainS4, data.testS4, classts) # <-- not properly working!!
  # #nnet <- net[[1]]
  # mlp <- net[[2]]
  # mlpML <- net[[3]]
  # #avNNet <- net[[4]]
  
  treerpart <- tree.rpart(data.trainS4, data.testS4, classts)
  rpart <- treerpart[[1]]
  genes_rpart_treeBased <- treerpart[[2]]
  genes_rpart_treeBased <- list(names(genes_rpart_treeBased[[1]]))
  
  treecforest <- tree.cforest(data.trainS4, data.testS4, classts)
  cforest <- treecforest[[1]]
  genes_cforest_treeBased <- treecforest[[2]]
  
  treectree <- tree.ctree(data.trainS4, data.testS4, classts)
  ctree <- treectree[[1]]
  genes_ctree_treeBased <- treectree[[2]]
  
  treerf <- tree.rf(data.trainS4, data.testS4, classts)
  rf <- treerf[[1]]
  genes_rf_treeBased <- treerf[[2]]
  
  bagAda <- bagg.AdaBag(data.trainS4, data.testS4, classts)
  AdaBag <- bagAda[[1]]
  genes_AdaBag_baggedBased <- bagAda[[2]]
  genes_AdaBag_baggedBased <- list(names(genes_AdaBag_baggedBased[[1]][genes_AdaBag_baggedBased[[1]]>0]))
  
  bagtree <- bagg.treebag(data.trainS4, data.testS4, classts)
  treebag <- bagtree[[1]]
  genes_treebag_baggedBased <- bagtree[[2]]
  
  baggfda <- bagg.bagFDA(data.trainS4, data.testS4, classts)
  bagFDA <- baggfda[[1]]

  bstgam <- boost.gamboost(data.trainS4, data.testS4, classts)
  gamboost <- bstgam[[1]]
  genes_gamboost_boostBased <- bstgam[[2]]
  genes_gamboost_boostBased[[1]] <- sub("bbs\\(([^,]+),.*", "\\1", genes_gamboost_boostBased[[1]])
  
  bstbstsm <- boost.bstSm(data.trainS4, data.testS4, classts)
  bstSm <- bstbstsm[[1]]
  genes_bstSm_boostBased <- bstbstsm[[2]]
  
  bsttree <- boost.bstTree(data.trainS4, data.testS4, classts)
  bstTree <- bsttree[[1]]
  genes_bstTree_boostBased <- bsttree[[2]]
 
  plsgpls <- pls.gpls(data.trainS4, data.testS4, classts)
  gpls <- plsgpls[[1]]
  genes_gpls_plsBased <- plsgpls[[2]]
  
  plspls <- pls.pls(data.trainS4, data.testS4, classts)
  pls <- plspls[[1]]
  genes_pls_plsBased <- plspls[[2]]
  
  plsspls <- pls.spls(data.trainS4, data.testS4, classts)
  spls <- plsspls[[1]]
  genes_spls_plsBased <- plsspls[[2]]
  
  # Build Performance-Metrics Table
  acc.df <- data.frame(svmRadial = c(svmRadial$overall, svmRadial$byClass), 
                       svmPoly = c(svmPoly$overall, svmPoly$byClass), 
                       svmLinear = c(svmLinear$overall, svmLinear$byClass),
                       voomDLDA = c(voomDLDA$overall, voomDLDA$byClass),
                       voomDQDA = c(voomDQDA$overall, voomDQDA$byClass),
                       voomNSC = c(voomNSC$overall, voomNSC$byClass),
                       PLDA = c(PLDA$overall, PLDA$byClass),
                       PLDA2 = c(PLDA2$overall, PLDA2$byClass),
                       sparseLDA = c(sparseLDA$overall, sparseLDA$byClass),
                       # nnet = c(nnet$overall, nnet$byClass),
                       # mlp = c(mlp$overall, mlp$byClass),
                       # mlpML = c(mlpML$overall, mlpML$byClass),
                       # avNNet = c(avNNet$overall, avNNet$byClass),
                       rpart = c(rpart$overall, rpart$byClass),
                       cforest = c(cforest$overall, cforest$byClass),
                       ctree = c(ctree$overall, ctree$byClass),
                       rf = c(rf$overall, rf$byClass),
                       AdaBag = c(AdaBag$overall, AdaBag$byClass),
                       treebag = c(treebag$overall, treebag$byClass),
                       bagFDA = c(bagFDA$overall, bagFDA$byClass),
                       gamboost = c(gamboost$overall, gamboost$byClass),
                       bstSm = c(bstSm$overall, bstSm$byClass),
                       bstTree = c(bstTree$overall, bstTree$byClass),
                       gpls = c(gpls$overall, gpls$byClass),
                       pls = c(pls$overall, pls$byClass),
                       spls = c(spls$overall, spls$byClass))
  
  # Save Metrics-Performance Table
  write.csv2(t(acc.df), paste0("../Results/provaAccuracyTable_",i,".csv"))
  
  # Build Selected genes list of list
  list_genes <- list(genes_svmRadial_SVMBased = genes_svmRadial_SVMBased[[1]],
                     genes_svmPoly_SVMBased = genes_svmPoly_SVMBased[[1]],
                     genes_svmLinear_SVMBased = genes_svmLinear_SVMBased[[1]],
                     genes_voomDLDA_voomBased = genes_voomDLDA_voomBased[[1]],
                     genes_voomDQDA_voomBased = genes_voomDQDA_voomBased[[1]],
                     genes_voomNSC_voomBased = genes_voomNSC_voomBased[[1]],
                     genes_PLDA_LDABased = genes_PLDA_LDABased[[1]],
                     genes_PLDA2_LDABased = genes_PLDA2_LDABased[[1]],
                     genes_sparseLDA_LDABased = genes_sparseLDA_LDABased[[1]],
                     genes_rpart_treeBased = genes_rpart_treeBased[[1]],
                     genes_cforest_treeBased = genes_cforest_treeBased[[1]],
                     genes_ctree_treeBased = genes_ctree_treeBased[[1]],
                     genes_rf_treeBased = genes_rf_treeBased[[1]],
                     genes_AdaBag_baggedBased = genes_AdaBag_baggedBased[[1]],
                     genes_treebag_baggedBased = genes_treebag_baggedBased[[1]],
                     genes_gamboost_boostBased = genes_gamboost_boostBased[[1]],
                     genes_bstSm_boostBased = genes_bstSm_boostBased[[1]],
                     genes_bstTree_boostBased = genes_bstTree_boostBased[[1]],
                     genes_gpls_plsBased = genes_gpls_plsBased[[1]],
                     genes_pls_plsBased = genes_pls_plsBased[[1]],
                     genes_spls_plsBased = genes_spls_plsBased[[1]])
  
  # Save list of list of genes as RDS object
  saveRDS(list_genes, paste0("../Results/provalist_genes_",i,".rds"))
  
  return(list(acc.df, list_genes))
}

i = 1
cv <- 10
for (i in c(1:cv)) {
  
  print(paste0("Performing Cross-Validation of ",i," layer"))
  crossVal.1layer(seed = i, i = i, mincorr = 0.42)
  
}

a <- readRDS(file = "../Results/provalist_genes_1.rds")


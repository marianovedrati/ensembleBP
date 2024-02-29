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


#' @description Tests several SVM-based classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list of Confusion Matrices one for each svm-based method
svm.based <- function(data.trainS4, data.testS4, classts, 
                      tL = 2, n = 5, r = 2){
  
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
  
  print("Fitting SVM-Poly")
  set.seed(1510)
  # Support vector machines with poly basis function kernel
  fit.svmPoly <- classify(data = data.trainS4, method = "svmPoly",
                          preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                          control = svmControl)
  
  #Predicted class labels
  pred.svmPoly <- predict(fit.svmPoly, data.testS4)
  pred.svmPoly <- relevel(pred.svmPoly, ref = "D")
  
  tblPoly <- table(Predicted = pred.svmPoly, Actual = actual)
  svmPoly.cm <- confusionMatrix(tblPoly, positive = "D")
  
  print("Fitting SVM-Linear")
  set.seed(1510)
  # Support vector machines with linear basis function kernel
  fit.svmLinear <- classify(data = data.trainS4, method = "svmLinear",
                            preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                            control = svmControl)
  
  #Predicted class labels
  pred.svmLinear <- predict(fit.svmLinear, data.testS4)
  pred.svmLinear <- relevel(pred.svmLinear, ref = "D")
  
  tblLinear <- table(Predicted = pred.svmLinear, Actual = actual)
  svmLinear.cm <- confusionMatrix(tblLinear, positive = "D")
  
  print("Successfully accomplished SVM-based methods")
  
  return(list(svmRadial.cm, svmPoly.cm, svmLinear.cm))

}



#' @description Tests several VOOM-based classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list of Confusion Matrices one for each VOOM-based method
voom.based <- function(data.trainS4, data.testS4, classts, 
                       tL = 2, n = 5, r = 2){
  
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
  
  print("Fitting voom-DQDA")
  set.seed(1510)
  # voomDQDA
  fit.voomDQDA <- classify(data = data.trainS4, method = "voomDQDA",
                           normalize = "deseq", ref = "D",
                           control = voomControl)
  
  #Predicted class labels
  pred.voomDQDA <- predict(fit.voomDQDA, data.testS4)
  pred.voomDQDA <- relevel(pred.voomDQDA, ref = "D")
  
  tblDQDA <- table(Predicted = pred.voomDQDA, Actual = actual)
  voomDQDA.cm <- confusionMatrix(tblDQDA, positive = "D")
  
  print("Fitting voom-NSC")
  set.seed(1510)
  # voomNSC
  fit.voomNSC <- classify(data = data.trainS4, method = "voomNSC",
                           normalize = "deseq", ref = "D",
                           control = voomControl)
  
  #Predicted class labels
  pred.voomNSC <- predict(fit.voomNSC, data.testS4)
  pred.voomNSC <- relevel(pred.voomNSC, ref = "D")
  
  tblNSC <- table(Predicted = pred.voomNSC, Actual = actual)
  voomNSC.cm <- confusionMatrix(tblNSC, positive = "D")
  
  print("Successfully accomplished VOOM-based methods")
  
  return(list(voomDLDA.cm, voomDQDA.cm, voomNSC.cm))
  
}



#' @description Tests several linear-based classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list of Confusion Matrices one for each linear-based method
linear.based <- function(data.trainS4, data.testS4, classts, 
                          tL = 2, n = 5, r = 2){
  
  # Define control function for all voom.based classifiers
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
  
  print("Fitting PLDA2")
  set.seed(1510)
  # PLDA2
  fit.PLDA2 <- classify(data = data.trainS4, method = "PLDA2",
                        normalize = "deseq", ref = "D",
                        control = linearControl)
  
  #Predicted class labels
  pred.PLDA2 <- predict(fit.PLDA2, data.testS4)
  pred.PLDA2 <- relevel(pred.PLDA2, ref = "D")
  
  tblPLDA2 <- table(Predicted = pred.PLDA2, Actual = actual)
  PLDA2.cm <- confusionMatrix(tblPLDA2, positive = "D")
  
  print("Fitting NBLDA")
  set.seed(1510)
  # NBLDA
  fit.NBLDA <- classify(data = data.trainS4, method = "NBLDA",
                        normalize = "deseq", ref = "D",
                        control = linearControl)
  
  #Predicted class labels
  pred.NBLDA <- predict(fit.NBLDA, data.testS4)
  pred.NBLDA <- relevel(pred.NBLDA, ref = "D")
  
  tblNBLDA <- table(Predicted = pred.NBLDA, Actual = actual)
  NBLDA.cm <- confusionMatrix(tblNBLDA, positive = "D")
  
  print("Successfully accomplished linear-based methods")
  
  return(list(PLDA.cm, PLDA2.cm, NBLDA.cm))
  
}



#' @description Tests several sparse-based classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list of Confusion Matrices one for each sparse-based method
sparse.based <- function(data.trainS4, data.testS4, classts, 
                         tL = 2, n = 2, r = 2){
  
  library(sdwd)
  library(sparseLDA)
  library(spls)
  
  # Define control function for all sparse.based classifiers
  sparseControl <- trainControl(method = "repeatedcv", number = n,
                                repeats = r, classProbs = TRUE)
  
  print("Fitting SDWD")
  # set.seed(1510)
  # # Sparse Distance Weighted Discrimination
  # fit.sdwd <- classify(data = data.trainS4, method = "sdwd",
  #                           preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
  #                           control = sparseControl)
  # 
  # #Predicted class labels
  # pred.sdwd <- predict(fit.sdwd, data.testS4)
  # pred.sdwd <- relevel(pred.sdwd, ref = "D")
  actual <- relevel(classts$condition, ref = "D")
  # 
  # tblsdwd <- table(Predicted = pred.sdwd, Actual = actual)
  # sdwd.cm <- confusionMatrix(tblsdwd, positive = "D")
  sdwd.cm <- 1
  
  print("Fitting sparseLDA")
  set.seed(1510)
  # Sparse linear discriminant analysis
  fit.sparseLDA <- classify(data = data.trainS4, method = "sparseLDA",
                          preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                          control = sparseControl)
  
  #Predicted class labels
  pred.sparseLDA <- predict(fit.sparseLDA, data.testS4)
  pred.sparseLDA <- relevel(pred.sparseLDA, ref = "D")
  
  tblsparseLDA <- table(Predicted = pred.sparseLDA, Actual = actual)
  sparseLDA.cm <- confusionMatrix(tblsparseLDA, positive = "D")
  
  print("Fitting SPLS")
  # set.seed(1510)
  # # Sparse partial least squares
  # fit.spls <- classify(data = data.trainS4, method = "spls",
  #                           preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
  #                           control = sparseControl)
  # 
  # #Predicted class labels
  # pred.spls <- predict(fit.spls, data.testS4)
  # pred.spls <- relevel(pred.spls, ref = "D")
  # 
  # tblspls <- table(Predicted = pred.spls, Actual = actual)
  # spls.cm <- confusionMatrix(tblspls, positive = "D")
  spls.cm <- 1
  
  print("Successfully accomplished sparse-based methods")
  
  return(list(sdwd.cm, sparseLDA.cm, spls.cm))
  
}



#' @description Tests several NNet-based classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list of Confusion Matrices one for each NNet-based method
nnet.based <- function(data.trainS4, data.testS4, classts, 
                         tL = 2, n = 2, r = 2){
  
  # Define control function for all NNet.based classifiers
  nnetControl <- trainControl(method = "repeatedcv", number = n,
                                repeats = r, classProbs = TRUE)
  
  print("Fitting nnet")
  # set.seed(1510)
  # # Sparse Distance Weighted Discrimination
  # fit.nnet <- classify(data = data.trainS4, method = "nnet",
  #                      preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
  #                      control = nnetControl)
  # 
  # #Predicted class labels
  # pred.nnet <- predict(fit.nnet, data.testS4)
  # pred.nnet <- relevel(pred.nnet, ref = "D")
  actual <- relevel(classts$condition, ref = "D")

  # tblnnet <- table(Predicted = pred.nnet, Actual = actual)
  # nnet.cm <- confusionMatrix(tblnnet, positive = "D")
  nnet.cm=1
  
  print("Fitting mlp")
  set.seed(1510)
  # Sparse linear discriminant analysis
  fit.mlp <- classify(data = data.trainS4, method = "mlp",
                            preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                            control = nnetControl)
  
  #Predicted class labels
  pred.mlp <- predict(fit.mlp, data.testS4)
  pred.mlp <- relevel(pred.mlp, ref = "D")
  
  tblmlp <- table(Predicted = pred.mlp, Actual = actual)
  mlp.cm <- confusionMatrix(tblmlp, positive = "D")
  
  print("Fitting mlpML")
  set.seed(1510)
  # Sparse partial least squares
  fit.mlpML <- classify(data = data.trainS4, method = "mlpML",
                       preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                       control = nnetControl)
  
  #Predicted class labels
  pred.mlpML <- predict(fit.mlpML, data.testS4)
  pred.mlpML <- relevel(pred.mlpML, ref = "D")
  
  tblmlpML <- table(Predicted = pred.mlpML, Actual = actual)
  mlpML.cm <- confusionMatrix(tblmlpML, positive = "D")
  
  # print("Fitting avNNet")
  # set.seed(1510)
  # # Sparse partial least squares
  # fit.avNNet <- classify(data = data.trainS4, method = "avNNet",
  #                       preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
  #                       control = nnetControl)
  # 
  # #Predicted class labels
  # pred.avNNet <- predict(fit.avNNet, data.testS4)
  # pred.avNNet <- relevel(pred.avNNet, ref = "D")
  # 
  # tblavNNet <- table(Predicted = pred.avNNet, Actual = actual)
  # avNNet.cm <- confusionMatrix(tblavNNet, positive = "D")
  avNNet.cm = 1
  
  print("Successfully accomplished NNet-based methods")
  
  return(list(nnet.cm, mlp.cm, mlpML.cm, avNNet.cm))
  
}


#' @description Tests several tree-based classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list of Confusion Matrices one for each tree-based method
tree.based <- function(data.trainS4, data.testS4, classts, 
                       tL = 2, n = 2, r = 2){
  
  # Define control function for all NNet.based classifiers
  treeControl <- trainControl(method = "repeatedcv", number = n,
                              repeats = r, classProbs = TRUE)
  
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
  
  print("Fitting cforest")
  set.seed(1510)
  # Sparse linear discriminant analysis
  fit.cforest <- classify(data = data.trainS4, method = "cforest",
                      preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                      control = treeControl)
  
  #Predicted class labels
  pred.cforest <- predict(fit.cforest, data.testS4)
  pred.cforest <- relevel(pred.cforest, ref = "D")
  
  tblcforest <- table(Predicted = pred.cforest, Actual = actual)
  cforest.cm <- confusionMatrix(tblcforest, positive = "D")
  
  print("Fitting ctree")
  set.seed(1510)
  # Sparse partial least squares
  fit.ctree <- classify(data = data.trainS4, method = "ctree",
                        preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                        control = treeControl)
  
  #Predicted class labels
  pred.ctree <- predict(fit.ctree, data.testS4)
  pred.ctree <- relevel(pred.ctree, ref = "D")
  
  tblctree <- table(Predicted = pred.ctree, Actual = actual)
  ctree.cm <- confusionMatrix(tblctree, positive = "D")
  
  print("Fitting rf")
  set.seed(1510)
  # Sparse partial least squares
  fit.rf <- classify(data = data.trainS4, method = "rf",
                        preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                        control = treeControl)

  #Predicted class labels
  pred.rf <- predict(fit.rf, data.testS4)
  pred.rf <- relevel(pred.rf, ref = "D")

  tblrf <- table(Predicted = pred.rf, Actual = actual)
  rf.cm <- confusionMatrix(tblrf, positive = "D")
  
  print("Successfully accomplished tree-based methods")
  
  return(list(rpart.cm, cforest.cm, ctree.cm, rf.cm))
  
}




#' @description Tests several bagging-based classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list of Confusion Matrices one for each bagging-based method
bagg.based <- function(data.trainS4, data.testS4, classts, 
                         tL = 2, n = 2, r = 2){
  
  library(adabag)
  library(earth)
  
  # Define control function for all bagg.based classifiers
  baggControl <- trainControl(method = "repeatedcv", number = n,
                                repeats = r, classProbs = TRUE)
  
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
  
  print("Fitting treebag")
  set.seed(1510)
  # Sparse linear discriminant analysis
  fit.treebag <- classify(data = data.trainS4, method = "treebag",
                            preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                            control = baggControl)
  
  #Predicted class labels
  pred.treebag <- predict(fit.treebag, data.testS4)
  pred.treebag <- relevel(pred.treebag, ref = "D")
  
  tbltreebag <- table(Predicted = pred.treebag, Actual = actual)
  treebag.cm <- confusionMatrix(tbltreebag, positive = "D")
  
  print("Fitting bagFDA")
  set.seed(1510)
  # Sparse partial least squares
  fit.bagFDA <- classify(data = data.trainS4, method = "bagFDA",
                       preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                       control = baggControl)
  
  #Predicted class labels
  pred.bagFDA <- predict(fit.bagFDA, data.testS4)
  pred.bagFDA <- relevel(pred.bagFDA, ref = "D")
  
  tblbagFDA <- table(Predicted = pred.bagFDA, Actual = actual)
  bagFDA.cm <- confusionMatrix(tblbagFDA, positive = "D")
  
  print("Successfully accomplished bagging-based methods")
  
  return(list(AdaBag.cm, treebag.cm, bagFDA.cm))
  
}



#' @description Tests several boost-based classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list of Confusion Matrices one for each boost-based method
boost.based <- function(data.trainS4, data.testS4, classts, 
                       tL = 2, n = 2, r = 2){
  
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
  
  print("Fitting bstSm")
  set.seed(1510)
  # Sparse linear discriminant analysis
  fit.bstSm <- classify(data = data.trainS4, method = "bstSm",
                          preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                          control = boostControl)
  
  #Predicted class labels
  pred.bstSm <- predict(fit.bstSm, data.testS4)
  pred.bstSm <- relevel(pred.bstSm, ref = "D")
  
  tblbstSm <- table(Predicted = pred.bstSm, Actual = actual)
  bstSm.cm <- confusionMatrix(tblbstSm, positive = "D")
  
  print("Fitting bstTree")
  set.seed(1510)
  # Sparse partial least squares
  fit.bstTree <- classify(data = data.trainS4, method = "bstTree",
                         preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                         control = boostControl)
  
  #Predicted class labels
  pred.bstTree <- predict(fit.bstTree, data.testS4)
  pred.bstTree <- relevel(pred.bstTree, ref = "D")
  
  tblbstTree <- table(Predicted = pred.bstTree, Actual = actual)
  bstTree.cm <- confusionMatrix(tblbstTree, positive = "D")
  
  print("Successfully accomplished boost-based methods")
  
  return(list(gamboost.cm, bstSm.cm, bstTree.cm))
  
}




#' @description Tests several pls-based classifiers
#' @param data.trainS4
#' @param data.testS4
#' @param classts
#' @param tL tune Length
#' @param n number of CV
#' @param r number of repeats for CV
#' @returns A list of Confusion Matrices one for each pls-based method
pls.based <- function(data.trainS4, data.testS4, classts, 
                         tL = 2, n = 2, r = 2){
  
  library(gpls)
  
  # Define control function for all sparse.based classifiers
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
  
  print("Fitting pls")
  set.seed(1510)
  # partial least squares
  fit.pls <- classify(data = data.trainS4, method = "pls",
                            preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                            control = plsControl)
  
  #Predicted class labels
  pred.pls <- predict(fit.pls, data.testS4)
  pred.pls <- relevel(pred.pls, ref = "D")
  
  tblpls <- table(Predicted = pred.pls, Actual = actual)
  pls.cm <- confusionMatrix(tblpls, positive = "D")
  
  print("Fitting SPLS")
  set.seed(1510)
  # Sparse partial least squares
  fit.spls <- classify(data = data.trainS4, method = "spls",
                       preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                       control = plsControl)
  
  #Predicted class labels
  pred.spls <- predict(fit.spls, data.testS4)
  pred.spls <- relevel(pred.spls, ref = "D")
  
  tblspls <- table(Predicted = pred.spls, Actual = actual)
  spls.cm <- confusionMatrix(tblspls, positive = "D")
  
  print("Successfully accomplished pls-based methods")
  
  return(list(gpls.cm, pls.cm, spls.cm))
  
}





dfsImport <- dfs.import()
df <- dfsImport[[1]]
#df <- df[1:1500, ]
class <- dfsImport[[2]]

keep <- rowSums(df > 10) > round(ncol(df)/3)
df <- df[keep, ]


crossVal.1layer <- function(seed, i, mincorr = 0.3){
  
  tts <- trainTest.split(df, class, mincorr = mincorr, seed = seed)
  data.trainS4 <- tts[[1]]
  data.testS4 <- tts[[2]]
  classts <- tts[[3]]
  # mini-check per vedere se i geni filtati sono gli stessi
  sum(rownames(assay(data.trainS4)) == rownames(assay(data.testS4)))

  svm <- svm.based(data.trainS4, data.testS4, classts)
  svmRadial <- svm[[1]] 
  svmPoly <- svm[[2]]
  svmLinear <- svm[[3]]

  voom <- voom.based(data.trainS4, data.testS4, classts)
  voomDLDA <- voom[[1]]
  voomDQDA <- voom[[2]]
  voomNSC <- voom[[3]]

  lin <- linear.based(data.trainS4, data.testS4, classts)
  PLDA <- lin[[1]]
  PLDA2 <- lin[[2]]
  NBLDA <- lin[[3]]

  sparse <- sparse.based(data.trainS4, data.testS4, classts) # <-- too slow!!
  #sdwd <- sparse[[1]]
  sparseLDA <- sparse[[2]]
  #spls <- sparse[[3]]

  net <- nnet.based(data.trainS4, data.testS4, classts) # <-- not properly working!!
  #nnet <- net[[1]]
  mlp <- net[[2]]
  mlpML <- net[[3]]
  #avNNet <- net[[4]]

  tree <- tree.based(data.trainS4, data.testS4, classts)
  rpart <- tree[[1]]
  cforest <- tree[[2]]
  ctree <- tree[[3]]
  rf <- tree[[4]]

  bag <- bagg.based(data.trainS4, data.testS4, classts)
  AdaBag <- bag[[1]]
  treebag <- bag[[2]]
  bagFDA <- bag[[3]]

  bst <- boost.based(data.trainS4, data.testS4, classts)
  gamboost <- bst[[1]]
  bstSm <- bst[[2]]
  bstTree <- bst[[3]]

  partls <- pls.based(data.trainS4, data.testS4, classts)
  gpls <- partls[[1]]
  pls <- partls[[2]]
  spls <- partls[[3]]


  acc.df <- data.frame(svmRadial = c(svmRadial$overall, svmRadial$byClass), 
                       svmPoly = c(svmPoly$overall, svmPoly$byClass), 
                       svmLinear = c(svmLinear$overall, svmLinear$byClass),
                       voomDLDA = c(voomDLDA$overall, voomDLDA$byClass),
                       voomDQDA = c(voomDQDA$overall, voomDQDA$byClass),
                       voomNSC = c(voomNSC$overall, voomNSC$byClass),
                       PLDA = c(PLDA$overall, PLDA$byClass),
                       PLDA2 = c(PLDA2$overall, PLDA2$byClass),
                       NBLDA = c(NBLDA$overall, NBLDA$byClass),
                       sparseLDA = c(sparseLDA$overall, sparseLDA$byClass),
                       # nnet = c(nnet$overall, nnet$byClass),
                       mlp = c(mlp$overall, mlp$byClass),
                       mlpML = c(mlpML$overall, mlpML$byClass),
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
  #i=1
  write.csv2(t(acc.df), paste0("../Results/AccuracyTable_",i,".csv"))
}


cv <- 5
for (i in c(1:cv)) {
  
  print(paste0("Performing Cross-Validation of ",i," layer"))
  crossVal.1layer(seed = i, i = i, mincorr = 0.3)
  
}



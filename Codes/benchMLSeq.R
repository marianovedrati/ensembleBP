setwd("/Users/giorgiomontesi/Desktop/Universita_di_Siena/A_PhD_Project/Biomarker_Prediction/ensembleBP/Codes")

library(MLSeq)
library(DESeq2)
library(edgeR)
library(VennDiagram)
library(pamr)
library(caret)


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
#' @returns data.trainS4
#' @returns data.testS4
#' @returns classts: real test labels
trainTest.split <- function(df, class, ratio = 0.3){

  set.seed(2128)
  data <- df
  nTest <- ceiling(ncol(data) * ratio)
  ind <- sample(ncol(data), nTest, FALSE)

  # Minimum count is set to 1 in order to prevent 0 division problem within
  # classification models.
  data.train <- as.matrix(data[ ,-ind] + 1)
  data.test <- as.matrix(data[ ,ind] + 1)
  classtr <- DataFrame(condition = class[-ind, ])
  classts <- DataFrame(condition = class[ind, ])


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
  
  print(selectedGenes(fit.svmPoly))
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




dfsImport <- dfs.import()
df <- dfsImport[[1]]
class <- dfsImport[[2]]

tts <- trainTest.split(df, class)
data.trainS4 <- tts[[1]]
data.testS4 <- tts[[2]]
classts <- tts[[3]]

svm <- svm.based(data.trainS4, data.testS4, classts)
voom <- voom.based(data.trainS4, data.testS4, classts)








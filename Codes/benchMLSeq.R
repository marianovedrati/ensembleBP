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
  
  svmControl <- trainControl(method = "repeatedcv", number = n,
                             repeats = r, classProbs = TRUE)
  set.seed(2128)
  # Support vector machines with radial basis function kernel
  fit.svmRadial <- classify(data = data.trainS4, method = "svmRadial",
                            preProcessing = "deseq-vst", ref = "D", tuneLength = tL,
                            control = svmControl)

  # show(fit.svm)
  # trained(fit.svm)
  # plot(fit.svm)
  #Predicted class labels
  pred.svmRadial <- predict(fit.svmRadial, data.testS4)
  pred.svmRadial <- relevel(pred.svmRadial, ref = "D")
  actual <- relevel(classts$condition, ref = "D")

  tblRadial <- table(Predicted = pred.svmRadial, Actual = actual)
  svmRadial.cm <- confusionMatrix(tblRadial, positive = "D")
  
  return(list(svmRadial.cm))

}




dfsImport <- dfs.import()
df <- dfsImport[[1]]
class <- dfsImport[[2]]

tts <- trainTest.split(df, class)
data.trainS4 <- tts[[1]]
data.testS4 <- tts[[2]]
classts <- tts[[3]]

svm <- svm.based(data.trainS4, data.testS4, classts)









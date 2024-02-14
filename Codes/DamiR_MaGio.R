setwd("C:/Users/Maria/Documents/TRANSVAC/TRANSVAC_LEISHMANIA/Analisi_Dati")
library(tidyverse)
#importo dati
Count <- read.csv2("Count_table_Leishmania.csv",row.names=1)
colnames(Count)
Count <- Count[,-61]
row.names(Count) <- gsub("-", "", rownames(Count))
Count <- Count[,-61]
Desc <- read.csv2("DESC_tab_PKDL.csv")



#match per avere count e Desc ordinati nello stesso modo
m <- match(colnames(Count),Desc$Sample_ID)
Desc <- Desc[m,]
colnames(Count) == Desc$Sample_ID


#filtering samples with less than 5* 10^6 counts
keep_sample<-colSums(Count)>4000000
Count<-Count[,keep_sample]
Desc<-Desc[keep_sample,]


Desc_d0 <- Desc[Desc$Timepoint == "0",]
Count_d0 <- Count[, colnames(Count) %in% Desc_d0$Sample_ID]
colnames(Count_d0) == Desc_d0$Sample_ID
Desc_d0$class <- character(length(Desc_d0$Coing))

for (i in 1:nrow(Desc)) {
  if (Desc_d0$Coing[i] %in% c("VR", "PR")) {
    Desc_d0$class[i] <- "R"
  } else if (Desc_d0$Coing[i] %in% c("VN", "PN")) {
    Desc_d0$class[i] <- "NR"
  }
}
row.names(Desc_d0) <- Desc_d0$Sample_ID
Desc_d0$Sample_ID<-as.factor(Desc_d0$Sample_ID)
Desc_d0$Timepoint<-as.factor(Desc_d0$Timepoint)
Desc_d0$age_group<-as.factor(Desc_d0$age_group)
Desc_d0$Randomization<-as.factor(Desc_d0$Randomization)
Desc_d0$Coing<-as.factor(Desc_d0$Coing)
Desc_d0$X.initia.disease.d90<-as.factor(Desc_d0$X.initia.disease.d90)
Desc_d0$class<-as.factor(Desc_d0$class)


library(DaMiRseq)
library(knitr)
SE<-DaMiR.makeSE(x=Count_d0, y=Desc_d0)


####SALTA da qui a...###############
# normalizzo con vst,controllo le amples ridondanti
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7,
                                 hyper = "yes",th.cv=3)

#This step introduces a sample quality checkpoint
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.9)


#Adjusting Data
Desc_d0 <- Desc_d0[!row.names(Desc_d0) %in% ("Y031_D0"),] 
Desc_d0$class <- as.factor(Desc_d0$class)


sv <- DaMiR.SV(data_filt)

pr <- colData(data_filt)
t <- DaMiR.corrplot(sv,pr,type="spearman",sig.level = 0.01)
class(Desc_d0) #non torna il corr plot

#aggiunto per sv
data_adjust<-DaMiR.SVadjust(data_filt, sv, n.sv=9)
DaMiR.Allplot(data_adjust, colData(data_adjust))
#######

# Dataset for prediction
set.seed(10101)
#Sampl_cl1 <- 5
#nSampl_cl2 <- 5
## May create unbalanced Learning and Test sets
idx_test <- sample(1:ncol(data_adjust), 10)
# Create balanced Learning and Test sets
#idx_test_cl1<-sample(1:(ncol(data_adjust)/2), nSampl_cl1)
#idx_test_cl2<-sample(1:(ncol(data_adjust)/2), nSampl_cl2) + ncol(data_adjust)/2
#idx_test <- c(idx_test_cl1, idx_test_cl2)
Test_set <- data_adjust[, idx_test, drop=FALSE]
Learning_set <- data_adjust[, -idx_test, drop=FALSE]

# Training and Test into a 'nfold' Cross Validation
nfold <- 9
cv_sample <- c(rep(seq_len(nfold), each=ncol(Learning_set)/(2*nfold)),
               rep(seq_len(nfold), each=ncol(Learning_set)/(2*nfold)))
# Variables initialization
cv_models <- list()
cv_predictors <- list()
res_df <- data.frame(matrix(nrow = nfold, ncol = 7))
colnames(res_df) <- c("Accuracy",
                      "N.predictors",
                      "MCC",
                      "sensitivity",
                      "Specificty",
                      "PPV",
                      "NPV")

cv_fold=1
for (cv_fold in seq_len(nfold)){
  # Create Training and Validation Sets
  set.seed(cv_fold)
  idx_cv<-sample(1:ncol(Learning_set), 30)
  
  #idx_cv <- which(cv_sample != cv_fold)
  TR_set <- Learning_set[,idx_cv, drop=FALSE]
  Val_set <- Learning_set[,-idx_cv, drop=FALSE]
  #### Feature selection
  set.seed(123)
  data_reduced <- DaMiR.FSelect(t(assay(TR_set)),
                                as.data.frame(colData(TR_set)),
                                th.corr=0.4)
  set.seed(123)
  data_reduced <- DaMiR.FReduct(data_reduced$data,th.corr = 0.9)
  set.seed(123)
  df_importance <- DaMiR.FSort(data_reduced,
                               as.data.frame(colData(TR_set)))
  set.seed(123)
  selected_features <- DaMiR.FBest(data_reduced,
                                   ranking=df_importance,
                                   autoselect = "yes")
  # update datasets
  TR_set <- TR_set[selected_features$predictors,, drop=FALSE]
  Val_set <- Val_set[selected_features$predictors,drop=FALSE]
  ### Model building
  ensl_model <- DaMiR.EnsL_Train(TR_set,fSample.tr.w = 0.5,
                                 cl_type = c("SVM", "LDA", "LR", "NB", "NN", "PLS"))#,cl_type = c("RF","LR"))
  # Store all trained models
  cv_models[[cv_fold]] <- ensl_model
  ### Model testing
  res_Val <- DaMiR.EnsL_Test(Val_set,
                             EnsL_model = ensl_model)
  # Store all ML results
  res_df[cv_fold,1] <- res_Val$accuracy[1] # Accuracy
  res_df[cv_fold,2] <- length(res_Val$predictors) # N. of predictors
  res_df[cv_fold,3] <- res_Val$MCC[1]
  res_df[cv_fold,4] <- res_Val$sensitivity[1]
  res_df[cv_fold,5] <- res_Val$Specificty[1]
  res_df[cv_fold,6] <- res_Val$PPV[1]
  res_df[cv_fold,7] <- res_Val$NPV[1]
  cv_predictors[[cv_fold]] <- res_Val$predictors
  
  print(paste("cv_fold TERMINATA numero:",cv_fold, sep=''))
}


# Model Selection
res_df[,1:5]

idx_best_model <- DaMiR.ModelSelect(res_df,
                                    type.sel = c("median"),
                                    npred.sel = c("min"))

# Prediction on the the independent test set
res_predict <- DaMiR.EnsL_Predict(Test_set,
                                  bestModel = cv_models[[idx_best_model]])
# Predictors
cv_predictors[[idx_best_model]]

# Prediction assessment for Ensemble learning
id_classifier <- 1 # Prendo la colonna dell' Ensemble Learning
table(colData(Test_set)$class, res_predict[,id_classifier])

# Prediction assessment for Logistic regression
id_classifier <- 2 # Logistic regression
table(colData(Test_set)$class, res_predict[,id_classifier])

#   NR R
#NR  6 1
#R   1 2



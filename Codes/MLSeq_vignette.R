## ----knitr_options, echo=FALSE, results="hide", warning=FALSE-----------------
library(knitr)
opts_chunk$set(tidy = FALSE, dev = "pdf", fig.show = "hide", message = FALSE, fig.align = "center", cache = FALSE)

## ----load_packages, echo=FALSE, results="hide", warning=FALSE-----------------
library(MLSeq)
library(DESeq2)
library(edgeR)
library(VennDiagram)
library(pamr)
library(caret)

## ----eval = FALSE-------------------------------------------------------------
#  library(MLSeq)

## ----file_path_cervical-------------------------------------------------------
filepath <- system.file("extdata/cervical.txt", package = "MLSeq")

## ----read_cervical_data-------------------------------------------------------
cervical <- read.table(filepath, header=TRUE)

## ----head_cervical------------------------------------------------------------
head(cervical[ ,1:10]) # Mapped counts for first 6 features of 10 subjects.

## ----define_class_labels------------------------------------------------------
class <- DataFrame(condition = factor(rep(c("N","T"), c(29, 29))))
class

## ----data_splitting-----------------------------------------------------------
library(DESeq2)

set.seed(2128)

# We do not perform a differential expression analysis to select differentially
# expressed genes. However, in practice, DE analysis might be performed before
# fitting classifiers. Here, we selected top 100 features having the highest
# gene-wise variances in order to decrease computational cost.
vars <- sort(apply(cervical, 1, var, na.rm = TRUE), decreasing = TRUE)
data <- cervical[names(vars)[1:100], ]
nTest <- ceiling(ncol(data) * 0.3)
ind <- sample(ncol(data), nTest, FALSE)

# Minimum count is set to 1 in order to prevent 0 division problem within
# classification models.
data.train <- as.matrix(data[ ,-ind] + 1)
data.test <- as.matrix(data[ ,ind] + 1)
classtr <- DataFrame(condition = class[-ind, ])
classts <- DataFrame(condition = class[ind, ])

## ----DESeqDataSets------------------------------------------------------------
data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,
                                      design = formula(~condition))
data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,
                                     design = formula(~condition))

## ----mwe_limitations_on_continuous_classifiers, eval = FALSE, message=FALSE----
#  # Support Vector Machines with Radial Kernel
#  fit <- classify(data = data.trainS4, method = "svmRadial",
#                   preProcessing = "deseq-rlog", ref = "T",
#                   control = trainControl(method = "repeatedcv", number = 2,
#                                          repeats = 2, classProbs = TRUE))
#  show(fit)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  set.seed(2128)
#  
#  # Voom based Nearest Shrunken Centroids.
#  fit <- classify(data = data.trainS4, method = "voomNSC",
#                   normalize = "deseq", ref = "T",
#                   control = voomControl(tuneLength = 20))
#  
#  trained(fit)  ## Trained model summary

## ----Optimizing_model_parameters_example, eval = TRUE, echo = TRUE------------
set.seed(2128)

# Support vector machines with radial basis function kernel
fit.svm <- classify(data = data.trainS4, method = "svmRadial",
                    preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                    control = trainControl(method = "repeatedcv", number = 5,
                                           repeats = 10, classProbs = TRUE))

show(fit.svm)

## ----fitted_model_svm---------------------------------------------------------
trained(fit.svm)

## ----eval = FALSE-------------------------------------------------------------
plot(fit.svm)

## ----fitted_model_svm_figure, echo = FALSE, results='hide'--------------------
cairo_pdf(filename = "fitted_model_svm_figure.pdf", height = 5.5)
plot(fit.svm)
dev.off()

## ----control_svm_model_example, eval = FALSE----------------------------------
#  # Define control list
#  ctrl.svm <- trainControl(method = "repeatedcv", number = 5, repeats = 1)
#  ctrl.plda <- discreteControl(method = "repeatedcv", number = 5, repeats = 1,
#                               tuneLength = 10)
#  ctrl.voomDLDA <- voomControl(method = "repeatedcv", number = 5, repeats = 1,
#                               tuneLength = 10)
#  
#  # Support vector machines with radial basis function kernel
#  fit.svm <- classify(data = data.trainS4, method = "svmRadial",
#                   preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
#                   control = ctrl.svm)
#  
#  # Poisson linear discriminant analysis
#  fit.plda <- classify(data = data.trainS4, method = "PLDA", normalize = "deseq",
#                       ref = "T", control = ctrl.plda)
#  
#  # Voom-based diagonal linear discriminant analysis
#  fit.voomDLDA <- classify(data = data.trainS4, method = "voomDLDA",
#                           normalize = "deseq", ref = "T", control = ctrl.voomDLDA)

## ----echo = FALSE-------------------------------------------------------------
# Define control list
ctrl.voomDLDA <- voomControl(method = "repeatedcv", number = 5, repeats = 1,
                             tuneLength = 10)

# Voom-based diagonal linear discriminant analysis
fit.voomDLDA <- classify(data = data.trainS4, method = "voomDLDA",
                         normalize = "deseq", ref = "T", control = ctrl.voomDLDA)

## -----------------------------------------------------------------------------
trained(fit.voomDLDA)

## -----------------------------------------------------------------------------
#Predicted class labels
pred.svm <- predict(fit.svm, data.testS4)
pred.svm

## -----------------------------------------------------------------------------
pred.svm <- relevel(pred.svm, ref = "T")
actual <- relevel(classts$condition, ref = "T")

tbl <- table(Predicted = pred.svm, Actual = actual)
confusionMatrix(tbl, positive = "T")

## ----results='hide', message=FALSE--------------------------------------------
set.seed(2128)

# Define control lists.
ctrl.continuous <- trainControl(method = "repeatedcv", number = 5, repeats = 10)
ctrl.discrete <- discreteControl(method = "repeatedcv", number = 5, repeats = 10,
                                 tuneLength = 10)
ctrl.voom <- voomControl(method = "repeatedcv", number = 5, repeats = 10,
                         tuneLength = 10)

# 1. Continuous classifiers, SVM and NSC
fit.svm <- classify(data = data.trainS4, method = "svmRadial",
                    preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                    control = ctrl.continuous)

fit.NSC <- classify(data = data.trainS4, method = "pam",
                    preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
                    control = ctrl.continuous)

# 2. Discrete classifiers
fit.plda <- classify(data = data.trainS4, method = "PLDA", normalize = "deseq",
                     ref = "T", control = ctrl.discrete)

fit.plda2 <- classify(data = data.trainS4, method = "PLDA2", normalize = "deseq",
                      ref = "T", control = ctrl.discrete)

fit.nblda <- classify(data = data.trainS4, method = "NBLDA", normalize = "deseq",
                      ref = "T", control = ctrl.discrete)

# 3. voom-based classifiers
fit.voomDLDA <- classify(data = data.trainS4, method = "voomDLDA",
                         normalize = "deseq", ref = "T", control = ctrl.voom)

fit.voomNSC <- classify(data = data.trainS4, method = "voomNSC",
                        normalize = "deseq", ref = "T", control = ctrl.voom)

# 4. Predictions
pred.svm <- predict(fit.svm, data.testS4)
pred.NSC <- predict(fit.NSC, data.testS4)
# ... truncated

## ----echo = FALSE, results='asis', message=FALSE------------------------------
library(xtable)

pred.svm <- predict(fit.svm, data.testS4)
pred.NSC <- predict(fit.NSC, data.testS4)
pred.plda <- predict(fit.plda, data.testS4)
pred.nblda <- predict(fit.nblda, data.testS4)
pred.voomDLDA <- predict(fit.voomDLDA, data.testS4)
pred.voomNSC <- predict(fit.voomNSC, data.testS4)

actual <- data.testS4$condition
nn <- length(actual)
diag.svm <- sum(diag(table(pred.svm, actual)))
diag.NSC <- sum(diag(table(pred.NSC, actual)))
diag.plda <- sum(diag(table(pred.plda, actual)))
diag.nblda <- sum(diag(table(pred.nblda, actual)))
diag.voomDLDA <- sum(diag(table(pred.voomDLDA, actual)))
diag.voomNSC <- sum(diag(table(pred.voomNSC, actual)))

acc <- c(diag.svm, diag.NSC, diag.plda, diag.nblda, diag.voomDLDA, diag.voomNSC) / nn
sparsity <- c(NA, trained(fit.NSC)$finalModel$nonzero/nrow(data.testS4),
              length(selectedGenes(fit.plda))/nrow(data.testS4), NA, NA,
              length(selectedGenes(fit.voomNSC))/nrow(data.testS4))

tbl <- data.frame(Classifier = c("SVM", "NSC", "PLDA (Transformed)", "NBLDA", "voomDLDA", "voomNSC"), Accuracy = acc, Sparsity = sparsity)

xtbl <- xtable(tbl, caption = "Classification results for cervical data.", label = "tbl:accRes", align = "lp{4cm}p{2cm}c")

digits(xtbl) <- c(0, 0, 3, 3)
print.xtable(xtbl, caption.placement = "top", include.rownames = FALSE, booktabs = TRUE)

## ----echo = FALSE-------------------------------------------------------------
best_in_accuracy <- as.character(tbl$Classifier[which(acc == max(acc, na.rm = TRUE))])
best_in_acc_text <- paste("\\textbf{", best_in_accuracy, "}", sep = "")

if (length(best_in_accuracy) >= 2){
  best_in_acc_text <- paste(paste(best_in_acc_text[-length(best_in_acc_text)], collapse = ", "), best_in_acc_text[length(best_in_acc_text)], sep = " and ")
}

best_in_sparsity <- as.character(tbl$Classifier[which(sparsity == min(sparsity, na.rm = TRUE))])
best_in_sparsity_text <- paste("\\textbf{", best_in_sparsity, "}", sep = "")

if (length(best_in_sparsity) >= 2){
  best_in_sparsity_text <- paste(paste(best_in_sparsity_text[-length(best_in_sparsity_text)], collapse = ", "), best_in_sparsity_text[length(best_in_sparsity_text)], sep = " and ")
}

## -----------------------------------------------------------------------------
selectedGenes(fit.voomNSC)

## ----all_common_features, echo = FALSE----------------------------------------
pam.final <- trained(fit.NSC)$finalModel   ## 'pamrtrained' object.
geneIdx <- pamr:::pamr.predict(pam.final, pam.final$xData, threshold = pam.final$threshold, type = "nonzero")

genes.pam <- colnames(pam.final$xData)[geneIdx]
genes.plda <- selectedGenes(fit.plda)
genes.plda2 <- selectedGenes(fit.plda2)
genes.vnsc <- selectedGenes(fit.voomNSC)

tmp.list <- list(genes.pam, genes.plda, genes.plda2, genes.vnsc)

nn <- c(length(genes.pam), length(genes.plda), length(genes.plda2), length(genes.vnsc))
ooo <- order(nn, decreasing = TRUE)

tmp.list <- tmp.list[ooo]

common <- tmp.list[[1]]
for (i in 2:(length(tmp.list))){
  tmp2 <- tmp.list[[i]]
  tmp <- common[common %in% tmp2]
  common <- tmp
}

## ----venn_diagram, echo = FALSE-----------------------------------------------
venn.plot <- venn.diagram(
  x = list(voomNSC = genes.vnsc, NSC = genes.pam, PLDA = genes.plda, PLDA2 = genes.plda2),
  height = 1200, width = 1200,
  resolution = 200,
  filename = "Selected_features.png", imagetype = "png",
  col = "black",
  fill = c("khaki1", "skyblue", "tomato3", "darkolivegreen3"),
  alpha = 0.50,
  cat.cex = 1.2,
  cex = 1.5,
  cat.fontface = "bold"
)

## -----------------------------------------------------------------------------
set.seed(2128)

ctrl <- discreteControl(method = "repeatedcv", number = 5, repeats = 2,
                        tuneLength = 10)

# PLDA without power transformation
fit <- classify(data = data.trainS4, method = "PLDA", normalize = "deseq",
                ref = "T", control = ctrl)
show(fit)

## -----------------------------------------------------------------------------
method(fit) <- "PLDA2"
show(fit)

## -----------------------------------------------------------------------------
ref(fit) <- "N"
normalization(fit) <- "TMM"
metaData(fit)

## -----------------------------------------------------------------------------
fit <- update(fit)
show(fit)

## ----echo = TRUE, message=FALSE, error=FALSE, eval = FALSE--------------------
#  method(fit) <- "rpart"
#  update(fit)

## ----echo = FALSE, message=FALSE, error=TRUE----------------------------------
method(fit) <- "rpart"
tmp <- try(update(fit))

## -----------------------------------------------------------------------------
control(fit) <- trainControl(method = "repeatedcv", number = 5, repeats = 2)

# 'normalize' is not valid for continuous classifiers. We use 'preProcessing'
# rather than 'normalize'.
preProcessing(fit) <- "tmm-logcpm"

fit <- update(fit)
show(fit)

## ----session_info-------------------------------------------------------------
sessionInfo()

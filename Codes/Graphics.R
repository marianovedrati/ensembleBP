setwd("/Users/giorgiomontesi/Desktop/Universita_di_Siena/A_PhD_Project/Biomarker_Prediction/ensembleBP/Codes")


library(ggplot2)



df1 <- read.csv2(paste0("../Results/AccuracyTable_1.csv"))
df2 <- read.csv2(paste0("../Results/AccuracyTable_2.csv"))
df3 <- read.csv2(paste0("../Results/AccuracyTable_3.csv"))
df4 <- read.csv2(paste0("../Results/AccuracyTable_4.csv"))
df5 <- read.csv2(paste0("../Results/AccuracyTable_5.csv"))

df1$Group <- c(rep("SVM-based",each=3), rep("voom-based",each=3), 
               rep("LDA-based",each=4), rep("NNet-based",each=2),
               rep("Tree-based",each=4), rep("Bagged",each=3),
               rep("Boosted",each=3), rep("PLS-based",each=3))

df2$Group <- c(rep("SVM-based",each=3), rep("voom-based",each=3), 
               rep("LDA-based",each=4), rep("NNet-based",each=2),
               rep("Tree-based",each=4), rep("Bagged",each=3),
               rep("Boosted",each=3), rep("PLS-based",each=3))

df3$Group <- c(rep("SVM-based",each=3), rep("voom-based",each=3), 
               rep("LDA-based",each=4), rep("NNet-based",each=2),
               rep("Tree-based",each=4), rep("Bagged",each=3),
               rep("Boosted",each=3), rep("PLS-based",each=3))

df4$Group <- c(rep("SVM-based",each=3), rep("voom-based",each=3), 
               rep("LDA-based",each=4), rep("NNet-based",each=2),
               rep("Tree-based",each=4), rep("Bagged",each=3),
               rep("Boosted",each=3), rep("PLS-based",each=3))

df5$Group <- c(rep("SVM-based",each=3), rep("voom-based",each=3), 
               rep("LDA-based",each=4), rep("NNet-based",each=2),
               rep("Tree-based",each=4), rep("Bagged",each=3),
               rep("Boosted",each=3), rep("PLS-based",each=3))




df <- do.call("rbind", list(df1, df2, df3, df4, df5))

ggplot(df, aes(x = as.factor(X), y = Accuracy, fill = as.factor(Group))) +
  geom_boxplot() +
  geom_point() +
  scale_fill_brewer(palette = c("Paired")) +
  theme_minimal()














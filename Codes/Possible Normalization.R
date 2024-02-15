library("DESeq2")
#For training dataset
ddsTrain <- DESeqDataSetFromMatrix(countData = trainData, colData = colDataTrain, design = ~condition)
ddsTrain <- DESeq(ddsTrain)
vstNormalizedExpressionDataForTrain <- varianceStabilizingTransformation(ddsTrain, blind = FALSE)

#For testing dataset
ddsTest <- DESeqDataSetFromMatrix(countData = testData, colData = colDataTest, design = ~condition)
ddsTest <- DESeq(ddsTest)
dispersionFunction(ddsTest) <- dispersionFunction(ddsTrain)
vstNormalizedExpressionDataForTest <- varianceStabilizingTransformation(ddsTest, blind = FALSE)
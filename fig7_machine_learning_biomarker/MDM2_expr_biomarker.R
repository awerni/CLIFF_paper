library(tidyverse)
library(extrafont)
library(XIFF)
library(CLIFF)

s <- getSettings()
print(s$db)
setDbOptions(s)

source("MDM2_analysis_functions.R")

cl_anno <- getCelllineAnno_4paper()
gene_anno <- getGeneAnno_4paper()

gene_anno_MDM2 <- gene_anno |> filter(symbol == "MDM2")

# ---------------- MDM2 CRISPR dependency data ----------------
MDM2_sens <- getMDM2SensitivityClassification(gene_anno_MDM2, cl_anno)
MDM2_avana <- MDM2_sens$MDM2_avana
ca_sens <- MDM2_sens$ca_sens

# --------- Machine Learning -----------
#sigSymbols <- c("MDM2", "CDKN1A", "ZMAT3", "DDB2", "FDXR", "RPS27L", "BAX", "RRM2B", "SESN1", "CCNG1", "XPC", "TNFRSF10B", "AEN")

sigSymbols <- c("MDM2", "ZMAT3", "DDB2", "CDKN1A", "RPS27L")
sig_genes <- gene_anno |> filter(symbol %in% sigSymbols)

sets <- XIFF::splitTrainingTestSets(ca_sens, 0.2)
trainingSet <- sets$training
testSet <- sets$test

defaultFit <- XIFF::buildMachineLearning(trainingSet, sig_genes$ensg, gene_anno, method = "rf", p_test = 0)

modelTestData <- getDataForModel(assignment = testSet, features = defaultFit$bestFeatures)


modelTestData |> select(celllinename, class) |>
  mutate(predicted = predict(defaultFit, newdata = modelTestData)) |>
  mutate(correct = ifelse(class == predicted, "correct", "incorrect")) |>
  View()
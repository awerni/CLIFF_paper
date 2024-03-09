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

# ----- Machine Learning 2 (finding MOA by hallmark prediction)------

sets_MDM2 <- XIFF::splitTrainingTestSets(ca_sens, 0.2)

trainingSet <- sets_MDM2$training
testSet <- sets_MDM2$test

msig_data <- getGSEAdata("human", "hallmark")

# ---- create models for every hallmark gene set
MDM2_model <- lapply(msig_data, function(m) {
  XIFF::buildMachineLearning(trainingSet, m, gene_anno, method = "glmnet", p_test = 0, maxFeatures = 25)
})

#save(ca_sens, sets_MDM2, MDM2_model, file = "MDM2_ML_hallmarks.Rdata")
load("MDM2_ML_gmlnet_max25_hallmarks.Rdata")

# ---- apply model ----
modelTestData <- getDataForModel(
  assignment = testSet,
  features = unlist(msig_data)
)

ml_result_MDM2 <- lapply(MDM2_model, function(m) {
  r <- modelTestData %>% select(celllinename, class) %>%
    mutate(class = factor(class, levels = c("sensitive", "resistant")),
           predicted = predict(m, newdata = modelTestData))
  XIFF::generateTestPerformanceData(table(r$predicted, r$class))
})

# ---- make statistics ------
bind_rows(ml_result_MDM2, .id = "hallmark") %>% View()
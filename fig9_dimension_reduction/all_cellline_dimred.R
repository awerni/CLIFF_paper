library(tidyverse)
library(extrafont)
library(XIFF)
library(CLIFF)

s <- getSettings()
print(s$db)
setDbOptions(s)

cl_anno <- getCellLineAnno("human")

cl_anno <- cl_anno |>
  mutate(tumortype_coarse = case_when(
    grepl("leukemia", tumortype) ~ "leukemia",
    grepl("lymphoma", tumortype) ~ "lymphoma",
    TRUE ~ tumortype
  ))

gene_expr_data <- getExpressionDimRedData(sampleClasses = cl_anno$celllinename, p = TRUE)

dimRedData <- getExpressionDimRed(gene_expr_data, cl_anno, method = "umap", numGenes = 300)
generateExpressionDimRedPlot(dimRedData, colorCol = "tumortype")

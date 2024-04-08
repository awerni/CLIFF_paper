library(tidyverse)
library(extrafont)
library(XIFF)
library(CLIFF)

s <- getSettings()
print(s$db)
setDbOptions(s)

cl_anno <- getCellLineAnno("human") |>
  filter(histology_type == "lung adenocarcinoma")

gene_anno <- getGeneAnno("human")

gene_expr_data <- getExpressionDimRedData(sampleClasses = cl_anno$celllinename, p = TRUE)

dimRedData <- getExpressionDimRed(gene_expr_data, cl_anno, method = "phate", numGenes = 300)
generateExpressionDimRedPlot(dimRedData, colorCol = "cluster", fontSize = 18)

lung_adeno_class <- unstack(dimRedData, celllinename ~ cluster) |>
  classAssignment()

diff_expr <- differentialGeneExpression_LimmaVoom(
  lung_adeno_class,
  gene_anno,
  p = TRUE
)

gseaData <- getGSEAdata("human", "hallmark")
gseaResult <- getGSEA(diffExResult = diff_expr, geneSets = gseaData, rankType = "p.valueDir")

gseaResult |> head()

top_geneset <- gseaResult[[1, "geneset"]]

p <- generateGSEA_plot(diff_expr, gseaData[[top_geneset]], "p.valueDir")


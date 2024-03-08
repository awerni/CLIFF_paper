library(tidyverse)
library(extrafont)
library(XIFF)
library(CLIFF)

s <- getSettings()
print(s$db)
setDbOptions(s)

cl_anno <- getCellLineAnno("human") |>
  filter(grepl("(lung|head|esopha)", tumortype))

gene_anno <- getGeneAnno("human") |>
  filter(grepl("^.{1,2}:", location))

EGFR_anno <- gene_anno |> filter(symbol == "EGFR")

expr_EGFR <- getCelllineDataGeneExpressionById(EGFR_anno$ensg, cl_anno$celllinename) |>
  mutate(log2tpm = log2(tpm)) |>
  left_join(cl_anno, by = "celllinename") |>
  reorderByScore(valueCol = "log2tpm")

generateExpressionWaterfallPlot(expr_EGFR)
ggsave("fig3_EGFR_expr_waterfall_4TT.pdf", width = 15, height = 5)

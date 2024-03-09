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
gene_anno_PPM1D <- gene_anno |> filter(symbol == "PPM1D")
gene_anno_CDKN1A <- gene_anno |> filter(symbol == "CDKN1A")
gene_anno_TP53  <- gene_anno |> filter(symbol == "TP53")

# ---------------- MDM2 CRISPR dependency data ----------------
MDM2_sens <- getMDM2SensitivityClassification(gene_anno_MDM2, cl_anno)
MDM2_avana <- MDM2_sens$MDM2_avana
ca_sens <- MDM2_sens$ca_sens

# ------------ differential Avana Dependency -----------
diff_Avana <- CLIFF::getDepletionAssociation(ca_sens, "Avana", "chronos", gene_anno, p = TRUE)
head(diff_Avana, n = 3)

getGeneTitle <- function(a) paste(a$symbol, "-", a$ensg)

TP53_Avana <- CLIFF::getCelllineDataDepletionById(gene_anno_TP53$ensg, ca_sens, "Avana", "chronos")
generateDepletionPlot(TP53_Avana, plotType = "point", dataCol = "chronos", title = getGeneTitle(gene_anno_TP53))
ggsave("fig6A_MDM2_sens_TP53_Avana.pdf", width = 6, height = 6)

CDKN1A_Avana <- CLIFF::getCelllineDataDepletionById(gene_anno_CDKN1A$ensg, ca_sens, "Avana", "chronos")
generateDepletionPlot(CDKN1A_Avana, plotType = "point", dataCol = "chronos", title = getGeneTitle(gene_anno_CDKN1A))
ggsave("fig6B_MDM2_sens_CDKN1A_Avana.pdf", width = 6, height = 6)

PPM1D_Avana <- CLIFF::getCelllineDataDepletionById(gene_anno_PPM1D$ensg, ca_sens, "Avana", "chronos")
generateDepletionPlot(PPM1D_Avana, plotType = "point", dataCol = "chronos", title = getGeneTitle(gene_anno_PPM1D))
ggsave("fig6C_MDM2_sens_PPM1D_Avana.pdf", width = 6, height = 6)

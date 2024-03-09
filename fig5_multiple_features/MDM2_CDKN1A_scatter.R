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
gene_anno_CDKN1A <- gene_anno |> filter(symbol == "CDKN1A")
gene_anno_TP53 <- gene_anno |> filter(symbol == "TP53")

# ---------------- MDM2 CRISPR dependency data ----------------
MDM2_sens <- getMDM2SensitivityClassification(gene_anno_MDM2, cl_anno)

MDM2_avana <- MDM2_sens$MDM2_avana |>
  rename(MDM2_chronos = chronos) |>
  select(-ensg)

# ---------------- CDKN1A gene expression data ----------------
CDKN1A_expr <- CLIFF::getCelllineDataGeneExpressionById(gene_anno_CDKN1A$ensg, cl_anno$celllinename) |>
  rename(CDKN1A_tpm = tpm) |>
  select(-ensg) 

# ---------------- TP53 mutation data ----------------
TP53_mut <- CLIFF::getCelllineDataMutationById(gene_anno_TP53$ensg, cl_anno$celllinename) |>
  rename(TP53_status = aamutated) |>
  select(-ensg)

# ---------------- data assembly ----------------
data_all <- cl_anno |>
  inner_join(MDM2_avana, by = "celllinename") |>
  inner_join(CDKN1A_expr, by = "celllinename") |>
  inner_join(TP53_mut, by = "celllinename")

# ---------------- plot data ----------------
ggplot(data_all, aes(x = CDKN1A_tpm, y = MDM2_chronos, color = TP53_status)) +
  geom_point() +
  scale_x_log10() +
  theme_light() +
  theme(text = element_text(family = "Roboto", size = 18), legend.position = "bottom") +
  xlab("CDKN1A gene expression [TPM]") +
  ylab("MDM2 CRISPR dependency [chronos]")

ggsave("fig5_CDKN1A_expr_vs_MDM2_chronos.pdf", width = 10, height = 6)
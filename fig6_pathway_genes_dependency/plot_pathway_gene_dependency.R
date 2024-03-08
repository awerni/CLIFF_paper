library(tidyverse)
library(extrafont)
library(XIFF)
library(CLIFF)

s <- getSettings()
print(s$db)
setDbOptions(s)

cl_anno <- getCellLineAnno("human") |>
  filter(tumortype %in% c("ovarian cancer", "melanoma", "renal cell carcinoma", "colorectal cancer", "esophagogastric cancer")) |>
  filter(!histology_type %in% c("esophageal adenocarcinoma", "esophageal squamous cell carcinoma")) |>
  mutate(tumortype = ifelse(tumortype == "esophagogastric cancer", "gastric cancer", as.character(tumortype))) |>
  mutate(tumortype = as.factor(tumortype))

gene_anno <- getGeneAnno("human") |>
  filter(grepl("^.{1,2}:", location))

gene_anno_MDM2 <- gene_anno |> filter(symbol == "MDM2")
gene_anno_PPM1D <- gene_anno |> filter(symbol == "PPM1D")
gene_anno_CDKN1A <- gene_anno |> filter(symbol == "CDKN1A")

# ---------------- MDM2 CRISPR dependency data ----------------
MDM2_avana <- getCelllineDataDepletionById(gene_anno_MDM2$ensg, celllineClasses = cl_anno$celllinename, study = "Avana", scores = "chronos") |>
  mutate(class = case_when(
    chronos < -1  ~ "sensitive",
    chronos > -0.5 ~ "resistant",
    TRUE  ~ "intermediate"
  )) |>
  left_join(cl_anno, by = "celllinename") |>
  reorderByScore(valueCol = "chronos")

# ---------------- MDM2 CRISPR dependency classification ----------------
ca_sens <- split(as.character(MDM2_avana$celllinename), MDM2_avana$class)[c("sensitive", "resistant")] |>
  classAssignment()

# ------------ differential Avana Dependency -----------

diff_Avana <- CLIFF::getDepletionAssociation(ca_sens, "Avana", "chronos", gene_anno, p = TRUE)

getGeneTitle <- function(a) paste(a$symbol, "-", a$ensg)

PPM1D_Avana <- CLIFF::getCelllineDataDepletionById(gene_anno_PPM1D$ensg, ca_sens, "Avana", "chronos")
generateDepletionPlot(PPM1D_Avana, plotType = "point", dataCol = "chronos", title = getGeneTitle(gene_anno_PPM1D))
ggsave("MDM2_sens_PPM1D_Avana.pdf", width = 6, height = 6)

TP53_Avana <- CLIFF::getCelllineDataDepletionById(gene_anno_TP53$ensg, ca_sens, "Avana", "chronos")
generateDepletionPlot(TP53_Avana, plotType = "point", dataCol = "chronos", title = getGeneTitle(gene_anno_TP53))
ggsave("MDM2_sens_TP53_Avana.pdf", width = 6, height = 6)

CDKN1A_Avana <- CLIFF::getCelllineDataDepletionById(gene_anno_CDKN1A$ensg, ca_sens, "Avana", "chronos")
generateDepletionPlot(CDKN1A_Avana, plotType = "point", dataCol = "chronos", title = getGeneTitle(gene_anno_CDKN1A))
ggsave("MDM2_sens_CDKN1A_Avana.pdf", width = 6, height = 6)

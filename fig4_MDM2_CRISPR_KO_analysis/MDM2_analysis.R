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

# ------------ MDM2 CRISPR dependency Waterfall plot with tumortype ------
generateWaterfallPlot(MDM2_avana, "chronos", ylabel = "chronos") + geom_hline(yintercept = -1, linetype = "dashed") +
  theme(text = element_text(family = "Roboto", size = 20))
ggsave("fig4X_MDM2_sens.pdf", width = 14, height = 5)

# ------------ MDM2 CRISPR dependency Waterfall plot with TP53 status ------
gene_anno_TP53 <- gene_anno |> filter(symbol == "TP53")

TP53_mut <- getCelllineDataMutationById(gene_anno_TP53$ensg, MDM2_avana$celllinename)

MDM2_avana |> 
  left_join(TP53_mut, by = "celllinename") |>
  rename(TP53_status = aamutated) |>
  XIFF::reorderByScore(valueCol = "chronos") |>
  generateWaterfallPlot("chronos", ylabel = "chronos", fill = "TP53_status") + 
  geom_hline(yintercept = -1, linetype = "dashed") +
  theme(text = element_text(family = "Roboto", size = 20))

ggsave("fig4A_MDM2_sens_TP53status.pdf", width = 14, height = 5)

# ---------------- mutation plot ----------------------------------
diff_mut <- CLIFF::getMutationAssociation(ca_sens, gene_anno, p = TRUE)
top_mut <- getCelllineDataMutationById(diff_mut[[1, "ensg"]], ca_sens)
gene_anno |> filter(ensg %in% unique(top_mut$ensg))

generateMutationPlot(top_mut, "bar")

generateMutationPlot(top_mut, "pie") + theme(text = element_text(size = 24))
ggsave("fig4I_MDM2_sens_TP53.pdf", width = 10, height = 5)

generateMutationPlot(top_mut, "coverage")

# ---------------- differential expression plot ----------------------------------
diff_expr <- differentialGeneExpression_LimmaVoom(ca_sens, gene_anno, p = TRUE)
XIFF::generateVolcanoPlot(diff_expr, minuslog10pval = 10, minuslog10adjpval = 1, 
                          log2FC = 1, classLabels = c("sensitive", "resistant"))
ggsave("Fig4E_MDM2_volcano.pdf", width = 15, height = 8)

# ----------- tumor type plot ----------------------------------
XIFF::generateClassSelectionPlot(
  sampleClasses = ca_sens,
  prop1 = "tumortype",
  prop2 = "none",
  n_classes = 20,
  annotationFocus = cl_anno
) + theme(text = element_text(size = 24))

ggsave("fig4B_MDM2_tumortype.pdf", width = 16, height = 3)

# ---------------- CDKN1A gene expression ----------

gene_anno_CDKN1A <- gene_anno |> filter(symbol == "CDKN1A")
CDKN1A_expr <- CLIFF::getCelllineDataGeneExpressionById(gene_anno_CDKN1A$ensg, ca_sens) 

generateExpressionPlot(CDKN1A_expr, plotType = "point", title = paste("CDKN1A -", gene_anno_CDKN1A$ensge))
ggsave("fig4C_MDM2_CDKN1A_expr_point.pdf", width = 6, height = 6)

generateExpressionPlot(CDKN1A_expr, plotType = "roc", title = paste("CDKN1A -", gene_anno_CDKN1A$ensg))
ggsave("fig4D_MDM2_CDKN1A_expr_ROC.pdf", width = 6, height = 6)

# ---------------- CDKN1A protein expression ----------
antibodyAnno <- getAntibodyInformation()
#massSpecExpr <- differentialMassSpecProteinExpression(ca_sens)

prot_title <- "CDN1A_HUMAN - P38936-0"
prot_id <- gsub(" ", "", prot_title)

CDKN1A_protexpr <- getCelllineDataMassSpecById(prot_id, ca_sens)

generateMassSpecProteinExpressionPlot(df = CDKN1A_protexpr, plotType = "point", title = prot_title)
ggsave("fig4F_MDM2_Protein_CDKN1A_Point.pdf", width = 6, height = 6)

generateMassSpecProteinExpressionPlot(df = CDKN1A_protexpr, plotType = "roc", title = prot_title)
ggsave("fig4G_MDM2_Protein_CDKN1A_AUC.pdf", width = 6, height = 6)

generateMassSpecProteinExpressionPlot(df = CDKN1A_protexpr, plotType = "coverage", title = prot_title)
ggsave("fig4H_MDM2_Protein_CDKN1A_coverage.pdf", width = 6, height = 6)

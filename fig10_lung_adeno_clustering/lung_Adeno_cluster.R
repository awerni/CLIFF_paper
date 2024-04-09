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

lung_adeno_class <- dimRedData |>
  mutate(cluster_label = case_when(
    cluster == "1" ~ "cluster 1 (mesenchymal)",
    cluster == "2" ~ "cluster 2 (epithelial)"
  )) |>
  unstack(celllinename ~ cluster_label) |>
  classAssignment()

gene_expr_data <- getExpressionDimRedData(sampleClasses = lung_adeno_class, p = TRUE)
dimRedData <- getExpressionDimRed(gene_expr_data, cl_anno, method = "phate", numGenes = 300)
generateExpressionDimRedPlot(dimRedData, colorCol = "class", fontSize = 18) + 
  ggtitle(paste0("PHATE plot\n", attr(dimRedData, "title")))

ggsave("fig10A_LUAD_dimRed_plot.pdf", width = 7, height = 7)

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
ggsave("fig10B_EMT_GSEA_plot.pdf", plot = p, width = 12, height = 5)

# ---------------- gene expression ----------------
lung_adeno_short_class <- lung_adeno_class
setClassLabel(lung_adeno_short_class, "class1") <- "mesenchymal"
setClassLabel(lung_adeno_short_class, "class2") <- "epithelial"

gene_anno_VIM <- gene_anno |> filter(symbol == "VIM")
VIM_expr <- CLIFF::getCelllineDataGeneExpressionById(gene_anno_VIM$ensg, lung_adeno_short_class) 
generateExpressionPlot(VIM_expr, plotType = "point", title = paste("VIM -", gene_anno_VIM$ensg))
#ggsave("fig10C_VIM_gene_expr_point.pdf", width = 6, height = 6)

gene_anno_CDH1 <- gene_anno |> filter(symbol == "CDH1" & grepl("ENSG", ensg))
CDH1_expr <- CLIFF::getCelllineDataGeneExpressionById(gene_anno_CDH1$ensg, lung_adeno_short_class) 
generateExpressionPlot(CDH1_expr, plotType = "point", title = paste("CDH1 -", gene_anno_VIM$ensg))
#ggsave("fig10C_CDH1_gene_expr_point.pdf", width = 6, height = 6)

# ----------------- protein expression -----------------
VIM_prot_expr <- getMassSpecProteinExpressionByGeneData(gene_anno_VIM$ensg) |>
  mutate(id = paste0(uniprotid, "-", accession, "-", isoform))

generateMassSpecProteinExpressionPlot(VIM_prot_expr, "point", gene_anno_VIM$name, ca = lung_adeno_short_class)
ggsave("fig10C_VIM_prot_expr_point.pdf", width = 6, height = 6)

CDH1_prot_expr <- getMassSpecProteinExpressionByGeneData(gene_anno_CDH1$ensg) |>
  mutate(id = paste0(uniprotid, "-", accession, "-", isoform))

generateMassSpecProteinExpressionPlot(CDH1_prot_expr, "point", gene_anno_CDH1$name, ca = lung_adeno_short_class)
ggsave("fig10C_CDH1_prot_expr_point.pdf", width = 6, height = 6)

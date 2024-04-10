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

all_hallmarks <- CLIFF::getAvailableHallmarkGeneSets()

performance2 <- lapply(all_hallmarks, function(h) {
  print(h)
  geneSet <- XIFF::getGSEAdata("human", "hallmark", h)

  defaultFit <- XIFF::buildMachineLearning(
    cs = ca_sens,
    geneSet = geneSet,
    geneAnno = gene_anno,
    method = "svmLinear2",
    p_test = 0,
    maxFeatures = 50
  )

  defaultFit$resample$Accuracy #|> summary()
})
save(MDM2_avana, ca_sens, performance2, all_hallmarks, file = "CLIFF_ML_hallmarks_SVM.Rdata")

performance_stat2 <- tibble(hallmark = all_hallmarks, accuracy = performance2) |> 
  unnest(cols = c(accuracy)) |>
  arrange(accuracy)

ggplot(performance_stat2, aes(x = accuracy, y = forcats::fct_inorder(hallmark))) + geom_point()

performance_stat2 |>
  group_by(hallmark) |>
  summarise(accuracy = mean(accuracy)) |>
  mutate(hallmark = gsub("HALLMARK_", "", hallmark)) |>
  mutate(hallmark = gsub("_", " ", hallmark)) |>
  arrange(desc(accuracy)) |> slice(1:10) |>
  XIFF::reorderByScore(orderCol = "hallmark", valueCol = "accuracy") |>
  ggplot(aes(y = hallmark, x = accuracy, fill = hallmark)) + geom_bar(stat = "identity") +
    theme(legend.pos = "none", text = element_text(size = 20)) +
    coord_cartesian(xlim = c(0.88, 0.93)) + ylab("") + xlab("Average accuracy")

ggsave("fig8_CLIFF_ML_hallmarks_top10.pdf", height = 3.28, width = 8.02)

performance_stat2 |>
  mutate(hallmark = gsub("HALLMARK_", "", hallmark)) |>
  mutate(hallmark = gsub("_", " ", hallmark)) |>
  arrange(desc(accuracy)) |>
  ggplot(aes(x = accuracy, y = forcats::fct_reorder(hallmark, accuracy, mean))) + geom_boxplot() +
    theme(text = element_text(size = 20)) + xlab("Accuracy") + ylab("Hallmark") +
    theme(legend.position = "none")

ggsave("fig8_CLIFF_ML_hallmarks_boxplot.pdf", height = 15, width = 8.02)

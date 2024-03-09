getCelllineAnno_4paper <- function() {
  cl_anno <- getCellLineAnno("human") |>
    filter(tumortype %in% c("ovarian cancer", "melanoma", "renal cell carcinoma", "colorectal cancer", "esophagogastric cancer")) |>
    filter(!histology_type %in% c("esophageal adenocarcinoma", "esophageal squamous cell carcinoma")) |>
    mutate(tumortype = ifelse(tumortype == "esophagogastric cancer", "gastric cancer", as.character(tumortype))) |>
    mutate(tumortype = as.factor(tumortype))
}

getGeneAnno_4paper <- function() {
  gene_anno <- getGeneAnno("human") |>
    filter(grepl("^.{1,2}:", location))
}

getMDM2SensitivityClassification <- function(gene_anno_MDM2, cl_anno) {
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

  list(MDM2_avana = MDM2_avana, ca_sens = ca_sens)
}
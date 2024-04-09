library(XIFF)
library(CLIFF)
library(tidyverse)

s <- getSettings()
s <- getCustomSettings(dbName = "bioinfo_23Q2.hg38", dbHost = "charlotte")
print(s$db)
setDbOptions(s)

gene_anno <- getGeneAnno("human") |>
  filter(grepl("ENSG", ensg))

gene_anno2 <- gene_anno |> 
  filter(ensg %in% getCelllineDataGeneExpression("127399_SOFT_TISSUE")$ensg)

cl_anno <- getCellLineAnno("human")
gene_expr <- getCelllineDataGeneExpressionById("ENSG00000000003", cl_anno$celllinename)
cl_anno2 <- cl_anno |>
  filter(celllinename %in% gene_expr$celllinename)

rnaseq_run <- getPostgresql("select * from cellline.rnaseqrun")

getTime <- function(sql, var) {
  con <- XIFF:::getPostgresqlConnection()

  d <- sapply(var, function(v) {
    sql2 <- gsub("%", v, sql)
    rs <- RPostgres::dbSendQuery(con, sql2)
    x <- list(rs = rs, con = con)
    data <- try(RPostgres::dbFetch(x$rs, n = -1))
    RPostgres::dbClearResult(x$rs)
    gsub("Execution Time: ", "", data[nrow(data),])
  })
  RPostgres::dbDisconnect(con)
  d
}

ensg <- sample(gene_anno2$ensg, 1000)
sql <- "explain analyze select * from cellline.processedrnaseq where ensg = '%'"
d1 <- getTime(sql, ensg)

d1a <- data.frame(t = gsub(" ms", "", d1) |> as.numeric())
ggplot(d1a, aes(x = t)) + geom_density() + ggtitle("filter per ENSG")

# -----------------

cl <- sample(cl_anno2$celllinename, 1000)
sql <- "explain analyze select * from cellline.processedrnaseqview where celllinename = '"
d2 <- getTime(sql, cl)
d2a <- data.frame(t = gsub(" ms", "", d2) |> as.numeric())

ggplot(d2a, aes(x = log10(t))) + geom_density() + ggtitle("filter per celllinename")

# ---------------------
rsr <- sample(rnaseq_run$rnaseqrunid, 1000)
sql <-  "explain analyze select * from cellline.processedrnaseq where rnaseqrunid = '"
d3 <- getTime(sql, rsr)

d3a <- data.frame(t = gsub(" ms", "", d3) |> as.numeric())
ggplot(d3a, aes(x = log10(t))) + geom_density() + ggtitle("filter per rnaseqrun")

#save(d1a, d2a, d3a, file = "C:/Users/andre/Boehringer/Broad/CLIFF_Paper/google_performance.Rdata")
#save(d1a, d2a, d3a, file = "C:/Users/andre/Boehringer/Broad/CLIFF_Paper/wsl_performance.Rdata")
#save(d1a, d2a, file = "C:/Users/andre/Boehringer/Broad/CLIFF_Paper/charlotte_performance.Rdata")
load("C:/Users/andre/Boehringer/Broad/CLIFF_Paper/google_performance.Rdata")

result_db <- data.frame(
  source = c("DB", "DB"),
  label = c("per gene", "per cell line"),
  mean = c(mean(d1a$t), mean(d2a$t)),
  sd = c(sd(d1a$t), sd(d2a$t))
)

# ---------------- from flat file ----------------
setwd("\\\\wsl.localhost\\Ubuntu\\home\\andreas\\CLIFF_performance_comparison")

getFileStat <- function(filename) {
  read_delim(filename, delim = " ", col_names = FALSE) |>
    filter(grepl("elapsed", X3)) |>
    mutate(X1 = 1000 * as.numeric(gsub("user", "", X1)),
           X2 = 1000 * as.numeric(gsub("system", "", X2)),
           X3 = gsub("elapsed", "", X3),
           X3 = 1000 * as.numeric(gsub("^0:", "", X3))) |>
    rename(user = X1, system = X2, elapsed = X3)
}

stat <- getFileStat("statFindENSG.txt")
stat2 <- getFileStat("statFindRNAseqrun.txt")

result_flatfile = data.frame(
  source = c("file", "file"), 
  label = c("per gene", "per cell line"),
  mean = c(mean(stat$elapsed), mean(stat2$elapsed)),
  sd = c(sd(stat$elapsed), sd(stat2$elapsed))
)

# ------ result + visualization -----------
result <- result_db |>
  bind_rows(result_flatfile)

p <- ggplot(result, aes(x=label, y = mean, fill = source)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .2,
                position=position_dodge(.9)) 

p + labs(x = "", y = "mean retrieval time [ms]") +
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00')) +
  scale_y_sqrt(breaks = c(10, 30, 100, 300, 1000, 2000, 3000)) +
  theme(text = element_text(size = 20))

ggsave("CLIFF_rnaseq_file_db_performance.pdf", width = 6, height = 7)

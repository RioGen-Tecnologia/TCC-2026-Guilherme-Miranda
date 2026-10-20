# ============== CARREGANDO PACOTES ==============

message("===========================================================================")
message("CARREGANDO PACOTES")
message("===========================================================================")

library(here)
source(here("scripts", "setup.r"))


# Carregar o modelo treinado
modelo_enet <- readRDS(here("results","machine_learning","modelo_elastic_net.rds"))
desempenho <- read.csv(here("results","machine_learning","desempenho_elastic_net.csv"))
genes_modelo <- read.csv(here("results","machine_learning","genes_selecionados_lambda_min_com_limma.csv"))

# ============== GSE166716 ==============

message("===========================================================================")
message("BAIXANDO GSE166716")
message("===========================================================================")

gse <- getGEO(
  "GSE166716",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

message("===========================================================================")
message("ORGANIZANDO GSE166716")
message("===========================================================================")

length(gse)
gse <- gse[[1]]

expr_gse <- exprs(gse)
pheno_gse <- pData(gse)

dim(expr_gse)
head(pheno_gse[, 1:5])

grupo_gse <- ifelse(
  grepl(
    "normal urothelium",
    pheno_gse$`tissue:ch1`,
    ignore.case = TRUE
  ),
  "non_tumor",
  "tumor"
)

table(grupo_gse)

message("===========================================================================")
message("ANOTANDO GSE166716")
message("===========================================================================")

annotation <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(expr_gse),
  columns = c("ENTREZID"),
  keytype = "PROBEID"
)

annotation <- annotation[
  !is.na(annotation$ENTREZID),
]

# fazendo média de probes duplicados

expr_gse <- expr_gse[
  rownames(expr_gse) %in% annotation$PROBEID,
]

annotation <- annotation[
  match(rownames(expr_gse), annotation$PROBEID),
]

rownames(expr_gse) <- annotation$ENTREZID

expr_gse_entrez <- limma::avereps(
  expr_gse,
  ID = rownames(expr_gse)
)

dim(expr_gse_entrez)
head(rownames(expr_gse_entrez))

# conferindo presença de genes

coeficientes <- as.matrix(
  coef(modelo_enet, s = "lambda.min")
)

coeficientes <- coeficientes[
  rownames(coeficientes) != "(Intercept)",
  ,
  drop = FALSE
]

genes_enet <- rownames(coeficientes)[
  coeficientes[, 1] != 0
]

genes_deg <- genes_modelo$EntrezID[
  genes_modelo$DEG_limma
]

cat(
  "Genes do modelo encontrados no GSE166716:",
  sum(genes_enet %in% rownames(expr_gse_entrez)),
  "de",
  length(genes_enet),
  "\n"
)

cat(
  "DEGs do modelo encontrados no GSE166716:",
  sum(genes_deg %in% rownames(expr_gse_entrez)),
  "de",
  length(genes_deg),
  "\n"
)









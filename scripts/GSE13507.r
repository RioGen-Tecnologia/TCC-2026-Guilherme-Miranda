# GSE13507.R
# Guilherme Moret Miranda - Riogen
# 27/04/2026
# Download: GEOquery
# Normalização normalizeBetweenArrays com pacote limma
# Anotação: anotação fornecida nos dados

# ============== PACOTES ==============

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(R.utils)
library(illuminaHumanv2.db)
library(dplyr)

# ============== EXTRAÇÃO DE DADOS ==============

id_projeto <- "GSE13507"

# carregando o master manifesto
metadata <- read.csv(metadata_path)
metadata <- metadata[metadata$study_ID == "GSE13507", c("sample_ID", "sample_type", "characteristics")]
rownames(metadata) <- metadata$sample_ID
metadata$sample_ID <- NULL

# --- Etapa de Download ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Baixando dados brutos para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Download (usando o geo_dir definido no mestre)
if(!dir.exists(file.path(geo_dir, id_projeto))) {
  getGEOSuppFiles(id_projeto, baseDir = geo_dir)
}

# definindo pasta do projeto
projeto_dir <- file.path(geo_dir, id_projeto)

# Descompactando
arquivos_tar <- list.files(
  path = projeto_dir,
  pattern = "\\.tar$",
  full.names = TRUE
)
untar(arquivos_tar, exdir = projeto_dir)

# ====== Leitura dos dados brutos illumina ======

# --- Etapa de leitura ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Lendo dados para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

arquivos <- list.files(
  path = projeto_dir,
  pattern = ".txt$",
  full.names = TRUE
)
cels.GSE13507 <- read.delim(arquivos, stringsAsFactors = FALSE)

# extrai apenas os valores como numeric
dados_brutos_GSE13507 <- matrix(as.numeric(as.matrix(cels.GSE13507[-1, -1])),
                       nrow = nrow(cels.GSE13507)-1,
                       ncol = ncol(cels.GSE13507)-1)

# define os rownames novamente
rownames(dados_brutos_GSE13507) <- cels.GSE13507[-1, 1]
colnames(dados_brutos_GSE13507) <- toupper(sub(".*(GSM[0-9]+).*", "\\1", colnames(cels.GSE13507[-1, -1])))

# ====== Normalização dos dados ======

# --- Etapa de Normalização ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Normalizando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Log2 transform
exprs_log <- log2(pmax(dados_brutos_GSE13507, 1))

# Quantile normalization
exprs_norm <- normalizeBetweenArrays(exprs_log, method = "quantile")

# removem genes que apresentem NA
norm_corrigido_GSE13507 <- exprs_norm[!is.na(rownames(exprs_norm)), ]
norm_corrigido_GSE13507 <- norm_corrigido_GSE13507[rowSums(is.na(norm_corrigido_GSE13507)) == 0, ]
norm_corrigido_GSE13507 <- t(norm_corrigido_GSE13507)

# ====== remove amostras fora do manifesto ======

ids <- intersect(rownames(norm_corrigido_GSE13507), rownames(metadata))

norm_corrigido_GSE13507 <- norm_corrigido_GSE13507[ids, ]
metadata <- metadata[ids, ]


# ====== CONTROLE DE QUALIDADE DOS ARRAYS ======

message("\n", paste(rep("=", 30), collapse = ""))
message("Realizando controle de qualidade de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Diretório para resultados de QC
qc_dir <- file.path(processed_dir, id_projeto, "QC_2")

if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
}

# ---------- RESUMO DA DISTRIBUIÇÃO DAS INTENSIDADES ----------

id <- "GSE13507"

qc_intensidade <- data.frame(
  Sample = rownames(norm_corrigido_GSE13507),
  Median = apply(
    norm_corrigido_GSE13507,
    1,
    median,
    na.rm = TRUE
  ),
  Q1 = apply(
    norm_corrigido_GSE13507,
    1,
    quantile,
    probs = 0.25,
    na.rm = TRUE
  ),
  Q3 = apply(
    norm_corrigido_GSE13507,
    1,
    quantile,
    probs = 0.75,
    na.rm = TRUE
  )
)

qc_intensidade$IQR <-
  qc_intensidade$Q3 - qc_intensidade$Q1

qc_intensidade <- qc_intensidade[
  order(qc_intensidade$Median),
]

print(qc_intensidade)

# ---------- PCA ----------

pca <- prcomp(
  norm_corrigido_GSE13507,
  center = TRUE,
  scale. = FALSE
)

var_exp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

pca_df$Group <- metadata[rownames(pca_df), "sample_type"]

png(
  file.path(qc_dir, "PCA_GSE13507_limpo.png"),
  width = 1800,
  height = 1400,
  res = 150
)

cores <- as.numeric(as.factor(pca_df$Group))

plot(
  pca_df$PC1,
  pca_df$PC2,
  pch = 19,
  col = cores,
  xlab = paste0("PC1 (", round(var_exp[1], 2), "%)"),
  ylab = paste0("PC2 (", round(var_exp[2], 2), "%)"),
  main = "PCA - GSE13507"
)

legend(
  "topright",
  legend = levels(as.factor(pca_df$Group)),
  col = seq_along(levels(as.factor(pca_df$Group))),
  pch = 19,
  bty = "n"
)

dev.off()

# ---------- CORRELAÇÃO ENTRE AMOSTRAS ----------

cor_mat <- cor(
  t(norm_corrigido_GSE13507),
  method = "pearson",
  use = "pairwise.complete.obs"
)

png(
  file.path(qc_dir, "Correlacao_GSE13507.png"),
  width = 1800,
  height = 1600,
  res = 150
)

heatmap(
  cor_mat,
  symm = TRUE,
  scale = "none",
  margins = c(10, 10),
  main = "Correlação de Pearson - GSE13507"
)

dev.off()

mean_cor <- sapply(
  1:nrow(cor_mat),
  function(i) mean(cor_mat[i, -i], na.rm = TRUE)
)

mean_cor <- sort(mean_cor)

mean_cor

# ---------- HIERARCHICAL CLUSTERING ----------

dist_mat <- dist(
  norm_corrigido_GSE13507,
  method = "euclidean"
)

hc <- hclust(
  dist_mat,
  method = "complete"
)

png(
  file.path(qc_dir, "Hierarchical_GSE13507.png"),
  width = 2200,
  height = 1400,
  res = 150
)

plot(
  hc,
  labels = rownames(norm_corrigido_GSE13507),
  main = "Hierarchical clustering - GSE13507",
  xlab = "",
  sub = "",
  cex = 0.55
)

dev.off()


# ====== anotação com EntrezID ======

# --- Etapa de anotação ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

arquivos_gunzip <- list.files(path = projeto_dir,pattern = "\\.bgx\\.gz$",full.names = TRUE)
gunzip(arquivos_gunzip,remove = FALSE,overwrite = TRUE)
arquivos_gunzip <- list.files(path = projeto_dir,pattern = "\\.bgx$",full.names = TRUE)

# Ler as primeiras linhas para inspecionar
linhas <- readLines(arquivos_gunzip, n = 100)
head(linhas, 20)  # ver onde começa "ProbeID" ou "ID"

# Encontrar a linha que contém o cabeçalho
linha_cabecalho <- grep("^Species", linhas)
linha_cabecalho

anotação_ref <- read.delim(arquivos_gunzip,
                           skip = linha_cabecalho - 1,
                           stringsAsFactors = FALSE)


# 1. criar tabela de anotação com apenas Probe_Id e EntrezID
anot <- anotação_ref %>%
  select(Probe_Id, Entrez_Gene_ID) %>%
  filter(!is.na(Entrez_Gene_ID))  # remove probes sem EntrezID

# 2. filtrar a matriz para conter apenas probes com EntrezID
expr_filtrada <- norm_corrigido_GSE13507[, colnames(norm_corrigido_GSE13507) %in% anot$Probe_Id]

# 3. garantir que a ordem dos probes da anotação bate com a matriz
anot <- anot[match(colnames(expr_filtrada), anot$Probe_Id), ]

# 4. substituir colnames pelos EntrezIDs
colnames(expr_filtrada) <- anot$Entrez_Gene_ID

# Remove genes NA
norm_corrigido_GSE13507 <- expr_filtrada[, !is.na(colnames(expr_filtrada))]
rm(expr_filtrada)

# Seleciona genes duplicados que tem maior variância
# Caso haja genes duplicados, é feito uma correlação de pearson.
# se a correlação for alta (0,7), é feita média dos sinais
# se a correlação for baixa, é selecionado o probe com maior variância
# se existirem probes triplicados ou mais, a comparação é realizada em clusters

# identificar grupos de probes (genes duplicados)
genes <- colnames(norm_corrigido_GSE13507)
grupos <- split(seq_along(genes), genes)

idx_final <- unlist(lapply(grupos, function(i) {
  
  # caso não haja duplicata
  if (length(i) == 1) return(i)
  
  submat <- norm_corrigido_GSE13507[, i, drop = FALSE]
  
  # correlação entre probes do mesmo gene
  cor_mat <- cor(submat, use = "pairwise.complete.obs")
  
  # se tudo NA → fallback variância
  if (all(is.na(cor_mat))) {
    vars <- apply(submat, 2, var, na.rm = TRUE)
    return(i[which.max(vars)])
  }
  
  # distância baseada em correlação
  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat, method = "average")
  
  clusters <- cutree(hc, h = 1 - 0.7)  # threshold 0.7
  
  # maior cluster
  tab <- table(clusters)
  main_cluster <- as.numeric(names(tab)[which.max(tab)])
  idx_cluster <- i[clusters == main_cluster]
  
  # se cluster confiável, usa média; senão variância
  if (length(idx_cluster) >= 2) {
    return(idx_cluster[1])  # ou poderia usar média depois
  } else {
    vars <- apply(submat, 2, var, na.rm = TRUE)
    return(i[which.max(vars)])
  }
}))

norm_corrigido_GSE13507 <- norm_corrigido_GSE13507[, idx_final]

# ====== Salvar arquivo ======

# confere se o pData e matriz estão realmente alinhados
all(rownames(norm_corrigido_GSE13507) == rownames(metadata))

#confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(processed_dir, id_projeto)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# arquivo de matriz de expressão
saveRDS(norm_corrigido_GSE13507,
        file = file.path(processed_dir,
                         id_projeto,
                         "exprs_GSE13507.rds"))
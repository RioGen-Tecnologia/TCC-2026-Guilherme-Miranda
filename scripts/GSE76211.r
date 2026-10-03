
# GSE76211.R
# Guilherme Moret Miranda - Riogen
# 27/04/2026
# Download: GEOquery
# Normalização RMA com pacote oligo
# Anotação: hta20transcriptcluster.db

# ============== PACOTES ==============

library(GEOquery)
library(oligo)
library(limma)
library(AnnotationDbi)
library(hta20transcriptcluster.db)
library(pd.hta.2.0)
library(RSQLite)
library(DBI)

# ============== EXTRAÇÃO DE DADOS ==============

id_projeto <- "GSE76211"

# carregando o master manifesto
metadata <- read.csv(metadata_path)
metadata <- metadata[metadata$study_ID == "GSE76211", c("sample_ID", "sample_type", "characteristics")]
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

# ====== Leitura dos dados brutos ======

# --- Etapa de leitura ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Lendo dados para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Cria uma variável com arquivos que possuem "CEL.gz" e os imprime (exclue os CHP.gz)
cels.GSE76211 <- list.files(
  path = projeto_dir,
  pattern = "[cC][eE][lL]\\.gz$",
  full.names = TRUE
)

# Lendo os dados brutos (raw data)
dados_brutos_GSE76211 <- oligo::read.celfiles(cels.GSE76211)


# ====== CONTROLE DE QUALIDADE DOS ARRAYS ======

# objeto para QC
dados_QC_GSE76211 <- dados_brutos_GSE76211
sampleNames(dados_QC_GSE76211) <- sub("_.*$","",sampleNames(dados_QC_GSE76211))

message("\n", paste(rep("=", 30), collapse = ""))
message("Realizando controle de qualidade de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

medianas_arrays <- apply(
  log2(Biobase::exprs(dados_QC_GSE76211)),
  2,
  median,
  na.rm = TRUE
)

mediana_global <- median(
  medianas_arrays,
  na.rm = TRUE
)

# Diretório para resultados de QC
qc_dir <- file.path(processed_dir, id_projeto, "QC_2")

if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
}

# ---------- DISTRIBUIÇÃO DAS INTENSIDADES BRUTAS ----------

png(
  filename = file.path(qc_dir, "Intensidades_brutas_GSE76211.png"),
  width = 1800,
  height = 1200,
  res = 150
)

boxplot(
  dados_QC_GSE76211,
  target = "core",
  main = paste("Distribuição das intensidades brutas -", id_projeto),
  ylab = "Log2(intensidade)",
  las = 2,
  outline = FALSE
)
abline(
  h = mediana_global,
  lty = 2,
  col = "black",
  lwd = 1.5
)

dev.off()


# ====== INVESTIGAÇÃO ADICIONAL DE QUALIDADE ======

message("\n", paste(rep("=", 30), collapse = ""))
message("Gerando análises adicionais de QC...")
message(paste(rep("=", 30), collapse = ""))

# ---------- PCA DAS INTENSIDADES BRUTAS ----------
exprs_log2_GSE76211 <- log2(Biobase::exprs(dados_QC_GSE76211))

pca_GSE76211 <- prcomp(
  t(exprs_log2_GSE76211),
  center = TRUE,
  scale. = FALSE
)

# Variância explicada
var_pca_GSE76211 <- pca_GSE76211$sdev^2
var_pca_GSE76211 <- var_pca_GSE76211 / sum(var_pca_GSE76211) * 100


png(filename = file.path(qc_dir,"PCA_GSE76211.png"),
    width = 1800,
    height = 1400,
    res = 150)
plot(
  pca_GSE76211$x[, 1],
  pca_GSE76211$x[, 2],
  pch = 19,
  xlab = paste0(
    "PC1 (",
    round(var_pca_GSE76211[1], 2),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(var_pca_GSE76211[2], 2),
    "%)"
  ),
  main = paste(
    "PCA das intensidades brutas -",
    id_projeto
  )
)

text(
  pca_GSE76211$x[, 1],
  pca_GSE76211$x[, 2],
  labels = colnames(dados_QC_GSE76211),
  pos = 3,
  cex = 0.55
)
dev.off()


# ---------- CORRELAÇÃO ENTRE ARRAYS ----------
cor_GSE76211 <- cor(
  exprs_log2_GSE76211,
  method = "pearson",
  use = "pairwise.complete.obs")

png(filename = file.path(qc_dir,"Correlacao_arrays_GSE76211.png"),width = 1800,height = 1600,res = 150)
heatmap(
  cor_GSE76211,
  Rowv = NA,
  Colv = NA,
  scale = "none",
  symm = TRUE,
  margins = c(10, 10),
  main = paste(
    "Correlação entre arrays -",
    id_projeto
  )
)
dev.off()

# ---------- DISTÂNCIA ENTRE ARRAYS ----------
dist_GSE76211 <- dist(t(exprs_log2_GSE76211))

png(filename = file.path(qc_dir,"Cluster_arrays_GSE76211.png"),width = 1800,height = 1400,res = 150)

plot(
  hclust(dist_GSE76211),
  main = paste(
    "Agrupamento hierárquico dos arrays -",
    id_projeto
  ),
  xlab = "",
  sub = "",
  cex = 0.6
)

dev.off()


# ====== Normalização dos dados ======

# --- Etapa de Normalização ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Normalizando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# normalização
dados_norm_GSE76211 <- oligo::rma(dados_brutos_GSE76211)

# extração de matriz de expressão normalizada
norm_GSE76211 <- exprs(dados_norm_GSE76211)

# remove o ".CEL.gz" ou outras interferências da tabela
colnames(norm_GSE76211) <- toupper(sub(".*(GSM[0-9]+).*", "\\1", colnames(norm_GSE76211)))

# removem genes que apresentem NA
norm_corrigido_GSE76211 <- norm_GSE76211[!is.na(rownames(norm_GSE76211)), ]
norm_corrigido_GSE76211 <- norm_corrigido_GSE76211[rowSums(is.na(norm_corrigido_GSE76211)) == 0, ]
norm_corrigido_GSE76211 <- t(norm_corrigido_GSE76211)

# ====== remove amostras fora do manifesto ======

ids <- intersect(rownames(norm_corrigido_GSE76211), rownames(metadata))

norm_corrigido_GSE76211 <- norm_corrigido_GSE76211[ids, ]
metadata <- metadata[ids, ]

# ====== anotação com EntrezID ======

# --- Etapa de anotação ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

entrez_ids <- mapIds(hta20transcriptcluster.db,
                     keys = colnames(norm_corrigido_GSE76211),
                     column = "ENTREZID",
                     keytype = "PROBEID",
                     multiVals = "first")

colnames(norm_corrigido_GSE76211) <- entrez_ids

# Remove genes NA
norm_corrigido_GSE76211 <- norm_corrigido_GSE76211[, !is.na(colnames(norm_corrigido_GSE76211))]



# Caso haja genes duplicados, é feito uma correlação de pearson.
# se a correlação for alta (0,7), é selecionado o primeiro probe.
# se a correlação for baixa, é selecionado o probe com maior variância
# se existirem probes triplicados ou mais, a comparação é realizada em clusters

# identificar grupos de probes (genes duplicados)
genes <- colnames(norm_corrigido_GSE76211)
grupos <- split(seq_along(genes), genes)

idx_final <- unlist(lapply(grupos, function(i) {
  
  # caso não haja duplicata
  if (length(i) == 1) return(i)
  
  submat <- norm_corrigido_GSE76211[, i, drop = FALSE]
  
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

# subset final
norm_corrigido_GSE76211 <- norm_corrigido_GSE76211[, idx_final]

# ====== Salvar arquivo ======

# confere se o pData e matriz estão realmente alinhados
all(rownames(norm_corrigido_GSE76211) == rownames(metadata))

#confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(processed_dir, id_projeto)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# arquivo de matriz de expressão
saveRDS(norm_corrigido_GSE76211,
        file = file.path(processed_dir,
                         id_projeto,
                         "exprs_GSE76211.rds"))

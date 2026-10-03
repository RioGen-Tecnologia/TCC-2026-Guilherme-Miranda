# GSE52519.R
# Guilherme Moret Miranda - Riogen
# 27/04/2026
# Download: GEOquery
# Normalização neqc com pacote limma
# Anotação: anotação fornecida nos dados

# ============== PACOTES ==============

library(GEOquery)
library(limma)
library(AnnotationDbi)

# ============== EXTRAÇÃO DE DADOS ==============

id_projeto <- "GSE52519"

# carregando o master manifesto
metadata <- read.csv(metadata_path)
metadata <- metadata[metadata$study_ID == "GSE52519", c("sample_ID", "sample_type", "characteristics")]
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

arquivos_gunzip <- list.files(
  path = projeto_dir,
  pattern = "\\.gz$",
  full.names = TRUE)

for (arquivo in arquivos_gunzip) {
  gunzip(
    arquivo,
    remove = FALSE,
    overwrite = TRUE)}

# ====== Leitura dos dados brutos illumina ======
# A normalização neqc() do pacote limma foi feito para trabalhar com objetos EListRaw contendo:

# --- Etapa de leitura ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Lendo dados para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

non_normalized <- list.files(
  path = projeto_dir,
  pattern = "non-normalized.txt$",
  full.names = TRUE
)

# Cria uma variável com arquivos .txt
cels.GSE52519 <- read.delim(non_normalized, check.names = FALSE)

# Remover coluna de IDs
dados <- cels.GSE52519[,-1]

# Colunas de expressão (1,3,5...)
expr <- dados[, seq(1, ncol(dados), by=2)]

# Colunas de detection p-value (2,4,6...)
detP <- dados[, seq(2, ncol(dados), by=2)]

#baixa o pData para renomear colunas
pData_GSE52519 <- pData(getGEO("GSE52519")[[1]])

colnames(expr) <- rownames(pData_GSE52519)
colnames(detP) <- rownames(pData_GSE52519)

rownames(expr) <- cels.GSE52519$ID_REF
rownames(detP) <- cels.GSE52519$ID_REF

elist <- new("EListRaw")
elist$E <- as.matrix(expr)
elist$other$Detection <- as.matrix(detP)

# ====== Normalização dos dados ======

# --- Etapa de Normalização ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Normalizando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# checagem
stopifnot(length(non_normalized) == 1)

dados <- read.delim(non_normalized, check.names = FALSE)

expr <- dados[, grep("Sample", colnames(dados))]
detP <- dados[, grep("Detection Pval", colnames(dados))]

pData_GSE52519 <- pData(getGEO("GSE52519")[[1]])

stopifnot(ncol(expr) == nrow(pData_GSE52519))

colnames(expr) <- rownames(pData_GSE52519)
colnames(detP) <- rownames(pData_GSE52519)

rownames(expr) <- dados$ID_REF
rownames(detP) <- dados$ID_REF

elist <- new("EListRaw")
elist$E <- as.matrix(expr)
elist$other$Detection <- as.matrix(detP)

# filtro importante
keep <- rowSums(detP < 0.05) >= (0.5 * ncol(detP))
elist <- elist[keep, ]

norm_GSE52519 <- neqc(elist)

norm_corrigido_GSE52519 <- t(norm_GSE52519$E)

# ====== remove amostras fora do manifesto ======

ids <- intersect(rownames(norm_corrigido_GSE52519), rownames(metadata))

norm_corrigido_GSE52519 <- norm_corrigido_GSE52519[ids, ]
metadata <- metadata[ids, ]

# ====== CONTROLE DE QUALIDADE ======

# criando diretório
qc_dir <- file.path(processed_dir, id_projeto, "QC_2")

if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
}

# ---------- DISTRIBUIÇÃO DAS INTENSIDADES BRUTAS ----------

medianas_arrays <- apply(log2(elist$E), 2, median, na.rm = TRUE)
mediana_global <- median(medianas_arrays, na.rm = TRUE)

png(
  file.path(qc_dir, "Intensidades_brutas_GSE52519.png"),
  1800, 1200, res = 150
)

boxplot(
  log2(elist$E),
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

# ---------- DISTRIBUIÇÃO APÓS NEQC ----------

medianas_norm <- apply(norm_GSE52519$E, 2, median, na.rm = TRUE)
mediana_global_norm <- median(medianas_norm, na.rm = TRUE)

png(
  file.path(qc_dir, "Intensidades_normalizadas_GSE52519.png"),
  1800, 1200, res = 150
)

boxplot(
  norm_GSE52519$E,
  main = paste("Distribuição das intensidades após neqc -", id_projeto),
  ylab = "Log2(intensidade)",
  las = 2,
  outline = FALSE
)

abline(
  h = mediana_global_norm,
  lty = 2,
  col = "black",
  lwd = 1.5
)

dev.off()

# ---------- PCA ----------

exprs_log2_GSE52519 <- norm_GSE52519$E

pca_GSE52519 <- prcomp(
  t(exprs_log2_GSE52519),
  center = TRUE,
  scale. = FALSE
)

var_pca_GSE52519 <- pca_GSE52519$sdev^2

var_pca_GSE52519 <- var_pca_GSE52519 /
  sum(var_pca_GSE52519) * 100

png(
  file.path(qc_dir, "PCA_GSE52519.png"),
  1800, 1400, res = 150
)

plot(
  pca_GSE52519$x[,1],
  pca_GSE52519$x[,2],
  pch = 19,
  xlab = paste0(
    "PC1 (",
    round(var_pca_GSE52519[1], 2),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(var_pca_GSE52519[2], 2),
    "%)"
  ),
  main = paste(
    "PCA das expressões normalizadas -",
    id_projeto
  )
)

text(
  pca_GSE52519$x[,1],
  pca_GSE52519$x[,2],
  labels = rownames(pca_GSE52519$x),
  pos = 3,
  cex = 0.6
)

dev.off()

# ---------- CORRELAÇÃO ENTRE ARRAYS ----------

cor_GSE52519 <- cor(
  exprs_log2_GSE52519,
  method = "pearson",
  use = "pairwise.complete.obs"
)

png(
  file.path(qc_dir, "Correlacao_arrays_GSE52519.png"),
  1800, 1600, res = 150
)

heatmap(
  cor_GSE52519,
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

# ---------- CLUSTERING HIERÁRQUICO ----------

dist_GSE52519 <- dist(
  t(exprs_log2_GSE52519)
)

png(
  file.path(qc_dir, "Cluster_arrays_GSE52519.png"),
  1800, 1400, res = 150
)

plot(
  hclust(dist_GSE52519),
  main = paste(
    "Agrupamento hierárquico dos arrays -",
    id_projeto
  ),
  xlab = "",
  sub = "",
  cex = 0.7
)

dev.off()

# ---------- DETECTION P-VALUE ----------

detP <- elist$other$Detection

proporcao_detectados <- colMeans(
  detP < 0.05,
  na.rm = TRUE
)

png(
  file.path(qc_dir, "Detection_P_GSE52519.png"),
  1800, 1200, res = 150
)

barplot(
  proporcao_detectados,
  names.arg = names(proporcao_detectados),
  las = 2,
  ylim = c(0, 1),
  ylab = "Proporção de probes detectados (p < 0,05)",
  main = paste(
    "Detecção de expressão -",
    id_projeto
  )
)

abline(
  h = mean(proporcao_detectados),
  lty = 2,
  col = "black",
  lwd = 1.5
)

dev.off()


# ====== anotação com EntrezID ======

# --- Etapa de anotação ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

bgx <- list.files(
  path = projeto_dir,
  pattern = "\\.bgx$",
  full.names = TRUE)


# Ler as primeiras linhas para inspecionar
linhas <- readLines(bgx, n = 100)
head(linhas, 20)  # ver onde começa "ProbeID" ou "ID"

# Encontrar a linha que contém o cabeçalho
linha_cabecalho <- grep("^Species", linhas)
linha_cabecalho

anotação_ref <- read.delim(bgx,
                           skip = linha_cabecalho - 1,
                           stringsAsFactors = FALSE)


# 1. criar tabela de anotação com apenas Probe_Id e EntrezID
anot <- anotação_ref %>%
  select(Probe_Id, Entrez_Gene_ID) %>%
  filter(!is.na(Entrez_Gene_ID))  # remove probes sem EntrezID

# 2. filtrar a matriz para conter apenas probes com EntrezID
expr_filtrada <- norm_corrigido_GSE52519[, colnames(norm_corrigido_GSE52519) %in% anot$Probe_Id]

# 3. garantir que a ordem dos probes da anotação bate com a matriz
anot <- anot[match(colnames(expr_filtrada), anot$Probe_Id), ]

# 4. substituir colnames pelos EntrezIDs
colnames(expr_filtrada) <- anot$Entrez_Gene_ID

# Remove genes NA
norm_corrigido_GSE52519 <- expr_filtrada[, !is.na(colnames(expr_filtrada))]
rm(expr_filtrada)

# Seleciona genes duplicados que tem maior variância
# Caso haja genes duplicados, é feito uma correlação de pearson.
# se a correlação for alta (0,7), é feita média dos sinais
# se a correlação for baixa, é selecionado o probe com maior variância
# se existirem probes triplicados ou mais, a comparação é realizada em clusters

# identificar grupos de probes (genes duplicados)
genes <- colnames(norm_corrigido_GSE52519)
grupos <- split(seq_along(genes), genes)

idx_final <- unlist(lapply(grupos, function(i) {
  
  # caso não haja duplicata
  if (length(i) == 1) return(i)
  
  submat <- norm_corrigido_GSE52519[, i, drop = FALSE]
  
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

norm_corrigido_GSE52519 <- norm_corrigido_GSE52519[, idx_final]

# ====== Salvar arquivo ======

# confere se o pData e matriz estão realmente alinhados
all(rownames(norm_corrigido_GSE52519) == rownames(metadata))

#confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(processed_dir, id_projeto)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# arquivo de matriz de expressão
saveRDS(norm_corrigido_GSE52519,
        file = file.path(processed_dir,
                         id_projeto,
                         "exprs_GSE52519.rds"))

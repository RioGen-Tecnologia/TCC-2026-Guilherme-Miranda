
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

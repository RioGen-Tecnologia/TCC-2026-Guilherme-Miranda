
# GSE37817.R
# Guilherme Moret Miranda - Riogen
# 27/04/2026
# Download: GEOquery
# Anotação: illuminaHumanv2.db

# ============== PACOTES ==============

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(R.utils)
library(illuminaHumanv2.db)

# ============== EXTRAÇÃO DE DADOS ==============

id_projeto <- "GSE37817"

# carregando o master manifesto
metadata <- read.csv(metadata_path)
metadata <- metadata[metadata$study_ID == "GSE37817", c("sample_ID", "sample_type", "characteristics")]
rownames(metadata) <- metadata$sample_ID
metadata$sample_ID <- NULL


# Criação da pasta se não existir
# Download (usando o geo_dir definido no mestre)
if(!dir.exists(file.path(geo_dir, id_projeto))) {
  getGEOSuppFiles(id_projeto, baseDir = geo_dir)
}

# ====== Leitura dos dados brutos illumina ======
# Este projeto não possui dados brutos disponíveis. A única opção é baixar a matriz normalizada do geo

# --- Etapa de Download ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Baixando dados brutos para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Baixa o GSE
gse <- getGEO('GSE37817', GSEMatrix = TRUE)

# --- Etapa de leitura ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Lendo dados para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Normalmente é a primeira plataforma
norm_GSE37817 <- exprs(gse[[1]])

# remove o ".CEL.gz" ou outras interferências da tabela
colnames(norm_GSE37817) <- toupper(sub(".*(GSM[0-9]+).*", "\\1", colnames(norm_GSE37817)))

# removem genes que apresentem NA
norm_corrigido_GSE37817 <- norm_GSE37817[!is.na(rownames(norm_GSE37817)), ]
norm_corrigido_GSE37817 <- norm_corrigido_GSE37817[rowSums(is.na(norm_corrigido_GSE37817)) == 0, ]
norm_corrigido_GSE37817 <- t(norm_corrigido_GSE37817)

# ====== remove amostras fora do manifesto ======

ids <- intersect(rownames(norm_corrigido_GSE37817), rownames(metadata))

norm_corrigido_GSE37817 <- norm_corrigido_GSE37817[ids, ]
metadata <- metadata[ids, ]

# ====== anotação com EntrezID ======

# --- Etapa de anotação ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))


entrez_ids <- mapIds(illuminaHumanv2.db,
                     keys = colnames(norm_corrigido_GSE37817),
                     column = "ENTREZID",
                     keytype = "PROBEID",
                     multiVals = "first")

colnames(norm_corrigido_GSE37817) <- entrez_ids

# Remove genes NA
norm_corrigido_GSE37817 <- norm_corrigido_GSE37817[, !is.na(colnames(norm_corrigido_GSE37817))]

# Seleciona genes duplicados que tem maior variância
# Caso haja genes duplicados, é feito uma correlação de pearson.
# se a correlação for alta (0,7), é feita média dos sinais
# se a correlação for baixa, é selecionado o probe com maior variância
# se existirem probes triplicados ou mais, a comparação é realizada em clusters

# identificar grupos de probes (genes duplicados)
genes <- colnames(norm_corrigido_GSE37817)
grupos <- split(seq_along(genes), genes)

idx_final <- unlist(lapply(grupos, function(i) {
  
  # caso não haja duplicata
  if (length(i) == 1) return(i)
  
  submat <- norm_corrigido_GSE37817[, i, drop = FALSE]
  
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

norm_corrigido_GSE37817 <- norm_corrigido_GSE37817[, idx_final]

# ====== Salvar arquivo ======

# confere se o pData e matriz estão realmente alinhados
all(rownames(norm_corrigido_GSE37817) == rownames(metadata))

#confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(processed_dir, id_projeto)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# arquivo de matriz de expressão
saveRDS(norm_corrigido_GSE37817,
        file = file.path(processed_dir,
                         id_projeto,
                         "exprs_GSE37817.rds"))

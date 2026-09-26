
# GSE65635.R
# Guilherme Moret Miranda - Riogen
# 27/04/2026
# Download: GEOquery
# Normalização neqc com pacote limma
# Anotação: illuminaHumanv4.db

# ============== PACOTES ==============

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(R.utils)
library(illuminaHumanv4.db)

# ============== EXTRAÇÃO DE DADOS ==============

id_projeto <- "GSE65635"

# carregando o master manifesto
metadata <- read.csv(metadata_path)
metadata <- metadata[metadata$study_ID == "GSE65635", c("sample_ID", "sample_type", "characteristics")]
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
# extraíndo
arquivo_zip <- list.files(
  path = projeto_dir,
  pattern = "\\.gz$",
  full.names = TRUE)
R.utils::gunzip(arquivo_zip, remove = FALSE, overwrite = TRUE)

# ====== Leitura dos dados brutos illumina ======

# --- Etapa de leitura ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Lendo dados para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Cria uma variável com arquivos .txt

arquivo_txt <- file.path(projeto_dir, "GSE65635_non-normalized.txt")
cels.GSE65635 <- read.delim(arquivo_txt,check.names = FALSE)

# Remover coluna de IDs
dados <- cels.GSE65635[,-1]

# Colunas de expressão (1,3,5...)
expr <- dados[, seq(1, ncol(dados), by=2)]

# Colunas de detection p-value (2,4,6...)
detP <- dados[, seq(2, ncol(dados), by=2)]

#baixa o pData para renomear colunas
pData_GSE65635 <- pData(getGEO("GSE65635")[[1]])

colnames(expr) <- rownames(pData_GSE65635)
colnames(detP) <- rownames(pData_GSE65635)

rownames(expr) <- cels.GSE65635$ID_REF
rownames(detP) <- cels.GSE65635$ID_REF

elist <- new("EListRaw")
elist$E <- as.matrix(expr)
elist$other$Detection <- as.matrix(detP)

# ====== Normalização dos dados ======

# --- Etapa de Normalização ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Normalizando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# normalização neqc()
dados_norm_GSE65635 <- neqc(elist)

# extração de matriz de expressão normalizada
norm_GSE65635 <- dados_norm_GSE65635$E

# remove o ".CEL.gz" ou outras interferências da tabela
colnames(norm_GSE65635) <- toupper(sub(".*(GSM[0-9]+).*", "\\1", colnames(norm_GSE65635)))

# removem genes que apresentem NA
norm_corrigido_GSE65635 <- norm_GSE65635[!is.na(rownames(norm_GSE65635)), ]
norm_corrigido_GSE65635 <- norm_corrigido_GSE65635[rowSums(is.na(norm_corrigido_GSE65635)) == 0, ]
norm_corrigido_GSE65635 <- t(norm_corrigido_GSE65635)

# ====== remove amostras fora do manifesto ======

ids <- intersect(rownames(norm_corrigido_GSE65635), rownames(metadata))

norm_corrigido_GSE65635 <- norm_corrigido_GSE65635[ids, ]
metadata <- metadata[ids, ]

# ====== anotação com EntrezID ======

# --- Etapa de anotação ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

entrez_ids <- mapIds(illuminaHumanv4.db,
                     keys = colnames(norm_corrigido_GSE65635),
                     column = "ENTREZID",
                     keytype = "PROBEID",
                     multiVals = "first")

colnames(norm_corrigido_GSE65635) <- entrez_ids

# Remove genes NA
norm_corrigido_GSE65635 <- norm_corrigido_GSE65635[, !is.na(colnames(norm_corrigido_GSE65635))]

# Seleciona genes duplicados que tem maior variância
# Caso haja genes duplicados, é feito uma correlação de pearson.
# se a correlação for alta (0,7), é feita média dos sinais
# se a correlação for baixa, é selecionado o probe com maior variância
# se existirem probes triplicados ou mais, a comparação é realizada em clusters

# identificar grupos de probes (genes duplicados)
genes <- colnames(norm_corrigido_GSE65635)
grupos <- split(seq_along(genes), genes)

idx_final <- unlist(lapply(grupos, function(i) {
  
  # caso não haja duplicata
  if (length(i) == 1) return(i)
  
  submat <- norm_corrigido_GSE65635[, i, drop = FALSE]
  
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
norm_corrigido_GSE65635 <- norm_corrigido_GSE65635[, idx_final]

# ====== Salvar arquivo ======

# confere se o pData e matriz estão realmente alinhados
all(rownames(norm_corrigido_GSE65635) == rownames(metadata))

#confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(processed_dir, id_projeto)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# arquivo de matriz de expressão
saveRDS(norm_corrigido_GSE65635,
        file = file.path(processed_dir,
                         id_projeto,
                         "exprs_GSE65635.rds"))
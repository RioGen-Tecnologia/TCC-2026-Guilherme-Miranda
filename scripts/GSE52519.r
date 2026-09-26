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

# ============== CARREGANDO PACOTES ==============

message("===========================================================================")
message("CARREGANDO PACOTES")
message("===========================================================================")

library(here)
source(here("scripts", "setup.r"))
library("hgu133plus2cdf")


# Carregar o modelo treinado
modelo_enet <- readRDS(here("results","machine_learning","modelo_elastic_net.rds"))
desempenho <- read.csv(here("results","machine_learning","desempenho_elastic_net.csv"))
genes_modelo <- read.csv(here("results","machine_learning","genes_selecionados_lambda_min_com_limma.csv"))

dir <- here("data", "temp")
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

# ============== GSE166716 ==============

message("===========================================================================")
message("BAIXANDO GSE166716")
message("===========================================================================")

id_projeto <- "GSE166716"

# baixando dados
getGEOSuppFiles(id_projeto, baseDir = dir)

# definindo pasta do projeto
projeto_dir <- file.path(dir, id_projeto)

# Descompactando
arquivos_tar <- list.files(
  path = projeto_dir,
  pattern = "\\.tar$",
  full.names = TRUE
)
untar(arquivos_tar, exdir = projeto_dir)

message("===========================================================================")
message("ORGANIZANDO GSE166716")
message("===========================================================================")

# Cria uma variável com arquivos que possuem "CEL.gz" e os imprime (exclue os CHP.gz)
cels.GSE166716 <- list.files(
  path = projeto_dir,
  pattern = "[cC][eE][lL]\\.gz$",
  full.names = TRUE
)

# Lendo os dados brutos (raw data)
dados_brutos_GSE166716 <- ReadAffy(filenames=cels.GSE166716)

# removendo arquivos
file.remove(list.files(projeto_dir, full.names = TRUE))

message("===========================================================================")
message("REALIZANDO QC DE GSE166716")
message("===========================================================================")

# Ajuste do modelo probe-level
plm_GSE166716 <- affyPLM::fitPLM(dados_brutos_GSE166716)

# RLE
Mbox(
  plm_GSE166716,
  main = paste("RLE -", id_projeto),
  las = 2
)

# NUSE
affyPLM::NUSE(
  plm_GSE166716,
  main = paste("NUSE -", id_projeto),
  las = 2
)

# DISTRIBUIÇÃO DAS INTENSIDADES BRUTAS
boxplot(
  dados_brutos_GSE166716,
  main = paste("Distribuição das intensidades brutas -", id_projeto),
  ylab = "Intensidade de fluorescência",
  las = 2,
  outline = FALSE
)

message("===========================================================================")
message("NORMALIZANDO GSE166716")
message("===========================================================================")

# normalização
dados_norm_GSE166716 <- affy::rma(dados_brutos_GSE166716)

# extração de matriz de expressão normalizada
norm_GSE166716 <- exprs(dados_norm_GSE166716)

# remove o ".CEL.gz" ou outras interferências da tabela
colnames(norm_GSE166716) <- toupper(sub(".*(GSM[0-9]+).*", "\\1", colnames(norm_GSE166716)))

# removem genes que apresentem NA
norm_corrigido_GSE166716 <- norm_GSE166716[!is.na(rownames(norm_GSE166716)), ]
norm_corrigido_GSE166716 <- norm_corrigido_GSE166716[rowSums(is.na(norm_corrigido_GSE166716)) == 0, ]
norm_corrigido_GSE166716 <- t(norm_corrigido_GSE166716)
norm_corrigido_GSE166716 <- norm_corrigido_GSE166716[,!grepl("^AFFX", colnames(norm_corrigido_GSE166716))] #remove probes de controle

message("===========================================================================")
message("ANOTANDO GSE166716")
message("===========================================================================")

# --- Etapa de anotação ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))


entrez_ids <- mapIds(hgu133plus2.db,
                     keys = colnames(norm_corrigido_GSE166716),
                     column = "ENTREZID",
                     keytype = "PROBEID",
                     multiVals = "first")

colnames(norm_corrigido_GSE166716) <- entrez_ids

# Remove genes NA
norm_corrigido_GSE166716 <- norm_corrigido_GSE166716[, !is.na(colnames(norm_corrigido_GSE166716))]


message("===========================================================================")
message("REMOVENDO DUPLICATAS DE GSE166716")
message("===========================================================================")

# Caso haja genes duplicados, é feito uma correlação de pearson.
# se a correlação for alta (0,7), é selecionado o primeiro probe.
# se a correlação for baixa, é selecionado o probe com maior variância
# se existirem probes triplicados ou mais, a comparação é realizada em clusters

# identificar grupos de probes (genes duplicados)
genes <- colnames(norm_corrigido_GSE166716)
grupos <- split(seq_along(genes), genes)

idx_final <- unlist(lapply(grupos, function(i) {
  
  # caso não haja duplicata
  if (length(i) == 1) return(i)
  
  submat <- norm_corrigido_GSE166716[, i, drop = FALSE]
  
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
norm_corrigido_GSE166716 <- norm_corrigido_GSE166716[, idx_final]

# ============== METADADOS ==============

# obtendo metadados
gse <- getGEO(id_projeto, GSEMatrix = TRUE)
pData <- pData(gse[[1]])
grupo_gse <- ifelse(grepl("normal urothelium", pData$`tissue:ch1`, ignore.case = TRUE), "non_tumor", "tumor")


# ============== VALIDAÇÃO TODOS GENES ==============
library(here)
library(pROC)
library(caret)
library(glmnet)

# Carregar o modelo e os arquivos salvos anteriormente (caso já não estejam no ambiente)
modelo_enet <- readRDS(here("results", "machine_learning", "modelo_elastic_net.rds"))
desempenho <- read.csv(here("results", "machine_learning", "desempenho_elastic_net.csv"))

# O objeto 'norm_corrigido_GSE166716' já deve estar carregado no seu ambiente R.


## 1. PREPARAÇÃO DA MATRIZ DE VALIDAÇÃO

# A matriz 'norm_corrigido_GSE166716' já possui linhas = amostras e colunas = genes (Entrez IDs).
X_gse <- as.matrix(norm_corrigido_GSE166716)


## 2. ALINHAMENTO DOS GENES (COLUNAS)

# Recuperar exatamente os nomes dos genes (EntrezIDs) que o modelo espera
genes_modelo_nomes <- rownames(coef(modelo_enet, s = "lambda.min"))
genes_modelo_nomes <- genes_modelo_nomes[genes_modelo_nomes != "(Intercept)"]

# Verificar se todos os genes do modelo estão presentes na base de validação
cat("Genes do modelo encontrados no GSE166716:", sum(genes_modelo_nomes %in% colnames(X_gse)), "de", length(genes_modelo_nomes), "\n")

# Garantir que a matriz de validação contenha exatamente as mesmas colunas (genes) 
# na mesma ordem exigida pelo modelo
X_gse_alinhado <- X_gse[, genes_modelo_nomes, drop = FALSE]


## 3. DEFINIÇÃO DO DESFECHO (Y_GSE)
# (Certifique-se de que o vetor 'grupo_gse' ou 'y_gse' foi recriado a partir do pheno_gse)
# Exemplo baseado no seu mapeamento anterior:
# grupo_gse <- ifelse(grepl("normal urothelium", pheno_gse$`tissue:ch1`, ignore.case = TRUE), "non_tumor", "tumor")
y_gse <- ifelse(grupo_gse == "tumor", 1, 0)


## 4. PREDIÇÃO DAS PROBABILIDADES
prob_gse <- predict(
  object = modelo_enet, 
  newx = X_gse_alinhado, 
  s = "lambda.min", 
  type = "response"
)

# Converter para vetor numérico
prob_gse <- as.numeric(prob_gse)

# Verificar o range das probabilidades (agora elas não devem estar mais travadas em 1.0)
print("Range das probabilidades previstas:")
print(range(prob_gse))


## 5. AVALIAÇÃO POR CURVA ROC E AUC
roc_gse <- roc(
  response = y_gse, 
  predictor = prob_gse, 
  direction = "<"
)

auc_gse <- auc(roc_gse)
cat("AUC no GSE166716 após pré-processamento correto:", auc_gse, "\n")

# Plotar a curva ROC externa
plot(roc_gse, print.auc = TRUE, main = "Curva ROC — Validação Externa (GSE166716)")


## 6. CLASSIFICAÇÃO E MATRIZ DE CONFUSÃO

# Recuperar o threshold ótimo salvo no arquivo de desempenho do treinamento
threshold <- desempenho$Threshold

# Classificar as amostras usando o threshold
classe_gse_predita <- ifelse(prob_gse >= threshold, 1, 0)

# Gerar Matriz de Confusão
matriz_conf_gse <- confusionMatrix(
  factor(classe_gse_predita, levels = c(0, 1)), 
  factor(y_gse, levels = c(0, 1)), 
  positive = "1"
)

print(matriz_conf_gse)



# 1. Recalcular o melhor threshold para o GSE166716 usando o índice de Youden
ponto_corte_gse <- coords(
  roc_gse, 
  x = "best", 
  best.method = "youden", 
  ret = c("threshold", "sensitivity", "specificity")
)

threshold_gse <- as.numeric(ponto_corte_gse["threshold"])
cat("Novo threshold otimizado para o GSE166716:", threshold_gse, "\n")

# 2. Classificar as amostras usando este novo threshold específico da validação
classe_gse_predita_otimizada <- ifelse(prob_gse >= threshold_gse, 1, 0)

# 3. Gerar a nova Matriz de Confusão
matriz_conf_gse_otimizada <- confusionMatrix(
  factor(classe_gse_predita_otimizada, levels = c(0, 1)), 
  factor(y_gse, levels = c(0, 1)), 
  positive = "1"
)

print(matriz_conf_gse_otimizada)

# Análise diferencial com dados obtidos do TCGA (tumoral) e GTex (não-tumoral)
# como validação de análise principal
# Guilherme Moret Miranda - Riogen
# 11/05/2026
# padronização do Recount3: monorail
# análise estatística: limma (voom)

message("\n", paste(rep("=", 30), collapse = ""))
message("Baixando dados do recount3")
message(paste(rep("=", 30), collapse = ""))

library(recount3)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(pROC)

# Listar projetos e filtrar TCGA (Bexiga) e GTEx (Bexiga)
projects <- available_projects()

proj_tcga <- subset(projects, project == "BLCA" & project_home == "data_sources/tcga")
proj_gtex <- subset(projects, project == "BLADDER" & project_home == "data_sources/gtex")

# apaga os projetos para liberar memória
rm(projects)
gc()

# Criar objetos RSE e converter para contagens brutas
rse_tcga <- create_rse(proj_tcga)
rse_gtex <- create_rse(proj_gtex)

message("\n", paste(rep("=", 30), collapse = ""))
message("Convertendo contagem bruta")
message(paste(rep("=", 30), collapse = ""))

assay(rse_tcga, "counts") <- transform_counts(rse_tcga)
assay(rse_gtex, "counts") <- transform_counts(rse_gtex)

# --- 2. HARMONIZAÇÃO DOS METADADOS ---
# Extrair tipos de amostra do TCGA usando a coluna correta identificada
# --- 2. HARMONIZAÇÃO DOS METADADOS ---

message("\n", paste(rep("=", 30), collapse = ""))
message("Lendo metadados")
message(paste(rep("=", 30), collapse = ""))

# Extrair metadados originais
tcga_coldata <- colData(rse_tcga)

# Tipo de amostra
tcga_types <- tcga_coldata$tcga.gdc_cases.samples.sample_type

# Estadiamento T
tcga_stage <- tcga_coldata$tcga.cgc_case_pathologic_t

# Limpar TX (opcional, recomendado)
tcga_stage[tcga_stage == "TX"] <- NA

# Criar metadados TCGA mantendo estágio
meta_tcga <- data.frame(
  sample_id = colnames(rse_tcga),
  condition = ifelse(tcga_types == "Primary Tumor", "Tumor", "Normal_Adj"),
  study = "TCGA",
  T_stage = tcga_stage,
  row.names = colnames(rse_tcga)
)

# GTEx (sem estágio)
meta_gtex <- data.frame(
  sample_id = colnames(rse_gtex),
  condition = "Healthy",
  study = "GTEx",
  T_stage = NA,
  row.names = colnames(rse_gtex)
)

# Substituir colData
colData(rse_tcga) <- as(meta_tcga, "DFrame")
colData(rse_gtex) <- as(meta_gtex, "DFrame")

# Combinar
rse_combined <- cbind(rse_tcga, rse_gtex)

message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando de ENSEMBL para ENTREZ")
message(paste(rep("=", 30), collapse = ""))

# --- 3. MAPEAMENTO DE IDs (ENSEMBL PARA ENTREZ) ---
# Limpar versões dos IDs (.1, .2, etc)
ensembl_ids_clean <- gsub("\\..*", "", rownames(rse_combined))

# Mapear IDs
mapping <- mapIds(org.Hs.eg.db, keys = ensembl_ids_clean, 
                  column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

message("\n", paste(rep("=", 30), collapse = ""))
message("Removendo valores nulos e duplicatas")
message(paste(rep("=", 30), collapse = ""))

# Filtrar NAs e resolver duplicatas (manter o gene com maior soma de contagens)
rse_filtered <- rse_combined[!is.na(mapping), ]
new_entrez <- mapping[!is.na(mapping)]

all_counts <- assay(rse_filtered, "counts")
row_sums <- rowSums(all_counts)

# Ordenar por Entrez e depois por expressão para remover duplicatas
best_idx <- order(new_entrez, row_sums, decreasing = TRUE)
rse_final <- rse_filtered[best_idx, ]
final_entrez <- new_entrez[best_idx]

rse_final <- rse_final[!duplicated(final_entrez), ]
rownames(rse_final) <- final_entrez[!duplicated(final_entrez)]

# --- 4. ANÁLISE ESTATÍSTICA (LIMMA-VOOM) ---

message("\n", paste(rep("=", 30), collapse = ""))
message("Realizando análise estatística (VOOM)")
message(paste(rep("=", 30), collapse = ""))

# Criar objeto DGEList e filtrar genes pouco expressos
dge <- DGEList(counts = assay(rse_final, "counts"), samples = colData(rse_final))
keep_exprs <- filterByExpr(dge, group = dge$samples$condition)
dge <- dge[keep_exprs, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Criar Modelo de Grupos (para evitar colinearidade entre Study e Condition)
dge$samples$group <- factor(paste(dge$samples$study, dge$samples$condition, sep = "_"))
design <- model.matrix(~ 0 + group, data = dge$samples)
colnames(design) <- levels(dge$samples$group)

# Transformação Voom
v <- voom(dge, design)

# Ajuste do Modelo Linear e Contrastes
fit <- lmFit(v, design)

cont.matrix <- makeContrasts(
  TumorVsHealthy = TCGA_Tumor - GTEx_Healthy,
  TumorVsNormalAdj = TCGA_Tumor - TCGA_Normal_Adj,
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

## arrumando dados de expressão para análise ROC/AUC
# criando anotação se amostra é tumoral ou não (amostras saudáveis GTex e adjacentes TCGA )
group <- dge$samples$condition
group_roc <- ifelse(group == "Tumor",
                    "tumor",
                    "non_tumor")

group_roc <- factor(group_roc, levels = c("non_tumor", "tumor"))

#conferindo se está alinhado
if (!all(colnames(v$E) == rownames(dge$samples))) {
  stop("Erro: Anotação não está alinhada!")
}

# --- 5. EXTRAÇÃO DE RESULTADOS PARA META-ANÁLISE ---
results <- topTable(fit2,
                    coef = "TumorVsHealthy",
                    number = Inf,
                    adjust.method = "BH",
                    sort.by = "none")

# --- PCA DOS DADOS DE EXPRESSÃO ---

message("\n", paste(rep("=", 30), collapse = ""))
message("Gerando PCA")
message(paste(rep("=", 30), collapse = ""))

# PCA sobre a matriz voom-normalizada
pca <- prcomp(
  t(v$E),
  center = TRUE,
  scale. = FALSE
)

# Variância explicada
percent_var <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

# Metadados
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = dge$samples$condition,
  study = dge$samples$study
)

# Criar rótulo combinando estudo e condição
pca_df$group <- factor(
  paste(pca_df$study, pca_df$condition, sep = "_")
)

# --- SALVAR FIGURA ---

png(
  filename = file.path(
    figures_dir,
    "PCA_Recount3.png"
  ),
  width = 1800,
  height = 1600,
  res = 250
)

plot(
  pca_df$PC1,
  pca_df$PC2,
  pch = 19,
  cex = 1.2,
  xlab = paste0("PC1 (", round(percent_var[1], 1), "%)"),
  ylab = paste0("PC2 (", round(percent_var[2], 1), "%)"),
  main = "PCA - TCGA vs GTEx",
  col = as.numeric(pca_df$group)
)

legend(
  "topright",
  legend = levels(pca_df$group),
  col = seq_along(levels(pca_df$group)),
  pch = 19,
  bty = "n"
)

dev.off()


message("\n", paste(rep("=", 30), collapse = ""))
message("Análise finalizada")
message(paste(rep("=", 30), collapse = ""))

# Análise integrada de projetos de microarray obtidas do Gene Expression Omnibus
# Guilherme Moret Miranda - Riogen
# 26/09/2026
# Esse script obtém, normaliza e análisa cada um dos projetos buscados previamente,
# combina os dados de cada projeto por comBat e realiza posteriores análises
# para determinar biomarcadores gênicos de câncer de bexiga.

# ============== CARREGANDO PACOTES ==============

message("===========================================================================")
message("CARREGANDO PACOTES")
message("===========================================================================")

library(here)
source(here("scripts", "setup.r"))

# ============== DEFINIÇÃO DE DIRETÓRIO E ARQUIVOS ==============
# Definir a pasta de trabalho e o master manifesto, com anotações das amostras

geo_dir      <- here("data", "raw", "GEO")
processed_dir <- here("data", "processed")
results_dir  <- here("results")
figures_dir  <- here("figures")
scripts_dir  <- here("scripts")
metadata_dir <- here("metadata")

dirs <- c(geo_dir,processed_dir,results_dir,figures_dir,scripts_dir)

# cria os diretórios se eles não existirem
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# Master Manifesto com anotações das amostras
metadata_path <- list.files(
  here("metadata"),
  pattern = "bladder_cancer_metadata_amostras_filtradas_QC_por_dataset.csv",
  full.names = TRUE
)

metadata <- read.csv(metadata_path)

# limpando
rm(d,dirs)
gc()


# ============== ANÁLISE DE PROJETOS ==============
# Esta sessão executa os scripts secundários que extrem os dados do GEO,
# normalizam, anotam os genes em entrezID e produz matriz de expressão

message("===========================================================================")
message("ANALISANDO DATASETS")
message("===========================================================================")

# Lista scripts que começam com "GSE"
scripts_projetos <- list.files(
  path = scripts_dir,
  pattern = "^GSE.*\\.[Rr]$",
  full.names = TRUE
)

# Loop de execução dos scripts secundários
for (script in scripts_projetos) {
  message("Processando: ", basename(script))
  env <- new.env()
  source(script, local = env)
  rm(env)
  gc()
}

# Lista das pastas GSE
gse_dirs <- list.dirs(processed_dir, recursive = FALSE)

# Nome dos datasets
gse_names <- basename(gse_dirs)

message("================================")
message("CARREGANDO DATASETS")
message("================================")

## Carregar todas as matrizes .rds
matrix_list <- lapply(gse_dirs, function(dir) {
  rds_file <- list.files(
    dir,
    pattern = "exprs_.*\\.rds$",
    full.names = TRUE
  )
  readRDS(rds_file)
})

# Nomear os elementos da lista
names(matrix_list) <- gse_names


#limpando
rm(arquivo,scripts_projetos,script,id)
gc()


# ============== INTEGRAÇÃO COMBAT ==============
# Integração das matrizes de expressão em uma super matriz através do pacote
# sva pelo método comBat, com os batches sendo os projetos individuais e as variações
# biológicas o estado de tumor e não-tumor

# genes comuns entre todos os estudos
genes_comuns <- Reduce(
  intersect,
  lapply(matrix_list, colnames)
)

message(
  "Número de genes comuns: ",
  length(genes_comuns)
)

# padronizando genes
matrix_list_common <- lapply(
  matrix_list,
  function(x) x[, genes_comuns, drop = FALSE]
)

message("================================")
message("INTEGRANDO DATASETS")
message("================================")

## integrando amostras
expr <- do.call(
  rbind,
  matrix_list_common
)

# manter somente amostras presentes no metadata
common_samples <- intersect(
  rownames(expr),
  metadata$sample_ID
)

expr <- expr[common_samples, ]

metadata <- metadata[
  match(common_samples, metadata$sample_ID),
]

## CORREÇÃO DE BATCHES

message("================================")
message("CORRIGINDO BATCH EFFECT")
message("================================")

# batch = estudo GEO
batch <- metadata$study_ID

# grupos biológicos
group <- metadata$sample_type

# características clínicas
characteristics <- metadata$characteristics

stopifnot(
  all(rownames(expr) == metadata$sample_ID)
)

# variável para visualização
sample_class <- characteristics

sample_class <- factor(
  sample_class,
  levels = c(
    "MIBC",
    "NMIBC",
    "non_cancer_individual"
  )
)


# modelo preservando grupos biológicos
mod <- model.matrix(~ group)


# ComBat
expr_combat <- ComBat(
  dat   = t(expr),
  batch = batch,
  mod   = mod
)

## Gráficos PCA

message("===========================================================================")
message("GERANDO E EXPORTANDO GRÁFICOS DE PCA (BEFORE VS AFTER COMBAT)")
message("===========================================================================")


# cálculo e preparação dos dados de pca

# A) Antes do ComBat (matriz expr: amostras nas linhas, genes nas colunas)
# Usamos prcomp diretamente em expr (genes com variação > 0)
expr_var <- expr[, apply(expr, 2, var) > 0]
pca_before <- prcomp(expr_var, scale. = TRUE)
var_explained_before <- round(100 * (pca_before$sdev^2 / sum(pca_before$sdev^2)), 1)

df_pca_before <- data.frame(
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  PC3 = pca_before$x[, 3],
  Study = metadata$study_ID,
  Group = metadata$sample_type,
  Class = sample_class
)

# B) Depois do ComBat (expr_combat: genes nas linhas, amostras nas colunas)
# Transpomos para amostras nas linhas
expr_combat_t <- t(expr_combat)
expr_combat_var <- expr_combat_t[, apply(expr_combat_t, 2, var) > 0]
pca_after <- prcomp(expr_combat_var, scale. = TRUE)
var_explained_after <- round(100 * (pca_after$sdev^2 / sum(pca_after$sdev^2)), 1)

df_pca_after <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  PC3 = pca_after$x[, 3],
  Study = metadata$study_ID,
  Group = metadata$sample_type,
  Class = sample_class
)


# Painés 2D (COLORIDOS POR ESTUDO E POR GRUPO BIOLÓGICO) ---

# A) Antes do ComBat - Por Estudo (Batch)
p1 <- ggplot(df_pca_before, aes(x = PC1, y = PC2, color = Study, shape = Group)) +
  geom_point(alpha = 0.7, size = 2.5) +
  labs(
    title = "Antes do ComBat (por Estudo)",
    x = paste0("PC1 (", var_explained_before[1], "%)"),
    y = paste0("PC2 (", var_explained_before[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.position = "right"
  )

# B) Depois do ComBat - Por Estudo (Batch)
p2 <- ggplot(df_pca_after, aes(x = PC1, y = PC2, color = Study, shape = Group)) +
  geom_point(alpha = 0.7, size = 2.5) +
  labs(
    title = "Depois do ComBat (por Estudo)",
    x = paste0("PC1 (", var_explained_after[1], "%)"),
    y = paste0("PC2 (", var_explained_after[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.position = "right"
  )

# C) Antes do ComBat - Por Classe Clínica / Grupo
p3 <- ggplot(df_pca_before, aes(x = PC1, y = PC2, color = Group, shape = Class)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("tumor" = "#D95F02", "non_tumor" = "#1B9E77")) +
  labs(
    title = "Antes do ComBat (por Grupo Biológico)",
    x = paste0("PC1 (", var_explained_before[1], "%)"),
    y = paste0("PC2 (", var_explained_before[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.position = "right"
  )

# D) Depois do ComBat - Por Classe Clínica / Grupo
p4 <- ggplot(df_pca_after, aes(x = PC1, y = PC2, color = Group, shape = Class)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("tumor" = "#D95F02", "non_tumor" = "#1B9E77")) +
  labs(
    title = "Depois do ComBat (por Grupo Biológico)",
    x = paste0("PC1 (", var_explained_after[1], "%)"),
    y = paste0("PC2 (", var_explained_after[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.position = "right"
  )


# montagem e exportação dos painéis lado a lado

# Painel 1: Avaliação do Efeito de Lote (Batch)
pca_batch_panel <- (p1 | p2) +
  plot_annotation(
    title    = "Remoção do Efeito de Lote via ComBat",
    subtitle = "Comparação da distribuição das amostras por Estudo GEO (Batch)",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  )

# Painel 2: Preservação da Biologia (Tumor vs Não-Tumor)
pca_biology_panel <- (p3 | p4) +
  plot_annotation(
    title    = "Preservação da Variabilidade Biológica",
    subtitle = "Separação dos grupos Tumor vs Não-Tumor antes e depois da correção",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  )


# Salvar imagens em alta resolução
ggsave(
  filename = file.path(figures_dir, "pca_batch_effect_combat_QC.png"),
  plot     = pca_batch_panel,
  width    = 14,
  height   = 6,
  dpi      = 300
)

ggsave(
  filename = file.path(figures_dir, "pca_biology_combat_QC.png"),
  plot     = pca_biology_panel,
  width    = 14,
  height   = 6,
  dpi      = 300
)


# gráficos pca 3d interativos (antes e depois do combat)


# Garantir ordem dos grupos para definir os símbolos
df_pca_before$Group <- factor(
  df_pca_before$Group,
  levels = c("non_tumor", "tumor")
)

df_pca_after$Group <- factor(
  df_pca_after$Group,
  levels = c("non_tumor", "tumor")
)


# A) PCA 3D Antes do ComBat
p_3d_before <- plot_ly(
  df_pca_before,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~Study,
  symbol = ~Group,
  symbols = c("circle", "diamond"),
  text = ~paste(
    "Amostra:", metadata$sample_ID,
    "<br>Estudo:", Study,
    "<br>Grupo:", Group,
    "<br>Classe:", Class
  ),
  hoverinfo = "text",
  marker = list(size = 5)
) %>%
  add_markers() %>%
  layout(
    title = "PCA 3D Antes do ComBat (Efeito de Lote Nítido)",
    scene = list(
      xaxis = list(
        title = paste0("PC1 (", var_explained_before[1], "%)")
      ),
      yaxis = list(
        title = paste0("PC2 (", var_explained_before[2], "%)")
      ),
      zaxis = list(
        title = paste0("PC3 (", var_explained_before[3], "%)")
      )
    )
  )


# B) PCA 3D Depois do ComBat
p_3d_after <- plot_ly(
  df_pca_after,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~Study,
  symbol = ~Group,
  symbols = c("circle", "diamond"),
  text = ~paste(
    "Amostra:", metadata$sample_ID,
    "<br>Estudo:", Study,
    "<br>Grupo:", Group,
    "<br>Classe:", Class
  ),
  hoverinfo = "text",
  marker = list(size = 5)
) %>%
  add_markers() %>%
  layout(
    title = "PCA 3D Pós-ComBat (Remoção do Efeito de Lote)",
    scene = list(
      xaxis = list(
        title = paste0("PC1 (", var_explained_after[1], "%)")
      ),
      yaxis = list(
        title = paste0("PC2 (", var_explained_after[2], "%)")
      ),
      zaxis = list(
        title = paste0("PC3 (", var_explained_after[3], "%)")
      )
    )
  )


# Salvar ambos os HTMLs interativos
saveWidget(p_3d_before, file = file.path(figures_dir, "pca_3d_pre_combat_QC.html"), selfcontained = TRUE)
saveWidget(p_3d_after,  file = file.path(figures_dir, "pca_3d_post_combat_QC.html"), selfcontained = TRUE)

# Remover pastas temporárias geradas pelos widgets
unlink(
  file.path(figures_dir, c(
    "pca_3d_pre_combat_QC_files",
    "pca_3d_post_combat_QC_files"
  )),
  recursive = TRUE,
  force = TRUE
)

# Limpeza de memória
rm(
  expr_var, pca_before, var_explained_before, df_pca_before,
  expr_combat_t, expr_combat_var, pca_after, var_explained_after, df_pca_after,
  p1, p2, p3, p4, pca_batch_panel, pca_biology_panel, p_3d_before, p_3d_after
)
gc()


# ============== ANÁLISE ESTATÍSTICA LIMMA ==============
# Analisando a matriz integrada pelo método empirical bayes do limma para
# identificar genes diferencialmente expressos entre os genes

# Reordena o metadata para ficar exatamente na mesma ordem das colunas da matriz de expressão
metadata_alinhado <- metadata[
  match(colnames(expr_combat), metadata$sample_ID),
]

# definindo o grupo experimental
grupo <- factor(
  metadata_alinhado$sample_type,
  levels = c("non_tumor", "tumor")
)

# contruindo modelo
design <- model.matrix(~ grupo)

message("===========================================================================")
message("REALIZANDO ANÁLISE LIMMA")
message("===========================================================================")

# rodando limma
fit <- lmFit(expr_combat, design)
fit <- eBayes(fit)

# extraindo resultados
resultados_limma <- topTable(
  fit,
  coef = "grupotumor",
  number = Inf,
  adjust.method = "BH"
)

# Adicionar Entrez ID como coluna
resultados_limma$EntrezID <- rownames(resultados_limma)

# Converter Entrez ID para Gene Symbol
resultados_limma$GeneSymbol <- mapIds(
  illuminaHumanv2.db,
  keys = resultados_limma$EntrezID,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Genes diferencialmente expressos
DEGs <- resultados_limma[
  abs(resultados_limma$logFC) > 1 &
    resultados_limma$adj.P.Val < 0.05,
]

# Genes diferencialmente expressos up-regulados
DEGs_up <- resultados_limma[
  resultados_limma$logFC > 1 &
    resultados_limma$adj.P.Val < 0.05,
]

# Genes diferencialmente expressos down-regulados
DEGs_down <- resultados_limma[
  resultados_limma$logFC < -1 &
    resultados_limma$adj.P.Val < 0.05,
]


# Exportar os resultados em CSV
write.csv(resultados_limma, file.path(results_dir, "limma/limma_tabela_completa.csv"), row.names = FALSE)
write.csv(DEGs, file.path(results_dir, "limma/limma_DEGs.csv"), row.names = FALSE)
write.csv(DEGs_up, file.path(results_dir, "limma/limma_DEGs_Upregulated.csv"), row.names = FALSE)
write.csv(DEGs_down, file.path(results_dir, "limma/limma_DEGs_Downregulated.csv"), row.names = FALSE)


# ============== GRÁFICOS DE EXPRESSÃO ==============
# Gerando um gráfico volcano plot e heatmap com os genes diferencialmente expressos

message("===========================================================================")
message("GERANDO GRÁFICOS DE EXPRESSÃO")
message("===========================================================================")

## ==== VOLCANO PLOT ====

# Adicionar coluna de status biológico para o plot
resultados_limma$diffexpression <- "Não significativo"
resultados_limma$diffexpression[resultados_limma$logFC > 1 & resultados_limma$adj.P.Val < 0.05] <- "Upregulated"
resultados_limma$diffexpression[resultados_limma$logFC < -1 & resultados_limma$adj.P.Val < 0.05] <- "Downregulated"

resultados_limma$diffexpression <- factor(resultados_limma$diffexpression, 
                                          levels = c("Downregulated", "Não significativo", "Upregulated"))

# Selecionar os 15 genes mais significativos por p-valor ajustado para rotular
top_genes <- head(resultados_limma[order(resultados_limma$adj.P.Val), ], 15)

# Construir o Volcano Plot (usando Gene_Symbol nos rótulos)
volcano_plot <- ggplot(resultados_limma, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpression)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(
    values = c("Downregulated" = "#377EB8", "Não significativo" = "grey70", "Upregulated" = "#E41A1C"),
    name = "Expressão Diferencial"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.7) +
  geom_text_repel(
    data = top_genes,
    aes(label = GeneSymbol),
    size = 3.5,
    max.overlaps = 15,
    show.legend = FALSE
  ) +
  labs(
    title = "Volcano Plot: Tumor vs Non-Tumor",
    subtitle = "FDR < 0.05 & |logFC| > 1",
    x = expression("log"[2]*" Fold Change"),
    y = expression("-log"[10]*" (FDR)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

# Salvar Volcano Plot em PNG
ggsave(
  filename = file.path(figures_dir, "volcano_plot.png"),
  plot = volcano_plot,
  width = 8,
  height = 7,
  dpi = 300
)


## ==== HEATMAP ====

# 1. Identificar qual coluna do DEGs contém os Entrez IDs
coluna_id <- intersect(c("Entrez_ID", "gene_id", "ID", "Entrez"), colnames(DEGs))
if (length(coluna_id) > 0) {
  degs_ids <- as.character(DEGs[[coluna_id[1]]])
} else {
  degs_ids <- as.character(rownames(DEGs))
}

# 2. Ordenar DEGs pelo valor de adj.P.Val
ordem_deg <- order(DEGs$adj.P.Val)
degs_ids_ordenados <- degs_ids[ordem_deg]
degs_ordenados_df <- DEGs[ordem_deg, ]

# 3. Filtrar apenas os genes que realmente existem nas linhas de expr_combat
genes_presentes <- intersect(degs_ids_ordenados, rownames(expr_combat))
top50_entrez <- head(genes_presentes, 50)

if (length(top50_entrez) < 2) {
  stop(sprintf("Apenas %d genes válidos foram encontrados na interseção.", length(top50_entrez)))
}

## ==== 4. EXTRAIR SUBMATRIZ E NORMALIZAR (Z-SCORE) ====
matriz_heatmap <- expr_combat[top50_entrez, , drop = FALSE]
matriz_heatmap_scaled <- t(scale(t(matriz_heatmap)))
matriz_heatmap_scaled <- matriz_heatmap_scaled[complete.cases(matriz_heatmap_scaled), , drop = FALSE]


## ==== 5. MAPEAR GENE SYMBOLS CORRESPONDENTES ====
coluna_symbol <- intersect(c("GeneSymbol", "symbol", "SYMBOL", "Gene_Symbol"), colnames(DEGs))
if (length(coluna_symbol) > 0) {
  symbols_top50 <- degs_ordenados_df[[coluna_symbol[1]]][match(rownames(matriz_heatmap_scaled), degs_ids_ordenados)]
} else {
  symbols_top50 <- rownames(matriz_heatmap_scaled)
}

symbols_top50[is.na(symbols_top50) | symbols_top50 == ""] <- rownames(matriz_heatmap_scaled)[is.na(symbols_top50) | symbols_top50 == ""]
rownames(matriz_heatmap_scaled) <- make.unique(as.character(symbols_top50))


## ==== ORDENAR E PREPARAR AMOSTRAS ====

# Detectar automaticamente o formato dos nomes das colunas na matriz
if (any(metadata_alinhado$sample_ID %in% colnames(matriz_heatmap_scaled))) {
  # As colunas usam apenas o Sample ID (ex: GSM123456)
  metadata_alinhado$matrix_id <- metadata_alinhado$sample_ID
} else {
  # As colunas usam StudyID_SampleID (ex: GSE13507_GSM123456)
  metadata_alinhado$matrix_id <- paste0(metadata_alinhado$study_ID, "_", metadata_alinhado$sample_ID)
}

# Ordenar o metadado pelo tipo biológico para manter os grupos juntos
ordem_amostras <- order(metadata_alinhado$sample_type)
metadata_alinhado_ordenado <- metadata_alinhado[ordem_amostras, ]

# Identificar as amostras comuns
amostras_comuns <- intersect(metadata_alinhado_ordenado$matrix_id, colnames(matriz_heatmap_scaled))

if (length(amostras_comuns) == 0) {
  # Imprime uma amostra dos nomes para ajudar na depuração se ainda falhar
  message("Colunas na matriz (head): ", paste(head(colnames(matriz_heatmap_scaled)), collapse=", "))
  message("IDs gerados no metadado (head): ", paste(head(metadata_alinhado$matrix_id), collapse=", "))
  stop("Nenhuma amostra em comum encontrada. Verifique as mensagens acima para alinhar os identificadores.")
}

# Filtrar o metadado mantendo a ordenação biológica
metadata_alinhado_ordenado <- metadata_alinhado_ordenado[metadata_alinhado_ordenado$matrix_id %in% amostras_comuns, ]

# Reordenar as colunas da matriz forçando a ordem do metadado
matriz_heatmap_scaled <- matriz_heatmap_scaled[, metadata_alinhado_ordenado$matrix_id, drop = FALSE]

# Atribuir nomes únicos definitivos
colnames(matriz_heatmap_scaled) <- make.unique(metadata_alinhado_ordenado$matrix_id)


## ==== ANOTAÇÕES E GERAÇÃO DO HEATMAP ====

# Criar dataframe de anotação para as colunas do heatmap
anotacoes_heatmap <- data.frame(
  Tipo = metadata_alinhado_ordenado$sample_type,
  Estudo = metadata_alinhado_ordenado$study_ID,
  row.names = colnames(matriz_heatmap_scaled)
)

heatmap_plot <- pheatmap(
  matriz_heatmap_scaled,
  annotation_col = anotacoes_heatmap,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  cluster_cols = FALSE, # Impede a reordenação, mantendo a separação por Tipo
  cluster_rows = TRUE,
  main = "Heatmap: Top 50 Genes Diferencialmente Expressos (Z-Score)"
)

ggsave(
  filename = file.path(figures_dir, "heatmap_top50_DEGs.png"),
  plot = heatmap_plot,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

# limpeza
rm(
  volcano_plot, top_genes, 
  coluna_id, degs_ids, ordem_deg, degs_ids_ordenados, degs_ordenados_df,
  genes_presentes, top50_entrez, matriz_heatmap, 
  coluna_symbol, symbols_top50, ordem_amostras, metadata_alinhado_ordenado,
  amostras_comuns, anotacoes_heatmap, heatmap_filepath
)
gc()

# ============== ANÁLISE DE ENRIQUECIMENTO FUNCIONAL ==============
# GO Biological Process, KEGG e Reactome

message("\n", paste(rep("=", 50), collapse = ""))
message("Análise de enriquecimento funcional")
message(paste(rep("=", 50), collapse = ""))

# --- Preparação dos conjuntos de genes ---

genes_up <- unique(as.character(DEGs_up$EntrezID))
genes_down <- unique(as.character(DEGs_down$EntrezID))

# Universo = todos os genes efetivamente testados pelo limma
universe <- unique(as.character(resultados_limma$EntrezID))

genes_up <- intersect(genes_up[!is.na(genes_up)], universe)
genes_down <- intersect(genes_down[!is.na(genes_down)], universe)
universe <- universe[!is.na(universe)]

message("Genes upregulated: ", length(genes_up))
message("Genes downregulated: ", length(genes_down))
message("Genes no universo: ", length(universe))


# --- Função para realizar o enriquecimento ---

realizar_enriquecimento <- function(gene_list, universe, nome_grupo) {
  
  message("\nGrupo: ", nome_grupo)
  
  # GO Biological Process
  ora_go <- clusterProfiler::enrichGO(
    gene = gene_list, universe = universe, OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
  )
  
  # KEGG
  ora_kegg <- clusterProfiler::enrichKEGG(
    gene = gene_list, universe = universe, organism = "hsa",
    pvalueCutoff = 0.05, pAdjustMethod = "BH"
  )
  
  ora_kegg <- clusterProfiler::setReadable(
    ora_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"
  )
  
  # Reactome
  ora_reactome <- ReactomePA::enrichPathway(
    gene = gene_list, universe = universe, organism = "human",
    pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE
  )
  
  list(GO = ora_go, KEGG = ora_kegg, Reactome = ora_reactome)
}


# --- Executar para genes upregulated e downregulated ---

enrichment_up <- realizar_enriquecimento(genes_up, universe, "Upregulated")
enrichment_down <- realizar_enriquecimento(genes_down, universe, "Downregulated")


# --- Exportação dos resultados ---

enrichment_dir <- file.path(results_dir, "enrichment")
dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(as.data.frame(enrichment_up$GO), file.path(enrichment_dir, "ORA_GO_BP_upregulated.csv"), row.names = FALSE)
write.csv(as.data.frame(enrichment_up$KEGG), file.path(enrichment_dir, "ORA_KEGG_upregulated.csv"), row.names = FALSE)
write.csv(as.data.frame(enrichment_up$Reactome), file.path(enrichment_dir, "ORA_Reactome_upregulated.csv"), row.names = FALSE)

write.csv(as.data.frame(enrichment_down$GO), file.path(enrichment_dir, "ORA_GO_BP_downregulated.csv"), row.names = FALSE)
write.csv(as.data.frame(enrichment_down$KEGG), file.path(enrichment_dir, "ORA_KEGG_downregulated.csv"), row.names = FALSE)
write.csv(as.data.frame(enrichment_down$Reactome), file.path(enrichment_dir, "ORA_Reactome_downregulated.csv"), row.names = FALSE)


# --- Função para gerar dotplots ---

make_dotplot <- function(obj, titulo, n_cat = 10) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste(titulo, "\nNenhum termo enriquecido")))
  }
  
  enrichplot::dotplot(obj, showCategory = n_cat, title = titulo) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          axis.text.y = element_text(size = 8),
          legend.position = "right")
}


# --- Dotplots: genes upregulated ---

p_go_up <- make_dotplot(enrichment_up$GO, "GO Biological Process — Upregulated")
p_kegg_up <- make_dotplot(enrichment_up$KEGG, "KEGG — Upregulated")
p_reactome_up <- make_dotplot(enrichment_up$Reactome, "Reactome — Upregulated")

combined_up <- (p_go_up | p_kegg_up | p_reactome_up) +
  patchwork::plot_annotation(
    title = "Enriquecimento funcional — genes upregulated",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(figures_dir, "ora_dotplots_upregulated.png"),
       combined_up, width = 20, height = 7, dpi = 300)


# --- Dotplots: genes downregulated ---

p_go_down <- make_dotplot(enrichment_down$GO, "GO Biological Process — Downregulated")
p_kegg_down <- make_dotplot(enrichment_down$KEGG, "KEGG — Downregulated")
p_reactome_down <- make_dotplot(enrichment_down$Reactome, "Reactome — Downregulated")

combined_down <- (p_go_down | p_kegg_down | p_reactome_down) +
  patchwork::plot_annotation(
    title = "Enriquecimento funcional — genes downregulated",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(figures_dir, "ora_dotplots_downregulated.png"),
       combined_down, width = 20, height = 7, dpi = 300)


# --- Limpeza ---

rm(genes_up, genes_down, universe, enrichment_up, enrichment_down, p_go_up, p_kegg_up, p_reactome_up, combined_up,
   p_go_down, p_kegg_down, p_reactome_down, combined_down, make_dotplot, realizar_enriquecimento)

gc()


# ============== MACHINE LEARNING ELASTIC NET ==============
#
# Nesta etapa, será utilizado o método Elastic Net para construir um modelo de
# classificação capaz de distinguir amostras de tecido tumoral e não tumoral a
# partir dos dados de expressão gênica.
#
# Primeiramente, a matriz de expressão integrada pelo ComBat será dividida em
# conjuntos de treinamento (80%) e teste (20%). O modelo Elastic Net será ajustado
# somente com as amostras de treinamento, utilizando validação cruzada de 10
# folds para determinar o parâmetro de regularização mais adequado.
#
# Em seguida, os genes com coeficientes diferentes de zero serão identificados
# como genes selecionados pelo modelo. O desempenho preditivo será avaliado nas
# amostras de teste por meio da curva ROC, área sob a curva (AUC), sensibilidade,
# especificidade e matriz de confusão. Por fim, os genes selecionados pelo Elastic
# Net serão comparados aos genes diferencialmente expressos identificados pela
# análise de limma.
#
# IMPORTANTE:
# Esta etapa representa uma validação interna por divisão treino/teste. Como a
# correção de batch pelo ComBat foi realizada anteriormente sobre o conjunto
# integrado completo, os resultados preditivos não devem ser interpretados como
# uma validação externa independente.


## 1. PREPARAÇÃO DA MATRIZ DE EXPRESSÃO=======================

# O glmnet espera:
#   linhas  = amostras
#   colunas = variáveis (genes)
#
# A matriz expr_combat possui:
#   linhas  = genes
#   colunas = amostras
#
# Portanto, é necessário transpô-la.

X <- t(expr_combat)


## DEFINIÇÃO DO DESFECHO

# Transformar o tipo de amostra em uma variável binária:
#
#   0 = não-tumor
#   1 = tumor
#
# Essa variável será utilizada como desfecho pelo modelo
# de regressão logística.

y <- ifelse(
  metadata_alinhado$sample_type == "tumor",
  1,
  0
)

## 3. VERIFICAÇÃO DO ALINHAMENTO DAS AMOSTRAS

# Os nomes das linhas de X devem corresponder aos IDs das amostras
# no metadata, e estar exatamente na mesma ordem.

rownames(X) <- metadata_alinhado$sample_ID

stopifnot(
  nrow(X) == length(y),
  identical(
    rownames(X),
    metadata_alinhado$sample_ID
  )
)


## 4. DIVISÃO EM TREINAMENTO E TESTE

# O conjunto completo será dividido em:
#
#   80% = treinamento
#   20% = teste
#
# A divisão é estratificada pelo desfecho para preservar,
# aproximadamente, a proporção de tumores e não-tumores
# nos dois conjuntos.

set.seed(123)

indices_treino <- createDataPartition(
  y = factor(y),
  p = 0.80,
  list = FALSE
)

# Criar os conjuntos de expressão

X_train <- X[
  indices_treino,
  ,
  drop = FALSE
]

X_test <- X[
  -indices_treino,
  ,
  drop = FALSE
]

# Criar os respectivos vetores de resposta

y_train <- y[indices_treino]

y_test <- y[-indices_treino]


# Verificar a distribuição das classes

table(y_train)

table(y_test)


## 5. AJUSTE DO MODELO ELASTIC NET

# O Elastic Net combina duas formas de regularização:
#
#   LASSO (L1)  -> favorece seleção de variáveis
#   Ridge (L2)  -> ajuda a lidar com variáveis correlacionadas
#
# alpha = 0.5 produz uma combinação equilibrada entre as duas.
#
# O modelo é ajustado SOMENTE no conjunto de treinamento.
#
# cv.glmnet realiza uma validação cruzada de 10 folds dentro
# do conjunto de treinamento para determinar o valor de lambda.
#
# type.measure = "auc" faz com que a seleção de lambda seja
# baseada na capacidade discriminativa medida pela AUC.
#
# standardize = TRUE padroniza as variáveis antes do ajuste,
# procedimento particularmente importante quando as variáveis
# apresentam escalas diferentes.

set.seed(123)

modelo_enet <- cv.glmnet(
  x = X_train,
  y = y_train,
  family = "binomial",
  
  # Elastic Net: combinação de LASSO e Ridge
  alpha = 0.5,
  
  # Validação cruzada interna
  nfolds = 10,
  
  # Critério utilizado para selecionar lambda
  type.measure = "auc",
  
  # Padronização das variáveis
  standardize = TRUE
)


## 6. EXTRAÇÃO DOS GENES SELECIONADOS

# O lambda.min corresponde ao valor de lambda que apresentou
# a maior AUC média durante a validação cruzada.
#
# Extraímos os coeficientes do modelo nesse ponto.

coef_enet <- coef(
  modelo_enet,
  s = "lambda.min"
)

coef_enet <- as.matrix(coef_enet)


# Criar uma tabela contendo os genes e seus coeficientes

genes_selecionados <- data.frame(
  EntrezID = rownames(coef_enet),
  Coeficiente = coef_enet[, 1],
  row.names = NULL
)


# Manter somente os genes com coeficiente diferente de zero.
#
# Genes com coeficiente = 0 não foram selecionados pelo Elastic Net
# para o modelo final.

# Remover o intercepto e manter somente os genes
# com coeficientes diferentes de zero.
genes_selecionados <- genes_selecionados[
  genes_selecionados$EntrezID != "(Intercept)" &
    genes_selecionados$Coeficiente != 0,
]

## 7. ANOTAÇÃO DOS GENES SELECIONADOS

# Converter os Entrez IDs para Gene Symbols.

genes_selecionados$GeneSymbol <- mapIds(
  illuminaHumanv2.db,
  keys = genes_selecionados$EntrezID,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)


# Ordenar os genes pelo valor absoluto do coeficiente.
#
# Valores maiores em módulo representam maior contribuição
# do gene para o modelo, embora o coeficiente não deva ser
# interpretado isoladamente como medida de importância biológica.

genes_selecionados <- genes_selecionados[
  order(
    abs(genes_selecionados$Coeficiente),
    decreasing = TRUE
  ),
]


# Criar diretório para os resultados do machine learning

dir.create(
  file.path(results_dir, "machine_learning"),
  recursive = TRUE,
  showWarnings = FALSE
)


# Exportar a lista inicial de genes selecionados

write.csv(
  genes_selecionados,
  file.path(
    results_dir,
    "machine_learning",
    "genes_selecionados_lambda_min.csv"
  ),
  row.names = FALSE
)


## 8. PREDIÇÃO DAS AMOSTRAS DO CONJUNTO DE TESTE

# Utilizar o modelo treinado para calcular a probabilidade
# de cada amostra do conjunto de teste pertencer à classe tumor.

prob_tumor <- predict(modelo_enet,newx = X_test,s = "lambda.min",type = "response")

# Converter o resultado para vetor numérico

prob_tumor <- as.numeric(prob_tumor)


## 9. AVALIAÇÃO POR CURVA ROC

# A curva ROC avalia a capacidade do modelo de distinguir
# tumores de amostras não tumorais considerando diferentes
# pontos de corte.
#
# AUC = área sob a curva ROC.
# Quanto mais próxima de 1, maior a capacidade discriminativa.

roc_modelo <- roc(response = y_test,predictor = prob_tumor,direction = "<")

auc_modelo <- auc(roc_modelo)

auc_modelo


# Gerar gráfico da curva ROC
plot(roc_modelo,print.auc = TRUE,legacy.axes = TRUE,main = "Curva ROC — Elastic Net")


## 10. DEFINIÇÃO DO PONTO DE CORTE

# O ponto de corte será determinado pelo índice de Youden.
#
# O índice de Youden procura um ponto da curva ROC que maximize
# simultaneamente a sensibilidade e a especificidade.
#
# O resultado será utilizado para transformar as probabilidades
# previstas pelo modelo em classes:
#
#   abaixo do ponto de corte = não-tumor
#   acima do ponto de corte  = tumor

ponto_corte <- coords(roc_modelo,x = "best",best.method = "youden",ret = c("threshold","sensitivity","specificity"))

ponto_corte


# Converter o threshold para um valor numérico simples
threshold <- as.numeric(ponto_corte["threshold"])


## 11. CLASSIFICAÇÃO DAS AMOSTRAS DO TESTE

# Transformar as probabilidades previstas em classificações
# binárias utilizando o ponto de corte definido anteriormente.

classe_predita <- ifelse(prob_tumor >= threshold,1,0)


# Verificar se o número de previsões corresponde ao número
# de amostras do conjunto de teste.

stopifnot(length(classe_predita) == length(y_test))


## 12. MATRIZ DE CONFUSÃO

# Comparar as classes previstas pelo modelo com as classes
# observadas nas amostras de teste.
#
# Classe positiva = tumor (1).

matriz_confusao <- confusionMatrix(factor(classe_predita,levels = c(0, 1)),factor(y_test,levels = c(0, 1)),positive = "1")

matriz_confusao


## 13. COMPARAÇÃO COM OS RESULTADOS DO LIMMA

# Verificar quais genes selecionados pelo Elastic Net também
# foram identificados como DEGs pelo limma.
#
# Para isso, recuperar os valores de logFC e FDR correspondentes
# na tabela de resultados do limma.

genes_selecionados$logFC <- resultados_limma[match(genes_selecionados$EntrezID,resultados_limma$EntrezID),"logFC"]
genes_selecionados$FDR <- resultados_limma[match(genes_selecionados$EntrezID,resultados_limma$EntrezID),"adj.P.Val"]


# Definir se cada gene selecionado pelo Elastic Net também
# satisfaz os critérios utilizados para definir DEGs:
#   FDR < 0.05
#   |logFC| > 1

genes_selecionados$DEG_limma <- (genes_selecionados$FDR < 0.05 &abs(genes_selecionados$logFC) > 1)


# Contar quantos genes selecionados pelo Elastic Net também
# são DEGs pelo limma.
table(genes_selecionados$DEG_limma)


# Exportar tabela final dos genes selecionados, agora contendo
# também os resultados correspondentes da análise de limma.
write.csv(genes_selecionados,file.path(results_dir,"machine_learning","genes_selecionados_lambda_min_com_limma.csv"),row.names = FALSE)


## 14. RESUMO DO DESEMPENHO DO MODELO
# Criar uma tabela contendo os principais resultados da análise.
# Os valores de acurácia e balanced accuracy são obtidos
# diretamente da matriz de confusão, evitando a necessidade
# de inserir os valores manualmente.

resultados_modelo <- data.frame(
  AUC = as.numeric(auc_modelo),
  Threshold = threshold,
  Sensibilidade = as.numeric(matriz_confusao$byClass["Sensitivity"]),
  Especificidade = as.numeric(matriz_confusao$byClass["Specificity"]),
  Acuracia = as.numeric(matriz_confusao$overall["Accuracy"]),
  Balanced_Accuracy = as.numeric(matriz_confusao$byClass["Balanced Accuracy"]),
  Genes_Selecionados = nrow(genes_selecionados),
  Genes_DEG_Limma = sum(genes_selecionados$DEG_limma)
)


# Visualizar o resumo no console
resultados_modelo


# Exportar os resultados
write.csv(resultados_modelo,file.path(results_dir,"machine_learning","desempenho_elastic_net.csv"),row.names = FALSE)

# liberar memória
rm(X,y,indices_treino,X_train,X_test,y_train,y_test,coef_enet,prob_tumor,
   roc_modelo,ponto_corte,threshold,classe_predita,matriz_confusao)
gc()



# ============== Validação no TCGA ==============
# aqui é feito uma análise diferencial de dados padronizados de RNA-seq da base de
# dados Recount3 que possui amostras BLCA tumorais do TCGA e amostras não-tumorais
# de bexiga do GTex padronizadas por monorail. Os dados pré-processados foram
# extraídos e análisados por limma (voom) para servirem como validação dos resultados.

message("===========================================================================")
message("INICIANDO ANÁLISE DE VALIDAÇÃO...")
message("BASE DE DADOS RECOUNT3 (TCGA / GTEX)")
message("===========================================================================")

# rodando o script em novo ambiente e extraíndo os resultados
recount <- new.env()
source(file.path(scripts_dir,"TCGA + GTex validation.r"),local = recount)
validation_results <- get("results", envir = recount)
v <- get("v", envir = recount)
group_roc <- get("group_roc", envir = recount)
group <- get("group", envir = recount)
validation_exprs <- v$E
# remover Normal_Adj
keep_samples <- group != "Normal_Adj"
validation_exprs <- validation_exprs[, keep_samples]
group_roc <- group_roc[keep_samples]
rm(v,group,keep_samples,recount)
gc()

# filtando os resultados de validação com os genes da análise principal
genes_selecionados_limma <- genes_selecionados[
  genes_selecionados$DEG_limma == TRUE,
]

validation_results_filtered <- validation_results[
  rownames(validation_results) %in%
    as.character(genes_selecionados_limma$EntrezID),
]

#salvando os dados
# confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(results_dir, "TCGA + GTex validation")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

write.csv(validation_results,
          file.path(results_dir, 'TCGA + GTex validation', "validation_analysis_results.csv"),row.names = TRUE)
write.csv(validation_results_filtered,
          file.path(results_dir,"TCGA + GTex validation","validation_analysis_comparison.csv"),row.names = TRUE)

# ============== ROC/AUC ==============
# calculado as curvas ROC para cada gene que refletem a sensibilidade e especificidade
# dos genes para identificar condição de tumor e não tumor e o AUC que é um valor
# resume e escala os dados ROC para serem usados como critério de qualidade
# do gene como candidato a biomarcador.

message("===========================================================================")
message("ANALISANDO CURVAS ROC E AUC...")
message("===========================================================================")

# Criar labels binários (tumor/não-tumor)
labels <- ifelse(group_roc == "tumor", 1, 0)

# estimando curvas AUC dos genes
roc_results <- lapply(rownames(validation_exprs), function(gene_id) {
  values <- as.numeric(validation_exprs[gene_id, ])
  r <- pROC::roc(labels, values, quiet = TRUE)
  data.frame(
    gene = gene_id,
    auc = as.numeric(pROC::auc(r))
  )
})

#resultado final
roc_df <- do.call(rbind, roc_results)

# filtrando para DEGs identificados
roc_filtered <- roc_df %>%
  filter(as.character(gene) %in% as.character(genes_selecionados_limma$EntrezID))

## checando se o diretório existe
out_dir <- file.path(results_dir, "auc")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

#salvando resultados
write.csv(roc_filtered,
          file.path(results_dir, 'auc',"auc.csv"))


## criando um plot de exemplo
gene <- roc_filtered[1,1]
values <- as.numeric(validation_exprs[gene,])
roc_curve <- pROC::roc(labels, values)
auc_value <- auc(roc_curve)
png(file.path("figures",paste("ROC_curve_gene",gene,".png")),height = 1800, width = 1800,res = 300)
ggroc(roc_curve) +
  geom_abline(
    intercept = 1,
    slope = 1,
    linetype = "dashed",
    color = "gray"
  ) +
  ggtitle(paste("ROC curve - gene", gene)) +
  annotate(
    "text",
    x = 0.65,
    y = 0.2,
    label = paste("AUC =", round(auc_value, 3))
  ) +
  theme_minimal()
dev.off()

# limpando
rm(roc_df,gene,values,roc_curve)
gc()


# ============================================================
# REDE DE INTERAÇÃO PROTEÍNA-PROTEÍNA (PPI)
# STRING + DEGs + destaque dos candidatos finais
# ============================================================

message("\n", paste(rep("=", 50), collapse = ""))
message("Construindo rede PPI")
message(paste(rep("=", 50), collapse = ""))


# ------------------------------------------------------------
# 1. Preparação dos genes
# ------------------------------------------------------------

# Todos os DEGs da análise principal
genes_ppi <- unique(as.character(DEGs$EntrezID))
genes_ppi <- genes_ppi[!is.na(genes_ppi)]

# 15 candidatos finais:
# genes selecionados pelo Elastic Net que também são DEGs
genes_candidatos <- unique(
  as.character(
    genes_selecionados$EntrezID[
      genes_selecionados$DEG_limma == TRUE
    ]
  )
)

genes_candidatos <- genes_candidatos[!is.na(genes_candidatos)]

message("DEGs para a rede: ", length(genes_ppi))
message("Candidatos finais: ", length(genes_candidatos))


# ------------------------------------------------------------
# 2. Conectar ao STRING
# ------------------------------------------------------------

string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 700,
  input_directory = ""
)


# ------------------------------------------------------------
# 3. Mapear Entrez IDs para STRING IDs
# ------------------------------------------------------------

genes_df <- data.frame(
  ENTREZID = genes_ppi,
  stringsAsFactors = FALSE
)

mapped_genes <- string_db$map(
  genes_df,
  "ENTREZID",
  removeUnmappedRows = TRUE
)

message("Genes mapeados pelo STRING: ", nrow(mapped_genes))

# Verificar quantos candidatos finais foram mapeados
mapped_candidates <- mapped_genes[
  mapped_genes$ENTREZID %in% genes_candidatos,
]

message(
  "Candidatos finais mapeados: ",
  nrow(mapped_candidates)
)


# ------------------------------------------------------------
# 4. Obter interações entre os DEGs
# ------------------------------------------------------------

interactions <- string_db$get_interactions(
  mapped_genes$STRING_id
)

# Manter apenas interações em que os dois parceiros
# pertencem ao conjunto de DEGs
string_ids <- unique(mapped_genes$STRING_id)

interactions <- interactions[
  interactions$from %in% string_ids &
    interactions$to %in% string_ids,
]

message("Interações encontradas: ", nrow(interactions))


# ------------------------------------------------------------
# 5. Converter STRING IDs novamente para Entrez IDs
# ------------------------------------------------------------

id_map <- mapped_genes[, c("STRING_id", "ENTREZID")]

interactions <- merge(
  interactions,
  id_map,
  by.x = "from",
  by.y = "STRING_id",
  all.x = TRUE
)

names(interactions)[names(interactions) == "ENTREZID"] <- "ENTREZID_from"

interactions <- merge(
  interactions,
  id_map,
  by.x = "to",
  by.y = "STRING_id",
  all.x = TRUE
)

names(interactions)[names(interactions) == "ENTREZID"] <- "ENTREZID_to"


# ------------------------------------------------------------
# 6. Criar tabela de nós
# ------------------------------------------------------------

nodes <- data.frame(
  name = mapped_genes$ENTREZID,
  STRING_id = mapped_genes$STRING_id,
  ENTREZID = mapped_genes$ENTREZID,
  stringsAsFactors = FALSE
)

nodes$GeneSymbol <- resultados_limma$GeneSymbol[
  match(nodes$ENTREZID, resultados_limma$EntrezID)
]

nodes$logFC <- resultados_limma$logFC[
  match(nodes$ENTREZID, resultados_limma$EntrezID)
]

nodes$FDR <- resultados_limma$adj.P.Val[
  match(nodes$ENTREZID, resultados_limma$EntrezID)
]

nodes$Candidato <- nodes$ENTREZID %in% genes_candidatos


# ------------------------------------------------------------
# 7. Criar tabela de arestas
# ------------------------------------------------------------

edges <- interactions[, c(
  "ENTREZID_from",
  "ENTREZID_to",
  "combined_score"
)]

names(edges) <- c("from", "to", "score")

edges <- unique(edges)


# ------------------------------------------------------------
# 8. Criar objeto igraph
# ------------------------------------------------------------

ppi_graph <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = FALSE
)

message("Nós no grafo: ", vcount(ppi_graph))
message("Arestas no grafo: ", ecount(ppi_graph))


# ------------------------------------------------------------
# 9. CONFIGURAÇÃO DA REDE
# ------------------------------------------------------------

# Grau de cada gene = número de interações na rede
V(ppi_graph)$degree <- degree(ppi_graph)

# Destacar os candidatos finais
V(ppi_graph)$color <- ifelse(
  V(ppi_graph)$Candidato,
  "red",
  "lightgray"
)

# Candidatos maiores que os demais genes
V(ppi_graph)$size <- ifelse(
  V(ppi_graph)$Candidato,
  14,
  5
)

# Somente os candidatos recebem rótulos
V(ppi_graph)$label <- ifelse(
  V(ppi_graph)$Candidato,
  V(ppi_graph)$GeneSymbol,
  NA
)

# Tamanho dos rótulos
V(ppi_graph)$label.cex <- ifelse(
  V(ppi_graph)$Candidato,
  0.8,
  0
)

# Arestas discretas
E(ppi_graph)$color <- "gray80"
E(ppi_graph)$width <- 0.5


# ------------------------------------------------------------
# 10. LAYOUT DA REDE
# ------------------------------------------------------------

set.seed(123)

layout_ppi <- layout_with_fr(
  ppi_graph,
  weights = E(ppi_graph)$score / 1000
)

# Normalizar escala do layout
layout_ppi <- norm_coords(
  layout_ppi,
  ymin = -1,
  ymax = 1,
  xmin = -1,
  xmax = 1
)


# ------------------------------------------------------------
# 11. EXPORTAÇÃO DA REDE
# ------------------------------------------------------------

png(
  filename = file.path(
    figures_dir,
    "PPI_DEGs_STRING.png"
  ),
  width = 5000,
  height = 5000,
  res = 400
)

plot(
  ppi_graph,
  layout = layout_ppi,
  vertex.frame.color = "gray30",
  vertex.label.color = "black",
  edge.curved = 0.05,
  main = "Rede PPI baseada no STRING"
)

dev.off()


# Versão vetorial para publicação
pdf(
  file = file.path(
    figures_dir,
    "PPI_DEGs_STRING.pdf"
  ),
  width = 12,
  height = 12
)

plot(
  ppi_graph,
  layout = layout_ppi,
  vertex.frame.color = "gray30",
  vertex.label.color = "black",
  edge.curved = 0.05,
  main = "Rede PPI baseada no STRING"
)

dev.off()


# ------------------------------------------------------------
# 12. Exportar tabelas para Cytoscape
# ------------------------------------------------------------

write.csv(
  nodes,
  file.path(results_dir,"PPI", "PPI_nodes.csv"),
  row.names = FALSE
)

write.csv(
  interactions,
  file.path(results_dir,"PPI", "PPI_interactions.csv"),
  row.names = FALSE
)

write.csv(
  nodes[nodes$Candidato, ],
  file.path(results_dir,"PPI", "PPI_candidates_15.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Resumo
# ------------------------------------------------------------

message("\nRede PPI concluída.")
message("Nós: ", vcount(ppi_graph))
message("Arestas: ", ecount(ppi_graph))
message(
  "Candidatos finais mapeados na rede: ",
  sum(V(ppi_graph)$Candidato)
)


# ------------------------------------------------------------
# 15. Limpeza
# ------------------------------------------------------------

rm(
  genes_ppi, genes_candidatos, genes_df, mapped_genes,
  mapped_candidates, string_ids, id_map, nodes, edges,
  interactions, layout_ppi
)

gc()



# ============================================================
# RESULTADO FINAL — CANDIDATOS A BIOMARCADORES
# ============================================================

message("\n", paste(rep("=", 50), collapse = ""))
message("Gerando tabela final de candidatos a biomarcadores")
message(paste(rep("=", 50), collapse = ""))


# ------------------------------------------------------------
# 1. Selecionar os 15 candidatos finais
# ------------------------------------------------------------

candidatos_finais <- genes_selecionados[
  genes_selecionados$DEG_limma == TRUE,
]

candidatos_finais <- candidatos_finais[
  order(candidatos_finais$FDR),
]


# ------------------------------------------------------------
# 2. Adicionar resultados da validação externa
# ------------------------------------------------------------

# Resultados de expressão diferencial no recount3
validacao_expr <- validation_results[
  rownames(validation_results) %in% candidatos_finais$EntrezID,
  c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"),
  drop = FALSE
]

validacao_expr$EntrezID <- rownames(validacao_expr)

names(validacao_expr)[names(validacao_expr) == "logFC"] <- "logFC_validacao"
names(validacao_expr)[names(validacao_expr) == "adj.P.Val"] <- "FDR_validacao"


# AUC individual na validação externa
validacao_auc <- do.call(rbind, roc_results)
rownames(validacao_auc) <- NULL

validacao_auc$EntrezID <- as.character(validacao_auc$gene)
validacao_auc <- validacao_auc[, c("EntrezID", "auc")]


# ------------------------------------------------------------
# 3. Preparar resultados da análise principal
# ------------------------------------------------------------

candidatos_finais$logFC_principal <- candidatos_finais$logFC
candidatos_finais$FDR_principal <- candidatos_finais$FDR

candidatos_finais$logFC <- NULL
candidatos_finais$FDR <- NULL


# ------------------------------------------------------------
# 4. Combinar os resultados
# ------------------------------------------------------------

resultado_final <- merge(
  candidatos_finais,
  validacao_expr,
  by = "EntrezID",
  all.x = TRUE
)

resultado_final <- merge(
  resultado_final,
  validacao_auc,
  by = "EntrezID",
  all.x = TRUE
)


# ------------------------------------------------------------
# 5. Verificar concordância da expressão & Verificar se os candidatos são DEGs na validação
# ------------------------------------------------------------

resultado_final$DEG_validacao <- (
  abs(resultado_final$logFC_validacao) > 1 &
    resultado_final$FDR_validacao < 0.05
)

resultado_final$direcao_concordante <- sign(
  resultado_final$logFC_principal
) == sign(
  resultado_final$logFC_validacao
)


# ------------------------------------------------------------
# 6. Organizar as colunas
# ------------------------------------------------------------

resultado_final <- resultado_final[
  , c(
    "EntrezID",
    "GeneSymbol",
    "Coeficiente",
    "logFC_principal",
    "FDR_principal",
    "DEG_limma",
    "logFC_validacao",
    "FDR_validacao",
    "DEG_validacao",
    "auc",
    "direcao_concordante"
  )
]

# ------------------------------------------------------------
# 7. Exportar resultado final
# ------------------------------------------------------------

final_dir <- file.path(results_dir, "final_candidates")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  resultado_final,
  file.path(final_dir, "candidatos_biomarcadores_finais.csv"),
  row.names = FALSE
)

message("\nTabela final gerada:")
message(file.path(final_dir, "candidatos_biomarcadores_finais.csv"))
message("Número de candidatos: ", nrow(resultado_final))

print(resultado_final)


# ------------------------------------------------------------
# 8. Resumo
# ------------------------------------------------------------

message("\nTabela final gerada:")
message(
  file.path(
    final_dir,
    "candidatos_biomarcadores_finais.csv"
  )
)

message("Número de candidatos: ", nrow(resultado_final))

message(
  "Candidatos com expressão diferencial reproduzida: ",
  sum(resultado_final$FDR_validacao < 0.05, na.rm = TRUE)
)

message(
  "Candidatos com direção concordante: ",
  sum(resultado_final$direcao_concordante, na.rm = TRUE)
)

print(resultado_final)


# ============== HEATMAP FINAL — 13 CANDIDATOS ==============

# organizando genes finais
genes_finais <- as.character(resultado_final$EntrezID)
symbols_finais <- as.character(resultado_final$GeneSymbol)

names(symbols_finais) <- genes_finais

print(data.frame(
  EntrezID = genes_finais,
  GeneSymbol = symbols_finais
))


## HEATMAP — DESCOBERTA

# Verificar quais genes estão presentes
genes_descoberta <- intersect(genes_finais,rownames(expr_combat))

if (length(genes_descoberta) < length(genes_finais)) {
  
  message("Genes ausentes na descoberta: ",
    paste(setdiff(genes_finais, genes_descoberta),collapse = ", "))}
if (length(genes_descoberta) < 2) {stop("Menos de 2 genes finais encontrados na matriz de descoberta.")}

# Manter a ordem original dos genes finais
genes_descoberta <- genes_finais[genes_finais %in% genes_descoberta]

# Extrair expressão
matriz_descoberta <- expr_combat[genes_descoberta,
  ,drop = FALSE]

# Z-score por gene
matriz_descoberta_z <- t(scale(t(matriz_descoberta)))

# Remover genes que eventualmente tenham variância zero
matriz_descoberta_z <- matriz_descoberta_z[complete.cases(matriz_descoberta_z),
  ,drop = FALSE]

# Substituir Entrez pelos símbolos
symbols_descoberta <- symbols_finais[rownames(matriz_descoberta_z)]

symbols_descoberta[
  is.na(symbols_descoberta) |
    symbols_descoberta == ""
] <- rownames(matriz_descoberta_z
)[is.na(symbols_descoberta) |
    symbols_descoberta == ""]
rownames(matriz_descoberta_z) <- make.unique(symbols_descoberta)


# alinhando metadados de descoberta
# Identificar os IDs das amostras na matriz
if (exists("metadata_alinhado")) {metadata_heatmap_desc <- metadata_alinhado
  if (
    "sample_ID" %in% colnames(metadata_heatmap_desc) &&
    all(metadata_heatmap_desc$sample_ID %in% colnames(matriz_descoberta_z))
  ) {
    metadata_heatmap_desc$matrix_id <- metadata_heatmap_desc$sample_ID
  } else {metadata_heatmap_desc$matrix_id <- paste0(metadata_heatmap_desc$study_ID,"_",metadata_heatmap_desc$sample_ID)
  }
} else {stop("O objeto 'metadata_alinhado' não está disponível.")}

# Manter somente amostras presentes na matriz
metadata_heatmap_desc <-metadata_heatmap_desc[metadata_heatmap_desc$matrix_id %in% colnames(matriz_descoberta_z),
    ,drop = FALSE]

# Ordenar biologicamente
ordem_desc <- order(metadata_heatmap_desc$sample_type)

metadata_heatmap_desc <- metadata_heatmap_desc[ordem_desc, , drop = FALSE]

# Manter somente as amostras presentes em ambos
amostras_desc <- intersect(metadata_heatmap_desc$matrix_id,colnames(matriz_descoberta_z))

metadata_heatmap_desc <- metadata_heatmap_desc[metadata_heatmap_desc$matrix_id %in% amostras_desc,
    ,drop = FALSE]

# Reordenar a matriz de acordo com o metadata
matriz_descoberta_z <-matriz_descoberta_z[,metadata_heatmap_desc$matrix_id,drop = FALSE]

# anotações de descoberta
anotacoes_desc <- data.frame(
  Tipo = metadata_heatmap_desc$sample_type,
  Estudo = metadata_heatmap_desc$study_ID,
  row.names = colnames(matriz_descoberta_z))


# salvando heatmap de descoberta

arquivo_descoberta <- file.path(
  figures_dir,
  "heatmap_13_candidatos_descoberta.png"
)

png(
  filename = arquivo_descoberta,
  width = 3000,
  height = 2200,
  res = 300
)

pheatmap(
  matriz_descoberta_z,
  annotation_col = anotacoes_desc,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 10,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "13 candidatos finais — Descoberta",
  border_color = NA
)

dev.off()


## HEATMAP — VALIDAÇÃO recount3

# Verificar genes presentes
genes_validacao <- intersect(
  genes_finais,
  rownames(validation_exprs)
)

if (length(genes_validacao) < length(genes_finais)) {
  
  message(
    "Genes ausentes na validação: ",
    paste(
      setdiff(genes_finais, genes_validacao),
      collapse = ", "
    )
  )
}

if (length(genes_validacao) < 2) {
  stop("Menos de 2 genes finais encontrados na validação.")
}

# Manter ordem original
genes_validacao <- genes_finais[
  genes_finais %in% genes_validacao]

# Extrair expressão
matriz_validacao <- validation_exprs[
  genes_validacao,
  ,drop = FALSE]

# Z-score independente dentro da validação
matriz_validacao_z <- t(
  scale(t(matriz_validacao)))

# Remover genes sem variância
matriz_validacao_z <- matriz_validacao_z[
  complete.cases(matriz_validacao_z),
  ,drop = FALSE]

# Usar símbolos
symbols_validacao <- symbols_finais[
  rownames(matriz_validacao_z)
]

symbols_validacao[
  is.na(symbols_validacao) |
    symbols_validacao == ""
] <- rownames(
  matriz_validacao_z
)[
  is.na(symbols_validacao) |
    symbols_validacao == ""
]

rownames(matriz_validacao_z) <- make.unique(symbols_validacao)


# organizando metadados de validação
# group_roc já está alinhado com validation_exprs
if (length(group_roc) != ncol(matriz_validacao_z)) {
  stop("group_roc não possui o mesmo número de amostras da validação.")}

# Criar metadata
metadata_heatmap_val <- data.frame(
  Tipo = group_roc,row.names = colnames(matriz_validacao_z))

# Ordenar: não-tumor primeiro, tumor depois
ordem_val <- order(metadata_heatmap_val$Tipo)

metadata_heatmap_val <-
  metadata_heatmap_val[
    ordem_val,
    ,
    drop = FALSE
  ]

# Reordenar matriz
matriz_validacao_z <- matriz_validacao_z[,rownames(metadata_heatmap_val),drop = FALSE]

# salvando heatmap de validação

arquivo_validacao <- file.path(
  figures_dir,
  "heatmap_13_candidatos_validacao.png"
)

png(
  filename = arquivo_validacao,
  width = 4000,
  height = 2200,
  res = 300
)

pheatmap(
  matriz_validacao_z,
  annotation_col = metadata_heatmap_val,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 10,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "13 candidatos finais — Validação",
  border_color = NA
)

dev.off()

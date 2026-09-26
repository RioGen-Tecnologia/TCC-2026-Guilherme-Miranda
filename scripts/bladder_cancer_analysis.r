
# Análise integrada de projetos de microarray obtidas do Gene Expression Omnibus
# Guilherme Moret Miranda - Riogen
# 25/04/2026
# Esse script obtém, normaliza e análisa cada um dos projetos buscados previamente,
# combina os dados de cada projeto por meta-análise e realiza posteriores análises
# para determinar biomarcadores gênicos de câncer de bexiga.

# ============== CARREGANDO PACOTES ==============

message("===========================================================================")
message("CARREGANDO PACOTES")
message("===========================================================================")

library(here)
source(here("scripts", "setup.r"))

# ============== DEFINIÇÃO DE DIRETÓRIO E ARQUIVOS ==============
# Definir as pastas de trabalho e o master manifesto, com anotações das amostras

geo_dir      <- here("data", "raw", "GEO")
processed_dir <- here("data", "processed")
results_dir  <- here("results")
figures_dir  <- here("figures")
scripts_dir  <- here("scripts")
metadata_dir <- here("metadata")

dirs <- c(
  geo_dir,
  processed_dir,
  results_dir,
  figures_dir,
  scripts_dir
)

# cria os diretórios se eles não existirem
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# Master Manifesto com anotações das amostras
metadata_path <- list.files(
  here("metadata"),
  pattern = "bladder_cancer_metadata.csv",
  full.names = TRUE
)

# Master Manifesto com anotações das amostras
metadata <- read.csv(file.path(metadata_dir,"bladder_cancer_metadata.csv"))

# limpando
rm(d,dirs)
gc()

# ============== ANÁLISE DE PROJETOS ==============
# Esta sessão executa os scripts secundários que extraem os dados do GEO, normalizam,
# anotam e exportam as matrizes de expressão para cada projeto.

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

# carregando as matrizes de expressão
arquivos_expressao <- list()
for (id in gsub("\\.r$", "", basename(scripts_projetos))) {
  arquivo <- file.path(
    processed_dir,
    id,
    paste0("exprs_", id, ".rds")
  )
  if (file.exists(arquivo)) {
    message("Carregando: ", arquivo)
    arquivos_expressao[[id]] <- readRDS(arquivo)
  } else {
    warning("Arquivo não encontrado: ", arquivo)
  }
}

#limpando
rm(arquivo,scripts_projetos,script,id)
gc()

# ============== INTEGRAÇÃO COMBAT ==============
# Integração das matrizes de expressão em uma super matriz através do pacote
# sva pelo método comBat, com os batches sendo os projetos individuais e as variações
# biológicas o estado de tumor e não-tumor

# Encontrando a interseção de genes (Entrez IDs) presentes em todos os GSEs
genes_comuns <- Reduce(intersect, lapply(arquivos_expressao, colnames))

# Filtrando as matrizes apenas para os genes comuns, transpor e converter para matriz
lista_matrizes <- lapply(arquivos_expressao, function(df) {
  mat <- as.matrix(df[, genes_comuns])
  return(t(mat)) # Linhas = genes, Colunas = amostras
})

# Combinar tudo em uma única super matriz
super_matrix <- do.call(cbind, lista_matrizes)

# Alinhar perfeitamente o metadata com as colunas da super matriz
metadata_alinhado <- metadata[match(colnames(super_matrix), metadata$sample_ID), ]

# Criar a matriz do modelo biológico (tumor x não-tumor)
mod_biologico <- model.matrix(~ as.factor(sample_type), data = metadata_alinhado)

# Definir a variável de lote garantindo o mesmo alinhamento
lote <- as.factor(metadata_alinhado$study_ID)

# Executar o ComBat preservando a variação biológica
super_matrix_corrigida <- ComBat(
  dat = super_matrix, 
  batch = lote, 
  mod = mod_biologico, 
  par.prior = TRUE
)

rm(genes_comuns,lista_matrizes)
gc()

# ============== PCA antes/depois ComBat ==============
# Gerando gráficos de PCA 2D e 3D sobre os dados antes e após a integração para
# analisar os efeitos de batch

# Função auxiliar para calcular PCA e formatar a tabela para plotagem
calcular_pca <- function(matriz_expressao, meta) {
  # prcomp espera amostras nas linhas e genes nas colunas, por isso a transposição t()
  pca <- prcomp(t(matriz_expressao), scale. = TRUE)
  # Percentual de variação explicada pelas 3 primeiras PCs
  var_expl <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 2)
  df_pca <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    study_ID = as.factor(meta$study_ID),
    sample_type = as.factor(meta$sample_type),
    sample_ID = meta$sample_ID
  )
  return(list(data = df_pca, var = var_expl))
}

#Calcular PCA para antes e depois do ComBat
pca_antes <- calcular_pca(super_matrix, metadata_alinhado)
pca_depois <- calcular_pca(super_matrix_corrigida, metadata_alinhado)

# -----------------------------------------------------------------------------
# VISUALIZAÇÃO 2D (ggplot2)
# -----------------------------------------------------------------------------

# Função para plotar 2D destacando o Lote (study_ID) e o Tipo de Amostra (sample_type)
plot_pca_2d <- function(pca_obj, titulo) {
  df <- pca_obj$data
  var <- pca_obj$var
  
  ggplot(df, aes(x = PC1, y = PC2, color = study_ID, shape = sample_type)) +
    geom_point(size = 2.5, alpha = 0.8) +
    labs(
      title = titulo,
      x = paste0("PC1 (", var[1], "%)"),
      y = paste0("PC2 (", var[2], "%)"),
      color = "Estudo (Lote)",
      shape = "Tipo de Amostra"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

p1 <- plot_pca_2d(pca_antes, "Antes do ComBat (Efeito de Lote)")
p2 <- plot_pca_2d(pca_depois, "Depois do ComBat (Harmonizado)")

# -----------------------------------------------------------------------------
# VISUALIZAÇÃO 3D (plotly)
# -----------------------------------------------------------------------------

plot_pca_3d <- function(pca_obj, titulo) {
  df <- pca_obj$data
  var <- pca_obj$var
  
  plot_ly(
    df, 
    x = ~PC1, y = ~PC2, z = ~PC3, 
    color = ~study_ID, 
    symbol = ~sample_type,
    symbols = c('circle', 'diamond'), # Forma para non_tumor e tumor
    text = ~paste("Amostra:", sample_ID, "<br>Estudo:", study_ID, "<br>Tipo:", sample_type),
    hoverinfo = 'text',
    type = 'scatter3d', 
    mode = 'markers',
    marker = list(size = 4)
  ) %>%
    layout(
      title = titulo,
      scene = list(
        xaxis = list(title = paste0("PC1 (", var[1], "%)")),
        yaxis = list(title = paste0("PC2 (", var[2], "%)")),
        zaxis = list(title = paste0("PC3 (", var[3], "%)"))
      )
    )
}

# Executar os plots 3D (abrem no viewer do RStudio)
plot_3d_antes  <- plot_pca_3d(pca_antes, "PCA 3D - Antes do ComBat")
plot_3d_depois <- plot_pca_3d(pca_depois, "PCA 3D - Depois do ComBat")


# Salvando o gráfico PCA 2D (Composição Lado a Lado)
p_2d_combinado <- p1 + p2

ggsave(
  filename = file.path(figures_dir, "pca_2d_antes_depois_comBat.png"),
  plot = p_2d_combinado,
  width = 12,
  height = 5,
  dpi = 300
)

# Salvando os gráficos PCA 3D Interativos (HTML)
saveWidget(
  widget = plot_3d_antes,
  file = file.path(figures_dir, "pca_3d_antes_comBat.html"),
  selfcontained = TRUE
)

saveWidget(
  widget = plot_3d_depois,
  file = file.path(figures_dir, "pca_3d_depois_comBat.html"),
  selfcontained = TRUE
)

# limpeza
rm(p1, p2, p_2d_combinado, plot_3d_antes, plot_3d_depois,
   pca_antes, pca_depois,calcular_pca, plot_pca_2d, plot_pca_3d)
gc()


# ============== ANÁLISE ESTATÍSTICA (limma) ==============
# realizando análise estatística sobre a super matriz com pacote limma com o método
# empirical bayes

# Definir os grupos biológicos e criar a matriz de desenho (design)
# Garantir que a coluna seja fator e definir "non_tumor" como nível de referência (controle)
grupo <- factor(metadata_alinhado$sample_type)
grupo <- relevel(grupo, ref = "non_tumor")

# Criar matriz de desenho sem intercepto para facilitar a definição de contrastes
design <- model.matrix(~ 0 + grupo)
colnames(design) <- levels(grupo)

# Ajustar o modelo linear aos dados de expressão
fit <- lmFit(super_matrix_corrigida, design)

# Construir a matriz de contraste (tumor vs non_tumor)
contrast_matrix <- makeContrasts(
  tumor_vs_non_tumor = tumor - non_tumor,
  levels = design
)

# Ajustar o modelo aos contrastes
fit_contrast <- contrasts.fit(fit, contrast_matrix)

# Aplicar o método Empírico de Bayes
fit_ebayes <- eBayes(fit_contrast)

# Extrair a tabela completa de resultados ---
# coef = 1 refere-se ao contraste tumor_vs_non_tumor
# number = Inf retorna todos os genes
# adjust.method = "BH" aplica a correção FDR (Benjamini-Hochberg)
tabela_completa <- topTable(fit_ebayes, coef = 1, number = Inf, adjust.method = "BH")

# Transformar o ID do gene (que está nos rownames) em uma coluna explícita
tabela_completa$Entrez_ID <- rownames(tabela_completa)

# Filtrar os DEGs com FDR < 0.05 e |logFC| > 1
degs_superam <- subset(
  tabela_completa, 
  adj.P.Val < 0.05 & abs(logFC) > 1
)

# Separar em genes Super-expressos (Upregulated) e Sub-expressos (Downregulated)
degs_up <- subset(degs_superam, logFC > 1)
degs_down <- subset(degs_superam, logFC < -1)

# Imprimir resumo no console
message("====================================================\n")
message("RESUMO DE GENES DIFERENCIALMENTE EXPRESSOS (DEGs):\n")
message("Total de DEGs significativos:", nrow(degs_superam), "\n")
message("Upregulated (Super-expressos no tumor):", nrow(degs_up), "\n")
message("Downregulated (Sub-expressos no tumor):", nrow(degs_down), "\n")
message("====================================================\n")

# --- 8. Exportar os resultados em CSV ---
write.csv(tabela_completa, file.path(results_dir, "limma/limma_tabela_completa.csv"), row.names = FALSE)
write.csv(degs_superam, file.path(results_dir, "limma/limma_DEGs_FDR005_logFC1.csv"), row.names = FALSE)
write.csv(degs_up, file.path(results_dir, "limma/limma_DEGs_Upregulated.csv"), row.names = FALSE)
write.csv(degs_down, file.path(results_dir, "limma/limma_DEGs_Downregulated.csv"), row.names = FALSE)

# --- 9. Limpeza de objetos intermediários da memória RAM ---
rm(design, fit, contrast_matrix, fit_contrast, fit_ebayes, grupo)
gc()


# ============== GRÁFICOS DE EXPRESSÃO ==============
# gerando gráficos de expressão, volcano plot e heatmap

message("===========================================================================")
message("GERANDO GRÁFICOS DE EXPRESSÃO")
message("===========================================================================")

# ============== GRÁFICOS DE EXPRESSÃO (Com Gene Symbols) ==============

message("===========================================================================")
message("GERANDO GRÁFICOS DE EXPRESSÃO (COM GENE SYMBOLS)")
message("===========================================================================")

# MAPEAMENTO DE ENTREZ ID PARA GENE SYMBOL
# Pegamos todos os Entrez IDs presentes na tabela e buscamos seus símbolos oficiais
all_entrezs <- as.character(tabela_completa$Entrez_ID)
symbol_mapping <- mapIds(org.Hs.eg.db,
                         keys = all_entrezs,
                         keytype = "ENTREZID",
                         column = "SYMBOL",
                         multiVals = "first")

# Adicionar a coluna Gene_Symbol na tabela completa
tabela_completa$Gene_Symbol <- symbol_mapping

# Caso algum gene não tenha símbolo mapeado, mantemos o Entrez ID como fallback
tabela_completa$Gene_Symbol[is.na(tabela_completa$Gene_Symbol)] <- tabela_completa$Entrez_ID[is.na(tabela_completa$Gene_Symbol)]


## ==== VOLCANO PLOT ====

# Adicionar coluna de status biológico para o plot
tabela_completa$diffexpression <- "Não significativo"
tabela_completa$diffexpression[tabela_completa$logFC > 1 & tabela_completa$adj.P.Val < 0.05] <- "Upregulated"
tabela_completa$diffexpression[tabela_completa$logFC < -1 & tabela_completa$adj.P.Val < 0.05] <- "Downregulated"

tabela_completa$diffexpression <- factor(tabela_completa$diffexpression, 
                                         levels = c("Downregulated", "Não significativo", "Upregulated"))

# Selecionar os 15 genes mais significativos por p-valor ajustado para rotular
top_genes <- head(tabela_completa[order(tabela_completa$adj.P.Val), ], 15)

# Construir o Volcano Plot (usando Gene_Symbol nos rótulos)
volcano_plot <- ggplot(tabela_completa, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpression)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(
    values = c("Downregulated" = "#377EB8", "Não significativo" = "grey70", "Upregulated" = "#E41A1C"),
    name = "Expressão Diferencial"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.7) +
  geom_text_repel(
    data = top_genes,
    aes(label = Gene_Symbol), # <--- AGORA USA O SÍMBOLO DO GENE
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
  filename = file.path(figures_dir, "volcano_plot_DEGs.png"),
  plot = volcano_plot,
  width = 8,
  height = 7,
  dpi = 300
)


## ==== HEATMAP ====

# Selecionar os 50 DEGs mais significativos (menor FDR)
top50_degs_df <- head(degs_superam[order(degs_superam$adj.P.Val), ], 50)
top50_entrez <- top50_degs_df$Entrez_ID

# Submeter a matriz corrigida apenas para os top 50 genes
matriz_heatmap <- super_matrix_corrigida[top50_entrez, ]

# Padronizar a expressão por gene (Z-score por linha)
matriz_heatmap_scaled <- t(scale(t(matriz_heatmap)))

# Substituir os nomes das linhas (Entrez IDs) pelos Gene Symbols correspondentes na matriz
# Usamos o mapeamento prévio garantindo a mesma ordem das linhas
symbols_top50 <- symbol_mapping[top50_entrez]
# Se houver duplicatas ou NAs de símbolos, garantimos unicidade limpa
symbols_top50[is.na(symbols_top50)] <- top50_entrez[is.na(symbols_top50)]
rownames(matriz_heatmap_scaled) <- make.unique(as.character(symbols_top50))

# Criar identificadores únicos combinando Estudo e Amostra
nomes_unicos <- paste0(metadata_alinhado$study_ID, "_", metadata_alinhado$sample_ID)
nomes_unicos_seguros <- make.unique(nomes_unicos)

# Atribuir os nomes únicos às colunas das matrizes
colnames(super_matrix_corrigida) <- nomes_unicos_seguros
colnames(matriz_heatmap_scaled) <- nomes_unicos_seguros

# Anotações laterais com chaves únicas garantidas
annotation_col <- data.frame(
  Tipo = metadata_alinhado$sample_type,
  Estudo = metadata_alinhado$study_ID,
  row.names = nomes_unicos_seguros
)

# Cores customizadas para as anotações
annotation_colors <- list(
  Tipo = c("non_tumor" = "#4DAF4A", "tumor" = "#984EA3")
)

# Caminho de saída do Heatmap
heatmap_filepath <- file.path(figures_dir, "heatmap_top50_DEGs.png")

# Gerar e salvar o Heatmap
png(filename = heatmap_filepath, width = 10, height = 8, units = "in", res = 300)
pheatmap(
  matriz_heatmap_scaled,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "Heatmap: Top 50 Genes Diferencialmente Expressos (Z-Score)"
)
dev.off()

# Limpeza de memória
rm(volcano_plot, top_genes, top50_degs_df, top50_entrez, matriz_heatmap, matriz_heatmap_scaled, 
   all_entrezs, symbol_mapping, symbols_top50,
   nomes_unicos, nomes_unicos_seguros, annotation_col, annotation_colors, heatmap_filepath)
gc()

# ============== ENRIQUECIMENTO FUNCIONAL ==============
# realiza enriquecimento funcional (Over-representation analysis) com genes
# diferencialmente expressos nas bases de dados KEGG, REACTOME e GO (biological
# process).

message("===========================================================================")
message("REALIZANDO ENRIQUECIMENTO FUNCIONAL E GERANDO RESULTADOS")
message("===========================================================================")

# PREPARAÇÃO & ANÁLISE ORA
gene_list <- as.character(degs_superam$Entrez_ID)

# definindo universo
universe <- as.character(rownames(super_matrix_corrigida))

# Execução dos testes de sobtrerrepresentação
ora_go       <- enrichGO(gene = gene_list, universe = universe, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
ora_kegg     <- setReadable(enrichKEGG(gene = gene_list, universe = universe, organism = 'hsa', pvalueCutoff = 0.05), OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ora_reactome <- enrichPathway(gene = gene_list, universe = universe, pvalueCutoff = 0.05, readable = TRUE)

# Garantir que as colunas textuais do resultado sejam vetores atômicos planos
ora_go@result$Description <- as.character(ora_go@result$Description)
ora_go@result$ID <- as.character(ora_go@result$ID)
ora_go@result$GeneRatio <- as.character(ora_go@result$GeneRatio)
ora_kegg@result$Description <- as.character(ora_kegg@result$Description)
ora_reactome@result$Description <- as.character(ora_reactome@result$Description)

# Exportação das tabelas
dir.create(file.path(results_dir, "enrichment"), showWarnings = FALSE, recursive = TRUE)
write.csv(as.data.frame(ora_go),       file.path(results_dir, "enrichment/ORA_GO_BP.csv"))
write.csv(as.data.frame(ora_kegg),     file.path(results_dir, "enrichment/ORA_KEGG.csv"))
write.csv(as.data.frame(ora_reactome), file.path(results_dir, "enrichment/ORA_Reactome.csv"))

# GERAÇÃO DOS PAINÉIS LADO A LADO
# Função auxiliar para padronizar o estilo visual
make_dotplot <- function(obj, titulo, n_cat = 10) {
  dotplot(obj, showCategory = n_cat, title = titulo) +
    theme_bw() +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.text.y  = element_text(size = 8),
      legend.position = "right"
    )
}

# Painéis individuais (exibindo as top 10 vias de cada)
p_go       <- make_dotplot(ora_go, "GO: Biological Process", n_cat = 10)
p_kegg     <- make_dotplot(ora_kegg, "KEGG Pathways", n_cat = 10)
p_reactome <- make_dotplot(ora_reactome, "Reactome Pathways", n_cat = 10)

# Junção e exportação via patchwork
combined_dotplots <- (p_go | p_kegg | p_reactome) + 
  plot_annotation(
    title    = "Análise de Enriquecimento Funcional (ORA)",
    subtitle = "Comparação entre GO Biological Process, KEGG e Reactome",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
                     plot.subtitle = element_text(size = 11, hjust = 0.5))
  )

ggsave(
  filename = file.path(figures_dir, "ora_dotplots_combinados.png"),
  plot     = combined_dotplots,
  width    = 20,
  height   = 7,
  dpi      = 300
)

# Limpeza do ambiente R
rm(ora_go, ora_kegg, ora_reactome, p_go, p_kegg, p_reactome, combined_dotplots, make_dotplot)
gc()


# ============== PPI NETWORK ==============
# análise de rede de interação proteína-proteína dos genes diferencialmente expressos
# com o objetivo de encontrar "hub genes", genes centrais na rede tumoral ou que
# parecem coordenar múltiplos processos.

message("===========================================================================")
message("ANÁLISE PROTEIN-PROTEIN INTERACTION")
message("INICIANDO STRINGdb")
message("===========================================================================")

# 1. Inicializar o STRINGdb
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700)

# 2. Mapear TODOS os Entrez IDs dos DEGs significativos
degs_mapped <- string_db$map(degs_superam, "Entrez_ID", removeUnmappedRows = TRUE)

# 3. Gerar e salvar o plot nativo do STRINGdb com TODOS os DEGs
# (Aviso: se forem muitos genes, a imagem pode ficar densa)
png(file.path(figures_dir, "ppi_network_all_degs.png"), width = 12, height = 12, units = "in", res = 300)
string_db$plot_network(degs_mapped$STRING_id)
dev.off()

# 4. Extrair a sub-rede completa para o igraph
ppi_igraph <- string_db$get_subnetwork(degs_mapped$STRING_id)

# 5. Verificação de segurança e cálculo dos Hubs (Via Matriz de Adjacência)
if (vcount(ppi_igraph) > 0 && ecount(ppi_igraph) > 0) {
  
  # Obter a matriz de adjacência da rede
  adj_matrix <- as.matrix(as_adjacency_matrix(ppi_igraph, sparse = FALSE))
  
  # O grau de cada nó é a soma das conexões (linhas) na matriz
  node_degrees <- rowSums(adj_matrix)
  
  # Garantir que é um vetor numérico atômico com os nomes dos genes
  node_degrees <- as.numeric(node_degrees)
  names(node_degrees) <- rownames(adj_matrix)
  
  # Ordenar do maior para o menor grau
  top_hubs <- sort(node_degrees, decreasing = TRUE)
  top_hubs_to_show <- head(top_hubs, min(10, length(top_hubs)))
  
  cat("====================================================\n")
  cat("Top Hub Genes (Proteínas com mais conexões na rede):\n")
  print(top_hubs_to_show)
  cat("====================================================\n")
  
} else {
  cat("Aviso: A rede gerada não possui conexões suficientes com o threshold de score > 700.\n")
}

# 6. Limpeza de RAM
rm(degs_mapped, ppi_igraph, adj_matrix, node_degrees)
gc()


# ============== MACHINE LEARNING ==============

message("===========================================================================")
message("MACHINE LEARNING: ELASTIC NET (REGRESSÃO LOGÍSTICA PENALIZADA)")
message("===========================================================================")

# --- 1. PREPARAÇÃO DA MATRIZ X E VETOR Y ---

# Extrair apenas os genes diferencialmente expressos (DEGs) para usar como features
features <- degs_superam$Entrez_ID

# Matriz X: Transpor para que Amostras = Linhas e Genes = Colunas
X <- t(super_matrix_corrigida[features, ])

# Vetor Y: Variável resposta binária (tumor vs. non_tumor) transformada em fator numérico para estabilidade
Y_factor <- as.factor(metadata_alinhado$sample_type)

# --- 2. DIVISÃO EM TREINO E TESTE (70/30) ---
set.seed(123) # Reprodutibilidade
train_index <- createDataPartition(Y_factor, p = 0.7, list = FALSE)

X_train <- X[train_index, ]
Y_train <- Y_factor[train_index]

X_test  <- X[-train_index, ]
Y_test  <- Y_factor[-train_index]

message(sprintf("Tamanho do Treino: %d amostras | Tamanho do Teste: %d amostras", nrow(X_train), nrow(X_test)))

# --- 3. TREINAMENTO DO MODELO (CROSS-VALIDATION) ---
# alpha = 0.5 (Elastic Net)
# family = "binomial" (Regressão Logística para 2 classes)
# type.measure = "auc" (Otimiza o modelo maximizando a área sob a curva ROC)

set.seed(123)
cv_elastic_net <- cv.glmnet(
  x = as.matrix(X_train), 
  y = Y_train, 
  family = "binomial", 
  alpha = 0.5, 
  type.measure = "auc",
  nfolds = 10
)

# Salvar gráfico da validação cruzada (mostra o Lambda ideal)
png(file.path(figures_dir, "elastic_net_cv_plot.png"), width = 8, height = 6, units = "in", res = 300)
plot(cv_elastic_net)
title("Cross-Validation: Elastic Net", line = 2.5)
dev.off()

# Obter o melhor Lambda (lambda.min minimiza o erro, lambda.1se é mais conservador/esparso)
best_lambda <- cv_elastic_net$lambda.min
message(sprintf("Melhor Lambda (lambda.min) encontrado: %f", best_lambda))

# --- 4. PREDIÇÃO E AVALIAÇÃO NO CONJUNTO DE TESTE ---

# Predizer as probabilidades (para ROC) e as classes (para Matriz de Confusão)
prob_preds <- predict(cv_elastic_net, newx = as.matrix(X_test), s = "lambda.min", type = "response")
class_preds <- predict(cv_elastic_net, newx = as.matrix(X_test), s = "lambda.min", type = "class")

# Converter previsões para fator com os mesmos níveis de Y_test
class_preds <- factor(class_preds, levels = levels(Y_test))

# Matriz de Confusão
conf_matrix <- confusionMatrix(class_preds, Y_test, positive = "tumor")
print(conf_matrix)

# Calcular e Plotar Curva ROC
roc_curve <- roc(response = Y_test, predictor = as.numeric(prob_preds))
auc_value <- auc(roc_curve)
message(sprintf("AUC no conjunto de teste: %.4f", auc_value))

png(file.path(figures_dir, "elastic_net_roc_curve.png"), width = 6, height = 6, units = "in", res = 300)
plot(roc_curve, col = "red", lwd = 2, main = paste0("Curva ROC (Elastic Net)\nAUC = ", round(auc_value, 4)))
dev.off()

# --- 5. EXTRAÇÃO DA ASSINATURA DE BIOMARCADORES ---

# Extrair os coeficientes do modelo otimizado (lambda.min)
coef_matrix <- coef(cv_elastic_net, s = "lambda.min")

# Converter para dataframe e remover os genes que foram zerados pela penalização
coef_df <- data.frame(
  Gene_Entrez = rownames(coef_matrix),
  Coefficient = as.numeric(coef_matrix)
)
assinatura_genes <- coef_df[coef_df$Coefficient != 0 & coef_df$Gene_Entrez != "(Intercept)", ]

# Ordenar por importância absoluta (módulo do coeficiente)
assinatura_genes <- assinatura_genes[order(abs(assinatura_genes$Coefficient), decreasing = TRUE), ]

# Anotar com Gene Symbols para interpretação biológica
assinatura_symbols <- mapIds(org.Hs.eg.db, 
                             keys = assinatura_genes$Gene_Entrez, 
                             keytype = "ENTREZID", 
                             column = "SYMBOL", 
                             multiVals = "first")

assinatura_genes$Gene_Symbol <- assinatura_symbols

# Salvar a assinatura de genes selecionada
write.csv(assinatura_genes, file.path(results_dir, "machine_learning", "ElasticNet_Biomarcadores.csv"), row.names = FALSE)

message("===========================================================================")
message(sprintf("Modelo treinado! O Elastic Net reduziu os %d DEGs originais para uma assinatura de %d biomarcadores preditivos.", length(features), nrow(assinatura_genes)))
message("Resultados salvos na pasta 'results/machine_learning/'.")

# Limpeza
rm(X, Y_factor, train_index, X_train, Y_train, X_test, Y_test, prob_preds, class_preds, coef_matrix, coef_df, assinatura_symbols)
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
validation_results_filtered <- validation_results[
  rownames(validation_results) %in% as.character(DEGs_filtered$Gene),
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
  filter(as.character(gene) %in% as.character(DEGs_filtered$Gene))

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
rm(roc_results,roc_df,gene,values,roc_curve)
gc()

# ============== COMPILAÇÃO DE RESULTADOS ==============
# Unindo todos os resultados para fazer a pontuação.

message("===========================================================================")
message("COMPILANDO RESULTADOS")
message("===========================================================================")

# criando dataframe dos resultados compilados
compiled_results <- tibble("Entrez"=DEGs_filtered$Gene,
                   "Gene_Symbol"=DEGs_filtered$Symbol,
                   "Hedges_g"=DEGs_filtered$HedgesG_pool,
                   "FDR"=DEGs_filtered$FDR,
                   "I2"=DEGs_filtered$I2)

# adicionando rede PPI
compiled_results <- compiled_results %>%
  left_join(
    hubs_table %>%
      select(
        Gene,
        degree
      ) %>%
      rename(
        Entrez = Gene,
        PPI_degree = degree
      ),
    by = "Entrez"
  )

# adicionando resultados de signficancia
compiled_results <- compiled_results %>%
  left_join(
    significance_results %>%
      select(gene,score) %>%
      rename(
        Entrez = gene,
        significance_score = score
      ),
    by = "Entrez"
  )

# dados de validação
compiled_results <- compiled_results %>%
  left_join(
    validation_results_filtered %>%
      rownames_to_column("Entrez") %>%
      select(Entrez, HedgesG_val, FDR_val),
    by = c("Entrez" = "Entrez")
  )

# dados AUC
compiled_results <- compiled_results %>%
  left_join(
    roc_filtered %>%
      select(gene, auc),
    by = c("Entrez" = "gene")
  )

#convertendo a dataframe formal
compiled_results <- as.data.frame(compiled_results)

#removendo os NA de PPI
compiled_results$PPI_degree[is.na(compiled_results$PPI_degree)] <- 0

cat(paste0(sum(is.na(compiled_results$HedgesG_val)), " genes não foram identificados na validação.\n"))

#removendo genes ausentes na validação
compiled_results <- compiled_results[!is.na(compiled_results$HedgesG_val), ]

# ============== PONTUAÇÃO ==============
# fase final da análise. Baseado nas análises feitas como critérios cada gene será
# pontuado num ranking de melhores candidatos a biomarcadores.

message("===========================================================================")
message("IDENTIFCANDO BIOMARCADORES...")
message("===========================================================================")

# Função para normalizar variáveis de cada critério para escala 0–1
safe_rescale <- function(x, to = c(0,1)) {
  # se todos forem NA
  if(all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  # intervalo dos dados
  r <- range(x, na.rm = TRUE)
  # evita divisão por zero
  # caso todos valores sejam iguais
  if(isTRUE(all.equal(r[1], r[2]))) {
    return(rep(mean(to), length(x)))
  }
  # normalização
  rescale(x, to = to, from = r)
}

# função para reduzir impacto de outliers extremos.
clip_quant <- function(x,
                       probs = c(0.05, 0.95)) {
  # calcula percentis
  q <- quantile(
    x,
    probs = probs,
    na.rm = TRUE,
    names = FALSE
  )
  # limita os valores
  pmin(pmax(x, q[1]), q[2])
}


## FUNÇÃO DE PONTUAÇÃO
score_biomarkers <- function(df) {
  df %>%
    mutate(
      
      # CRITÉRIO 1: mangnitude biológica (g de hedges)
      # priorização de up-regulados com pmax(LogFC, 0)
      # Genes negativos viram 0
      Hedges_g_pos = pmax(Hedges_g, 0),
      HedgesG_val_pos = pmax(HedgesG_val, 0),
      
      # CRITÉRIO 2: valor p ajustado por FDR
      # valores menores refletem melhor evidência estatística
      # usou-se -log10(FDR)
      # 0.01   -> 2
      # 0.001  -> 3
      # 1e-20  -> 20
      meta_fdr_score = -log10(pmax(FDR, 1e-300)),
      val_fdr_score = -log10(pmax(FDR_val, 1e-300)),
      
      # NORMALIZAÇÃO DAS VARIÁVEIS
      # aplicação da função em cada critério
      # Reduz outliers e converte para escala 0-1

      # magnitude biológica GEO
      meta_effect_n = safe_rescale( clip_quant(Hedges_g_pos)),
      
      # magnitude biológica TCGA
      val_effect_n = safe_rescale(clip_quant(HedgesG_val_pos)),
      
      # significância estatística GEO
      meta_fdr_n = safe_rescale(clip_quant(meta_fdr_score)),
      
      # significância estatística TCGA
      val_fdr_n = safe_rescale(clip_quant(val_fdr_score)),
      
      # AUC diagnóstica
      auc_n = safe_rescale(clip_quant(auc)),
      
      # consistência entre datasets
      significance_n = safe_rescale(clip_quant(significance_score)),
      
      # conectividade em PPI
      ppi_n = safe_rescale(clip_quant(PPI_degree)),
      
      # heterogeneidade
      i2_n = safe_rescale(clip_quant(I2)),
      
      # CRITÉRIO 3: robustez entre estudos na meta-análise
      # Definido por I² que reflete a heterogeneidade entre datasets
      # valores menores são mais homogêneros
      # I² alto -> score baixo
      # I² baixo -> score alto
      robustness_n = 1 - i2_n,
      
      
      # CRITÉRIO 4: Concordância direcional
      # Analisa se os genes são up ou down regulados na análise e validação ao
      # mesmo tempo
      # +1 = mesma direção
      #  0 = direção oposta
      direction_bonus = ifelse(sign(Hedges_g) == sign(HedgesG_val),1,0),
      
      
      # SCORE FINAL
      # Os pesos podem ser ajustados
      biomarker_score =
        100 * (
            # magnitude GEO
            0.14 * meta_effect_n +
            # magnitude TCGA
            0.14 * val_effect_n +
            # significância GEO
            0.12 * meta_fdr_n +
            # significância TCGA
            0.12 * val_fdr_n +
            # consistência entre datasets
            0.12 * significance_n +
            # AUC diagnóstica
            0.20 * auc_n +
            # baixa heterogeneidade
            0.10 * robustness_n +
            # contexto biológico PPI
            0.03 * ppi_n +
            # concordância GEO ↔ TCGA
            0.03 * direction_bonus
          
        )) %>%
    
  # ordenba os resultados pelo score
  arrange(desc(biomarker_score))}

# execução da pontuação
ranked_results <- score_biomarkers(compiled_results)

# VISUALIZAÇÃO DOS TOP GENES
head(
  ranked_results[, c(
    "Gene_Symbol",
    "biomarker_score",
    "Hedges_g",
    "FDR",
    "I2",
    "significance_score",
    "HedgesG_val",
    "FDR_val",
    "auc",
    "PPI_degree"
  )]
)

## checando se o diretório existe
out_dir <- file.path(results_dir, "biomarker_results")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

ranked_results_clean <- ranked_results %>%
  select(
    Gene_Symbol,
    Entrez,
    meta_effect_n,
    val_effect_n,
    meta_fdr_n,
    val_fdr_n,
    auc_n,
    significance_n,
    ppi_n,
    robustness_n,
    direction_bonus,
    biomarker_score
  )

message("=======================================")
message("EXPORTANDO RESULTADOS DE BIOMARCADOERES")
message("=======================================")

write.csv(ranked_results_clean,file.path(results_dir,"biomarker_results","Bladder_cancer_biomarker_rank.csv"))

# ===== SESSION INFO =====
writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "sessionInfo.txt")
)

gc()


# ============== Gráficos de resultados ==============
# Gráficos dos resultados obtidos na análise

message("===========================================================================")
message("GERANDO GRÁFICOS DOS RESULTADOS!")
message("===========================================================================")

## heatmap de resultados
# representação gráfica de pontuação dos genes

# compilando a pontuação
ranked_results_points <- ranked_results %>%
  select(
    Gene_Symbol,
    meta_effect_n,
    val_effect_n,
    meta_fdr_n,
    val_fdr_n,
    auc_n,
    significance_n,
    ppi_n,
    robustness_n,
    direction_bonus,
    biomarker_score
  )

# organizando os dados
# renomeando colunas
ranked_results_points <- ranked_results_points %>%
  rename(
    Hedges_G = meta_effect_n,
    HedgesG_val = val_effect_n,
    FDR = meta_fdr_n,
    FDR_val = val_fdr_n,
    AUC = auc_n,
    consistencia = significance_n,
    ppi = ppi_n,
    I2 = robustness_n,
    concordancia_direcional = direction_bonus
  )

# reordenando colunas
ranked_results_points <- ranked_results_points %>%
  select(Gene_Symbol,Hedges_G,FDR,I2,consistencia,ppi,HedgesG_val,
         FDR_val,concordancia_direcional,AUC,biomarker_score)

# dividindo os pontos de biomarcador por 100 para ficarem na escala 0-1
ranked_results_points$biomarker_score <- ranked_results_points$biomarker_score / 100

#filtrando o gráfico para melhor visualização
ranked_results_points_filtered <- ranked_results_points[(1:100),]

#criando heatmap de visualização de dados
heatmap_results <- ranked_results_points_filtered %>%
  rename(
    "g de Hedges" = Hedges_G,
    "g de Hedges validação" = HedgesG_val,
    "FDR validação" = FDR_val,
    "I²" = I2,
    "concordãncia" = concordancia_direcional,
    "pontuação final" = biomarker_score
  ) %>%
  column_to_rownames("Gene_Symbol") %>%
  as.matrix()
heatmap_results <- t(heatmap_results)

png(
  file.path("figures","biomarker_score_heatmap.png"),
  width = 4200,
  height = 1400,
  res = 300
)
Heatmap(
  heatmap_results,
  name = "Pontuação",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9),
  col = colorRamp2(
    c(0, 0.5, 1),
    c("#D73027", "#FFFFBF", "#1A9850")
  )
)
dev.off()


## Heatmap de expressão
# foi feito um gráfico de heatmap para os 89 DEGs e top 20 genes candidatos a
# biomarcador por pontuação.

## organizando amostras
# ordenar agrupando non-tumor e tumor
sample_order <- order(group_roc)

# aplicar ordenação
exprs_ordered <- validation_exprs[, sample_order]
group_ordered <- group_roc[sample_order]


## HEATMAP DOS GENES FILTRADOS HOMOGÊNEOS

# genes em ENTREZ
genes_filtered <- as.character(DEGs_filtered$Gene)

# filtrar matriz
heatmap_filtered <- exprs_ordered[rownames(exprs_ordered) %in% genes_filtered,]

# ordenar linhas igual DEGs_filtered
heatmap_filtered <- heatmap_filtered[match(genes_filtered, rownames(heatmap_filtered)),]

# substituir ENTREZ por símbolo gênico
rownames(heatmap_filtered) <- DEGs_filtered$Symbol

# Z-score por gene
heatmap_filtered_scaled <- t(scale(t(heatmap_filtered)))

# remover possíveis NAs
heatmap_filtered_scaled <- heatmap_filtered_scaled[complete.cases(heatmap_filtered_scaled),]

# anotação de grupos
ha_filtered <- HeatmapAnnotation(
  Group = group_ordered,
  col = list(Group = c("non_tumor" = "#3A7D44", "tumor" = "#7B2CBF")))

# salvar figura
png(file.path("figures","heatmap_filtered_genes.png"),
    width = 3200,height = 4200,res = 400)

Heatmap(
  heatmap_filtered_scaled,
  name = "Z-score",
  top_annotation = ha_filtered,
  cluster_rows = TRUE,
  # NÃO clusterizar amostras
  cluster_columns = FALSE,
  # separar grupos visualmente
  column_split = group_ordered,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "TCGA Tumor vs GTEx Healthy",
  heatmap_legend_param = list(
    title = "Expression"
  ),
  col = colorRamp2(
    c(-2, 0, 2),
    c("#2C7BB6", "#F0F0F0", "#D41159")
  )
)

dev.off()

## Heatmap dos top 20 biomarcadores

# top 20 genes
top20 <- ranked_results_clean$Entrez[1:20] %>%as.character()

# filtrar matriz
heatmap_top20 <- exprs_ordered[rownames(exprs_ordered) %in% top20,]

# ordenar linhas
heatmap_top20 <- heatmap_top20[match(top20, rownames(heatmap_top20)),]

# símbolos gênicos
rownames(heatmap_top20) <- ranked_results_clean$Gene_Symbol[1:20]

# Z-score
heatmap_top20_scaled <- t(
  scale(t(heatmap_top20)))

# remover NAs
heatmap_top20_scaled <- heatmap_top20_scaled[complete.cases(heatmap_top20_scaled),]

# anotação
ha_top20 <- HeatmapAnnotation(
  Group = group_ordered,
  col = list(Group = c("non_tumor" = "#3A7D44", "tumor" = "#7B2CBF")))

# salvar figura
png(file.path("figures","heatmap_top20_genes.png"),
    width = 2600,height = 2200,res = 400)

Heatmap(
  heatmap_top20_scaled,
  name = "Z-score",
  top_annotation = ha_top20,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = group_ordered,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_title = "Top 20 Biomarkers",
  heatmap_legend_param = list(
    title = "Expression"
  ),
  col = colorRamp2(
    c(-2, 0, 2),
    c("#2C7BB6", "#F0F0F0", "#D41159")))

dev.off()

## limpando
rm(
  sample_order,
  exprs_ordered,
  genes_filtered,
  heatmap_filtered,
  top20,
  heatmap_top20
)
gc()


message("===========================================================================")
message("ANÁLISE FINALIZADA!")
message("===========================================================================")

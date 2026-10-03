# ============== CARREGANDO PACOTES ==============

library(here) #pacote "aqui" auxilia com diretórios
library(R.utils) #ferramentas do R, permite usar gunzip
library(GEOquery) #busca e download de datasets do Gene Omnibus
library(affy) #pacote de normalização affymetrix
library(oligo) #pacote de normalização oligo
library(limma) #análise estatística e normalização
library(AnnotationDbi) # pacote de execução de anotação
library(hgu133plus2.db) #pacote de base de anotação
library(hta20transcriptcluster.db) #pacote de base de anotação
library(hgu133a.db) #pacote de base de anotação
library(hgu133acdf) #pacote de base de anotação
library(hgu133acdf) #pacote de base de anotação
library(illuminaHumanv2.db) #pacote de base de anotação
library(illuminaHumanv4.db) #pacote de base de anotação
library(affyPLM) # controle de qualidade de arrays affy
library(ggplot2) #pacote de gráficos de expressão
library(dplyr) #gerenciamento de dataframes
library(tibble) #gerenciamento de dataframes
library(tidyr) #gerenciamento de dataframes
library(sva) #comBat
library(patchwork) # Para organizar os plots 2D lado a lado
library(plotly)    # Para gráficos 3D interativos
library(htmlwidgets) #para salvar gráficos 3D
library(clusterProfiler) #enriquecimento funcional (GO e KEGG)
library(ReactomePA) #enriquecimento funcional (Reactome)
library(enrichplot) #pacote de gráficos de enriquecimento
library(STRINGdb) #pacote de rede PPi
library(igraph) #complemento do grafico de rede PPI
library(recount3) #pacote da base de dados recount3 para validação
library(edgeR) #pacote para validação
library(ggrepel)
library(pheatmap) # heatmap
library(circlize) #complemento de heatmap
library(glmnet) # Machine leaning
library(caret) # Para divisão dos dados e matriz de confusão
library(pROC) # Para curva ROC e AUC

# complementos a outros pacotes
library(Biobase)
library(BiocGenerics)
library(generics)
library(stats4)
library(IRanges)
library(S4Vectors)
library(R.oo)
library(R.methodsS3)
library(oligoClasses)
library(Biostrings)
library(DBI)
library(RSQLite)

# ==================================================
# Opções
# ==================================================

options(stringsAsFactors = FALSE)

set.seed(123)
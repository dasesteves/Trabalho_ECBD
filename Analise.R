# explicação dos dados na pasta "dataset" ou "https://www.cbioportal.org/study/summary?id=prad_tcga_pan_can_atlas_2018", sua origem e relevância (alguém pode fazer este?)
# - tarefas de preparação e de pré-processamento dos dados (a executar)
# - sumarização dos dados (estatística descritiva, exploração com recurso a gráficos)
# - análise estatística univariada; análise de expressão diferencial e de enriquecimento

# Carregar os dados
dados <- read.delim("dataset/prad_tcga_pan_can_atlas_2018_clinical_data.tsv")
summary(dados)

# Sumarização dos dados
summary(dados_limpos) # Estatísticas descritivas
hist(dados_limpos$Idade) # Histograma da idade dos pacientes

#'dados_extracao' é para criar o dataframe já com apenas as colunas de interesse 
dados_extracao <- dados_limpos[, colunas_geneticas] #'colunas_geneticas' são para substituir apenas pelas colunas de interesse, 
# podemos no entanto trabalhar diretamente com o dataset se especificarmos as colunas ou as criarmos para as diferentes analises

# Imputação pela média para variáveis contínuas
dados_limpos <- dados
for(i in 1:ncol(dados_limpos)) {
  if(is.numeric(dados_limpos[,i])) {
    dados_limpos[is.na(dados_limpos[,i]), i] <- mean(dados_limpos[,i], na.rm = TRUE)
  }
}

# Visualização com heatmap (exemplo)
pheatmap::pheatmap(dados_extracao)

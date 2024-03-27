# explicação dos dados na pasta "dataset" ou "https://www.cbioportal.org/study/summary?id=prad_tcga_pan_can_atlas_2018", sua origem e relevância (alguém pode fazer este?)
# - tarefas de preparação e de pré-processamento dos dados (a executar)
# - sumarização dos dados (estatística descritiva, exploração com recurso a gráficos)
# - análise estatística univariada; análise de expressão diferencial e de enriquecimento


# Carregar os dados
dados <- read.delim("caminho/para/o/arquivo.tsv")

# Imputação pela média para variáveis contínuas
dados_limpos <- dados
for(i in 1:ncol(dados_limpos)) {
  if(is.numeric(dados_limpos[,i])) {
    dados_limpos[is.na(dados_limpos[,i]), i] <- mean(dados_limpos[,i], na.rm = TRUE)
  }
}

# Sumarização dos dados
summary(dados_limpos) # Estatísticas descritivas
hist(dados_limpos$Idade) # Histograma da idade dos pacientes

# Supondo que 'dados_limpos' é o seu dataframe e 'colunas_geneticas' são as colunas de interesse
dados_geneticos <- dados_limpos[, colunas_geneticas]

# Visualização com heatmap (exemplo)
pheatmap::pheatmap(dados_geneticos)
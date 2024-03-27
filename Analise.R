# explicação dos dados na pasta "dataset" ou "https://www.cbioportal.org/study/summary?id=prad_tcga_pan_can_atlas_2018", sua origem e relevância
# - tarefas de preparação e de pré-processamento dos dados
# - sumarização dos dados (estatística descritiva, exploração com recurso a gráficos)
# - análise estatística univariada; análise de expressão diferencial e de enriquecimento

# Carregar os dados
dados <- read.delim("caminho/para/o/arquivo.tsv")

# Sumarização dos dados
summary(dados) # Estatísticas descritivas
hist(dados$Idade) # Histograma da idade dos pacientes

# Supondo que 'dados_limpos' é o seu dataframe e 'colunas_geneticas' são as colunas de interesse
dados_genetica <- dados_limpos[, colunas_geneticas]

# Visualização com heatmap (exemplo)
pheatmap::pheatmap(dados_geneticos)
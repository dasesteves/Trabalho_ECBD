# explicação dos dados na pasta "dataset" ou "https://www.cbioportal.org/study/summary?id=prad_tcga_pan_can_atlas_2018", sua origem e relevância
# - tarefas de preparação e de pré-processamento dos dados
# - sumarização dos dados (estatística descritiva, exploração com recurso a gráficos)
# - análise estatística univariada; análise de expressão diferencial e de enriquecimento


# Carregar os dados
dados <- read.delim("caminho/para/o/arquivo.tsv")

# Limpeza dos dados
dados_limpos <- na.omit(dados) # Remove linhas com valores NA

# Sumarização dos dados
summary(dados_limpos) # Estatísticas descritivas
hist(dados_limpos$Idade) # Histograma da idade dos pacientes

# Análise estatística univariada
t.test(Idade ~ Status, data = dados_limpos) # Teste t para comparar a idade entre dois grupos de status

# Análise de expressão diferencial (exemplo genérico)
# Supondo que 'expressao' é um dataframe com dados de expressão gênica
resultados_de <- DESeq2::DESeqDataSetFromMatrix(countData = expressao, ...)
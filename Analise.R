# explicação dos dados na pasta "dataset" ou "https://www.cbioportal.org/study/summary?id=prad_tcga_pan_can_atlas_2018", sua origem e relevância (alguém pode fazer este?)

Os dados que vão ser analisados neste trabalho provêm do projeto Pan-Cancer Atlas (TCGA), este projeto foi uma colaboração em larga escala de diversos investigadores de todo o mundo apoiado pelo National Cancer Institute (NCI) e pelo National Human Genome Research Institute (NHGRI).
Neste trabalho serão analisadas 494 amostras obtidas através do perfil molecular de mais de 11.000 tumores de 33 diferentes tipos de cancro, as obtenções destes dados foram realizadas através da utilização de diversas técnicas que permitiram examinar as alterações moleculares em diferentes níveis, quer no DNA, RNA, na proteína e na epigenética.
A relevância desses dados consiste na compreensão da origem, de onde, e de como os tumores se desenvolvem em humanos, a análise dos perfis moleculares de diferentes tipos de cancro terá como objetivo identificar padrões comuns, diferenças e temas emergentes entre eles, estes dados serão fundamentais para o desenvolvimento de tratamentos mais eficazes e personalizados ,para além disso, os dados servem como uma importante referência para futuras pesquisas sobre o cancro, fornecendo uma base sólida para futuros estudos e análises.

# - tarefas de preparação e de pré-processamento dos dados (a executar)
# - sumarização dos dados (estatística descritiva, exploração com recurso a gráficos)
# - análise estatística univariada; análise de expressão diferencial e de enriquecimento

# Carregar os dados
dados <- read.delim("dataset/prad_tcga_pan_can_atlas_2018_clinical_data.tsv")
summary(dados)

# Sumarização dos dados
summary(dados_limpos) # Estatísticas descritivas
hist(dados_limpos$Idade) # Histograma da idade dos pacientes

#'dados_extracao' é para criar o dataframe já com apenas as colunas de interesse, podemos usar o pandas para criar esse dataset (envolve outro código e import)
dados_extracao <- dados_limpos[, colunas_geneticas] #'colunas_geneticas' são para substituir apenas pelas colunas de interesse, 
# podemos no entanto trabalhar diretamente com o dataset se especificarmos as colunas ou as criarmos para as diferentes analises
# Study ID	Patient ID	Sample ID	Diagnosis Age	Disease Stage Cancer Type	Cancer Type Detailed	Ethnicity Category	Fraction Genome Altered	Genetic Ancestry Label	ICD-10 Classification	In PanCan Pathway Analysis	MSI MANTIS Score	MSIsensor Score	Mutation Count	Oncotree Code	Progress Free Survival (Months)	Number of Samples Per Patient	Sample Type	Sex	Somatic Status	Subtype	Tissue Prospective Collection Indicator	Tissue Retrospective Collection Indicator Tissue Source Site Tissue Source Site Code	TMB (nonsynonymous)	Tumor Disease Anatomic Site	Tumor Type	Patient Weight	Winter Hypoxia Score


# Imputação pela média para variáveis contínuas
dados_limpos <- dados
for(i in 1:ncol(dados_limpos)) {
  if(is.numeric(dados_limpos[,i])) {
    dados_limpos[is.na(dados_limpos[,i]), i] <- mean(dados_limpos[,i], na.rm = TRUE)
  }
}

# Visualização com heatmap (exemplo)
pheatmap::pheatmap(dados_extracao)

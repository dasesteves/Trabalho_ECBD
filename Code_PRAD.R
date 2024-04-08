# Este script orienta uma pipeline de análise de dados transcriptómicos do projeto TCGA-PRAD, incluindo recuperação de dados, pré-processamento,
# análise exploratória e análise de expressão diferencial usando o método DESeq2.

# Mudar o diretório de trabalho conforme necessário
#setwd("C:/Users/dases/Desktop/Trabalho_ECDB-main/Trabalho_ECDB-main") # Alterar 

# Instalar o pacote BiocManager caso ele não esteja já instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar o pacote TCGAbiolinks caso ele não esteja já instalado
if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")

BiocManager::install("SummarizedExperiment", dependencies = TRUE)

# Carregar os pacotes na sessão atual do R
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)

# Criar uma consulta ao Genomic Data Commons (GDC) e obter dados transcriptómicos do projeto TCGA sobre câncer da próstata (PRAD)
# Dados: https://portal.gdc.cancer.gov/projects/TCGA-PRAD
projeto <- "TCGA-PRAD"
query <- GDCquery(
  project = projeto,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Baixar os dados do Genomic Data Commons (GDC) com base nas especificações de uma consulta criada anteriormente
GDCdownload(query = query)
# Preparar os dados do Genomic Data Commons (GDC) para análise
rna_seq_PRAD <- GDCprepare(query = query, save = TRUE, save.filename = "mRNA_TCGA-PRAD.rda")

# Correção para carregar os dados corretamente
rna_seq_PRAD <- get(load("mRNA_TCGA-PRAD.rda"))

# Visualizar informações básicas do conjunto de dados
metadata(rna_seq_PRAD)
class(rna_seq_PRAD)

# Retornar as dimensões do objeto (matriz ou um dataframe; objeto mais complexo)
dim(rna_seq_PRAD)
rownames(rna_seq_PRAD)
colnames(rna_seq_PRAD)

# Extrair metadados relacionados com as linhas
row_metadados <- SummarizedExperiment::rowData(rna_seq_PRAD)[1:10, c('type', 'gene_type', 'gene_name')]

# Atribuir a um novo objeto chamado meta_PRAD, os metadados associados aos dados de expressão
amostras_metadados <- colData(rna_seq_PRAD)
class(amostras_metadados)
dim(amostras_metadados)
colnames(amostras_metadados)

# Código utilizado para extrair a informação relacionada à contagem da expressão dos genes do objeto rna_seq_PRAD
geneExp <- SummarizedExperiment::assay(rna_seq_PRAD)
# Para ver os outliers - no nosso caso é zero
pre <- TCGAanalyze_Preprocessing(rna_seq_PRAD) # Faz correlação e tenta identificar outliers
dim(geneExp)
dim(pre)

# Nota: A normalização pode ser realizada com o pacote DESeq2 e com a função de normalização
# Nota: o edgeR pode fazer outro processo de normalização através da função RPKM 
# (tem que se fazer a normalização ou antes ou depois da análise diferencial)
meta_PRAD = colData(rna_seq_PRAD)
meta_PRAD = as.data.frame(meta_PRAD)
#retorna as dimensões dos metadados 
dim(meta_PRAD)
#Ver nomes das colunas
colnames(meta_PRAD)
#nomes das linhas
names(meta_PRAD)
#extrair componentes de um objeto por nome (através de colunas)
meta_PRAD$patient

# Criar vetor lógico em que diz que linhas possuem grade diferente de NA
coluna <- !is.na(meta_PRAD$definition)

# Selecionar apenas as linhas que têm grade
meta <- meta_PRAD[coluna, c("definition", "age_at_diagnosis", "disease_type")]
meta$disease_type
table(meta_PRAD[rownames(meta), "paper_tumor_grade"]) # Confirmar que os pacientes selecionados têm um valor de grade atribuído na tabela de metadados original
head(meta) # Algo de errado se passa

summary(meta_PRAD$age_at_diagnosis) # Idade está em dias
boxplot(meta_PRAD$age_at_diagnosis, horizontal = TRUE)

anova <- aov(meta_PRAD$age_at_diagnosis ~ meta_PRAD$definition)
summary(anova)

TukeyHSD(anova)

boxplot(meta_PRAD$age_at_diagnosis ~ meta_PRAD$definition, horizontal = TRUE)
# Verificação da dimensão do dataframe meta após seleção de apenas as linhas com informação relativa ao grade
sum(table(meta_PRAD$definition))
dim(meta_PRAD)

# Criação do dataframe de expressão apenas para os pacientes com informação relativa ao grade
exp_grade <- geneExp[, rownames(meta)]
dim(exp_grade)

# Nomes das linhas
names(meta_PRAD)
# Extrair componentes de um objeto por nome (através de colunas)
meta_PRAD$patient
meta_PRAD$initial_weight

# Processamento de Dados Clínicos com TCGAbiolinks
# Obtém dados clínicos do projeto TCGA-PRAD (câncer da próstata),
# especificamente suplementos clínicos, através de GDCquery
query_clin <- GDCquery(project = "TCGA-PRAD", 
                       data.category = "Clinical",
                       data.type = "Clinical Supplement", 
                       data.format = "BCR Biotab")

# Baixa e prepara os dados clínicos
GDCdownload(query = query_clin)
clinical.PRAD <- GDCprepare(query = query_clin, save = TRUE, save.filename = "clinical_data_PRAD.rda")
# Listar os nomes dos componentes do objeto
names(clinical.PRAD)

head(clinical.PRAD$clinical_drug_prad)
# Converte os dados contidos em clinical.PRAD$clinical_patient_PRAD para um dataframe
df <- as.data.frame(clinical.PRAD$clinical_patient_prad)
# Utilize View(df) no RStudio para visualizar o dataframe

# Análise exploratória
# Pré-processamento e filtragem
# Retirar as colunas dos metadados onde havia mais de 50 elementos como: 
# “not reported/Not Reported/Reported” e/ou “NA”

cols_with_not_reported <- which(sapply(meta_PRAD, function(x) sum(x == "not reported", na.rm = TRUE)) > 50)
cols_with_Not_Reported <- which(sapply(meta_PRAD, function(x) sum(x == "Not Reported", na.rm = TRUE)) > 50)
cols_with_NA <- which(sapply(meta_PRAD, function(x) sum(is.na(x))) > 60)
# Remover as colunas baseadas nos critérios específicos de cima
metadata_matriz_clean <- meta_PRAD[, -c(cols_with_not_reported, cols_with_Not_Reported, cols_with_NA)]
dim(metadata_matriz_clean)

# Código utilizado para extrair a informação relacionada à contagem da expressão dos genes do objeto rna_seq_PRAD
geneExp <- SummarizedExperiment::assay(rna_seq_PRAD)

summary(as.data.frame(metadata_matriz_clean$age_at_diagnosis))
summary(as.data.frame(metadata_matriz_clean$age_at_index))
boxplot(as.data.frame(metadata_matriz_clean$age_at_index),horizontal=T)
hist(metadata_matriz_clean$age_at_index)

##########

# Análise de Expressão Diferencial
library(DESeq2) # DESeq2, uma ferramenta para análise de expressão diferencial de dados de contagem de sequenciamento de RNA (RNA-Seq)

# Filtra os dados de RNA para incluir apenas amostras não nulas
data_de <- rna_seq_PRAD[, !is.na(rna_seq_PRAD$colData$vital_status)]
table(as.data.frame(metadata_matriz_clean$vital_status))

countData <- assays(data_de)
# Cria um objeto DESeqDataSet para análise, especificando um design experimental que compara o status vital
ddsSE <- DESeqDataSetFromMatrix(countData = countData, 
                                colData = colData(data_de), 
                                design = ~ vital_status)
# Filtragem de Genes: Remove genes com contagens baixas (menos de 10) para melhorar a confiabilidade da análise de expressão diferencial.
keep <- rowSums(counts(ddsSE)) >= 10
ddsSE <- ddsSE[keep,]

# Executa a análise de expressão diferencial com a função DESeq
ddsSE <- DESeq(ddsSE)

# Visualizar nomes dos resultados disponíveis
resultsNames(ddsSE)

# Obter e visualizar resultados da comparação de interesse
res <- results(ddsSE)
dea <- as.data.frame(res)

# Resumir os resultados para obter uma visão geral dos achados estatísticos, como o número de genes significativamente diferencialmente expressos
summary(res)

#Este script orienta uma pipeline de análise de dados transcriptómicos do projeto TCGA-PRAD, incluindo recuperação de dados, pré-processamento,
#análise exploratória e análise de expressão diferencial usadno o método DESeq2.

#set wd (alterar para o nosso caso)


#instalar o pacote BiocManager caso ele não esteja já instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#instalar o pacote TCGAbiolinks caso ele não esteja já instalado
if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")

BiocManager::install("SummarizedExperiment", dependencies = TRUE)


#carregar os pacotes na sessão atual do R
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)


#criar uma consulta ao Genomic Data Commons (GDC) e obter dados transcriptómicos do projeto TCGA sobre dados do cancro da prostata (PRAD)
#dados : https://portal.gdc.cancer.gov/v1/projects/TCGA-PRAD
projeto <- "TCGA-PRAD"
query <- GDCquery(
  project = projeto,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

#baixar os dados do Genomic Data Commons (GDC) com base nas especificações 
#de uma consulta criada anteriormente
GDCdownload(query=query)
#preparar os dados baixados do Genomic Data Commons (GDC) para análise
rna_seq_PRAD  <- GDCprepare(query = query, save = TRUE, save.filename = "mRNA_TCGA-PRAD.rda")

#retorna a classe do objeto (tipo de dados ou a estrutura de dados que o objeto 
#representa, o que por sua vez determina quais funções podem ser aplicadas nele)
class(rna_seq_PRAD)
#retorna as dimensões do objeto(matriz ou um dataframe;objeto mais complexo)
dim(rna_seq_PRAD)

#código utilizado para extrair a informação relacionada à contagem da expressão dos genes do objeto rna_seq_PRAD
geneExp <- SummarizedExperiment::assay(rna_seq_PRAD)
#para ver os outliers - no nosso caso é zero
pre = TCGAanalyze_Preprocessing(rna_seq_PRAD) #faz correlação e tenta identificar outliers
dim(geneExp)
dim(pre)

#NOTA: a normalização pode ser realizada com o package deseq2 e com a função normalization
#NOTA: o edger pode fazer outro processo de normalização através da função rpkm 
#(tem que se fazer a normalização ou antes ou depois da analise diferencial)

#Ver nomes das colunas
colnames(rna_seq_PRAD)
#nomes das linhas
names(rna_seq_PRAD)


#atribui a um novo objeto chamado meta_PRAD, os metadados associados ao conjunto 
#de dados rna_seq_PRAD
meta_PRAD = colData(rna_seq_PRAD)
#retorna as dimensões dos metadados 
dim(meta_PRAD)
#Ver nomes das colunas
colnames(meta_PRAD)
#nomes das linhas
names(meta_PRAD)
#extrair componentes de um objeto por nome (através de colunas)
meta_PRAD$patient
meta_PRAD$paper_vital_status


##################################
#Processamento de Dados Clínicos com TCGAbiolinks
## dados clinicos

#buscar dados clínicos do projeto TCGA-PRAD (cancro da prostata), 
#especificamente suplementos clínicos, através de GDCquery
query_clin <- GDCquery(project = "TCGA-PRAD", 
                       data.category = "Clinical",
                       data.type = "Clinical Supplement", 
                       data.format = "BCR Biotab")
#Baixa e prepara os dados clínicos
GDCdownload(query = query_clin)
clinical.PRAD <- GDCprepare(query = query_clin, save = TRUE, save.filename = "clinical_data_PRAD.rda")
#listar os nomes dos componentes do objeto
names(clinical.PRAD)

#exibir as primeiras linhas da coluna ou lista clinical_drug_PRAD contida dentro do objeto clinical.PRAD
head (clinical.PRAD$clinical_drug_PRAD)
#converte os dados contidos em clinical.PRAD$clinical_patient_PRAD para um 
#dataframe e atribui esse novo dataframe à variável df
df = as.data.frame(clinical.PRAD$clinical_patient_PRAD)
View(df)

#######################################
#Análise exploratória
#Pré-processamento e filtragem
#retirar as colunas dos metadados onde havia mais de 60 elementos como: 
#“not/Not reported/Reported” e/ou “NA”

cols_with_not_reported <- which(sapply(meta_PRAD,function(x) sum(x == "not reported", na.rm = TRUE)) > 60)
cols_with_Not_Reported <- which(sapply(meta_PRAD,function(x) sum(x == "Not Reported", na.rm = TRUE)) > 60)
cols_with_NA <- which(sapply(meta_PRAD, function(x) sum(is.na(x))) > 60)
# remover as colunas baseadas nos critérios específicos de cima
metadata_matriz_clean <- meta_PRAD[, -c(cols_with_not_reported, cols_with_Not_Reported, cols_with_NA)] 
dim(metadata_matriz_clean)


#código utilizado para extrair a informação relacionada à contagem da expressão dos genes do objeto rna_seq_PRAD
geneExp <- SummarizedExperiment::assay(rna_seq_PRAD)



##########

#Análise de Expressão Diferencial com DESeq2
library(DESeq2)  #Nota: pacote DESeq2, uma ferramenta para análise de expressão diferencial de dados de contagem de sequenciamento de RNA (RNA-Seq)

#Filtra os dados de RNA para incluir apenas amostras não nulo
data_de <- rna_seq_PRAD[,!is.na(rna_seq_PRAD$paper_vital_status)]

#Cria um objeto DESeqDataSet para análise, especificando um design experimental 
#que compara o status de IDH
ddsSE <- DESeqDataSet(data_de, design = ~ paper_vital_status)

#Filtragem de Genes: Remove genes com contagens baixas (menos de 10) 
#para melhorar a confiabilidade da análise de expressão diferencial.
keep <- rowSums(counts(ddsSE)) >= 10

#Executa a análise de expressão diferencial com a função DESeq
ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)

resultsNames(ddsSE)
#comparação "WT vs Mutant" para o status do IDH, e converte os resultados para um dataframe
res <- results(ddsSE, name = "paper_IDH.status_WT_vs_Mutant")
dea <- as.data.frame(res)
#resume os resultados para obter uma visão geral dos achados estatísticos, 
#como o número de genes significativamente diferencialmente expressos
summary(res)

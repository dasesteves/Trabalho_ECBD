#Este script orienta uma pipeline de análise de dados transcriptómicos do projeto TCGA-PRAD, incluindo recuperação de dados, pré-processamento,
#análise exploratória e análise de expressão diferencial usadno o método DESeq2.

#setwd("C:/Users/dases/Desktop/Trabalho_ECDB-main/Trabalho_ECDB-main") (alterar para o nosso caso)


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

#baixa os dados do Genomic Data Commons (GDC) com base nas especificações de uma consulta criada anteriorm
GDCdownload(query=query)
#prepara os dados do Genomic Data Commons (GDC) para análise
rna_seq_PRAD  <- GDCprepare(query = query, save = TRUE, save.filename = "mRNA_TCGA-PRAD.rda")
metadata(rna_seq_PRAD)
rna_seq_PRAD= get(load("mRNA_TCGA-Prad.rda"))
class(rna_seq_PRAD)

#retorna as dimensões do objeto(matriz ou um dataframe;objeto mais complexo)
dim(rna_seq_PRAD)
rownames(rna_seq_PRAD)
colnames(rna_seq_PRAD)

row_metadados=SummarizedExperiment::rowData(rna_seq_PRAD);row_metadados[1:10,c('type','gene_type','gene_name')] # extarct metados relacionados com as linhas

#atribui a um novo objeto chamado meta_PRAD, os metadados associados aos dados de expressão
amostras_metadados = colData(rna_seq_PRAD)
class(amostras_metadados)
dim(amostras_metadados)
colnames(amostras_metadados)


#código utilizado para extrair a informação relacionada à contagem da expressão dos genes do objeto rna_seq_PRAD
geneExp <- SummarizedExperiment::assay(rna_seq_PRAD)
#para ver os outliers - no nosso caso é zero
pre = TCGAanalyze_Preprocessing(rna_seq_PRAD) #faz correlação e tenta identificar outliers
dim(geneExp)
dim(pre)


#NOTA: a normalização pode ser realizada com o package deseq2 e com a função normalization
#NOTA: o edger pode fazer outro processo de normalização através da função rpkm 
#(tem que se fazer a normalização ou antes ou depois da analise diferencial)

#atribui a um novo objeto chamado meta_PRAD, os metadados associados ao conjunto de dados rna_seq_PRAD
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



# criar vetor logica em que diz que linhas é que possuem grade diferente de na
coluna=!is.na(meta_PRAD$definition) 


# seleciona apenas as linhas que tem grade
meta=meta_PRAD[coluna,c("definition","age_at_diagnosis","disease_type")] 
meta$disease_type
table(meta_PRAD[rownames(meta),"paper_tumor_grade"]) # confirmar que os pacientes selecionados tem um valor de grade atribuido na tabela de metadados original
head(meta) #algo de errado se passa

summary(meta_PRAD$age_at_diagnosis) #idade está em dias
boxplot(meta_PRAD$age_at_diagnosis,horizontal=T)

anova=aov(meta_PRAD$age_at_diagnosis~meta_PRAD$definition)
summary(anova)

TukeyHSD(anova)

boxplot(meta_PRAD$age_at_diagnosis~meta_PRAD$definition,horizontal=T)
#verificaçáo da dimensão do dataframe meta após selação de apenas as linhas com informação relativa ao grade
sum(table(meta_PRAD$definition))
dim(meta_PRAD)


#criação do data frame de expressão apenas para os pacientes com informação relativa ao grade
exp_grade=geneExp[,rownames(meta)]
dim(exp_grade)



row.names(meta_PRAD[,"definition"])
select=meta_PRAD$definition
row_names=rownames(select)

#nomes das linhas
names(meta_PRAD)
#extrair componentes de um objeto por nome (através de colunas)
meta_PRAD$patient
meta_PRAD$initial_weight



#Processamento de Dados Clínicos com TCGAbiolinks 
#obtem dados clínicos do projeto TCGA-PRAD (cancro da prostata), 
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

head (clinical.PRAD$clinical_drug_prad)
#converte os dados contidos em clinical.PRAD$clinical_patient_PRAD para um 
#dataframe e atribui esse novo dataframe à variável df
df = as.data.frame(clinical.PRAD$clinical_patient_prad)
View(df)



#Análise exploratória
#Pré-processamento e filtragem
#retirar as colunas dos metadados onde havia mais de 50 elementos como: 
#“not/Not reported/Reported” e/ou “NA”

cols_with_not_reported <- which(sapply(meta_PRAD,function(x) sum(x == "not reported", na.rm = TRUE)) > 50)
cols_with_Not_Reported <- which(sapply(meta_PRAD,function(x) sum(x == "Not Reported", na.rm = TRUE)) > 50)
cols_with_NA <- which(sapply(meta_PRAD, function(x) sum(is.na(x))) > 60)
# remover as colunas baseadas nos critérios específicos de cima
metadata_matriz_clean <- meta_PRAD[, -c(cols_with_not_reported, cols_with_Not_Reported, cols_with_NA)] 
dim(metadata_matriz_clean)


#código utilizado para extrair a informação relacionada à contagem da expressão dos genes do objeto rna_seq_PRAD
geneExp <- SummarizedExperiment::assay(rna_seq_PRAD)



##########

#Análise de Expressão Diferencial
library(DESeq2) #DESeq2, uma ferramenta para análise de expressão diferencial de dados de contagem de sequenciamento de RNA (RNA-Seq)

#Filtra os dados de RNA para incluir apenas amostras não nulo
data_de <- rna_seq_PRAD[,!is.na(rna_seq_PRAD@colData@listData[["vital_status"]])]

#Cria um objeto DESeqDataSet para análise, especificando um design experimental 
#que compara o status de IDH
ddsSE <- DESeqDataSet(data_de, design = ~ vital_status)

#Filtragem de Genes: Remove genes com contagens baixas (menos de 10) 
#para melhorar a confiabilidade da análise de expressão diferencial.
keep <- rowSums(counts(ddsSE)) >= 10

#Executa a análise de expressão diferencial com a função DESeq
ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)

resultsNames(ddsSE)
#comparação "WT vs Mutant" para o status do IDH, e converte os resultados para um dataframe
res <- results(ddsSE, name = "ddsSE@colData@listData[["paper_IDH1_mut"]].ddsSE@colData@listData[["paper_Mutations"]]")
dea <- as.data.frame(res)
#resume os resultados para obter uma visão geral dos achados estatísticos, 
#como o número de genes significativamente diferencialmente expressos
summary(res)

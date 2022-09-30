install.packages("BiocManager")
BiocManager::install('GEOquery')
BiocManager::install('oligo')
BiocManager::install('oligoClasses')
BiocManager::install('tidyverse')
BiocManager::install('limma')
BiocManager::install('dplyr',force = TRUE)
BiocManager::install("mta10transcriptcluster.db")

library(GEOquery)
gse81580 <- getGEO('GSE81580',GSEMatrix=TRUE,AnnotGPL=TRUE,destdir="./")
gse81580 <- gse81580[[1]]
getGEOSuppFiles('GSE81580')

library(oligo)
library(oligoClasses)
gse81580_celdata <- read.celfiles(list.celfiles('GSE81580_RAW',
                                                full.names=TRUE,listGzipped=TRUE))

library(tidyverse)
varLabels(gse81580)
gse81580$supplementary_file
pd <- pData(gse81580)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
gse81580_celdata <- read.celfiles(paste0('GSE81580_RAW/',pd$cel_file),
                                  phenoData=phenoData(gse81580))
colnames(pData(gse81580_celdata))
pData(gse81580_celdata)[,c("geo_accession","sex chromosome complement:ch1","surgery:ch1")]
names(pData(gse81580_celdata))[45] <- 'sex_chromosome_complement'
names(pData(gse81580_celdata))[46] <- 'surgery'

# Normalization
gse81580_eset <- rma(gse81580_celdata)
exprs(gse81580_eset)

library(limma)
library(dplyr)
pd <- pData(gse81580_eset)
pd$surgery <- as.factor(pd$surgery)
levels(pd$surgery) <- c("OVX","sham")
pd$surgery 
pd$group <- as.factor(paste0(pd$sex_chromosome_complement,pd$surgery))
levels(pd$group) <- c("XX.OVX","XX.sham","XY.OVX","XY.sham")
# Factorial design
# OVX stands for ovariectomy

design <- model.matrix(~ 0 + pd$group)
colnames(design) <- levels(pd$group)
design

contrast_matrix <- makeContrasts(genotype = (XX.OVX+XX.sham)-(XY.OVX+XY.sham),
                                 surgery = (XX.sham+XY.sham) - (XX.OVX+XY.OVX),
                                 interaction=(XX.sham - XY.OVX) - (XX.OVX - XY.sham),
                                 levels=design)
contrast_matrix 

gse81580_fit <- lmFit(gse81580_eset,design)
gse81580_fit2 <- contrasts.fit(gse81580_fit,contrasts = contrast_matrix)
gse81580_fit2 <- eBayes(gse81580_fit2)
summary(decideTests(gse81580_fit2,lfc=1))

# Annotation
ps <- rownames(topTable(gse81580_fit2,number=Inf,p.value = 0.05,lfc=1))
ps 
#probesets

BiocManager::install("mta10transcriptcluster.db")
library(mta10transcriptcluster.db)
# The annotation the microarray used
columns(mta10transcriptcluster.db)
keytypes(mta10transcriptcluster.db)
head(keys(mta10transcriptcluster.db,keytype="PROBEID"))
df <- AnnotationDbi::select(mta10transcriptcluster.db,ps,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
df
df2 <- dplyr::mutate(df,GENENAME=stringr::str_trunc(GENENAME,30))
df2

# Volcanoplot
interest_genes <- topTable(gse81580_fit2,number=Inf,p.value = 0.05,lfc=1)
volcanoplot(gse81580_fit2, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interest_genes)))

# Heatmap
eset_of_interest <- gse81580_eset[rownames(interest_genes),]
heatmap(exprs(eset_of_interest), labRow = df2$SYMBOL, 
        labCol=pd$group,
        margins = c(6,6))

vennDiagram(decideTests(gse81580_fit2,lfc=1),
            circle.col = c("red2", "blue3", "orange2"))



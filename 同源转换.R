

library(Seurat)
library(data.table)
library(tidyverse)
setwd("C:/Users/jk/Desktop/COF/单细胞/细胞通讯/")
counts <- fread(input = "excellType1/counts.txt") %>% tibble::column_to_rownames("ENSEMBL")
counts[1:4,1:5]

library(homologene)

library(clusterProfiler)

symbol <- bitr(rownames(counts),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Mm.eg.db")
symbol <- symbol[!duplicated(symbol$SYMBOL),]
symbol <- symbol[!duplicated(symbol$ENSEMBL),]
symbol <- filter(symbol,ENSEMBL%in%rownames(counts))
symbol <- tibble::column_to_rownames(symbol,var = "ENSEMBL")
counts <- counts[rownames(symbol),]
rownames(counts) <- symbol$SYMBOL
homo <- homologene::homologene(symbol$SYMBOL,inTax =10090,outTax = 9606)
homo <- homo[!duplicated(homo$`9606`),]
counts <- counts[homo$`10090`,]
rownames(counts) <- homo$`9606`


write.table(counts,file = "human_counts.txt",sep = "\t")
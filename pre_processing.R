##### Pre Processing Data #####
args = commandArgs(trailingOnly = T)

GSE = args[1]
disease = args[2]
outname = args[3]

##### Loading libraries #####
library(limma)
library(matrixStats)
library(GEOquery)
library(affy)
library(biomaRt)
library(NOISeq)
library(edgeR)
library(DESeq2)
library(minfi)
library(DMRcate)
library(wateRmelon)
library(lumi)
library(hipathia)
library(KEGG.db)
library(KEGGprofile)
library(rowr)
library(jsonlite)

source("lib.R")

hipathia_pathways <- load_pathways(species="hsa")
geneMap <- read.delim('GCF_000001405.38_GRCh38.p12_map.csv',stringsAsFactors = F,header = F)
colnames(geneMap) <- c('entrez','gensymb')
geneMap <- unique(geneMap)
kegg_universe <- unique(read.delim("KEGG_genes.tsv",stringsAsFactors = F)[,'entrez'])

results <- pre_processing(GSE,disease,hipathia_pathways,geneMap,kegg_universe)
save(results,file=paste0(outname,'.RData'))

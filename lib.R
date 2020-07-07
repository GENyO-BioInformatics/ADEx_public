doExtractInfo <- function(GSE,disease){
  GSE_pheno = read.delim(paste0(GSE,'.txt'),stringsAsFactors = F)
  disease_samples = GSE_pheno[GSE_pheno$Condition==disease,"GSM"]
  control_samples = GSE_pheno[GSE_pheno$Condition=="Healthy","GSM"]
  GSE_pheno <- GSE_pheno[which(GSE_pheno$GSM %in% disease_samples | GSE_pheno$GSM %in% control_samples), ]
  gpl <- unique(GSE_pheno$GPL)
  GPL_info <- getGEO(gpl)
  protocol <- GPL_info@header$manufacturer
  if (length(protocol)==0){
    protocol = 'RNA-Seq'
  }
  return_list <- list(protocol=protocol,gpl=gpl,pheno=GSE_pheno,
                      disease_samples = disease_samples,
                      control_samples = control_samples)
  return(return_list)
}

convert_table <- function(expressionMatrix,GPL_Table){
  GPL_Table <- GPL_Table[GPL_Table$probe %in% rownames(expressionMatrix),]
  GeneSymbols <- unique(GPL_Table[,'gene'])
  annotationMatrix <- matrix(NA,
                             ncol = ncol(expressionMatrix), 
                             nrow = (length(GeneSymbols)))
  rownames(annotationMatrix) <- GeneSymbols
  colnames(annotationMatrix) <- colnames(expressionMatrix)
  for (gene in GeneSymbols){
    probes <- GPL_Table[GPL_Table$gene==gene,'probe']
    if (length(probes) == 1){
      values <- expressionMatrix[probes,]
    } else{
      values <- colMedians(expressionMatrix[probes,])
    }
    annotationMatrix[gene,] <- values
  }
  return(annotationMatrix)
}

doAffymetrix <- function(pheno,GSE,gpl){
  GPL_Table <- read.csv(paste0('gpl/',gpl,'.csv'),stringsAsFactors = F)
  filenames <- pheno$GSM
  table.of.names <- data.frame(matrix(NA,ncol = 1,nrow = length(filenames)))
  rownames(table.of.names) <- filenames
  colnames(table.of.names) <- 'names.col'
  for (filename in filenames){
    gsm_info <- getGEO(filename)
    gsm_file <- basename(gsm_info@header$supplementary_file)[1]
    filepath <- list.files(GSE,recursive = F,full.names = T)
    idx <- grep(gsm_file,filepath)
    table.of.names[filename,"names.col"] <- basename(filepath[idx])
  }
  
  list.of.files <- list.files(GSE,full.names = T)
  list.of.files <- list.of.files[grep('.gz',list.of.files)]
  expressionMatrix_affy <- ReadAffy(filenames = list.of.files)
  expressionMatrix <- rma(expressionMatrix_affy,destructive = F,normalize = F,background = F)
  expressionMatrix <- exprs(expressionMatrix)
  x <- log2(100)
  p <- round(ncol(expressionMatrix)*0.1)
  keep <- c()
  for (i in rownames(expressionMatrix)){
    values <- expressionMatrix[i,]
    values <- values[values >= x]
    if (length(values)>=p){
      keep <- c(keep,i)
    }
  }
  expressionMatrix <- rma(expressionMatrix_affy)
  expressionMatrix <- exprs(expressionMatrix)
  expressionMatrix <- na.omit(expressionMatrix[keep,])
  
  expressionMatrix <- convert_table(expressionMatrix,GPL_Table)
  
  for (i in 1:ncol(expressionMatrix)){
    select_row <- which(table.of.names$names.col == colnames(expressionMatrix)[i])
    colnames(expressionMatrix)[i] <- rownames(table.of.names)[select_row]
  }
  return(expressionMatrix)
}

doIllumina <- function(pheno,GSE,gpl){
  GPL_Table <- read.csv(paste0('gpl/',gpl,'.csv'),stringsAsFactors = F)
  file_parsing <- read.delim(paste0(GSE,'/',GSE,'.tsv'),stringsAsFactors = F)
  rownames(file_parsing) <- file_parsing[,1]
  file_parsing <- file_parsing[,-1]
  detectionMatrix <- file_parsing[,grep("Detection.Pval",colnames(file_parsing))]
  expressionMatrix <- file_parsing[,-(grep("Detection.Pval",colnames(file_parsing)))]
  colnames(detectionMatrix) <- colnames(expressionMatrix) <- substr(colnames(detectionMatrix),0,regexpr("\\.",colnames(detectionMatrix))-1)
  expressionMatrix <- expressionMatrix[,pheno$GSM]
  detectionMatrix <- detectionMatrix[,pheno$GSM]
  expressionMatrix <- neqc(expressionMatrix,detection.p = detectionMatrix)
  detection.threshold <- 0.05
  samples.threshold.fraction <- 0.1
  number <- ncol(detectionMatrix)*samples.threshold.fraction
  keep <- c()
  for (probe in rownames(detectionMatrix)){
    values <- detectionMatrix[probe,]
    values <- values[values < detection.threshold]
    if (length(values)>number){
      keep <- c(keep,probe)
    }
  }
  
  expressionMatrix <- na.omit(expressionMatrix[keep,])
  expressionMatrix <- convert_table(expressionMatrix,GPL_Table)
  return(expressionMatrix)
  
}

doRNASeq <- function(GSE,pheno,disease){
  ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  geneInfo <- getBM(attributes=c('ensembl_gene_id','gene_biotype',
                                 'chromosome_name','start_position','end_position', "description",
                                 "external_gene_name", "percentage_gene_gc_content"),
                    mart = ensembl)
  rsem_files = list.files(full.names=T, 
                          pattern = "*genes.results",
                          path = GSE)
  verified_rsem <- c()
  for (rsem_file in rsem_files){
    gsm <- basename(rsem_file)
    gsm <- substr(gsm,1,regexpr('\\.',gsm)-1)
    if (gsm %in% pheno$GSM){
      verified_rsem <- c(verified_rsem,rsem_file)
    }
  }
  
  rsem_files <- verified_rsem
  targets = read.delim(rsem_files[1], 
                       as.is = TRUE)[,1]
  rsem_expcounts = data.frame(targets, row.names=1)
  rsem_FPKM = data.frame(targets, row.names=1)
  for (k in 1:length(rsem_files)){
    column <-  gsub(".genes.results", "", rsem_files[k])
    column <- gsub(paste0(GSE,"/"),"",column)
    rsem_expcounts[column] = read.delim(rsem_files[k], header=TRUE, as.is = TRUE)[,"expected_count"]
  }
  
  samples = colnames(rsem_expcounts)
  l <- length(samples)
  countData <- rsem_expcounts
  clinical <- data.frame(Condition=as.factor(pheno$Condition),row.names = colnames(countData))
  mygc = geneInfo$percentage_gene_gc_content
  names(mygc) = geneInfo$ensembl_gene_id
  mybiotypes = geneInfo$gene_biotype
  names(mybiotypes) = geneInfo$ensembl_gene_id
  mychroms = geneInfo[,c("chromosome_name", 'start_position','end_position')]
  rownames(mychroms) = geneInfo$ensembl_gene_id
  mylength = read.delim(rsem_files[1], header=TRUE, as.is = TRUE)[,"length"]
  names(mylength) = targets
  mydata = readData(data=countData, 
                    factors=clinical, 
                    biotype=mybiotypes, 
                    chromosome=mychroms, 
                    gc=mygc, 
                    length=mylength)
  countDataFiltered = filtered.data(countData, 
                                    factor = clinical$Condition, 
                                    norm=FALSE, 
                                    method=1, cpm=0.5)
  countData <- countData[rownames(countDataFiltered),]
  genes <- geneInfo[geneInfo$ensembl_gene_id %in% rownames(countData),]
  symbol_names <- unique(genes$external_gene_name)
  
  n_countData <- data.frame(matrix(0,nrow = length(symbol_names),ncol = ncol(countData)))
  rownames(n_countData) <- symbol_names
  colnames(n_countData) <- colnames(countData)
  for (gene in symbol_names){
    idx <- genes[genes$external_gene_name==gene,]
    ensembl_ids <- idx$ensembl_gene_id
    values <- countData[ensembl_ids,]
    if (nrow(idx)>1){
      n_countData[gene,] <- colMedians(as.matrix(values))
    } else{
      n_countData[gene,] <- values
    }
  }
  
  countData <- round(n_countData,0)
  normFactors = calcNormFactors(countData, method="TMM", refColumn = 1, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = T, Acutoff = -1e+10)
  normFactors = normFactors*(colSums(countData)/mean(colSums(countData)))
  names(normFactors) = colnames(countData)
  countData <- round(countData,0)
  ddsFromMatrix <- DESeqDataSetFromMatrix(countData = countData,
                                          colData = clinical,
                                          design = ~Condition)
  sizeFactors(ddsFromMatrix) = normFactors
  dds <- DESeq(ddsFromMatrix)
  res <- as.data.frame(results(dds, c("Condition",disease,"Healthy")))
  res <- na.omit(res)
  expressionMatrix <- as.data.frame(counts(dds, normalized=T))
  expressionMatrix <- log2(expressionMatrix+1)
  res <- res[,c("log2FoldChange","padj")]
  colnames(res) <- c('logFC','adj.P.Val')
  result <- list(diffExpression=res,expressionMatrix=expressionMatrix)
  return(result)
}


dolimma = function(data, input1, input2) {
  data = data[,c(input1, input2)]
  input1Name =  deparse(substitute(input1)) 
  input2Name =  deparse(substitute(input2)) 
  design = matrix(data=0,ncol=2,nrow=ncol(data),dimnames = list(c(input1, input2), c(input1Name, input2Name)))
  design[1:length(input1),1] = 1
  design[(length(input1)+1):(length(input1) + length(input2)),2] = 1
  fit <- lmFit(data, design)
  contrast = paste0(input1Name, "-", input2Name)
  cont.matrix <- makeContrasts(contrasts = contrast, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.05)
  limma <- topTable(fit2, adjust="fdr", sort.by="B", resort.by = "p", number=dim(data)[1])
  return(limma)
}


getInfoDataset <- function(GSE, gpl,detP = TRUE) {
  signalFile = read.delim(paste0(GSE, "/", GSE, ".tsv"), row.names = 1, sep = "\t",  stringsAsFactors = F)
  detP = signalFile[seq(3, ncol(signalFile), by = 3)]
  if (gpl == "GPL8490"){
    platform = "IlluminaHumanMethylation27k"
  } else{
    platform = "IlluminaHumanMethylation450k"
  }
  return_list=list(detP=detP,platform=platform)
  return(return_list)
}

process_meth <- function(GSE,gpl,Uname = ".Unmethylated.signal", Mname = ".Methylated.signal"){
  res <- getInfoDataset(GSE,gpl)
  detP <- res$detP
  platform <- res$platform
  signalIntensities <- readGEORawFile(paste0(GSE, "/", GSE, ".tsv"), sep="\t",
                                      array = platform,
                                      annotation = "ilmn12.hg19",
                                      Uname = Uname,
                                      Mname = Mname)
  
  
  # Minfi normalization
  mset <- ratioConvert(signalIntensities)
  mset <- mapToGenome(mset)
  
  # Removing poorly performing probes (p>0.05 in >10% of samples)
  if (!is.null(detP)) { # some datasets do not provide detection P values
    detP = na.omit(detP[rownames(detP) %in% rownames(mset),])
    discard <- rownames(detP)[rowSums(detP > 0.05) > ncol(mset)*0.1]
    mset <- mset[!rownames(mset) %in% discard,]
  }
  # Remove probes with SNPs
  mset <- dropLociWithSnps(mset)
  
  # Remove cross hibridizing and XY probes
  beta <- getBeta(mset)
  beta <- rmSNPandCH(beta, rmXY = TRUE)
  
  # Quantile normalization
  beta <- lumiN(beta, method = "quantile")
  
  # BMIQ normalization
  if (platform == "IlluminaHumanMethylation450k") {
    probeType = as.data.frame(getAnnotation(signalIntensities)[rownames(beta),c("Name","Type")])
    probeType = ifelse(probeType$Type %in% "I",1,2)
    
    beta <- apply(beta, 2, function(x) {BMIQ(x, design.v = probeType)$nbeta})
  }
  
  # Save data
  return(beta)
}

hipathia_function <- function(pheno,expressionMatrix,hipathia_pathways,disease){
  # Get ENTREZ IDs
  trans_data <- translate_data(expressionMatrix, "hsa")
  # scale data between 0 and 1
  exp_data <- normalize_data(trans_data)
  # Transduction signal calculation
  results = hipathia(exp_data, hipathia_pathways)
  # Comparing between conditions (Status column)
  path_vals <- get_paths_data(results)
  sample_group <- pheno[pheno$GSM %in% colnames(expressionMatrix),"Condition"]
  # From here, we use the same code with each comparison.
  # If the dataset have two comparison (SLE vs HC and RA vs HC)
  # the following code must be repeated twice, with the corresponding changes.
  comp <- do_wilcoxon(path_vals, sample_group, g1 = disease, g2 = "Healthy")
  pathways_summary <- get_pathways_summary(comp, hipathia_pathways) # Get the results
  return(pathways_summary)
}

getKeggRes <- function(geneSelected,geneMap,kegg_universe){
  geneSelected <- geneMap[geneMap$gensymb %in% geneSelected,"entrez"]
  keggResults <- find_enriched_pathway(species = "hsa",gene = geneSelected,returned_pvalue = 1,returned_genenumber = 1,returned_adjpvalue = 1,refGene = kegg_universe)
  if(nrow(keggResults$stastic) == 0){
    cat("NO significative keggResults") 
    return(NA)
  }
  pathways <- row.names(keggResults$stastic)
  pathways <- pathways[pathways != "01100"] # 01100 the broadest pathway "metabolism" not informative
  genesOfPathways <- c()
  for (pathway in pathways){
    entrezsOfpath <- keggResults$detail[[pathway]]  
    geneSymbsOfpath <- geneMap$gensymb[geneMap$entrez %in% entrezsOfpath]
    geneSymbsOfpathStr <- paste(geneSymbsOfpath,collapse = ", ")
    genesOfPathways <- c(genesOfPathways,geneSymbsOfpathStr)
  }
  keggResultsFull <- cbind(keggResults$stastic[pathways,],genes=genesOfPathways)
  keggResultsFull <- keggResultsFull[order(keggResultsFull$pvalueAdj),]
  return(keggResultsFull)
}

pre_processing <- function(GSE,disease,hipathia_pathways,geneMap,kegg_universe){
  GSE_info = doExtractInfo(GSE,disease)
  protocol = GSE_info$protocol
  gpl = GSE_info$gpl
  pheno = GSE_info$pheno
  disease_samples = GSE_info$disease_samples
  control_samples = GSE_info$control_samples
  if (protocol == "Affymetrix"){
    expressionMatrix <- doAffymetrix(pheno,GSE,gpl)
    diffExpression <- dolimma(expressionMatrix,disease_samples,control_samples)
    hipathia_res <- hipathia_function(pheno,expressionMatrix,hipathia_pathways,disease)
    geneSelected <- rownames(diffExpression[diffExpression$adj.P.Val<0.05,])
    kegg_res <- getKeggRes(geneSelected,geneMap,kegg_universe)
    return_list=list(expressionMatrix=expressionMatrix,
                     diffExpression=diffExpression,
                     hipathia_res=hipathia_res,
                     kegg_res=kegg_res)
    return(return_list)
  } else if (protocol == "Illumina Inc."){
    expressionMatrix <- doIllumina(pheno,GSE,gpl)
    diffExpression <- dolimma(expressionMatrix,disease_samples,control_samples)
    hipathia_res <- hipathia_function(pheno,expressionMatrix,hipathia_pathways,disease)
    geneSelected <- rownames(diffExpression[diffExpression$adj.P.Val<0.05,])
    kegg_res <- getKeggRes(geneSelected,geneMap,kegg_universe)
    return_list=list(expressionMatrix=expressionMatrix,
                     diffExpression=diffExpression,
                     hipathia_res=hipathia_res,
                     kegg_res=kegg_res)
    return(return_list)
  } else if (protocol == "RNA-Seq"){
    res <- doRNASeq(GSE,pheno,disease)
    expressionMatrix <- res$expressionMatrix
    diffExpression <- res$diffExpression
    hipathia_res <- hipathia_function(pheno,expressionMatrix,hipathia_pathways,disease)
    geneSelected <- rownames(diffExpression[diffExpression$adj.P.Val<0.05,])
    kegg_res <- getKeggRes(geneSelected,geneMap,kegg_universe)
    return_list=list(expressionMatrix=expressionMatrix,
                     diffExpression=diffExpression,
                     hipathia_res=hipathia_res,
                     kegg_res=kegg_res)
    return(return_list)
  } else{
    beta_matrix <- process_meth(GSE,gpl)
    return(beta_matrix)
  }
}
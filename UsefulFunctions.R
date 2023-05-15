dds_to_volcano <- function(dds) {
  library(gsubfn)
  require(ggplot2)
  require(ggrepel)
  keep <- rowSums(counts(dds)) >= ncol(dds)
  dds <- dds[keep,]
  dds = DESeq(dds)
  res = results(dds)
  resOrdered <- res[order(res$padj), ]
  genesInOrder <- as.data.frame(row.names(resOrdered))
  genesInOrder$order = 1:nrow(genesInOrder)
  genesInOrder = merge(genesInOrder,G_list[!duplicated(G_list$ensembl_gene_id),],
                       by.x = 1,by.y="ensembl_gene_id",all.x=TRUE)
  genesInOrder = genesInOrder[order(genesInOrder$order),]
  0 %in% match(genesInOrder$`row.names(resOrdered)`,row.names(resOrdered),nomatch=0)
  genesInOrder_2 = genesInOrder
  genesInOrder_2$ensembl_id = as.character(genesInOrder_2$`row.names(resOrdered)`)
  genesInOrder_2$hgnc_symbol_character = as.character(genesInOrder_2$hgnc_symbol)
  genesInOrder_2$hgnc_symbol_character[1:3]
  genesInOrder_2$hgnc_symbol_character = ifelse(genesInOrder_2$hgnc_symbol_character == "",
                                                genesInOrder_2$ensembl_id,print(genesInOrder_2$hgnc_symbol_character))
  genesInOrder_2 = genesInOrder_2[,c(-1,-3)]
  resO_2 = resOrdered
  resO_2$geneName = NA_character_
  resO_2$geneName = genesInOrder_2$hgnc_symbol
  grep("^H2", genesInOrder_2$hgnc_symbol_character)
  mhcIIgenes_VP = c("CD74",genesInOrder_2[c(grep("^HLA-D", genesInOrder_2$hgnc_symbol_character)),"hgnc_symbol_character"])
  housekeeping_genes = c("ACTB","GAPDH")
  relevantTranscripts = c(mhcIIgenes_VP,housekeeping_genes)
  keyvals <- rep('lightsteelblue1', nrow(resO_2))
  names(keyvals) <- rep('Other_Genes', nrow(resO_2))
  c(grep("^HLA-D", resO_2$geneName),grep("CD74",resO_2$geneName))
  keyvals[c(grep("^HLA-D", resO_2$geneName),grep("CD74",resO_2$geneName))] = "black"
  names(keyvals)[c(grep("^HLA-D", resO_2$geneName),grep("CD74",resO_2$geneName))] = "MHCII_Genes"
  keyvals[match(housekeeping_genes,resO_2$geneName)] = "red"
  names(keyvals)[match(housekeeping_genes,resO_2$geneName)] = "Housekeeping_Genes"
  unique(names(keyvals))
  unique(keyvals)
  keyvals.shape <- rep(3, nrow(resO_2))
  names(keyvals.shape) <- rep('PBC', nrow(resO_2))
  keyvals.shape[which(rownames(resO_2) %in% mhcIIgenes_VP)]
  names(keyvals.shape)[which(rownames(resO_2) %in% mhcIIgenes_VP)] <- 'MHCII Genes'
  df = resO_2
  df$pvalueOffAxis = ifelse(df$pvalue<10**-3,TRUE,FALSE)
  df$pvalueWithOffAxisValues = ifelse(df$pvalueOffAxis,10**-3,df$pvalue)
  df$log2FCOffAxis = ifelse(abs(df$log2FoldChange)>3,TRUE,FALSE)
  df$log2FCWithOffAxisValues = ifelse(df$log2FCOffAxis,
                                      sign(df$log2FoldChange)*3,
                                      df$log2FoldChange)
  df$offAxis = ifelse(df$pvalueOffAxis,TRUE,
                      ifelse(df$log2FCOffAxis,TRUE,FALSE))
  df$ht_mhcII = ifelse(grepl("^HLA-D",df$geneName),"HLA_gene","Other_gene")
  df$ht_mhcII = ifelse(df$geneName == "CD74","HLA_gene",df$ht_mhcII)
  df$ht_mhcII = ifelse(df$geneName == "GAPDH","Housekeeping_gene",df$ht_mhcII)
  df$ht_mhcII = ifelse(df$geneName == "ACTB","Housekeeping_gene",df$ht_mhcII)
  df$ht_mhcII = ifelse(df$geneName == "AICDA","Interesting_gene",df$ht_mhcII)
  df$ht_mhcII = ifelse(df$geneName == "CIITA","Interesting_gene",df$ht_mhcII)
  df$ht_mhcII = factor(df$ht_mhcII, levels = c("Other_gene","Housekeeping_gene","HLA_gene","Interesting_gene"))
  table(df$ht_mhcII)
  mhc_genes = c("B2M","CD58",genesInOrder_2[c(grep("^HLA-", genesInOrder_2$hgnc_symbol_character)),"hgnc_symbol_character"])
  mhcIgenes = mhc_genes[is.na(match(mhc_genes,mhcIIgenes_VP))]
  df$ht_mhcI = ifelse(df$geneName %in% mhcIgenes,"HLA_gene","Other_gene")
  df$ht_mhcI = ifelse(df$geneName %in% housekeeping_genes,"Housekeeping_gene",df$ht_mhcI)
  interesting_genes = c("AICDA","CIITA")
  df$ht_mhcI = ifelse(df$geneName %in% interesting_genes,"Interesting_gene",df$ht_mhcI)
  df$ht_mhcI = factor(df$ht_mhcI, levels = c("Other_gene","Housekeeping_gene","HLA_gene","Interesting_gene"))
  df = df[order(df$ht_mhcII),]
  df = as.data.frame(df)
  df[order(-abs(df$log2FoldChange)),]
  mhcIIfullDataGG = ggplot(df, aes(x = (log2FCWithOffAxisValues), y = -log10(pvalueWithOffAxisValues),
                                   color = ht_mhcII, shape = offAxis, label=df$geneName)) + 
    geom_point() + 
    scale_x_continuous(limits = c(-3, 3)) + 
    scale_y_continuous(limits = c(0,3)) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    theme_bw() +
    scale_color_manual(values=c("lightsteelblue1","indianred3","seagreen","darkorchid4")) + 
    geom_text_repel(aes(label = ifelse(as.character(df$ht_mhcII) != "Other_gene",
                                       as.character(df$geneName),'')),
                    nudge_x = 0.2, nudge_y = 0.1,
                    force = 5,
                    segment.size = 0.1, segment.alpha = 1, 
                    point.padding = NA, box.padding = 1,
                    size = 5)
  return(list(df,mhcIIfullDataGG))
}
removeENSGDecimal <- function(df){
  library(stringr)
  rownames(df) = str_split_fixed(rownames(df),"\\.",Inf)[,1]
  return(df)
}
addLengths <- function(df, which_col_Ensembl){
  gLG_2 = read.csv("~/Documents/Programming/reference_files/mm10_star_reference_files/mm10_gene_lengths.txt",row.names =  1)
  df$order = 1:nrow(df)
  df = merge(df,gLG_2,by.x = which_col_Ensembl,by.y=0,all.x=TRUE)
  df = df[order(df$order),]
  return(df)
 }
add_HGNC <- function(cts){
  G_list = read.csv("~/Documents/Programming/UsefulTables/Ensembl_To_Hgnc_Symbol.csv",row.names =  1)
  genesInOrder <- as.data.frame(row.names(cts))
  genesInOrder$order = 1:nrow(genesInOrder)
  genesInOrder = merge(genesInOrder,G_list[!duplicated(G_list$ensembl_gene_id),],
                       by.x = 1,by.y="ensembl_gene_id",all.x=TRUE)
  genesInOrder = genesInOrder[order(genesInOrder$order),]
  0 %in% match(genesInOrder$`row.names(cts)`,row.names(cts),nomatch=0)
  genesInOrder_2 = genesInOrder
  genesInOrder_2$ensembl_id = as.character(genesInOrder_2$`row.names(cts)`)
  genesInOrder_2$hgnc_symbol_character = as.character(genesInOrder_2$hgnc_symbol)
  genesInOrder_2$hgnc_symbol_character = ifelse(genesInOrder_2$hgnc_symbol_character == "",
                                                genesInOrder_2$ensembl_id,genesInOrder_2$hgnc_symbol_character)
  genesInOrder_2$hgnc_symbol_character = ifelse(is.na(genesInOrder_2$hgnc_symbol_character),
                                                genesInOrder_2$ensembl_id,genesInOrder_2$hgnc_symbol_character)
  genesInOrder_2 = genesInOrder_2[,c(-1,-3)]
  rownames(cts) = make.unique(genesInOrder_2$hgnc_symbol)
  return(cts)
}
scatterGenes <- function(df_cts,hasTumorType,tumorTypes,includeLineOfBestFit) {
  require(reshape2)
  require(ggplot2)
  require(ggrepel)
  require(ggpmisc)
  mhcIIgenes_VP = c("CD74","CIITA",rownames(df_cts)[c(grep("^HLA-D", rownames(df_cts)))])
  housekeeping_genes = c("ACTB","GAPDH")
  head(mhcIIgenes_VP)
  
  relevantTranscripts_temp = c(mhcIIgenes_VP,housekeeping_genes,"AICDA")
  relevantTranscripts = c("AICDA","ACTB","GAPDH","CD74","CIITA",rownames(df_cts)[c(grep("^HLA-D", rownames(df_cts)))])
  relMatrix = df_cts[row.names(df_cts) %in% relevantTranscripts,]
  relMatrix = relMatrix[c("AICDA","ACTB","GAPDH","CIITA","CD74",rownames(relMatrix)[c(grep("^HLA-D", rownames(relMatrix)))]),]
  if(hasTumorType){
    trm = as.data.frame(t(relMatrix))
    rmwc = cbind(trm,tumorTypes)
    rmwcm = melt(rmwc, id.vars = c("AICDA","tumorTypes"))
    rmwcm$AICDA = as.numeric(rmwcm$AICDA)
    rmwcm$value = as.numeric(rmwcm$value)
    rmwcm$color = ifelse(rmwcm$tumorTypes == "GCB","red","blue")
    scatterplotData = ggplot(rmwcm, aes(x = AICDA, y = value)) + 
      geom_point(color = rmwcm$color) + facet_wrap(~variable, scales = "free_y")
  }else{
    if (includeLineOfBestFit) {
      relMatrixMelt = melt(as.data.frame(t(relMatrix)), id.vars = "AICDA")
      scatterplotData = ggplot(relMatrixMelt, aes(x = AICDA, y = value)) +
        geom_point() + geom_smooth(method = 'lm', formula = y ~ x) +
        stat_poly_eq(formula = y ~ x,
                       aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                       parse = TRUE, size = 3) + 
        facet_wrap( ~ variable, scales = "free_y")
    } else {
      relMatrixMelt = melt(as.data.frame(t(relMatrix)), id.vars = "AICDA")
      scatterplotData = ggplot(relMatrixMelt, aes(x = AICDA, y = value)) +
        geom_point() +
        facet_wrap( ~ variable, scales = "free_y",
                    nrow = 8, ncol = 3) + theme(aspect.ratio = 1)
    }
    }
  return(scatterplotData)
}
countsToTPM <- function(cts,lengthMatrix) {
  rpk = (cts/lengthMatrix$Length)
  tpm_cts = sweep(rpk,2,colSums(rpk),`/`)*1e6
  
  return(tpm_cts)
}
countsToCPM <- function(cts) {
  return(sweep(cts,2,colSums(cts),`/`)*1e6)
}
geneCPMvsTPM <- function(cts,cML,geneName,color,colors = NULL){
  cpm_cts = countsToCPM(cts)
  tpm_cts = countsToTPM(cts,cML)
  gene_cpm_tpm = as.data.frame(cbind(t(cpm_cts[geneName,]),t(tpm_cts[geneName,])))
  colnames(gene_cpm_tpm) = c("CPM","TPM")
  gctc = merge(gene_cpm_tpm,colors,by.x = 0, by.y = 0)
  rownames(gctc) = gctc$Row.names
  gctc = gctc[,-1]
  if(!color){
    gene_cpm_vs_tpm = ggplot(gene_cpm_tpm, aes(x = CPM, y = TPM,)) + 
    geom_point() +
    scale_x_log10() + scale_y_log10()
  } else {
    gene_cpm_vs_tpm = ggplot(gctc, aes(x = CPM, y = TPM,)) + 
    geom_point(color = gctc$color) +
    scale_x_log10() + scale_y_log10()
  }
}
createSampleTable <- function(counts,gene,splitall,low,high){
  countsTransposed = as.data.frame(t(counts))
  if(splitall){
    countsTransposed$quantile = with(countsTransposed, factor(
    findInterval(AICDA, c(-Inf,
                         quantile(AICDA, probs=low), Inf)), 
    labels=c("AID_Low","AID_High")
  ))
  } else {
    countsTransposed$quantile = with(countsTransposed, factor(
      findInterval(AICDA, c(-Inf,
                            quantile(AICDA, probs=c(low,high)), Inf)), 
      labels=c("AID_Low","Middle_Fraction","AID_High")
    ))
  }
  sampleTable = as.data.frame(countsTransposed[,"quantile"])
  sampleTable$samples = row.names(countsTransposed)
  sampleTable = sampleTable[,c(2,1)]
  colnames(sampleTable)[2] = "quantile"
  return(sampleTable)
}
createDESeqObject = function(directory,sampleTable){
  require(DESeq2)
  stfDES = sampleTable[!sampleTable$quantile %in% "Middle_Fraction",]
  stfDES$quantile = factor(stfDES$quantile)
  stfDES$samples = ifelse(substr(stfDES$samples,1,1) == 'X', sub('.', '', stfDES$samples),stfDES$samples)
  stfDES$fileName = paste0(substr(stfDES$samples,1,8),"-",substr(stfDES$samples,10,13),"-",
                           substr(stfDES$samples,15,18),"-",substr(stfDES$samples,20,23),"-",
                           substr(stfDES$samples,25,53))
  stfDES = stfDES[,c("samples","fileName","quantile")]
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = stfDES,
                                         directory = directory,
                                         design= ~ quantile)
  rownames(ddsHTSeq) = rownames(removeENSGDecimal(as.data.frame(rownames(ddsHTSeq),
                                                                row.names = rownames(ddsHTSeq))))
  return(ddsHTSeq)
}
#Make tables for each LymphGen predicted subset from table containing a melted list
#See input file
makeSampleCSVs = function(dlbcl_subject_table) {
  require(rlist)
  sampleSubsets = (dlbcl_subject_table %>%
                     group_by(LymphGen.Subtype.prediction.final) %>%
                     group_split())
  write.list(file = "~/Documents/Programming/dbGAP/DLBCL/samples_subsets",sampleSubsets)
  setwd("~/Documents/Programming/dbGAP/DLBCL/Subjects/LymphGen/")
  sampleClassifications = as.character(as.data.frame(t(as.data.frame(lapply(sampleSubsets, function(l) l[[3]][1]))))$V1)
  names(sampleSubsets) = make.names(sampleClassifications)
  
  justSubsetList = list()
    
  # read each file as array element of DT and rename the last 2 cols
  # we created a list of single sample tables
  # Removed the step that removes the htseq_counts.txt from filename
  for (i in 1:length(sampleSubsets)) {
    write.table(as.character(sampleSubsets[[i]]$dbGAP.subject.ID),
                file = paste0("~/Documents/Programming/dbGAP/DLBCL/Subjects/LymphGen/",
                              make.names(as.character(sampleSubsets[[i]]$LymphGen.Subtype.prediction.final)[1]),
                              ".txt"), 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  write.table(as.character(sampleSubsets[[9]]$dbGAP.subject.ID),
              file = paste0("~/Documents/Programming/dbGAP/DLBCL/Subjects/LymphGen/",
                            "EZB-MYC-added",
                            ".txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}
#Make a pair of plots (1 tiff for viewing and 1 PDF for illustrator) from a ggplot
#Place in specified directory with specified plotname
doubleplot = function(ggplot_dp,write_directory_dp,plotname){
  ggsave(paste0(write_directory_dp,"/",plotname,".tiff"),
         plot = ggplot_dp, device = "tiff",
         width=10, height=10)
  ggsave(paste0(write_directory_dp,"/",plotname,".eps"),
         plot = ggplot_dp, device = "eps",
         width=10, height=10)
  pdf(file = paste0(write_directory_dp,"/",plotname,".pdf"),
    width = 10, height = 10)
  print(ggplot_dp)
  dev.off()
}
geneCombinationForGGPlot = function(written_directory,gene){
  require(RColorBrewer)
  require(ggplot2)
  setwd(written_directory)
  tpm_name = paste0(gene,"_TPM")
  gene_level_complete = data.frame(Type = character(),gene_TPM=double()) 
  
  for(sampleType in 
      list.dirs(recursive=F)[which(lengths(lapply(list.dirs(recursive=F), list.files)) > 0)]
      ){
    gene_level = as.data.frame(t(read.csv(paste0(sampleType,"/tpm_cts.csv"),row.names = 1)[gene,]))
    gl2 = cbind(rep.int(sampleType,ncol(gene_level)),gene_level)
    colnames(gl2) = c("type","gene_TPM")
    rownames(gl2) = NULL
    gene_level_complete = rbind(gene_level_complete,gl2)
  }
  gene_level_complete$type = as.factor(str_split_fixed(gene_level_complete$type,"/",Inf)[,2])
  return(gene_level_complete)
}
add_MGI = function(cts, ensembl_column=0){
  require(stringr)
  G_list = read.csv("~/Documents/Programming/Rahul_Sequencing/Project_10863_B/R_files/Ensembl_To_Mgi_Symbol.csv",row.names =  1)
  genesInOrder <- as.data.frame(row.names(cts))
  genesInOrder$order = 1:nrow(genesInOrder)
  genesInOrder$EnsemblID = str_split_fixed(genesInOrder$`row.names(cts)`,'\\.',2)[,1]
  genesInOrder = merge(genesInOrder,G_list[!duplicated(G_list$ensembl_gene_id),],
                       by.x = "EnsemblID",by.y="ensembl_gene_id",all.x=TRUE)
  genesInOrder = genesInOrder[order(genesInOrder$order),]
  0 %in% match(genesInOrder$`row.names(cts)`,row.names(cts),nomatch=0)
  genesInOrder_2 = genesInOrder
  genesInOrder_2$EnsemblID = as.character(genesInOrder_2$EnsemblID)
  genesInOrder_2$mgi_symbol_character = as.character(genesInOrder_2$external_gene_name)
  genesInOrder_2$mgi_symbol_character = ifelse(genesInOrder_2$mgi_symbol_character == "",
                                                genesInOrder_2$ensembl_id,genesInOrder_2$mgi_symbol_character)
  genesInOrder_2$mgi_symbol_character = ifelse(is.na(genesInOrder_2$mgi_symbol_character),
                                                genesInOrder_2$EnsemblID,genesInOrder_2$mgi_symbol_character)
  genesInOrder_3 = genesInOrder_2[,c("order","mgi_symbol_character")]
  rownames(cts) = make.unique(genesInOrder_2$mgi_symbol_character)
  return(cts)
}

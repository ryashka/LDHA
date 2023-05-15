make_volcano <- function(df,DESeqOrEdgeR,includeEdgeRUnknownGenes = T){
  require(dplyr)
  require(ggplot2)
  require(ggrepel)
  if(DESeqOrEdgeR == "D"){
    top_genes = rownames(slice_min(df[df$log2FoldChange > 0,], n=10, order_by=padj ))
    bot_genes = rownames(slice_min(df[df$log2FoldChange < 0,], n=10, order_by=padj ))
    df$relevant_genes = ifelse(rownames(df) %in% top_genes,"Top_10","Other_gene")
    df$relevant_genes = ifelse(rownames(df) %in% bot_genes,"Bot_10",df$relevant_genes)
    df$relevant_genes = ifelse(rownames(df) %in% c("Ldha"),"Interesting_Gene",df$relevant_genes)
    df$log2FCOffAxis = ifelse(abs(df$log2FoldChange)>3,TRUE,FALSE)
    df$log2FCWithOffAxisValues = ifelse(df$log2FCOffAxis,
                                        sign(df$log2FCOffAxis)*3,
                                        df$log2FoldChange)
    ggplot(df, aes(x = log2FCWithOffAxisValues, y = -log10(padj), 
                   color = relevant_genes, shape = log2FCOffAxis)) + 
      geom_point() + theme_bw() +
      geom_vline(xintercept = -2, linetype="dotted", color = "grey28", size = 1.5) +
      geom_vline(xintercept = 2, linetype="dotted", color = "grey28", size = 1.5) +
      scale_color_manual(values=c("seagreen","indianred3","grey50","darkorchid4")) + 
      geom_text_repel(aes(label = ifelse(as.character(`relevant_genes`) != "Other_gene",
                                         as.character(rownames(df)),'')),
                      nudge_x = 0.2, nudge_y = 0.1,
                      force = 20, force_pull = 3,
                      segment.size = 0.1, segment.alpha = 1, 
                      point.padding = NA, box.padding = 1,
                      size = 4, max.overlaps = Inf)
  } else {
    if(DESeqOrEdgeR == "E"){
      if(!includeEdgeRUnknownGenes){
        rowsToDrop = c(c(grep("^Gm", row.names(df))),
                       c(grep("^ENSMUSG", row.names(df))),
                       c(grep("-ps$", row.names(df))))
        df = df[-rowsToDrop,]
      }
      
      top_genes = rownames(slice_min(df[df$logFC > 0,], n=10, order_by=FDR ))
      bot_genes = rownames(slice_min(df[df$logFC < 0,], n=10, order_by=FDR ))
      df$relevant_genes = ifelse(rownames(df) %in% top_genes,"Top_10","Other_gene")
      df$relevant_genes = ifelse(rownames(df) %in% bot_genes,"Bot_10",df$relevant_genes)
      df$relevant_genes = ifelse(rownames(df) %in% c("Ldha"),"Interesting_Gene",df$relevant_genes)
      df$logFCOffAxis = ifelse(abs(df$logFC)>3,TRUE,FALSE)
      df$log2FCWithOffAxisValues = ifelse(df$logFCOffAxis,
                                          sign(df$logFC)*3,
                                          df$logFC)
      ggplot(df, aes(x = log2FCWithOffAxisValues, y = -log10(FDR),
                     color = relevant_genes, shape = logFCOffAxis)) + 
        geom_point() + theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_vline(xintercept = -2, linetype="dotted", color = "grey28", size = 1.5) +
        geom_vline(xintercept = 2, linetype="dotted", color = "grey28", size = 1.5) +
        scale_color_manual(values=c("seagreen","indianred3","grey60","darkorchid4")) + 
        geom_text_repel(aes(label = ifelse(as.character(`relevant_genes`) != "Other_gene",
                                           as.character(rownames(df)),'')),
                        nudge_x = 0.2, nudge_y = 0.1,
                        force = 50,
                        segment.size = 0.1, segment.alpha = 1, 
                        point.padding = NA, box.padding = 1,
                        size = 4, max.overlaps = Inf)
    } else {
      return("DESeqOrEdgeR variable must match D or E")
    }
  }
}
deseq_GSE128710 = function(){
  fc = read.delim("~/Documents/Programming/Rahul_Sequencing/GSE128710/combinedFeatureCounts.txt", row.names = 1, skip = 1)
  colnames(fc)[6:29] = str_split_fixed(str_split_fixed(colnames(fc),"\\.",Inf)[6:29,8],"_",Inf)[,1]
  head(fc)
  countMatrix = fc[,c(9:11,15:17)]
  sampleTypes = c(rep("Peak_GC",3),rep("Activated",3))
  library(reshape2)
  sampleTable = melt(data.frame(colnames(countMatrix),sampleTypes))
  colnames(sampleTable) <- c("sampleName","sampleType")
  sampleTable$sampleName = factor(sampleTable$sampleName)
  sampleTable$sampleType = factor(sampleTable$sampleType)
  library("DESeq2")
  rownames(sampleTable) = sampleTable$sampleName
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = sampleTable,
                                design = ~ sampleType)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds = DESeq(dds)
  res = results(dds)
  resOrdered <- res[order(res$padj), ]
  View(resOrdered)
  write.csv(resOrdered,"~/Documents/Programming/Rahul_Sequencing/GSE128710/DESeq2_Ensembl.csv")
  source("~/Documents/Programming/Cornell_DLBCL/DDS_To_Volcano.R")
  GSE128710_mgi = add_MGI(resOrdered)
  write.csv(GSE128710_mgi,"~/Documents/Programming/Rahul_Sequencing/GSE128710/DESeq2_MGI.csv")
  df = as.data.frame(GSE128710_mgi)
  GSE128710_volcano = make_volcano(df,"D",includeEdgeRUnknownGenes = TRUE)
}
deseq_ldha = function(){
  fc = read.delim("~/Documents/Programming/Rahul_Sequencing/Project_10863_B/featureCounts_files/combinedFeatureCounts.txt", row.names = 1, skip = 1)
  head(fc[6:11])
  countMatrix = fc[6:11]
  colnames(countMatrix) = paste0("Sample_",seq(1:6))
  type_1 = rep("Samples_1to3",3)
  type_2 = rep("Samples_4to6",3)
  sampleTypes = c(type_1,type_2)
  
  library(reshape2)
  sampleTable = melt(data.frame(paste0("Sample_",seq(1:6)),sampleTypes))
  colnames(sampleTable) <- c("sampleName","sampleType")
  sampleTable$sampleName = factor(sampleTable$sampleName)
  sampleTable$sampleType = factor(sampleTable$sampleType)
  library("DESeq2")
  rownames(sampleTable) = sampleTable$sampleName
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = sampleTable,
                                design = ~ sampleType)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds = DESeq(dds)
  res = results(dds)
  resOrdered <- res[order(res$padj), ]
  View(resOrdered)
  write.csv(resOrdered,"~/Documents/Programming/Rahul_Sequencing/LDHA_DESeq_Ensembl.csv")
  return(resOrdered)
}
make_gsea = function(ldha_mgi){
  library(fgsea)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  ldha_for_gsea = ldha_mgi %>%
    rownames_to_column(var = "SYMBOL") %>%
    dplyr::select(SYMBOL,stat) %>%
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(stat)) %>%
    deframe()
  library(reactome.db)
  pathways.baderlab <- gmtPathways("~/Documents/Programming/UsefulTables/Mouse_GOBP_AllPathways_with_GO_iea_October_01_2020_symbol.gmt")
  ranks_ldha = ldha_for_gsea$test_statistic
  names(ranks_ldha) = ldha_for_gsea$symbol
  fgseaRes_2 <- fgsea(pathways=pathways.baderlab, stats=ranks_ldha, nperm = 1000, minSize = 15, maxSize = 500)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>% 
    DT::datatable()
  ggplot(fgseaResTidy[abs(fgseaResTidy$NES)>2,], aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways NES from GSEA") + 
    theme_minimal()
  ldha_for_gsea = as.data.frame(ldha_mgi$stat)
  ldha_for_gsea$symbol = rownames(ldha_mgi)
  colnames(ldha_for_gsea)[1] = "test_statistic"
  ldha_for_gsea = ldha_for_gsea[,c('symbol','test_statistic')]
  fgsea(pathways = , 
        stats    = exampleRanks,
        minSize  = 15,
        maxSize  = 500)
}
runGo = function(DESeq2_results,fdr_threshold){
  library(goseq)
  library(org.Mm.eg.db)
  library(fgsea)
  assayed_genes = rownames(DESeq2_results)
  de_genes = rownames(DESeq2_results)[which(DESeq2_results$padj < fdr_threshold)]
  gene_vector = as.integer(assayed_genes %in% de_genes)
  names(gene_vector) = assayed_genes
  supportedOrganisms()[supportedOrganisms()$Genome == "mm9",]
  pwf = nullp(gene_vector, "mm9", "ensGene")
  GO.wall = goseq(pwf, "mm9", "ensGene")
  GO.samp=goseq(pwf,"mm9","ensGene",method="Sampling",repcnt=1000)
  head(Go.wall)
  head(Go.samp)
  plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
       xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
       xlim=c(-3,0))
  abline(0,1,col=3,lty=2)
  GO.BP.wall = goseq(pwf,"mm9","ensGene",test.cats=c("GO:BP"))
  GO.BP.samp = goseq(pwf,"mm9","ensGene",method="Sampling",repcnt=1000,test.cats=c("GO:BP"))
  plot(log10(GO.BP.wall[,2]), log10(GO.BP.samp[match(GO.BP.wall[,1],GO.BP.samp[,1]),2]),
       xlab="log10(BP Wallenius p-values)",ylab="log10(BP Sampling p-values)",
       xlim=c(-3,0))
  abline(0,1,col=3,lty=2)
  GO.wall.enriched = enrichedGO(GO.wall, "GO_Full_List_Wallenius")
  GO.BP.wall.enriched = enrichedGO(GO.BP.wall, "GO_Biological_Processes_Wallenius")
  GO.samp.enriched = enrichedGO(GO.samp, "GO_Full_List_Sampling")
  GO.BP.samp.enriched = enrichedGO(GO.BP.samp,"GO_Biological_Processes_Sampling")
}
seqcp = function(table_directory){
  require(clusterProfiler)
  require(enrichplot)
  require(ggplot2)
  organism = "org.Mm.eg.db"
  #BiocManager::install(organism, character.only = TRUE)
  deseq_output = read.csv(paste0(table_directory,"/deseq_output.csv"), row.names = 1)
  library(organism, character.only = TRUE)
  original_gene_list <- deseq_output$log2FoldChange
  names(original_gene_list) <- deseq_output$ENSEMBL
  gene_list = na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  #Maybe go back here and switch gseaMultilevel by removing nperm
  gse = gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "ENSEMBL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mm.eg.db, 
               pAdjustMethod = "none")
  
  require(DOSE)
  require(ggnewscale)
  gsea_dotplot = dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  eplot(ggplot_dp = gsea_dotplot,
             write_directory_dp = paste0(table_directory,"/Plots/"),
             plotname = "ldha_gsea_dotplot")
  #emap = emapplot(gse, showCategory = 10)
  cnet = cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
  eplot(ggplot_dp = cnet,
             write_directory_dp = paste0(table_directory,"/Plots/"),
             plotname = "ldha_gsea_cnetplot")
  ridge = ridgeplot(gse) + labs(x = "enrichment distribution")
  eplot(ggplot_dp = ridge,
             write_directory_dp = paste0(table_directory,"/Plots/"),
             plotname = "ldha_gsea_ridgeplot")
  for(i in 1:20){
    ldha_gsea = gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)
    doubleplot(ggplot_dp = ldha_gsea,
               write_directory_dp = paste0(table_directory,"/Plots/ENSEMBL/GSEA/"),
               plotname = paste0("ldha_gsea_set_",i))
  }
}
keggcp = function(original_gene_list_seqcp){
  organism = "org.Mm.eg.db"
  original_gene_list = original_gene_list_seqcp
  ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  
  
  
  kegg_organism = "mmu"
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid")
}
  

  
enrichedGO = function(GO.list,nameOfFile){
  library(GO.db)
  library(magrittr)

  enriched.GO=GO.list[p.adjust(GO.list$over_represented_pvalue,
                                        method="BH")<.05,]
  write.csv(file = paste0("~/Documents/Programming/Rahul_Sequencing/enriched_",nameOfFile,".csv"), enriched.GO)
  
  GO.list %>% 
    top_n(10, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc, 
               y=term, 
               colour=over_represented_pvalue, 
               size=numDEInCat)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="GO term", colour="p value", size="Count")
  return(enriched.GO)
}
plotLDHA = function(resOrdered){
  ldha_mgi = add_MGI(resOrdered)
  View(ldha_mgi)
  write.csv(ldha_mgi,"~/Documents/Programming/Rahul_Sequencing/LDHA_DESeq.csv")
  write.csv(add_MGI(counts(dds)), "~/Documents/Programming/Rahul_Sequencing/LDHA_DESeq_counts.csv")
  df = as.data.frame(ldha_mgi)
  with_GM = make_volcano(df)
  ggsave(file = "~/Documents/Programming/Rahul_Sequencing/LDHA_Volcano_with_GM.pdf", 
         plot = with_GM,
         device = "pdf")
  df_without_mtGM = df[-c(grep("^mt",rownames(df)),grep("^Gm",rownames(df))),]
  without_GM = make_volcano(df_without_mtGM)
  ggsave(file = "~/Documents/Programming/Rahul_Sequencing/LDHA_Volcano_without_Gm_or_mt_genes.pdf", 
         plot = without_GM,
         device = "pdf")
  
  ldha_mgi = read.csv("~/Documents/Programming/Rahul_Sequencing/LDHA_DESeq.csv", row.names = 1)
  
}
edgeR_ldha = function(for_camera){
  fc = read.delim("~/Documents/Programming/Rahul_Sequencing/Project_10863_B/featureCounts_files/combinedFeatureCounts.txt", row.names = 1, skip = 1)
  head(fc[6:11])
  countMatrix = fc[6:11]
  colnames(countMatrix) = paste0("Sample_",seq(1:6))
  type_1 = rep("Samples_1to3",3)
  type_2 = rep("Samples_4to6",3)
  sampleTypes = c(type_1,type_2)
  
  library(reshape2)
  sampleTable = melt(data.frame(paste0("Sample_",seq(1:6)),sampleTypes))
  colnames(sampleTable) <- c("sampleName","sampleType")
  sampleTable$sampleName = factor(sampleTable$sampleName)
  sampleTable$sampleType = factor(sampleTable$sampleType)
  library("edgeR")
  library("DESeq2")
  library("DEFormats")
  # rownames(sampleTable) = sampleTable$sampleName
  # 
  # 
  # dds <- DESeqDataSetFromMatrix(countData = countMatrix,
  #                               colData = sampleTable,
  #                               design = ~ sampleType)
  # dge = as.DGEList(dds)
  # voom(counts = dge)
  groups = c(rep("Group1",3),rep("Group2",3))
  if(for_camera){
    countMatrix = add_Entrez(countMatrix,FALSE)
    countMatrix = countMatrix[row.names(countMatrix!="No_Entrez_ID"),]
  }
  d0 = DGEList(counts = countMatrix, group = groups)
  keep = filterByExpr(d0)
  d0 = d0[keep,,keep.lib.sizes=FALSE]
  d0 = calcNormFactors(d0)
  source("~/Documents/Programming/Rahul_Sequencing/Project_10863_B/R_files/UsefulFunctions.R")
  if(!for_camera){
    d = add_MGI(d0)
  } else {
    d = d0
  }
  plotMDS(d, col = as.numeric(groups))
  mm <- model.matrix(~groups)
  d <- estimateDisp(d,mm)
  fit <- glmQLFit(d,mm)
  qlf = glmQLFTest(fit,coef=2)
  if(!for_camera){
    top.table = topTags(qlf, n = length(d$counts))
    write.table(top.table, file = "~/Documents/Programming/Rahul_Sequencing/edgeR_MGI.txt",
                row.names = T, sep = "\t", quote = F)
  } else {
    top.table = topTags(qlf, n = length(d$counts))
    write.table(top.table, file = "~/Documents/Programming/Rahul_Sequencing/edgeR_MGI_No_Entrez.txt",
                row.names = T, sep = "\t", quote = F)
  }
  return(qlf)
}
edgeR_volcano = function(){
  top.table_MGI = read.table("~/Documents/Programming/Rahul_Sequencing/edgeR_MGI.txt")
  edgeR_vp = make_volcano(top.table_MGI, "E",includeEdgeRUnknownGenes = T)
  ggsave("~/Documents/Programming/Rahul_Sequencing/EdgeRVolcano.pdf",
         plot = edgeR_vp, device = "pdf",
         width=10, height=10)
  ggsave("~/Documents/Programming/Rahul_Sequencing/EdgeRVolcano.tiff",
         plot = edgeR_vp, device = "tiff",
         width=10, height=10)
  return(edgeR_vp)
}

deseq_clean_results = function(output_directory){
  source("~/Documents/Programming/Rahul_Sequencing/Project_10863_B/R_files/UsefulFunctions.R")
  deseq_results = deseq_ldha()
  derows = as.data.frame(row.names(deseq_results))
  derows$ENSEMBL = str_split_fixed(derows$`row.names(deseq_results)`,'\\.',Inf)[,1]
  row.names(derows) = derows$ENSEMBL
  conversion_table = add_MGI(derows)
  row.names(deseq_results) = row.names(conversion_table)
  deseq_results$ENSEMBL = conversion_table$ENSEMBL
  write.csv(x = deseq_results, file = paste0(output_directory,"/deseq_output.csv"))
  return(deseq_results)
}
deseq_volcano = function(table_directory){
  deseq_output = read.csv(paste0(table_directory,"/deseq_output.csv"), row.names = 1)
  deseq_vp = make_volcano(deseq_output, "D")
  ggsave(paste0(table_directory,"/Plots/deseq_volcano.pdf"),
         plot = deseq_vp, device = "pdf",
         width=10, height=10)
  ggsave(paste0(table_directory,"/Plots/deseq_volcano.tiff"),
         plot = deseq_vp, device = "tiff",
         width=10, height=10)
  return(deseq_vp)
}

gc_tpm = function(output_file = "~/Documents/Programming/Rahul_Sequencing/GC_Bcell/220715/GC_TPM.xlsx"){
  source("~/Documents/Programming/programs/ScriptsAndFunctions/UsefulFunctions.R")
  fc = read.delim("~/Documents/Programming/Rahul_Sequencing/GC_Bcell/Project_10863_B/featureCounts_files/10863B_combinedFeatureCounts.txt", row.names = 1, skip = 1)
  cts = convertRowsToMGI(fc)
  colnames(cts)[6:11] = paste0("Sample_",seq(1:6))
  tpmvalues = countsToTPM(cts = cts[,6:11], countMatrix = cts)
  xlsx::write.xlsx(x = tpmvalues,output_file)
  return(tpmvalues)
}

ldhb_dotplot = function(){
  source("~/Documents/Programming/Rahul_Sequencing/GC_Bcell/Project_10863_B/R_files/LDHA_Sequencing.R")
  library(ggplot2)
  library(reshape2)
  df = gc_tpm()
  gc_tpm_rel_g = as.data.frame(t(df[c("Aicda","Cd19","Ldha","Ldhb"),]))
  gc_tpm_rel_g$SampleType = c(rep("WT",3), rep("KO",3))
  gc_tall = melt(gc_tpm_rel_g, id.vars = "SampleType")
  gc_tall$SampleType = as.factor(gc_tall$SampleType)
  gc_tall$SampleType = factor(gc_tall$SampleType, levels = c("WT","KO"))
  cd19p = ggplot(gc_tall[gc_tall$variable %in% c("Cd19","Ldhb"),], 
                 aes(variable,value, fill = SampleType)) + 
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1.2)) +
    geom_boxplot(position=position_dodge(1.2)) +
    ylab("Relative Expression (TPM)") +
    scale_fill_manual(values=c("lightgrey", "seagreen")) +
    theme_classic()
  aicdap = ggplot(gc_tall[gc_tall$variable %in% c("Aicda","Ldhb"),], 
                  aes(variable,value, fill = SampleType)) + 
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1.2)) +
    geom_boxplot(position=position_dodge(1.2)) +
    ylab("Relative Expression (TPM)") +
    scale_fill_manual(values=c("lightgrey", "seagreen")) +
    theme_classic()
  ggsave(filename = "~/Documents/Programming/Rahul_Sequencing/GC_Bcell/Cd19andLdhb.pdf",plot = cd19p,device = "pdf")
  ggsave(filename = "~/Documents/Programming/Rahul_Sequencing/GC_Bcell/AicdaandLdhb.pdf",plot = aicdap,device = "pdf")

}


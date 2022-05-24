## Over-Representation Analysis with ClusterProfiler
#  https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####   
  #### Install and load packages ####
    # BiocManager::install("clusterProfiler")
    # BiocManager::install("pathview")
    # BiocManager::install("enrichplot")
    # BiocManager::install("ggupset")
    # BiocManager::install(organism, character.only = TRUE)

    # install.versions(c('RColorBrewer'), c('1.1.2'))

  #### Load Packages ####
    library(clusterProfiler)
    library(enrichplot)
    library(ggupset)
    library(tidyverse)
  
    # Annotations
    organism = "org.Dm.eg.db" ## Genome wide annotation for Fly 
    library(organism, character.only = TRUE)
  
    
  ##### Function setting  ##### 
    ## Call function
    source("FUN_Beautify_ggplot.R")

  ##### Current path and new folder setting  ##### 
    Version = paste0(Sys.Date(),"_","ORA")
    Save.Path = paste0(getwd(),"/",Version)
    dir.create(Save.Path)
  

##### Prepare Input #####
  ## reading in input from deseq2
  df = read.csv("drosphila_example_de.txt", header=TRUE)
  
 ## For the universe in clusterProfiler
  # we want the log2 fold change 
  original_gene_list <- df$log2FoldChange
  
  # name the vector
  names(original_gene_list) <- df$X
  
  # omit any NA values 
  gene_list <- na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  ##---------------------------------------------##
    ## Use pipeline
    gene_list2 <- df %>% drop_na(.,3) %>% select(log2FoldChange) %>% 
                  unlist() %>% as.numeric()
    names(gene_list2) <- df %>% drop_na(.,3) %>% select(X)
    gene_list2 = sort(gene_list2, decreasing = TRUE)
    # Check
    sum(gene_list==gene_list2)
  ##---------------------------------------------##
  
 ## Gene list
  # Exctract significant results (padj < 0.05)
  sig_genes_df = subset(df, padj < 0.05)
  
  # From significant results, we want to filter on log2fold change
  genes <- sig_genes_df$log2FoldChange
  
  # Name the vector
  names(genes) <- sig_genes_df$X
  
  # omit NA values
  genes <- na.omit(genes)
  
  # filter on min log2fold change (log2FoldChange > 2)
  genes <- names(genes)[abs(genes) > 2]

##### Create enrichGO object #####  
  ## Create the object
  go_enrich <- enrichGO(gene = genes,
                        universe = names(gene_list),
                        OrgDb = organism, 
                        keyType = 'ENSEMBL',
                        readable = T,
                        ont = "BP",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
  
##### Outcome #####
  ## Upset Plot
  Upsetplot_GO <- upsetplot(go_enrich)
  Upsetplot_GO
  # https://alanlee.fun/2022/01/08/introducing-upsetplot/
  
  ## Barplot
  Barplot_GO <- barplot(go_enrich, 
                        drop = TRUE, 
                        showCategory = 10, 
                        title = "GO Biological Pathways",
                        font.size = 8)
  Barplot_GO 
  
  Barplot_GO <- Barplot_GO %>% BeautifyggPlot()
  Barplot_GO
  
  ## Dotplot
  Dotplot_GO <- dotplot(go_enrich)
  Dotplot_GO
  Dotplot_GO <- Dotplot_GO %>% BeautifyggPlot(LegPos = c(0.15, 0.65))
  Dotplot_GO 
  
  ## Encrichment map:
  emapplot(go_enrich)

  
  #-----------------------------------------------------------------------------------------------#
  ## error
  ## https://github.com/YuLab-SMU/enrichplot/issues/79
  ## Solution 1
  go_enrich2 <- pairwise_termsim(go_enrich) 
  Emapplot_GO <- emapplot(go_enrich2)
  Emapplot_GO
  
  # ## Solution 2
  # d <- GOSemSim::godata(organism, ont = "BP")    
  # compare_cluster_GO_emap <- enrichplot::pairwise_termsim(go_enrich, semData = d,  method="Wang")
  # emapplot(compare_cluster_GO_emap)
  #-----------------------------------------------------------------------------------------------#
  
  ## Enriched GO induced graph:
  Goplot_GO <- goplot(go_enrich, showCategory = 10)
  Goplot_GO
  
  ## Category Netplot
  # categorySize can be either 'pvalue' or 'geneNum'
  Cnetplot_GO <- cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)
  Cnetplot_GO

  
##### Export #####
  
  Plot.lt <- list(Upsetplot=Upsetplot_GO, Barplot=Barplot_GO, Dotplot=Dotplot_GO,
                  Emapplot=Emapplot_GO, Goplot=Goplot_GO, Cnetplot=Cnetplot_GO)
  
  
  ## Export PDF file
  pdf(
    file = paste0(Save.Path,"/ORA.pdf"),
    width = 10,  height = 8
  )
    Plot.lt
  
  dev.off()
  
  ## Export TIFF file
  for (i in 1:length(Plot.lt)) {
    try({
      tiff(file = paste0(Save.Path,"/",names(Plot.lt)[i],"_ORA.tif"), 
           width = 27, height = 27, units = "cm", res = 200)
      
        print(Plot.lt[i])
      
      graphics.off()
    })
  }
  

# ##### Error part (to be corrected) #####
#   
#   ##### Wordcloud #####
#   # install.packages("wordcloud")
#   library(wordcloud)
#   
#   ## Wordcloud
#   wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
#   wcdf$term<-go_enrich[,2]
#   wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)
#   
#       
## Over-Representation Analysis with ClusterProfiler
#  https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/

##### Load Packages  #####
source("FUN_Package_InstLoad.R")
PKG_Basic.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","DT")
PKG_BiocManager.set <- c("clusterProfiler","enrichplot","ggupset","limma")

FUN_Package_InstLoad(Basic.set = PKG_Basic.set, BiocManager.set = PKG_BiocManager.set)

  # Annotations
  # organism = "org.Dm.eg.db" ## Genome wide annotation for Fly
  organism = "org.Hs.eg.db" ## Genome wide annotation for Human
  # organism = "org.Mm.eg.db" ## Genome wide annotation for Mouse

  library(organism, character.only = TRUE)


##### Function setting  #####
  ## Call function
  source("FUN_Beautify_ggplot.R")

##### Prepare Input #####
  head(DE_Extract.df)
 ## For the universe in clusterProfiler
  # we want the log2 fold change
  original_gene_list <- DE_Extract.df$logFC

  # name the vector
  names(original_gene_list) <- DE_Extract.df$Gene

  # omit any NA values
  gene_list <- na.omit(original_gene_list)

  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)

  # ##---------------------------------------------##
  #   ## Use pipeline
  #   gene_list2 <- DE_Extract.df %>% drop_na(.,3) %>% select(log2FoldChange) %>%
  #                 unlist() %>% as.numeric()
  #   names(gene_list2) <- DE_Extract.df %>% drop_na(.,3) %>% select(X)
  #   gene_list2 = sort(gene_list2, decreasing = TRUE)
  #   # Check
  #   sum(gene_list==gene_list2)
  # ##---------------------------------------------##

 ## Gene list
  # Exctract significant results (padj < 0.05)
  # sig_genes_df = subset(DE_Extract.df, FDR < 0.05)
  sig_genes_df = subset(DE_Extract.df, PValue < 0.05)

  # From significant results, we want to filter on log2fold change
  genes <- sig_genes_df$logFC

  # Name the vector
  names(genes) <- sig_genes_df$Gene

  # omit NA values
  genes <- na.omit(genes)

  # filter on min log2fold change (log2FoldChange > 2)
  genes <- names(genes)[abs(genes) > 1]

##### Create enrichGO object #####
  ## Create the object
  go_enrich <- enrichGO(gene = genes,
                        universe = names(gene_list),
                        OrgDb = organism,
                        keyType = "SYMBOL", #'SYMBOL', #'ENSEMBL'
                        # http://bioconductor.org/help/course-materials/2014/useR2014/Integration.html
                        #readable = T,
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

  Barplot_GO <- Barplot_GO %>% FUN_BeautifyggPlot()
  Barplot_GO

  ## Dotplot
  Dotplot_GO <- dotplot(go_enrich)
  Dotplot_GO
  Dotplot_GO <- Dotplot_GO %>% FUN_BeautifyggPlot(LegPos = c(0.15, 0.65))
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

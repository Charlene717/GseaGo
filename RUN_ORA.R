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
  head(DEG_Extract.df)
 ## For the universe in clusterProfiler
  # we want the log2 fold change
  original_gene_list <- DEG_Extract.df$logFC

  # name the vector
  names(original_gene_list) <- DEG_Extract.df$Gene

  # omit any NA values
  ORA_GeneDiff_All <- na.omit(original_gene_list)

  # sort the list in decreasing order (required for clusterProfiler)
  ORA_GeneDiff_All = sort(ORA_GeneDiff_All, decreasing = TRUE)
  rm(original_gene_list)

  # ##---------------------------------------------##
  #   ## Use pipeline
  #   ORA_GeneDiff_All2 <- DEG_Extract.df %>% drop_na(.,3) %>% select(log2FoldChange) %>%
  #                 unlist() %>% as.numeric()
  #   names(ORA_GeneDiff_All2) <- DEG_Extract.df %>% drop_na(.,3) %>% select(X)
  #   ORA_GeneDiff_All2 = sort(ORA_GeneDiff_All2, decreasing = TRUE)
  #   # Check
  #   sum(ORA_GeneDiff_All==ORA_GeneDiff_All2)
  # ##---------------------------------------------##

 ## Gene list
  # Exctract significant results (padj < 0.05)
  # Sig_GeneExp.df = subset(DEG_Extract.df, FDR < 0.05)
  Sig_GeneExp.df = subset(DEG_Extract.df, PValue < 0.05)

  # From significant results, we want to filter on log2fold change
  ORA_GeneList_Sig <- Sig_GeneExp.df$logFC

  # Name the vector
  names(ORA_GeneList_Sig) <- Sig_GeneExp.df$Gene

  # omit NA values
  ORA_GeneList_Sig <- na.omit(ORA_GeneList_Sig)

  # filter on min log2fold change (log2FoldChange > 2)
  ORA_GeneList_Sig <- names(ORA_GeneList_Sig)[abs(ORA_GeneList_Sig) > 1]

##### Create enrichGO object #####
  ## Create the object
  ORA_GO_Result <- enrichGO(gene = ORA_GeneList_Sig,
                        universe = names(ORA_GeneDiff_All),
                        OrgDb = organism,
                        keyType = "SYMBOL", #'SYMBOL', #'ENSEMBL'
                        # http://bioconductor.org/help/course-materials/2014/useR2014/Integration.html
                        #readable = T,
                        ont = "BP",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.10)

##### Outcome #####
  ## Upset Plot
  Upsetplot_GO <- upsetplot(ORA_GO_Result)
  Upsetplot_GO
  # https://alanlee.fun/2022/01/08/introducing-upsetplot/

  ## Barplot
  Barplot_GO <- barplot(ORA_GO_Result,
                        drop = TRUE,
                        showCategory = 10,
                        title = "GO Biological Pathways",
                        font.size = 8)
  Barplot_GO

  Barplot_GO <- Barplot_GO %>% FUN_BeautifyggPlot()
  Barplot_GO

  ## Dotplot
  Dotplot_GO <- dotplot(ORA_GO_Result)
  Dotplot_GO
  Dotplot_GO <- Dotplot_GO %>% FUN_BeautifyggPlot(LegPos = c(0.15, 0.65))
  Dotplot_GO

  ## Encrichment map:
  try({emapplot(ORA_GO_Result)})

  #-----------------------------------------------------------------------------------------------#
  ## error
  ## https://github.com/YuLab-SMU/enrichplot/issues/79
  ## Solution 1
  ORA_GO_Result <- pairwise_termsim(ORA_GO_Result)
  Emapplot_GO <- emapplot(ORA_GO_Result)
  Emapplot_GO

  # ## Solution 2
  # d <- GOSemSim::godata(organism, ont = "BP")
  # compare_cluster_GO_emap <- enrichplot::pairwise_termsim(ORA_GO_Result, semData = d,  method="Wang")
  # emapplot(compare_cluster_GO_emap)
  #-----------------------------------------------------------------------------------------------#

  ## Enriched GO induced graph:
  Goplot_GO <- goplot(ORA_GO_Result, showCategory = 10)
  Goplot_GO

  ## Category Netplot
  # categorySize can be either 'pvalue' or 'geneNum'
  Cnetplot_GO <- cnetplot(ORA_GO_Result, categorySize = "pvalue", foldChange = ORA_GeneDiff_All)
  Cnetplot_GO


##### Export #####

  ORA_Plot.lt <- list(Upsetplot=Upsetplot_GO, Barplot=Barplot_GO, Dotplot=Dotplot_GO,
                  Emapplot=Emapplot_GO, Goplot=Goplot_GO, Cnetplot=Cnetplot_GO)


  ## Export PDF file
  pdf(
    file = paste0(Save_Path,"/ORA.pdf"),
    width = 10,  height = 8
  )

  try({print(ORA_Plot.lt)})

  dev.off()

  ## Export TIFF file
  for (i in 1:length(ORA_Plot.lt)) {
    try({
      tiff(file = paste0(Save_Path,"/",names(ORA_Plot.lt)[i],"_ORA.tif"),
           width = 27, height = 27, units = "cm", res = 200)

        print(ORA_Plot.lt[i])

      graphics.off()
    })
  }
  rm(i)

  rm(organism)

  rm(Upsetplot_GO, Barplot_GO, Dotplot_GO, Emapplot_GO, Goplot_GO,Cnetplot_GO,
     Sig_GeneExp.df)

# ##### Error part (to be corrected) #####
#
#   ##### Wordcloud #####
#   # install.packages("wordcloud")
#   library(wordcloud)
#
#   ## Wordcloud
#   wcdf<-read.table(text=ORA_GO_Result$GeneRatio, sep = "/")[1]
#   wcdf$term<-ORA_GO_Result[,2]
#   wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)
#
#

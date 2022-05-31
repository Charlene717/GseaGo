## RNA-seq analysis in R
## Ref: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

## Gene Set Enrichment Analysis with ClusterProfiler
## Ref: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

## Visualzation of GSEA results
## Ref: https://rpubs.com/shbrief/gsea_263

## GSEA Chard Liu
## Ref: http://rstudio-pubs-static.s3.amazonaws.com/514990_9690f31b5ef7488bb4f0bb6c10ac4da8.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Current path and new folder setting* #####
  ProjectName = "Example"
  Sampletype = "Mm"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

  ## Import information
  InputFolder = "Input_files_10x"
  InputAnno = "PBMC_Ano.csv"

  InputGSEA = "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"

##### Parameter setting* #####
  # Set the desired organism
  organism = "org.Mm.eg.db" # organism = "org.Hs.eg.db", "org.Mm.eg.db", "org.Dm.eg.db"

##### Load Packages  #####
  #### Basic installation ####
  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","ggplot2","msigdbr","forcats","ggupset")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c(organism, "org.Hs.eg.db", "org.Mm.eg.db", "org.Dm.eg.db",
                   "fgsea","clusterProfiler","enrichplot","pathview","enrichplot")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  options(stringsAsFactors = FALSE)

  #### Github installation ####
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
  detach("package:clusterProfiler.dplyr", unload = TRUE)
  devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
  devtools::install_github("lionel-/ggstance")
  library(clusterProfiler.dplyr)
  library(ggstance)


##### Load Files #####
  ## Load data
  setwd("../")
  load("Demo_data/Robjects/Annotated_Results_LvV.RData")
  ## Load pathways
  load("Demo_data/Robjects/mouse_H_v5.RData")
  pathwaysH <- Mm.H

  setwd("GSEA_Analysis")

##### GSEA analysis (fgsea) #####
  library(fgsea)

  #### Create ranks ####
    gseaDat <- filter(shrinkLvV, !is.na(Entrez))
    ranks <- gseaDat$logFC
    names(ranks) <- gseaDat$Entrez
    head(ranks)

    # Plot the ranked fold changes.
    barplot(sort(ranks, decreasing = T))

    pdf(file = paste0(Save.Path,"/",ProjectName,"_Rankbarplot.pdf"),
        width = 17,  height = 7
    )
    barplot(sort(ranks, decreasing = T))
    dev.off()


  #### Conduct analysis ####
    # fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
    # ## fgsea: What does fgseaMultilevel argument sampleSize mean/when to change it?
    # ## https://www.biostars.org/p/479821/

    fgseaRes <- fgseaMultilevel(pathwaysH, ranks, minSize=15, maxSize = 500, nPermSimple = 1000)
    ## Error when running parallelized process: Warning in serialize... package:stats may not be available when loading
    ## https://community.rstudio.com/t/error-when-running-parallelized-process-warning-in-serialize-package-stats-may-not-be-available-when-loading/110573
    ## https://stackoverflow.com/questions/27623901/r-warning-packagestats-may-not-be-available-when-loading

    head(fgseaRes[order(padj, -abs(NES)), ], n=10)

  #### Enrichment score plot ####
    plotEnrichment(pathwaysH[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]], ranks)
    dev.off()

  #### GSEA table plot ####
    topUp <- fgseaRes %>%
      filter(ES > 0) %>%
      top_n(10, wt=-padj)
    topDown <- fgseaRes %>%
      filter(ES < 0) %>%
      top_n(10, wt=-padj)
    topPathways <- bind_rows(topUp, topDown) %>%
      arrange(-ES)
    plotGseaTable(pathwaysH[topPathways$pathway],
                  ranks,
                  fgseaRes,
                  gseaParam = 0.5)
    dev.off()

    pdf(file = paste0(Save.Path,"/",ProjectName,"_tableplot.pdf"),
        width = 12,  height = 7
    )
    plotGseaTable(pathwaysH[topPathways$pathway],
                  ranks,
                  fgseaRes,
                  gseaParam = 0.5)
    dev.off()

##### GSEA analysis (clusterProfiler) #####
    library(clusterProfiler)
    #### Conduct analysis2 ####

    geneList <- sort(ranks, decreasing = T)
    ## MSigDB_C2
    library(msigdbr)
    msigdbr_species()
    m_c2 <- msigdbr(species = "Mus musculus", category = "C2") %>%
      dplyr::select(gs_name, entrez_gene)
    msC2_2 <- GSEA(geneList, TERM2GENE = m_c2)

    #### Visualization ####

    ## 2.1 Barplot
    library(clusterProfiler.dplyr) # devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
    y <- mutate(msC2_2, ordering = abs(NES)) %>%
      arrange(desc(ordering))

    library(ggstance) # devtools::install_github("lionel-/ggstance")
    library(enrichplot) # BiocManager::install("enrichplot")
    library(forcats)
    library(ggplot2)

    n <- 10
    y_bar <- group_by(y, sign(NES)) %>% slice(1:n)

    ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalues), showCategory=(n*2)) +
           geom_barh(stat='identity') +
           scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
           theme_minimal() + ylab(NULL) -> P2_1.Barplot
    P2_1.Barplot

    ## 2.2 Dotplot
    P2_2.Dotplot <- dotplot(msC2_2, showCategory = 10, font.size = 8,
                            x = "GeneRatio",   # option -> c("GeneRatio", "Count")
                            color = "p.adjust")   # option -> c("pvalue", "p.adjust", "qvalue")
    P2_2.Dotplot

    ## 2.3 Gene-Concept Network
    n <- 3
    p1 <- cnetplot(msC2_2, showCategory = (n*2), colorEdge = TRUE, node_label = "category")
    P2_3.GeneNet <- cowplot::plot_grid(p1, ncol=1, labels=LETTERS[1], rel_widths=c(1))
    P2_3.GeneNet

    ## 2.4 Heatmap-like functional classification
    ## Ref:https://bioc.ism.ac.jp/packages/3.7/bioc/vignettes/enrichplot/inst/doc/enrichplot.html#:~:text=upsetplot(ego)-,Heatmap%2Dlike%20functional%20classification,easy%20to%20identify%20expression%20patterns.
    msC2_3 <- msC2_2[msC2_2@geneList %in% msC2_2@geneList[1:50],]
    heatplot(msC2_3, showCategory = 30, foldChange=geneList)

    ## 2.5 Enrichment Map
    p2 <- emapplot(pairwise_termsim(y), showCategory = 20)
    p2
    P2_5.EnrichMap <- cowplot::plot_grid(p2, ncol = 1, lables = LETTERS[1])
    P2_5.EnrichMap


    ## 2.6 UpSet Plot
    library(ggupset) # https://github.com/const-ae/ggupset
    P2_6.UpSetPlot <- upsetplot(msC2_2)
    P2_6.UpSetPlot


    ## 2.7 ridgeline plot for expressiong distribution
    P2_7.RidgelinePlot <- ridgeplot(msC2_2)
    P2_7.RidgelinePlot


    ## 2.8 gseaplot
    y2 <- arrange(msC2_2, desc(NES))

    p1 <- gseaplot(y2, geneSetID = 1, title = y2$Description[1])   # max NES
    n <- nrow(y2)
    p2 <- gseaplot(y2, geneSetID = n, title = y2$Description[n])   # min NES
    P2_8.gseaplot <- cowplot::plot_grid(p1, p2, ncol = 1, labels = LETTERS[1:2])
    P2_8.gseaplot


      ## 2.8.1 gseaplot2
      p3 <- gseaplot2(y2, geneSetID = 1, title = y2$Description[1])
      p4 <- gseaplot2(y2, geneSetID = n, title = y2$Description[n])   # min NES
      cowplot::plot_grid(p3, p4, ncol = 1, labels = LETTERS[1:2])
      gseaplot2(y2, geneSetID = 10, title = y2$ID[10],color="red",pvalue_table=T)

      ## Use keyword (Overlay graphics)
      keyword <- "breast"
      ind <- grep(keyword, msC2_2$Description, ignore.case = TRUE)
      P2_8_1.gseaplot_Keyword <- gseaplot2(msC2_2, geneSetID = ind, pvalue_table=T,
                                           base_size = 16,title = keyword,subplots = 1:3,
                                           rel_heights = c(3, 0.5, 1))
      P2_8_1.gseaplot_Keyword

      ## Overlay graphics by ID
      P2_8_1.gseaplot_IDTop <- gseaplot2(y2, geneSetID = 1:10, pvalue_table=T,
                                           base_size = 16, title = "TOP 10",subplots = 1:3,
                                           rel_heights = c(3, 0.5, 1))
      P2_8_1.gseaplot_IDTop

      P2_8_1.gseaplot_IDBot <- gseaplot2(y2, geneSetID = (length(y2@result[["ID"]])-4):length(y2@result[["ID"]]),
                                         pvalue_table=T,
                                        base_size = 16, title = "Bottom 10",subplots = 1:3,
                                        rel_heights = c(3, 0.5, 1))
      P2_8_1.gseaplot_IDBot

      gseaplot2(y2, geneSetID = c(1:5,(length(y2@result[["ID"]])-4):
                length(y2@result[["ID"]])), pvalue_table=T,
                base_size = 16, title = "TOP 5 and Bottom 5",subplots = 1:3,
                rel_heights = c(3, 0.5, 1)) -> P2_8_1.gseaplot_IDTopBot
      P2_8_1.gseaplot_IDTopBot


      # ## Modify the function
      # ## https://www.biostars.org/p/9470087/
      # ## trace("gseaplot2", edit = TRUE)
      # ##> pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")] -> pd <- x[geneSetID, c("Description","NES", "pvalue", "p.adjust")]
      # gseaplot2(y2, geneSetID = 1:10, pvalue_table=T)

      ## 2.8.2 gsearank
      P2_8_2.gsearank <- gsearank(y2, geneSetID = 1, title = y2$Description[1])


    ## 2.9 PubMed trend of enriched terms
    terms <- msC2_2$Description[1:3]
    P2_9.PubMedTrend <- pmcplot(terms, 2010:2017, proportion=FALSE)
    P2_9.PubMedTrend


  #### Export PDF files ####
    pdf(file = paste0(Save.Path,"/",ProjectName,"_clusterProfiler.pdf"),
        width = 20,  height = 15
    )
      P2_1.Barplot
      P2_2.Dotplot
      P2_3.GeneNet
      P2_5.EnrichMap
      P2_6.UpSetPlot
      P2_7.RidgelinePlot
      P2_8.gseaplot
      P2_8_1.gseaplot_Keyword
      P2_8_1.gseaplot_IDTop
      P2_8_1.gseaplot_IDBot
      P2_8_1.gseaplot_IDTopBot
      P2_8_2.gsearank
      P2_9.PubMedTrend
    dev.off()




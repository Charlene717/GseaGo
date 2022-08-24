## Run GSEA analysis in R


## RNA-seq analysis in R
## Ref: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

## Gene Set Enrichment Analysis with ClusterProfiler
## Ref: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

## Visualzation of GSEA results
## Ref: https://rpubs.com/shbrief/gsea_263

## GSEA Chard Liu
## Ref: http://rstudio-pubs-static.s3.amazonaws.com/514990_9690f31b5ef7488bb4f0bb6c10ac4da8.html

FUN_GSEA_Analysis = function(DE_Extract.df, pathwayGeneSet = Pathway.all,
                             TarGeneName = TarGene_name, GroupMode = Mode_Group,
                             Species = "Homo sapiens", # Speices type can check by msigdbr_species()
                             Save.Path = Save.Path, SampleName = SampleName, AnnoName = "C2"
){


  ##### Load Packages  #####
    #### Basic installation ####
    ## Check whether the installation of those packages is required from basic
    Package.set <- c("tidyverse","ggplot2")
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
    Package.set <- c(organism,"fgsea","clusterProfiler","enrichplot","pathview")
    for (i in 1:length(Package.set)) {
      if (!requireNamespace(Package.set[i], quietly = TRUE)){
        BiocManager::install(Package.set[i])
      }
    }
    ## Load Packages
    lapply(Package.set, library, character.only = TRUE)
    rm(Package.set,i)

    options(stringsAsFactors = FALSE)


# #************************************************************************************************************************#
#     ##### Parameter setting* #####
#     ## Set the desired organism
#     organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db")   ##  c("org.Dm.eg.db")
#
#   ##### GSEA analysis (fgsea) #####
#
#
#     #### Create ranks ####
#       # gseaDat2 <- filter(shrinkLvV, !is.na(shrinkLvV[,"Entrez"]))
#       # gseaDat <- filter(shrinkLvV, !is.na(Entrez))
#       # ranks <- gseaDat$logFC
#       # names(ranks) <- gseaDat$Entrez
#       # head(ranks)
#       #
#       # # Plot the ranked fold changes.
#       # barplot(sort(ranks, decreasing = T))
#
#       ranks <- DE_Extract.df[,ThrSet[["LogFC"]][1]]
#       names(ranks) <- row.names(DE_Extract.df)
#       head(ranks)
#
#       # Plot the ranked fold changes.
#       barplot(sort(ranks, decreasing = T))
#
#     #### Transform pathwayGeneSet to list ####
#       # pathwayGeneSet.lt <- as.list(pathwayGeneSet[,c(-1,-2)])
#
#       # pathwayGeneSet.lt <- split(pathwayGeneSet[,c(-1,-2)], 1:nrow(pathwayGeneSet))
#       # names(pathwayGeneSet.lt) <- pathwayGeneSet[,1]
#
#       ## No use
#       # pathwayGeneSet.lt2 <- pathwayGeneSet.lt[!is.na(pathwayGeneSet.lt)]
#       # pathwayGeneSet.lt2 <- lapply(pathwayGeneSet.lt, function(x) unique(x))
#
#       ## EXist ""
#       # pathwayGeneSet.lt <- apply(pathwayGeneSet[,c(-1,-2)], 1, function(a)unique(a))
#       # names(pathwayGeneSet.lt) <- pathwayGeneSet[,1]
#
#       pathwayGeneSet.lt <- apply(pathwayGeneSet[,c(-1,-2)], 1, function(a)unique(a[a != ""]))
#       names(pathwayGeneSet.lt) <- pathwayGeneSet[,1]
#
#     #### Conduct analysis ####
#       fgseaRes <- fgsea(pathwayGeneSet.lt, ranks, minSize=15, maxSize = 500, nperm=1000)
#
#       ## fgsea: What does fgseaMultilevel argument sampleSize mean/when to change it?
#       ## https://www.biostars.org/p/479821/
#
#       fgseaRes <- fgseaMultilevel(pathwayGeneSet.lt, ranks, minSize=15, maxSize = 500)
#       ## Error when running parallelized process: Warning in serialize... package:stats may not be available when loading
#       ## https://community.rstudio.com/t/error-when-running-parallelized-process-warning-in-serialize-package-stats-may-not-be-available-when-loading/110573
#       ## https://stackoverflow.com/questions/27623901/r-warning-packagestats-may-not-be-available-when-loading
#
#       head(fgseaRes[order(padj, -abs(NES)), ], n=10)
#
#     #### Enrichment score plot ####
#       # plotEnrichment(pathwayGeneSet[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]], ranks)
#
#       plotEnrichment(pathwayGeneSet[pathwayGeneSet[,1]=="HALLMARK_ESTROGEN_RESPONSE_EARLY",], ranks)
#       dev.off()
#
#     #### GSEA table plot ####
#       # topUp <- fgseaRes %>%
#       #   filter(ES > 0) %>%
#       #   top_n(10, wt=-padj)
#       # topDown <- fgseaRes %>%
#       #   filter(ES < 0) %>%
#       #   top_n(10, wt=-padj)
#       # topPathways <- bind_rows(topUp, topDown) %>%
#       #   arrange(-ES)
#       # plotGseaTable(pathwayGeneSet[topPathways$pathway],
#       #               ranks,
#       #               fgseaRes,
#       #               gseaParam = 0.5)
#       # dev.off()
#
#       topUp <- fgseaRes %>%
#                filter(NES > 0) %>%
#                top_n(10, wt=NES)
#       topDown <- fgseaRes %>%
#                  filter(NES < 0) %>%
#                  top_n(10, wt=-NES)
#       topPathways <- bind_rows(topUp, topDown) %>%
#                      arrange(-NES)
#       plotGseaTable(pathwayGeneSet.lt[topPathways$pathway],
#                     ranks,
#                     fgseaRes,
#                     gseaParam = 0.5)
#       dev.off()
#
#       ## 2.1 Barplot
#       library(ggstance)
#       ggplot(topPathways, aes(NES, fct_reorder(pathway, NES), fill = pval), showCategory=(n*2)) +
#              geom_barh(stat='identity') +
#              scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
#              theme_minimal() + ylab(NULL)


#************************************************************************************************************************#
  ##### GSEA analysis #####
      #### Conduct analysis2 ####
      library(DOSE)
      # data(geneList)
      # x <- gseDO(geneList)
      # gseaplot(x, geneSetID=1)

      # geneList <- sort(ranks, decreasing = T)
      # #geneList2 <- names(ranks) %>% as.numeric()
      # x <- gseDO(geneList)
      # gseaplot(x, geneSetID=1)

      geneList <- sort(ranks, decreasing = T)

      #### Use online pathwayGeneSet ####
      ## MSigDB_C2
      library(msigdbr)
      msigdbr_species()
      m_c2 <- msigdbr(species = Species, category = "C2") %>%
              dplyr::select(gs_name, gene_symbol, entrez_gene)

      #### Transform Customized pathwayGeneSet ####
      Temp <- pathwayGeneSet[,-2]
      m_c2 <- Temp %>% pivot_longer(cols=2:ncol(.),names_to = "Temp", values_to = "Gene") %>%
                       select(cols=c(1,3))

      m_c2 <- m_c2[!m_c2$cols2 == "",]

      rm(Temp)

      #### RUN GSEA ####
      GSEA_Result <- GSEA(geneList, TERM2GENE = m_c2) #       GSEA_Result <- GSEA(geneList, TERM2GENE = m_c2)

      #### Visualization ####

      ## 2.1 Barplot
      library(clusterProfiler.dplyr) # devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
      y <- mutate(GSEA_Result, ordering = abs(NES)) %>%
           arrange(desc(ordering))

      library(ggstance)
      library(enrichplot)
      library(forcats)
      library(ggplot2)

      n <- 10
      y_bar <- group_by(y, sign(NES)) %>% slice(1:n)

      ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalues), showCategory=(n*2)) +
            geom_barh(stat='identity') +
            scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
            theme_minimal() + ylab(NULL)

      ## 2.2 Dotplot
      dotplot(GSEA_Result, showCategory = 10, font.size = 8,
              x = "GeneRatio",   # option -> c("GeneRatio", "Count")
              color = "p.adjust")   # option -> c("pvalue", "p.adjust", "qvalue")

      ## 2.3 Gene-Concept Network
      n <- 3
      p1 <- cnetplot(GSEA_Result, showCategory = (n*2), colorEdge = TRUE, node_label = "category")
      cowplot::plot_grid(p1, ncol=1, labels=LETTERS[1], rel_widths=c(1))

      ## 2.4 Heatmap-like functional classification



      ## 2.5 Enrichment Map
      p2 <- emapplot(pairwise_termsim(y), showCategory = 10)
      cowplot::plot_grid(p2, ncol = 1, lables = LETTERS[1])


      ## 2.6 UpSet Plot
      library(ggupset)
      upsetplot(GSEA_Result)

      ## 2.7 ridgeline plot for expressiong distribution
      ridgeplot(GSEA_Result)

      ## 2.8 gseaplot
      y2 <- arrange(GSEA_Result, desc(NES))

      p1 <- gseaplot(y2, geneSetID = 1, title = y2$Description[1])   # max NES
      n <- nrow(y2)
      p2 <- gseaplot(y2, geneSetID = n, title = y2$Description[n])   # min NES
      cowplot::plot_grid(p1, p2, ncol = 1, labels = LETTERS[1:2])

      ## 2.8.1 gseaplot2
      p3 <- gseaplot2(y2, geneSetID = 1, title = y2$Description[1])   # max NES
      p4 <- gseaplot2(y2, geneSetID = n, title = y2$Description[n])   # min NES
      cowplot::plot_grid(p3, p4, ncol = 1, labels = LETTERS[1:2])

      ## Use keyword (Overlay graphics)
      try({
        keyword <- "breast"
        ind <- grep(keyword, GSEA_Result$Description, ignore.case = TRUE)
        gseaplot2(GSEA_Result, geneSetID = ind, title = GSEA_Result$Description[ind])
      })

      ## Overlay graphics by ID
      gseaplot2(y2, geneSetID = 1:10)

      ## 2.8.2 gsearank
      gsearank(y2, geneSetID = 1, title = y2$Description[1])


      ## 2.9 PubMed trend of enriched terms
      terms <- GSEA_Result$Description[1:3]
      pmcplot(terms, 2010:2017, proportion=FALSE)
      # dev.off()

    # ##### Export Result #####
    #   pdf(
    #     file = paste0(Save.Path,"/",SampleName,"_",TarGeneName,"_DensityPlot.pdf"),
    #     width = 10,  height = 8
    #   )
    #   print(TGeneDen_SD.p)
    #   print(TGeneDen_Q.p)
    #   print(TGeneDen_SD_Q.p)
    #
    #   dev.off()

    #### Output ####
      Output <- list()
      Output[["GSEA_Result"]] <- GSEA_Result
      Output[["geneList_ranks"]] <- geneList

    return(Output)


return(Output)

}

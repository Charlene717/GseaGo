## RNA-seq analysis in R
## Ref: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

## Gene Set Enrichment Analysis with ClusterProfiler
## Ref: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

## Visualzation of GSEA results
## Ref: https://rpubs.com/shbrief/gsea_263

## GSEA Chard Liu
## Ref: http://rstudio-pubs-static.s3.amazonaws.com/514990_9690f31b5ef7488bb4f0bb6c10ac4da8.html

FUN_GSEA_ANAL = function(DE_Extract.df, CMGeneSet = Pathway.all, pathwayGeneSet = "C2",
                         NumGenesetsPlt = 10,
                         TarGeneName = TarGene_name,
                         ThrSet = DEGThr.lt, Species = "Homo sapiens", # Speices type can check by msigdbr_species()
                         Save.Path = Save.Path, SampleName = SampleName, AnnoName = "C2"
){

  ##### Load Packages  #####
    source("FUN_Package_InstLoad.R")
    Basic.set <- c("tidyverse","ggplot2","eoffice")
    BiocManager.set <- c("clusterProfiler","enrichplot","pathview")
    Package_InstLoad(Basic.set = Basic.set, BiocManager.set = BiocManager.set)

#************************************************************************************************************************#
  ##### GSEA analysis #####
    #### Conduct analysis2 ####
    library(DOSE)

    ranks <- DE_Extract.df[,ThrSet[["LogFC"]][1]]
    names(ranks) <- row.names(DE_Extract.df)
    head(ranks)

    geneList <- sort(ranks, decreasing = T)

    #### Use online GeneSet ####
    ## MSigDB_C2
    library(msigdbr)
    # msigdbr_species()
    # m_c2 <- msigdbr(species = Species, category = "C2") %>%
    #         dplyr::select(gs_name, gene_symbol, entrez_gene)

    #### Transform Customized CMGeneSet ####
    # CMGeneSet <- data.frame(Geneset = row.names(CMGeneSet),CMGeneSet)
    Temp <- CMGeneSet[,-2]
    m_c2 <- Temp %>% pivot_longer(cols=2:ncol(.),names_to = "Temp", values_to = "Gene") %>%
                     dplyr::select(cols=c(1,3))

    m_c2 <- m_c2[!m_c2$cols2 == "",]

    rm(Temp)

    #### RUN GSEA ####
    GSEA_Result <- GSEA(geneList, TERM2GENE = m_c2) #       GSEA_Result <- GSEA(geneList, TERM2GENE = m_c2)
    # nPerm=1000
    # https://support.bioconductor.org/p/99810/
    # https://bioinformatics.stackexchange.com/questions/149/are-fgsea-and-broad-institute-gsea-equivalent

    #### Visualization ####
    ## 2.1 Barplot
    library(clusterProfiler.dplyr) # devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
    y <- mutate(GSEA_Result, ordering = abs(NES)) %>%
         arrange(desc(ordering))

    library(ggstance)
    library(enrichplot)
    library(forcats)
    library(ggplot2)


    y_bar <- group_by(y, sign(NES)) %>% slice(1:NumGenesetsPlt)

    p1 <- ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalues), showCategory=(NumGenesetsPlt*2)) +
                 geom_barh(stat='identity') +
                 scale_fill_continuous(low='#d45772', high='#3b74bf', guide = guide_colorbar(reverse=TRUE)) +
                 theme_minimal() + ylab(NULL)


    ## 2.2 Dotplot
    p2 <- dotplot(GSEA_Result, showCategory = NumGenesetsPlt, font.size = 8,
                  x = "GeneRatio",   # option -> c("GeneRatio", "Count")
                  color = "p.adjust")   # option -> c("pvalue", "p.adjust", "qvalue")


    ## 2.3 Gene-Concept Network
    # n <- 3
    p3 <- cnetplot(GSEA_Result, showCategory = (NumGenesetsPlt*2), colorEdge = TRUE, node_label = "category")
    cowplot::plot_grid(p3, ncol=1, labels=LETTERS[1], rel_widths=c(1))

    ## 2.4 Heatmap-like functional classification



    ## 2.5 Enrichment Map
    p5 <- emapplot(pairwise_termsim(y), showCategory = NumGenesetsPlt/2)
    cowplot::plot_grid(p5, ncol = 1, lables = LETTERS[1])


    ## 2.6 UpSet Plot
    library(ggupset)
    p6 <- upsetplot(GSEA_Result, n = NumGenesetsPlt)

    ## 2.7 ridgeline plot for expressiong distribution
    p7 <- ridgeplot(GSEA_Result) +
          theme(axis.text.y = element_text(size = 10))

    ## 2.8 gseaplot
    y2 <- arrange(GSEA_Result, desc(NES))

    # p1 <- gseaplot(y2, geneSetID = 1, title = y2$Description[1])   # max NES
    # n <- nrow(y2)
    # p2 <- gseaplot(y2, geneSetID = n, title = y2$Description[n])   # min NES
    # cowplot::plot_grid(p1, p2, ncol = 1, labels = LETTERS[1:2])

    ## 2.8.1 gseaplot2
    n <- nrow(y2)
    p8_1 <- gseaplot2(y2, geneSetID = 1, title = y2$Description[1])   # max NES
    p8_1_2 <- gseaplot2(y2, geneSetID = 2, title = y2$Description[2])   # Sec NES

    p8_2 <- gseaplot2(y2, geneSetID = n, title = y2$Description[n])   # min NES
    p8_2_2 <- gseaplot2(y2, geneSetID = (n-1), title = y2$Description[(n-1)])   # min NES

    p8A <- cowplot::plot_grid(p8_1,p8_1_2, p8_2,p8_2_2, ncol = 2, labels = LETTERS[1:4])

    # ## Use keyword (Overlay graphics)
    # try({
    #   keyword <- "breast"
    #   ind <- grep(keyword, GSEA_Result$Description, ignore.case = TRUE)
    #   gseaplot2(GSEA_Result, geneSetID = ind, title = GSEA_Result$Description[ind])
    # })

    ## Overlay graphics by ID
    # p8B <- gseaplot2(y2, geneSetID = 1:10)
    p8B <- gseaplot2(y2, geneSetID = 1:NumGenesetsPlt)

    # ## 2.8.2 gsearank
    # gsearank(y2, geneSetID = 1, title = y2$Description[1])


    ## 2.9 PubMed trend of enriched terms
    terms <- GSEA_Result$Description[1:3]
    p9 <- pmcplot(terms, 2010:2017, proportion=FALSE)
    # dev.off()

  ##### Export Result #####
    pdf(
      file = paste0(Save.Path,"/GSEAResult_",AnnoName,"_",SampleName,"_",TarGeneName,".pdf"),
      width = 15,  height = 10
    )

      print(p1 %>% BeautifyggPlot(LegPos = c(0.9, 0.2),AxisTitleSize=1.7,YtextSize=14) )
      print(p2 %>% BeautifyggPlot(LegPos = c(0.9, 0.3),AxisTitleSize=1.7))
      print(p3)
      # print(p4)
      print(p5)
      print(p6)
      print(p7)
      print(p8A)
      print(p8B)
      print(p9)

    dev.off()

  #### Output ####
    Output <- list()
    Output[["GSEA_Result"]] <- GSEA_Result
    Output[["geneList_ranks"]] <- geneList
    Output[["GSEABar_Plot"]] <- p1
    Output[["GSEADot_Plot"]] <- p2
    Output[["p3"]] <- p3
    # Output[["p4"]] <- p4
    Output[["p5"]] <- p5
    Output[["UpSet_Plot"]] <- p6
    Output[["p7"]] <- p7
    Output[["Gsea_Plot"]] <- p8A
    Output[["OverlayGsea_Plot"]] <- p8B
    Output[["p9"]] <- p9

  return(Output)


}


# ##### Load Packages  #####
#   #### Basic installation ####
#   ## Check whether the installation of those packages is required from basic
#   Package.set <- c("tidyverse","ggplot2")
#   for (i in 1:length(Package.set)) {
#     if (!requireNamespace(Package.set[i], quietly = TRUE)){
#       install.packages(Package.set[i])
#     }
#   }
#   ## Load Packages
#   lapply(Package.set, library, character.only = TRUE)
#   rm(Package.set,i)
#
#
#   #### BiocManager installation ####
#   ## Set the desired organism
#   # organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db")   ##  c("org.Dm.eg.db")
#
#   ## Set the desired organism
#   organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db","org.Dm.eg.db")
#   c(organism,"fgsea")
#
#   ## Check whether the installation of those packages is required from BiocManager
#   if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   Package.set <- c("clusterProfiler","enrichplot","pathview") # c(organism,"fgsea","clusterProfiler","enrichplot","pathview")
#   for (i in 1:length(Package.set)) {
#     if (!requireNamespace(Package.set[i], quietly = TRUE)){
#       BiocManager::install(Package.set[i])
#     }
#   }
#   ## Load Packages
#   lapply(Package.set, library, character.only = TRUE)
#   rm(Package.set,i)
#
#   options(stringsAsFactors = FALSE)


# #************************************************************************************************************************#
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
#     #### Transform CMGeneSet to list ####
#       # CMGeneSet.lt <- as.list(CMGeneSet[,c(-1,-2)])
#
#       # CMGeneSet.lt <- split(CMGeneSet[,c(-1,-2)], 1:nrow(CMGeneSet))
#       # names(CMGeneSet.lt) <- CMGeneSet[,1]
#
#       ## No use
#       # CMGeneSet.lt2 <- CMGeneSet.lt[!is.na(CMGeneSet.lt)]
#       # CMGeneSet.lt2 <- lapply(CMGeneSet.lt, function(x) unique(x))
#
#       ## EXist ""
#       # CMGeneSet.lt <- apply(CMGeneSet[,c(-1,-2)], 1, function(a)unique(a))
#       # names(CMGeneSet.lt) <- CMGeneSet[,1]
#
#       CMGeneSet.lt <- apply(CMGeneSet[,c(-1,-2)], 1, function(a)unique(a[a != ""]))
#       names(CMGeneSet.lt) <- CMGeneSet[,1]
#
#     #### Conduct analysis ####
#       fgseaRes <- fgsea(CMGeneSet.lt, ranks, minSize=15, maxSize = 500, nperm=1000)
#
#       ## fgsea: What does fgseaMultilevel argument sampleSize mean/when to change it?
#       ## https://www.biostars.org/p/479821/
#
#       fgseaRes <- fgseaMultilevel(CMGeneSet.lt, ranks, minSize=15, maxSize = 500)
#       ## Error when running parallelized process: Warning in serialize... package:stats may not be available when loading
#       ## https://community.rstudio.com/t/error-when-running-parallelized-process-warning-in-serialize-package-stats-may-not-be-available-when-loading/110573
#       ## https://stackoverflow.com/questions/27623901/r-warning-packagestats-may-not-be-available-when-loading
#
#       head(fgseaRes[order(padj, -abs(NES)), ], n=10)
#
#     #### Enrichment score plot ####
#       # plotEnrichment(CMGeneSet[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]], ranks)
#
#       plotEnrichment(CMGeneSet[CMGeneSet[,1]=="HALLMARK_ESTROGEN_RESPONSE_EARLY",], ranks)
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
#       # plotGseaTable(CMGeneSet[topPathways$pathway],
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
#       plotGseaTable(CMGeneSet.lt[topPathways$pathway],
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

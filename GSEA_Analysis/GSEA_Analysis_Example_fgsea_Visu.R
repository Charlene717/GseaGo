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

##### Parameter setting* #####
  # Set the desired organism
  organism = "org.Dm.eg.db"

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
  Package.set <- c(organism,"fgsea","clusterProfiler","enrichplot","pathview","enrichplot")
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


  #### Conduct analysis ####
    fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
    ## fgsea: What does fgseaMultilevel argument sampleSize mean/when to change it?
    ## https://www.biostars.org/p/479821/

    fgseaRes <- fgseaMultilevel(pathwaysH, ranks, minSize=15, maxSize = 500)
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
    y_bar <- group_by(y, sign(NES)) %>%
      slice(1:n)

    ggplot(y_bar, aes(NES, fct_reorder(Description, NES), fill = qvalues), showCategory=(n*2)) +
      geom_barh(stat='identity') +
      scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
      theme_minimal() + ylab(NULL)

    ## 2.2 Dotplot
    dotplot(msC2_2, showCategory = 10, font.size = 8,
            x = "GeneRatio",   # option -> c("GeneRatio", "Count")
            color = "p.adjust")   # option -> c("pvalue", "p.adjust", "qvalue")

    ## 2.3 Gene-Concept Network
    n <- 3
    p1 <- cnetplot(msC2_2, showCategory = (n*2), colorEdge = TRUE, node_label = "category")
    cowplot::plot_grid(p1, ncol=1, labels=LETTERS[1], rel_widths=c(1))

    ## 2.4 Heatmap-like functional classification

    ## 2.5 Enrichment Map
    p2 <- emapplot(pairwise_termsim(y), showCategory = 10)
    cowplot::plot_grid(p2, ncol = 1, lables = LETTERS[1])


    ## 2.6 UpSet Plot
    library(ggupset) # https://github.com/const-ae/ggupset
    upsetplot(msC2_2)

    ## 2.7 ridgeline plot for expressiong distribution
    ridgeplot(msC2_2)

    ## 2.8 gseaplot
    y2 <- arrange(msC2_2, desc(NES))

    p1 <- gseaplot(y2, geneSetID = 1, title = y2$Description[1])   # max NES
    n <- nrow(y2)
    p2 <- gseaplot(y2, geneSetID = n, title = y2$Description[n])   # min NES
    cowplot::plot_grid(p1, p2, ncol = 1, labels = LETTERS[1:2])

      ## 2.8.1 gseaplot2
      p3 <- gseaplot2(y2, geneSetID = 1, title = y2$Description[1])
      p4 <- gseaplot2(y2, geneSetID = n, title = y2$Description[n])   # min NES
      cowplot::plot_grid(p3, p4, ncol = 1, labels = LETTERS[1:2])

      ## Use keyword (Overlay graphics)
      keyword <- "breast"
      ind <- grep(keyword, msC2_2$Description, ignore.case = TRUE)
      gseaplot2(msC2_2, geneSetID = ind, title = msC2_2$Description[ind])

      ## Overlay graphics by ID
      gseaplot2(y2, geneSetID = 1:10)

      ## 2.8.2 gsearank
      gsearank(y2, geneSetID = 1, title = y2$Description[1])


    ## 2.9 PubMed trend of enriched terms
    terms <- msC2_2$Description[1:3]
    pmcplot(terms, 2010:2017, proportion=FALSE)



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
  Package.set <- c(organism,"fgsea")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  options(stringsAsFactors = FALSE)


##### Load Files #####
  ## Load data
  setwd("../")
  load("Demo_data/Robjects/Annotated_Results_LvV.RData")
  ## Load pathways
  load("Demo_data/Robjects/mouse_H_v5.RData")
  pathwaysH <- Mm.H

  setwd("GSEA_Analysis")

##### GSEA analysis #####


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



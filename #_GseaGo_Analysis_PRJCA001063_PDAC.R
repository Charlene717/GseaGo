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
  ProjectName = "TOP2A"
  Sampletype = "PDAC"
  # ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

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
  Package.set <- c(organism,"fgsea","clusterProfiler","enrichplot","pathview","edgeR")
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
  detach("package:ggstance", unload = TRUE)
  devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
  devtools::install_github("lionel-/ggstance")
  library(clusterProfiler.dplyr)
  library(ggstance)


##### Load Files #####
  ## Load data
  load("PRJCA001063_PDAC_Raw.RData")

  ## Load pathways
  setwd("../")
  load("Demo_data/Robjects/mouse_H_v5.RData")
  pathwaysH <- Mm.H
  setwd("GSEA_Analysis")

##### Data preprocessing #####
# https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### Select cell type ####
  # row.names(matrix.df) <- matrix.df[,1]
  # matrix.df <- matrix.df[,-1]



  #### Differential Expression ####
  ## Ref: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

  library(edgeR)
  Anno_Ints.df <- Anno.df[Anno.df$ReCluster %in% c("AD","AC"),]  # c("CoreCD00","CDOri")
  matrix_Ints.df <- matrix.df[,colnames(matrix.df) %in% Anno_Ints.df$CELL]
  row.names(matrix_Ints.df) <- matrix.df[,1]


  Anno_Ints.df <- left_join(data.frame(CELL=colnames(matrix_Ints.df)),Anno_Ints.df)

  DGE_Ints.lt <- DGEList(counts=matrix_Ints.df, group=Anno_Ints.df$ReCluster, lib.size=rep(1000,ncol(matrix_Ints.df)))
  DE_Extract.lt <- exactTest(DGE_Ints.lt, dispersion=0.2)
  DE_Extract.df <- topTags(DE_Extract.lt,n = nrow(matrix_Ints.df)) %>%
                   as.data.frame() %>%
                   data.frame(Gene=row.names(.),.)


  # #### Differential Expression ####
  # library(edgeR)
  # # generate raw counts from NB, create list object
  # y <- matrix(rnbinom(80,size=1/0.2,mu=10),nrow=20,ncol=4)
  # d <- DGEList(counts=y, group=c(1,1,2,2), lib.size=rep(1000,4))
  #
  # de <- exactTest(d, dispersion=0.2)
  # topTags(de)
  #
  # # same p-values using low-level function directly
  # p.value <- exactTestDoubleTail(y[,1:2], y[,3:4], dispersion=0.2)
  # sort(p.value)[1:10]






##### GSEA analysis (fgsea) #####
  library(fgsea)

  #### Create ranks ####
    gseaDat <- filter(DE_Extract.df, !is.na(Gene))
    ranks <- gseaDat$logFC
    names(ranks) <- gseaDat$Gene
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
    m_c2 <- msigdbr(species = "Homo sapiens") %>% # category = "C2"
            dplyr::select(gs_name, gene_symbol)
    msC2_2 <- GSEA(geneList, TERM2GENE = m_c2,pvalueCutoff =1)

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
    dotplot(msC2_2, showCategory = 20, font.size = 8,
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
      gseaplot2(y2, geneSetID = 10, title = y2$ID[10],color="red",pvalue_table=T)

      ## Use keyword (Overlay graphics)
      keyword <- "breast"
      ind <- grep(keyword, msC2_2$Description, ignore.case = TRUE)
      gseaplot2(msC2_2, geneSetID = ind, title = msC2_2$Description[ind])

      ## Overlay graphics by ID
      gseaplot2(y2, geneSetID = 1:10)
      gseaplot2(y2, geneSetID = 1:10, pvalue_table=T)

      Int_Path.set <- c("GRUETZMANN_PANCREATIC_CANCER_UP",
                        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                        "HP_ABNORMALITY_OF_CHROMOSOME_STABILITY",
                        "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION",
                        "NAKAMURA_METASTASIS_MODEL_UP",
                        "REACTOME_G0_AND_EARLY_G1",
                        "REACTOME_INITIATION_OF_NUCLEAR_ENVELOPE_NE_REFORMATION",
                        # "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
                        # "REACTOME_SEMA4D_MEDIATED_INHIBITION_OF_CELL_ATTACHMENT_AND_MIGRATION",
                        "REACTOME_NUCLEAR_PORE_COMPLEX_NPC_DISASSEMBLY"
                        )
      gseaplot2(y2, geneSetID = Int_Path.set)

      ## Modify the function
      ## https://www.biostars.org/p/9470087/
      ## trace("gseaplot2", edit = TRUE)
      ##> pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")] -> pd <- x[geneSetID, c("Description","NES, "pvalue", "p.adjust")]

      P2 <- gseaplot2(y2, geneSetID = Int_Path.set, pvalue_table=T) #+ theme_classic()  # White background
      P2

      pdf(
        file = paste0(getwd(), "/",Version,"/", Sys.Date(), "_Gseaplot_Int.pdf"),
        width = 17, height = 12
      )

      P2
      graphics.off()
      ## 2.8.2 gsearank
      gsearank(y2, geneSetID = 1, title = y2$Description[1])


    ## 2.9 PubMed trend of enriched terms
    terms <- msC2_2$Description[1:3]
    pmcplot(terms, 2010:2017, proportion=FALSE)



## Run GSEA analysis in R


## RNA-seq analysis in R
## Ref: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

## Gene Set Enrichment Analysis with ClusterProfiler
## Ref: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

## Visualzation of GSEA results
## Ref: https://rpubs.com/shbrief/gsea_263

## GSEA Chard Liu
## Ref: http://rstudio-pubs-static.s3.amazonaws.com/514990_9690f31b5ef7488bb4f0bb6c10ac4da8.html

FUN_ORA_Analysis = function(DE_Extract.df, pathwayGeneSet = Pathway.all,
                             TarGeneName = TarGene_name, GroupMode = Mode_Group,
                             ThrSet = Thr.lt,Species = "Homo sapiens",
                             Save.Path = Save.Path, SampleName = SampleName
){

  ##### Parameter setting* #####
    ## Set the desired organism
    organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db")   ##  c("org.Dm.eg.db")

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



  ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    library(clusterProfiler)
    library(enrichplot)
    # we use ggplot2 to add x axis labels (ex: ridgeplot)
    library(ggplot2)


    ##### Prepare Input #####
      # reading in data from deseq2
      df = read.csv("drosphila_example_de.txt", header=TRUE)

      # we want the log2 fold change
      original_gene_list <- df$log2FoldChange

      # name the vector
      names(original_gene_list) <- df$X

      # omit any NA values
      gene_list<-na.omit(original_gene_list)

      # sort the list in decreasing order (required for clusterProfiler)
      gene_list = sort(gene_list, decreasing = TRUE)

    ##### Gene Set Enrichment #####
      gse <- gseGO(geneList=gene_list,
                   ont ="ALL",
                   keyType = "ENSEMBL",
                   nPerm = 10000,
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = 0.05,
                   verbose = TRUE,
                   OrgDb = organism,
                   #OrgDb = org.Dm.eg.db, # https://github.com/YuLab-SMU/clusterProfiler/issues/279
                   pAdjustMethod = "none")
        # Error # Error in (function (classes, fdef, mtable)  :
        #             unable to find an inherited method for function ‘species’ for signature ‘"character"’
        #           In addition: There were 14 warnings (use warnings() to see them)
        ## https://stackoverflow.com/questions/62147679/unable-to-find-an-inherited-method-for-function-species-for-signature-charac
        # https://github.com/YuLab-SMU/clusterProfiler/issues/279

    ##### Dotplot #####
      require(DOSE)
      dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


    ##### Encrichment Map #####
      emapplot(gse, showCategory = 10)
      # Error # Term similarity matrix not available. Please use pairwise_termsim function to deal with the results of enrichment analysis.
      # Solved: https://github.com/YuLab-SMU/enrichplot/issues/79
      d <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")
      compare_cluster_GO_emap <- enrichplot::pairwise_termsim(gse, semData = d,  method="Wang")
      emapplot(compare_cluster_GO_emap)
      # Error # Error in loadNamespace(x) : 不存在叫 ‘ggnewscale’ 這個名稱的套件
      ## install.packages("ggnewscale")


    ##### Category Netplot #####
      # categorySize can be either 'pvalue' or 'geneNum'
      cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 5)

    ##### Ridgeplot #####
      ridgeplot(gse) + labs(x = "enrichment distribution")

    ##### GSEA Plot #####
      # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
      gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

    ##### PubMed trend of enriched terms #####
      terms <- gse$Description[1:3]
      pmcplot(terms, 2010:2018, proportion=FALSE)
      # Error # Error in loadNamespace(x) : 不存在叫 ‘europepmc’ 這個名稱的套件
      ## install.packages("europepmc")

  # ##### KEGG Gene Set Enrichment Analysis #####
    ## Prepare Input
      # Convert gene IDs for gseKEGG function
      # We will lose some genes here because not all IDs will be converted
      ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
      # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

      # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
      df2 = df[df$X %in% dedup_ids$ENSEMBL,]

      # Create a new column in df2 with the corresponding ENTREZ IDs
      df2$Y = dedup_ids$ENTREZID

      # Create a vector of the gene unuiverse
      kegg_gene_list <- df2$log2FoldChange

      # Name vector with ENTREZ ids
      names(kegg_gene_list) <- df2$Y

      # omit any NA values
      kegg_gene_list<-na.omit(kegg_gene_list)

      # sort the list in decreasing order (required for clusterProfiler)
      kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

    ## Create gseKEGG object
      kegg_organism = "dme"
      kk2 <- gseKEGG(geneList     = kegg_gene_list,
                     organism     = kegg_organism,
                     nPerm        = 10000,
                     minGSSize    = 3,
                     maxGSSize    = 800,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "none",
                     keyType       = "ncbi-geneid")
    ## Dotplot
      dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

    ## Encrichment map:
      emapplot(kk2)
      # Error # # Error in has_pairsim(x) :
      # Term similarity matrix not available. Please use pairwise_termsim function to deal with the results of enrichment analysis.
      # Solved: https://github.com/YuLab-SMU/enrichplot/issues/79
      d <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")
      compare_cluster_GO_emap_kk2 <- enrichplot::pairwise_termsim(kk2, semData = d)
      emapplot(compare_cluster_GO_emap_kk2)

    ## Category Netplot:
      ## categorySize can be either 'pvalue' or 'geneNum'
      cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
      ## Ridgeplot
      ridgeplot(kk2) + labs(x = "enrichment distribution")

    ## GSEA Plot
      # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
      gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

      ## Ref: https://rpubs.com/shbrief/gsea_263
      gseaplot2(kk2, geneSetID = 1:10)
      p3 <- gseaplot2(kk2, geneSetID = 1, title = kk2$Description[1])
      p4 <- gseaplot2(kk2, geneSetID = 2, title = kk2$Description[2])   # min NES
      cowplot::plot_grid(p3, p4, ncol = 1, labels = LETTERS[1:2])

    ## Pathview
      library(pathview)

      # Produce the native KEGG plot (PNG)
      dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme04130", species = kegg_organism)

      # Produce a different plot (PDF) (not displayed here)
      dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme04130", species = kegg_organism, kegg.native = F)
      knitr::include_graphics("dme04130.pathview.png")

return(Output)

}

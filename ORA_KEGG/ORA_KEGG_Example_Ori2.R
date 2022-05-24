## Over-Representation Analysis with ClusterProfiler
#  https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####     
  library("clusterProfiler")
  library(pathview)

##### Import files #####    
  # reading in input from deseq2
  df = read.csv("drosphila_example_de.txt", header=TRUE)

##### Prepare Input #####
  # we want the log2 fold change 
  original_gene_list <- df$log2FoldChange
  
  # name the vector
  names(original_gene_list) <- df$X

  # omit any NA values
  gene_list <- na.omit(original_gene_list)

  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)

  # Convert gene IDs for enrichKEGG function
  # We will lose some genes here because not all IDs will be converted
  ids <- bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Dm.eg.db") # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
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
  
  # Exctract significant results from df2
  kegg_sig_genes_df = subset(df2, padj < 0.05)
  
  # From significant results, we want to filter on log2fold change
  kegg_genes <- kegg_sig_genes_df$log2FoldChange
  
  # Name the vector with the CONVERTED ID!
  names(kegg_genes) <- kegg_sig_genes_df$Y
  
  # omit NA values
  kegg_genes <- na.omit(kegg_genes)
  
  # filter on log2fold change (PARAMETER)
  kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 2]

##### Create enrichKEGG object #####
  library("clusterProfiler")
  kegg_organism = "dme"
  kk <- enrichKEGG(gene=kegg_genes, universe=names(kegg_gene_list),organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")

##### Plot Chart #####
  ## Bar Chart
  barplot(kk, 
          showCategory = 10, 
          title = "Enriched Pathways",
          font.size = 8)
  ## Dot Chart
  dotplot(kk, 
          showCategory = 10, 
          title = "Enriched Pathways",
          font.size = 8)
  ## Category Netplot:
  # categorySize can be either 'pvalue' or 'geneNum'
  cnetplot(kk, categorySize="pvalue", foldChange=gene_list)


  library(pathview)
  
  ## KEGG plot
  # Produce the native KEGG plot (PNG)
  dme <- pathview(gene.data=gene_list, pathway.id="dme04080", species = kegg_organism, gene.idtype=gene.idtype.list[3])
  
  # Produce a different plot (PDF) (not displayed here)
  dme <- pathview(gene.data=gene_list, pathway.id="dme04080", species = kegg_organism, gene.idtype=gene.idtype.list[3], kegg.native = F)
  
  knitr::include_graphics("dme04080.pathview.png")

##### Presetting ######
  ## Clear variables
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Version information ######
  # platform       x86_64-w64-mingw32
  # arch           x86_64
  # os             mingw32
  # system         x86_64, mingw32
  # status
  # major          4
  # minor          1.2
  # year           2021
  # month          11
  # day            01
  # svn rev        81115
  # language       R
  # version.string R version 4.1.2 (2021-11-01)
  # nickname       Bird Hippie

##### Load library #####
  library(data.table)
  library(plyr)
  library(dplyr) # Basic data manupilation tools

##### Current path and new folder setting* #####
  ## File setting*
  FileName <- "Xena_TCGA_LGG_GE"
  Target_gene_name <- "TP53"
  Result_Folder_Name <- paste0(Target_gene_name,"_",Sys.Date()) ## Generate output folder automatically
  dir.create(Result_Folder_Name)

  ## Import genetic data file
  GeneExp.df <- read.table(FileName, header=T, row.names = 1, sep="\t")

##### Extract Target gene and Statistics ####
  # Extract data with Target_gene_name
  Target_gene_Mean <- GeneExp.df[Target_gene_name,] %>%
    as.numeric() %>%
    mean()

  #rowMeans(data.matrix(Target_gene))
  Target_gene_SD <- GeneExp.df[Target_gene_name,] %>%
    as.numeric() %>%
    sd()

##### Group the expression matrix according to the expression level of Target gene ####
  GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD]
  GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] <= Target_gene_Mean-Target_gene_SD]
  rm(Target_gene_Mean, Target_gene_SD)

##### Build Expression matrix for GSEA #####
  GeneExp_GSEA.df <- cbind(
    NAME=row.names(GeneExp.df),
    Description = rep("na", nrow(GeneExp.df)),
    GeneExp.df[,c(GeneExp_high.set, GeneExp_low.set)]
  )

  GSEA_SampleCol.df <- data.frame(t(colnames(GeneExp_GSEA.df)), stringsAsFactors=FALSE)
  colnames(GSEA_SampleCol.df) <- GSEA_SampleCol.df

  GeneExp_GSEA.df <- rbind(GSEA_SampleCol.df,GeneExp_GSEA.df)

  GeneExp_GSEA.df <- data.frame(
    "NAME" = c("#1.2",nrow(GeneExp.df)),
    "Description" = c('',length(c(GeneExp_high.set, GeneExp_low.set)))
  ) %>%
    rbind.fill(GeneExp_GSEA.df)
  rm(GSEA_SampleCol.df)

##### Build Group Files #####
  ## Set the group array
  Pheno_sum.df <- c(ncol(GeneExp_GSEA.df)-2,2,1) %>% t() %>% data.frame() %>%
    rbind.fill(c(paste0("#",Target_gene_name,"_high"),paste0(Target_gene_name,"_Low")) %>% t() %>% data.frame(stringsAsFactors=FALSE)) %>%
    rbind.fill(c(rep(0,length(GeneExp_high.set)),rep(1,length(GeneExp_low.set))) %>% t() %>% data.frame())


##### Export Result #####

  write.table(
    GeneExp_GSEA.df,
    file=paste0(Result_Folder_Name,"/",FileName,"_",Target_gene_name,"_collapsed.gct"),
    quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t'
  )
  write.table(
    Pheno_sum.df,
    file=paste0(Result_Folder_Name,"/",FileName,"_",Target_gene_name,".cls"),
    quote = FALSE,row.names = FALSE, na = "",col.names = FALSE
  )

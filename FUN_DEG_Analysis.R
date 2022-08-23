#### Differential Expression Gene Analysis ####
## Ref: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
## Ref: https://www.jianshu.com/p/b6912d318de5

FUN_DEG_Analysis = function(GeneExp.df,
                            TarGeneName = TarGene_name, GroupMode = Mode_Group,
                            Save.Path = Save.Path, SampleName = SampleName
){
  ##### Parameter setting* #####
  # Set the desired organism
  organism = "org.Dm.eg.db"


  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  # BiocManager::install()
  Package.set <- c(organism,"edgeR","baySeq")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### Differential Expression Gene Analysis ####

  #### Example ####
  library(baySeq)
  load("mobData.RData")

  con <- url("https://bios221.stanford.edu/data/mobData.RData") # Create connexion
  load(con) #Load the data

  head(mobData)

  #help(mobData)
  mobDataGroups <- c("MM", "MM", "WM", "WM", "WW", "WW")
  # MM="triple mutatnt shoot grafted onto triple mutant root"
  # WM="wild-type shoot grafted onto triple mutant root"
  # WW="wild-type shoot grafted onto wild-type root"
  data(mobAnnotation)
  #?mobAnnotation
  head(mobAnnotation)

  d <- DGEList(counts=mobData,group=factor(mobDataGroups))
  d

  ## Filtering the data
  dim(d)

  d.full <- d # keep the old one in case we mess up
  head(d$counts)
  head(cpm(d))

  apply(d$counts, 2, sum) # total gene counts per sample

  keep <- rowSums(cpm(d)>100) >= 2
  d <- d[keep,]
  dim(d)

  d$samples$lib.size <- colSums(d$counts)
  d$samples

  ## Normalizing the data
  d <- calcNormFactors(d)
  d

  ## Data Exploration
  plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
  legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

  ## Estimating the Dispersion
  d1 <- estimateCommonDisp(d, verbose=T)
  names(d1)

  d1 <- estimateTagwiseDisp(d1)
  names(d1)

  plotBCV(d1)

  ## GLM estimates of dispersion
  design.mat <- model.matrix(~ 0 + d$samples$group)
  colnames(design.mat) <- levels(d$samples$group)
  d2 <- estimateGLMCommonDisp(d,design.mat)
  d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
  # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
  # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
  d2 <- estimateGLMTagwiseDisp(d2,design.mat)
  plotBCV(d2)

  ## Comparing the models in DESeq and edgeR

  require(DESeq2)
  cds <- newCountDataSet( data.frame(d$counts), d$samples$group )
  cds <- estimateSizeFactors( cds )
  sizeFactors( cds )

  cds <- estimateDispersions( cds , method="blind")
  plotDispEsts(cds)




  ###############################################
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

  #### Export file ####


  #### Output ####


  return(Output)

}

#### Differential Expression Gene Analysis ####
## Ref: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

FUN_DEG_Analysis = function(GeneExp.df,
                            TarGeneName = TarGene_name, GroupMode = Mode_Group,
                            Save.Path = Save.Path, SampleName = SampleName
){

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


  return(Output)

}

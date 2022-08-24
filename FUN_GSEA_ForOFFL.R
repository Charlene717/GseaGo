## Build files for GSEA official input

FUN_GSEA_ForOFFL = function(GeneExp.df, Group1 = GeneExp_high.set, Group2 = GeneExp_low.set,
                            TarGeneName = TarGene_name, GroupMode = Mode_Group,
                            Save.Path = Save.Path, SampleName = SampleName
){

##### Build Expression matrix for GSEA #####
  GeneExp_GSEA.df <- cbind(
    NAME=row.names(GeneExp.df),
    Description = rep("na", nrow(GeneExp.df)),
    GeneExp.df[,c(Group1, Group2)]
  )

  GSEA_SampleCol.df <- data.frame(t(colnames(GeneExp_GSEA.df)), stringsAsFactors=FALSE)
  colnames(GSEA_SampleCol.df) <- GSEA_SampleCol.df

  GeneExp_GSEA.df <- rbind(GSEA_SampleCol.df,GeneExp_GSEA.df)

  GeneExp_GSEA.df <- data.frame(
    "NAME" = c("#1.2",nrow(GeneExp.df)),
    "Description" = c('',length(c(Group1, Group2)))
  ) %>%
    rbind.fill(GeneExp_GSEA.df)
  rm(GSEA_SampleCol.df)

##### Build Group Files #####
  ## Set the group array
  Pheno_sum.df <- c(ncol(GeneExp_GSEA.df)-2,2,1) %>% t() %>% data.frame() %>%
    rbind.fill(c(paste0("#",TarGeneName,"_high"),paste0(TarGeneName,"_Low")) %>% t() %>% data.frame(stringsAsFactors=FALSE)) %>%
    rbind.fill(c(rep(0,length(Group1)),rep(1,length(Group2))) %>% t() %>% data.frame())


##### Export Result #####
  if(GroupMode$Mode == "Mean"){
    write.table(
      GeneExp_GSEA.df,
      file=paste0(Save.Path,"/OFFL_",SampleName,"_",
                  GroupMode$Mode,GroupMode$SD,"SD_",
                  TarGeneName,"_collapsed.gct"),
      quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t'
    )
    write.table(
      Pheno_sum.df,
      file=paste0(Save.Path,"/OFFL_",SampleName,"_",
                  GroupMode$Mode,GroupMode$SD,"SD_",
                  TarGeneName,".cls"),
      quote = FALSE,row.names = FALSE, na = "",col.names = FALSE
    )
  }else{
    write.table(
      GeneExp_GSEA.df,
      file=paste0(Save.Path,"/OFFL_",SampleName,"_",
                  GroupMode$Mode,"Q2",GroupMode$Q2,"_",
                  TarGeneName,"_collapsed.gct"),
      quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t'
    )
    write.table(
      Pheno_sum.df,
      file=paste0(Save.Path,"/OFFL_",SampleName,"_",
                  GroupMode$Mode,"Q2",GroupMode$Q2,"_",
                  TarGeneName,".cls"),
      quote = FALSE,row.names = FALSE, na = "",col.names = FALSE
    )

  }


}

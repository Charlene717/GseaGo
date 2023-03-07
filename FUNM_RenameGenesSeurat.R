# Based on https://github.com/satijalab/seurat/issues/1049#issuecomment-453602745
# and on https://github.com/satijalab/seurat/issues/978#issuecomment-444526949
# using https://waldronlab.io/HGNChelper/articles/index.html

RenameGenesSeurat <- function(SeuObj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]) {
  print("Run this before integration. It only changes SeuObj@assays$RNA@counts, @data and @scale.data")
  RNA <- SeuObj@assays$RNA

  if (nrow(RNA) == nrow(newnames)) {
    if(l(RNA@counts)) RNA@counts@Dimnames[[1]] <-         newnames$Suggested.Symbol
    if(l(RNA@data)) RNA@data@Dimnames[[1]] <-             newnames$Suggested.Symbol
    if(l(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames$Suggested.Symbol
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  SeuObj@assays$RNA <- RNA
  return(SeuObj)
}

plot.UpdateStats <- function(genes = HGNC.updated[[i]]) {
  (MarkedAsUpdated <- genes[genes$Approved == FALSE, ])
  (AcutallyUpdated <- sum(MarkedAsUpdated[,1] != MarkedAsUpdated[,3]))
  (UpdateStats = c((AcutallyUpdated / nrow(genes)), AcutallyUpdated, nrow(genes)))
  return(UpdateStats)
}


# ------------------------------------------------------------------------------------------------------------------------------------------------

if (TRUE) {
  library(HGNChelper); getCurrentHumanMap(); # install.packages("HGNChelper")
  # new.hgnc.table <- getCurrentHumanMap()
  # write.simple.tsv(new.hgnc.table)
  HGNC.updated <- list.fromNames(samples)
  UpdateStatMat = matrix.fromNames(rowname_vec =samples, colname_vec = c( "percent","changed", "total"))

  tic(); for (i in 1:n.datasets ) {
    HGNC.updated[[i]]      <- checkGeneSymbols(rownames(ls.Seurat[[i]]), unmapped.as.na = FALSE, map = NULL, species = "human")
    ls.Seurat[[i]] <- RenameGenesSeurat(SeuObj = ls.Seurat[[i]], newnames = HGNC.updated[[i]])

    UpdateStatMat[i,] <- plot.UpdateStats(HGNC.updated[[i]])
  }; toc();

  HGNC.PercentUpdated <- cbind(100*UpdateStatMat[,1], UpdateStatMat[,2])
  wplot(HGNC.PercentUpdated)

  # Check
  if (FALSE) {
    goi = "FAM132A"
    goi = "ARMH1" # new C1orf228
    goi = "C1orf228" # old

    rownames(ls.Seurat[[i]]@assays$RNA)
    rownames(ls.Seurat_bac[[i]]@assays$RNA)
    grep(x = rownames(ls.Seurat[[i]]@assays$RNA), pattern = goi)
    grep(x = HGNC.updated[[i]]$Suggested.Symbol, pattern = goi)
    grep(x = HGNC.updated[[i]]$x, pattern = goi)
  }
}



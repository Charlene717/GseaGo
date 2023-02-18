####-------------------- Load Packages --------------------###

#### For few packages ####
# if(!require("tidyverse")) install.packages("tidyverse")
# library(tidyverse)

#### Set the desired organism ####
# organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db","org.Dm.eg.db")
# c(organism,"fgsea")

#### For Many packages  ####
FUN_Package_InstLoad = function( Basic.set = Basic.set,
                                 BiocManager.set = BiocManager.set)
{

  #### Basic installation ####
  ## Check whether the installation of those packages is required from basic
  if(exists("Basic.set")==TRUE){

    for (i in 1:length(Basic.set)) {
      if (!requireNamespace(Basic.set[i], quietly = TRUE)){
        install.packages(Basic.set[i])
      }
    }
    ## Load Packages
    lapply(Basic.set, library, character.only = TRUE)
    rm(Basic.set,i)

  }


  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if(exists("BiocManager.set")==TRUE){

    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    for (i in 1:length(BiocManager.set)) {
      if (!requireNamespace(BiocManager.set[i], quietly = TRUE)){
        BiocManager::install(BiocManager.set[i])
      }
    }
    ## Load Packages
    lapply(BiocManager.set, library, character.only = TRUE)
    rm(BiocManager.set,i)

  }

}

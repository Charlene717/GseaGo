Package_InstLoad = function( Basic.set = Basic.set,
                             BiocManager.set = BiocManager.set)
{

  ##### Load Packages #####
  # if(!require("tidyverse")) install.packages("tidyverse")
  # library(tidyverse)

  #### Basic installation ####
  ## Check whether the installation of those packages is required from basic
  for (i in 1:length(Basic.set)) {
    if (!requireNamespace(Basic.set[i], quietly = TRUE)){
      install.packages(Basic.set[i])
    }
  }
  ## Load Packages
  lapply(Basic.set, library, character.only = TRUE)
  rm(Basic.set,i)


  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager.set <- BiocManager.set
  for (i in 1:length(BiocManager.set)) {
    if (!requireNamespace(BiocManager.set[i], quietly = TRUE)){
      BiocManager::install(BiocManager.set[i])
    }
  }
  ## Load Packages
  lapply(BiocManager.set, library, character.only = TRUE)
  rm(BiocManager.set,i)

  # options(stringsAsFactors = FALSE)
  # Sys.setlocale(category = "LC_ALL", locale = "UTF-8")

}

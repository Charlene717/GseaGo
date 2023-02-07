## MSigDB Collections from GSEA
## http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####
  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)
  if(!require("gtools")) install.packages("gtools")
  library(gtools)
  # if(!require("XML")) install.packages("XML")
  # library(XML)

##### Current path and new folder setting* #####
  OutputFileName <- "Default"

##### Import files & Combine df #####
  #### Import XML file of GSEA GeneSet ####
  # # Ref: https://www.educative.io/answers/how-to-read-xml-files-in-r
  # GSEAGeneSet.df <- xmlToDataFrame(paste0("Input_Genesets/Gsea_Genesets_Hs/msigdb_v2022.1.Hs.xml"))
  # GSEAGeneSet.df <- read.delim2(paste0("Input_Genesets/Gsea_Genesets_Hs/msigdb_v2022.1.Hs.xml"),sep = "\t")
  #
  GSEAGeneSet_Hs_XML.df <- read.delim2(paste0("Input_Genesets/Gsea_Genesets_Hs/msigdb_v2022.1.Hs.txt"),sep = "\t")
  GSEAGeneSet_Mm_XML.df <- read.delim2(paste0("Input_Genesets/Gsea_Genesets_Mm/msigdb_v2022.1.Mm.txt"),sep = "\t")

  ## Clean up the data




  #### Import gmt file of GSEA GeneSet ####
  # Function for read and combine multiple files
  FUN_ImportGmt <- function(FilesList.set) {

    Nfiles = length(FilesList.set)

    for(i in 1:Nfiles){
      if(i==1){
        # Deal with different number of columns
        merge.df <- read.delim2(FilesList.set[1],
                                col.names = 1:max(count.fields(FilesList.set[1])),
                                header = F,sep = "\t")
      }else{
        new_1 <- read.delim2(paste0(FilesList.set[i]),
                             col.names = 1:max(count.fields(FilesList.set[i])),
                             header = F,sep = "\t")
        merge.df <- smartbind(merge.df,new_1)
      }

    }

    #### Clean up df ####
    ## Remove duplicated
    merge.df <- merge.df[!duplicated(merge.df[,2]), ]
    # ## Remove NA (Have set in the write.table) # Ref: https://www.delftstack.com/zh-tw/howto/r/replace-na-with-0-in-r/
    #   merge.df[is.na(merge.df)] <- ""

    return(merge.df)

  }

  ## read gmt file of Hs.symbols
  # InputFolder = "Gsea_Genesets_Hs"
  # # target.dir <- list.dirs( paste0("Input_Genesets/", InputFolder) )[-1]
  # FilesList.set <- list.files(paste0("Input_Genesets/", InputFolder),full.names = T)
  # FilesList.set <- str_subset(FilesList.set, pattern = "\\.symbols.gmt$")
  FilesList.set <- list.files(paste0("Input_Genesets/", "Gsea_Genesets_Hs"),full.names = T) %>%
    str_subset(., pattern = "\\.symbols.gmt$")

  GSEAGeneSet_Hs_symb_gm.df <- FUN_ImportGmt(FilesList.set)

  ## read gmt file of Hs.entrez
  FilesList.set <- list.files(paste0("Input_Genesets/", "Gsea_Genesets_Hs"),full.names = T) %>%
                str_subset(., pattern = "\\.entrez.gmt$")

  GSEAGeneSet_Hs_entrez_gmt.df <- FUN_ImportGmt(FilesList.set)

  ## read gmt file of Mm.symbols
  FilesList.set <- list.files(paste0("Input_Genesets/", "Gsea_Genesets_Mm"),full.names = T) %>%
    str_subset(., pattern = "\\.symbols.gmt$")

  GSEAGeneSet_Mm_symb_gmt.df <- FUN_ImportGmt(FilesList.set)

  ## read gmt file of Mm.entrez
  FilesList.set <- list.files(paste0("Input_Genesets/", "Gsea_Genesets_Mm"),full.names = T) %>%
    str_subset(., pattern = "\\.entrez.gmt$")

  GSEAGeneSet_Mm_entrez_gmt.df <- FUN_ImportGmt(FilesList.set)
  rm(FilesList.set)

##### Update gene name ####



#################################################################################################################################


#### Save RData ####
  save.image(paste0("Input_Genesets/Genesets_", OutputFileName,".RData"))

#################################################################################

  ## XML檔設定
  ## 加入其他篩選條件
  ## Geneset by 自己的實驗或線上數據(文字型,matrix型)




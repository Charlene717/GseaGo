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

##### Condition setting* #####




##### Current path and new folder setting* #####
  OutputFileName <- "Default"
  InputFolder = "Gsea_Genesets_Hs"
  OutputFolder <- paste0("Input_Genesets/", InputFolder, "_", OutputFileName)
  dir.create(OutputFolder) ## Generate output folder

##### Import files & Combine df #####
  #### Import XML file of GSEA GeneSet ####
  # # Ref: https://www.educative.io/answers/how-to-read-xml-files-in-r
  # GSEAGeneSet.df <- xmlToDataFrame(paste0("Input_Genesets/Gsea_Genesets_Hs/msigdb_v2022.1.Hs.xml"))
  # GSEAGeneSet.df <- read.delim2(paste0("Input_Genesets/Gsea_Genesets_Hs/msigdb_v2022.1.Hs.xml"),sep = "\t")
  #
  GSEAGeneSet_Hs_XML.df <- read.delim2(paste0("Input_Genesets/Gsea_Genesets_Hs/msigdb_v2022.1.Hs.txt"),sep = "\t")
  GSEAGeneSet_Mm_XML.df <- read.delim2(paste0("Input_Genesets/Gsea_Genesets_Mm/msigdb_v2022.1.Mm.txt"),sep = "\t")

  ## Import Customization
  # target.dir <- list.dirs( paste0("Input_Genesets/", InputFolder) )[-1]
  list.files <- list.files(paste0("Input_Genesets/", InputFolder),full.names = T)
  list.files <- str_subset(list.files, pattern = "\\.symbols.gmt$")

  Nfiles = length(list.files)

  for(i in 1:Nfiles){
    if(i==1){
      # Deal with different number of columns
      merge.df <- read.delim2(list.files[1],
                              col.names = 1:max(count.fields(list.files[1])),
                              header = F,sep = "\t")
    }else{
    new_1 <- read.delim2(paste0(list.files[i]),
                         col.names = 1:max(count.fields(list.files[i])),
                         header = F,sep = "\t")
    merge.df <- smartbind(merge.df,new_1)
    }

  }
  rm(new_1,i)

  #### Clean up df ####
  ## Remove duplicated
    merge.df <- merge.df[!duplicated(merge.df[,2]), ]

  # ## Remove NA (Have set in the write.table)
  # # Ref: https://www.delftstack.com/zh-tw/howto/r/replace-na-with-0-in-r/
  #   merge.df[is.na(merge.df)] <- ""

##### Update gene name ####



#################################################################################################################################


#### Save RData ####
  save.image(paste0("Input_Genesets/", InputFolder,".RData"))

#################################################################################

  ## XML檔設定
  ## 加入其他篩選條件
  ## Geneset by 自己的實驗或線上數據(文字型,matrix型)




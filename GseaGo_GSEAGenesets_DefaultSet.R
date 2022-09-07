## MSigDB Collections from GSEA
## Human
## http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
## Mouse
## https://www.gsea-msigdb.org/gsea/msigdb/mouse_geneset_resources.jsp

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####
  #library(dplyr)
  #library(plyr) #merge.df<- rbind.fill(merge.df,new_1)

  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)
  if(!require("gtools")) install.packages("gtools")
  library(gtools)


##### Current path and new folder setting* #####
  InputFolder = "Customized_GSEAGenesets_Pathway3D_Hm"
  OutputFolder <- paste0(InputFolder,"_",Sys.Date(),"_CF") ## Generate output folder automatically
  dir.create(OutputFolder)

  ExpFilName <- "ComB" #"ComB" # Combine


##### Import files & Combine df #####
  # target.dir <- list.dirs(InputFolder)[-1]
  list.files <- list.files(InputFolder,full.names = T)
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

  ##### Export Result of Combine #####
  ## Note ## Need to remove the quote
    write.table(merge.df,paste0(OutputFolder,"/",InputFolder,'_',ExpFilName ,'.txt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
    write.table(merge.df,paste0(OutputFolder,"/",InputFolder,'_',ExpFilName ,'.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")

#################################################################################


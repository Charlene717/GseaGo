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

  #library(dplyr)
  #library(plyr) #merge.df<- rbind.fill(merge.df,new_1)

##### Current path and new folder setting* #####
  OutputFileName <- "Convert"
  InputFolder = "Cust_GSEA_Genesets_Test"
  OutputFolder <- paste0("Input_Genesets/", InputFolder, "_", OutputFileName)
  dir.create(OutputFolder) ## Generate output folder

##### Import files & Combine df #####
  # target.dir <- list.dirs( paste0("Input_Genesets/", InputFolder) )[-1]
  list.files <- list.files(paste0("Input_Genesets/", InputFolder),full.names = T)
  list.files <- str_subset(list.files, pattern = "\\.gmt$")

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
    write.table(merge.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName ,'.txt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")
    write.table(merge.df,paste0(OutputFolder,"/",InputFolder,'_',OutputFileName ,'.gmt'),
                row.names = FALSE,col.names= FALSE,quote = FALSE, sep = '\t', na="")


#################################################################################

##### Convert GeneName #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  source("FUN_HSsymbol2MMsymbol.R")

  # # Inport Geneset
  # Pathway.all <- read.delim2(paste0(getwd(),"/GSEA_Geneset/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"),
  #                            col.names = 1:max(count.fields(paste0(getwd(),"/GSEA_Geneset/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"))),
  #                            header = F,sep = "\t")

  Pathway.all <- merge.df

  ##### Converting the Human gene name to Mouse gene name #####
  # #  Need to be optimized
  # # (Method1) bind the different length of column (Cannot use rbind)
  # # (Method2) Save the data as list first and than use do.call() to unlist to have dataframe
  #
  # ## (Ori method)
  # Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*2))
  # for (i in 1:nrow(Pathway.all)) {
  #   #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
  #   PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()
  #   colnames(PathwayN)="Temp"
  #   PathwayN <- HSsymbol2MMsymbol(PathwayN,"Temp")
  #   Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
  # }
  #
  # Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
  # colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
  #
  # rm(PathwayN)
  #
  # # assign(paste0("marrow_sub_DucT2_TOP2ACenter_T", i),marrow_sub_DucT2_TOP2ACenter_Tn)
  # # assign(colnames(Pathway.all)[i],Pathway.all[,i])

  ## (Method1)
  # Refer # https://stackoverflow.com/questions/3699405/how-to-cbind-or-rbind-different-lengths-vectors-without-repeating-the-elements-o
  # How to cbind or rbind different lengths vectors without repeating the elements of the shorter vectors?
  ## Modify by Charlene: Can use in MultRow
  bind_diff <- function(x, y){
    if(ncol(x) > ncol(y)){
      len_diff <- ncol(x) - ncol(y)
      y <- data.frame(y, rep(NA, len_diff) %>% t() %>% as.data.frame())
      colnames(x) <- seq(1:ncol(x))
      colnames(y) <- seq(1:ncol(y))
    }else if(ncol(x) < ncol(y)){
      len_diff <- ncol(y) - ncol(x)
      x <- data.frame(x, rep(NA, len_diff) %>% t() %>% as.data.frame())
      colnames(x) <- seq(1:ncol(x))
      colnames(y) <- seq(1:ncol(y))
    }
    rbind(x, y)
  }

  ## Converting
  for (i in 1:nrow(Pathway.all)) {
    PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()  %>% as.data.frame()
    colnames(PathwayN)="Temp"
    PathwayN <- HSsymbol2MMsymbol(PathwayN,"Temp")
    PathwayN <- PathwayN[!PathwayN$MM.symbol  == 0,]
    PathwayN <- PathwayN[!is.na(PathwayN$MM.symbol),]
    if(i==1){
      Pathway.all.MM <- unique(PathwayN$MM.symbol) %>% t()  %>% as.data.frame()
    }else{
      Pathway.all.MM <- bind_diff(Pathway.all.MM,unique(PathwayN$MM.symbol) %>% t()  %>% as.data.frame())
      # Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
    }
  }

  Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
  colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
  #Pathway.all.MM[Pathway.all.MM==0] <-NA

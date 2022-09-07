##### 09_0 GSEA Analysis (Geneset Prepare) #####
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
# install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
library(fgsea)
source("FUN_GSEA_LargeGeneSet.R")
source("FUN_HSsymbol2MMsymbol.R")
source("FUN_GSEA_ggplot.R")

# Geneset from GSEA
# Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
Pathway.all <- read.delim2(paste0(getwd(),"/GSEA_Geneset/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"),
                           col.names = 1:max(count.fields(paste0(getwd(),"/GSEA_Geneset/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"))),
                           header = F,sep = "\t")

# Convert Human gene to mouse
Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*1.5))
for (i in 1:nrow(Pathway.all)) {
  #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
  PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()
  colnames(PathwayN)="Test"
  PathwayN <- HSsymbol2MMsymbol(PathwayN,"Test")
  Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)

}

Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
Pathway.all.MM[is.na(Pathway.all.MM)] <- ""

#### Save RData ####
save.image(paste0(Save.Path,"/09_0_GSEA_Analysis(Geneset_Prepare).RData"))






#
# ##### Other #####
#
# ## Function setting
#
# ## Call function
# filePath <- ""
# # Import R files in the same folder
# getFilePath <- function(fileName) {
#   # Absolute path of project folder
#   # path <- setwd("~")
#   path <- setwd(getwd())
#   # Combine strings without gaps
#   # <<- Assigning values to global variable
#   filePath <<- paste0(path ,"/" , fileName)
#   # Load file
#   sourceObj <- source(filePath)
#   return(sourceObj)
# }
# getFilePath("HSsymbol2MMsymbol.R")
#
# ## Geneset from GSEA
# H.all <- read.delim(paste0(PathName,"/h.all.v7.4.symbols.gmt"),header = F)
#
# H.all.list <- list()
# for (i in c(1:length(H.all[,1]))) {
#   H.all.list.ori <- as.data.frame(t(H.all[i,3:length(H.all[i,])]))
#   colnames(H.all.list.ori)[[1]] <- c("Gene")
#   H.all.list.ori <- HSsymbol2MMsymbol(H.all.list.ori,"Gene")
#
#   # Delete NA(or 0)
#   H.all.list.ori <- H.all.list.ori[H.all.list.ori$MM.symbol!=0,]
#
#   H.all.list.ori <- unique(H.all.list.ori$MM.symbol)
#   H.all.list[[i]] <- as.character(H.all.list.ori)
#
#   rm(H.all.list.ori)
#   names(H.all.list)[[i]] <- H.all[i,1]
# }
#
#

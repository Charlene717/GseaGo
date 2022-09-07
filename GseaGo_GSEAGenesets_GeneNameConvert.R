

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

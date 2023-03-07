##### To do list #####
## -[] How to reduce computing time?

##### Update the genename ####
## Ref: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
library(limma)

## Set Species
if(Set_Species == "Homo sapiens"){
  Specie = "Hs"
}else if(Set_Species == "Mus musculus"){
  Specie = "Mm"
}else{
  Specie = "Hs"
}

## Update the gene names to the latest version
FUN_UpdateGene <- function(GeneName_Ori, Species = Specie, AvoidMult = TRUE) {
  UpdateGeneName <- alias2Symbol(GeneName_Ori, species = Species, expand.symbols = FALSE)
  if( length(UpdateGeneName) == 0 ){
    GeneName <- GeneName_Ori
  }else{
    GeneName <- UpdateGeneName
  }

  if(AvoidMult == TRUE){ ## Aovid one-to-many
    GeneName <- GeneName[1]
  }
  return(GeneName)
}

# ## Test FUN_UpdateGene function
# TestGene <- "SEPT1" # TestGene <- "HBII-52-46"
# TestGene <- FUN_UpdateGene(TestGene)

GeneNameUpdate.df <- lapply(row.names(GeneExp.df), function(x)FUN_UpdateGene(x))  %>% unlist() %>% as.data.frame()

UpdateGeneName_Compare.df <- cbind(row.names(GeneExp.df),GeneNameUpdate.df[,1]) %>% as.data.frame()


#### Find Duplicate name Aovid many-to-one ####
## Ref: http://guangzheng.name/2017/10/07/%E5%A6%82%E4%BD%95%E6%9F%A5%E6%89%BE%E6%95%B0%E6%8D%AE%E6%A1%86%E4%B8%AD%E9%87%8D%E5%A4%8D%E7%9A%84%E6%95%B0%E6%8D%AE/
library(dplyr)
## Extract duplicate name
UpdateGeneName_Compare.df %>% group_by(V2) %>%
                   dplyr::mutate(index = n()) %>%
                   filter(index > 1) %>%
                   dplyr::select(2) %>%
                   ungroup() %>%
                   unique() %>%
                   unlist() -> UpdateGeneName_Dup.Set

## -[] How to deal with duplicate name (many-to-one)?
# (1)Small number which can be ignored
# (2)Find the most varied gene as the representative gene (standard deviation, quartile, range)
# (3)The gene that shows the most difference between the two populations is used as the representative gene

## Deal with duplicate name (many-to-one)
FUN_DWMany2One <- function(df,x) {
  if( (df[x,2] %in% UpdateGeneName_Dup.Set)== TRUE ){
    df[x,1] = df[x,1] # (No change if encounter duplicate names)
  }else{
    df[x,1] = df[x,2]
  }
  return(df[x,1])
}

GeneNameDWM2O.df <- lapply(1:nrow(UpdateGeneName_Compare.df), function(x)FUN_DWMany2One(UpdateGeneName_Compare.df,x))  %>% as.data.frame() %>% t
row.names(GeneExp.df) <- GeneNameDWM2O.df


## Difference before and after statistics
UpdateGeneName_Compare.df <- cbind(UpdateGeneName_Compare.df, GeneNameDWM2O.df)
colnames(UpdateGeneName_Compare.df) <- c("Ori","UpGeneName","DUPEGene")
row.names(UpdateGeneName_Compare.df) <- seq(1:nrow(UpdateGeneName_Compare.df))

UpdateGeneName_Sumdf <- data.frame(
  OriToUpdate_Same = sum(UpdateGeneName_Compare.df[,1] == UpdateGeneName_Compare.df[,2]),
  OriToUpdate_Diff = sum(UpdateGeneName_Compare.df[,1] != UpdateGeneName_Compare.df[,2]),

  DWM2OToUpdate_Same = sum(UpdateGeneName_Compare.df[,2] == UpdateGeneName_Compare.df[,3]),
  DWM2OToUpdate_Diff = sum(UpdateGeneName_Compare.df[,2] != UpdateGeneName_Compare.df[,3])
)

rm(GeneNameUpdate.df, GeneNameDWM2O.df, Specie)

#************************************************************************************************************************#
# #### Old version 1 ####
#
# # ## Update the genename ##* Take very long time
# # UpdateGene <- "No"  # UpdateGene <- c("Yes","No")
# # if(UpdateGene == "Yes"){
# #   row.names(GeneExp.df) <- UpdateSymbolList(row.names(GeneExp.df))
# # }
#
#
# # #### Test ####
# # UpdateSymbolList("SEPT1")
# # A <- UpdateSymbolList(row.names(GeneExp.df))
# # B <- row.names(GeneExp.df)
# # # sum(c("a","c")==c("a","b"))
# # sum(A==B)
# # summary(A==B)
# #
# ## Error: Timeout was reached: [rest.genenames.org] Operation timed out after 10005 milliseconds with 0 bytes received
#
# # ## https://rdrr.io/github/vertesy/Seurat.utils/src/Development/Functions/Seurat.update.gene.symbols.HGNC.R
# # HGNC.EnforceUniquet("SEPT1")


#************************************************************************************************************************#
# #### Test2 ####
# ## Ref: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
# library(limma)
# alias2Symbol("SEPT1", species = "Hs", expand.symbols = FALSE)


#************************************************************************************************************************#
#### Backup ####
## https://www.nature.com/articles/s41588-020-0669-3#Sec18
## Genenames.org: the HGNC and VGNC resources in 2021
## https://academic.oup.com/nar/article/49/D1/D939/5957168

# HGNC
# https://www.genenames.org/

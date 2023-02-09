##### Update the genename ####
## Ref: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
library(limma)

if(Set_Species == "Homo sapiens"){
  Specie = "Hs"
}else if(Set_Species == "Mus musculus"){
  Specie = "Mm"
}else{
  Specie = "Hs"
}

UpdateGene <- function(GeneName_Ori, Species = Specie, AvoidMult = TRUE) {
  UpdateGeneName <- alias2Symbol(GeneName_Ori, species = Species, expand.symbols = FALSE)
  if( length(UpdateGeneName) == 0 ){
    GeneName <- GeneName_Ori
  }else{
    GeneName <- UpdateGeneName
  }

  if(AvoidMult == TRUE){
    GeneName <- GeneName[1]
  }
  return(GeneName)
}

# ## Test UpdateGene function
# TestGene <- "SEPT1" # TestGene <- "HBII-52-46"
# TestGene <- UpdateGene(TestGene)

GeneExp_Temp.df <- GeneExp.df
UpGeneName.df <- lapply(row.names(GeneExp_Temp.df), function(x)UpdateGene(x))  %>% unlist() %>% as.data.frame()

CompareGene.df <- cbind(row.names(GeneExp_Temp.df),UpGeneName.df[,1]) %>% as.data.frame()


## Find Duplicate name
## Ref: http://guangzheng.name/2017/10/07/%E5%A6%82%E4%BD%95%E6%9F%A5%E6%89%BE%E6%95%B0%E6%8D%AE%E6%A1%86%E4%B8%AD%E9%87%8D%E5%A4%8D%E7%9A%84%E6%95%B0%E6%8D%AE/
library(dplyr)
CompareGene.df %>% group_by(V2) %>%
                   mutate(index = n()) %>%
                   filter(index > 1) %>%
                   select(2) %>%
                   ungroup() %>%
                   unique() %>%
                   unlist() -> Dup.set

## Deal with Duplicate name (No change if encounter duplicate names)
UpdateGeneDUPE <- function(df,x) {
  if( (df[x,2] %in% Dup.set)== TRUE ){
    df[x,1] = df[x,1]
  }else{
    df[x,1] = df[x,2]
  }
  return(df[x,1])
}

UpGeneName2.df <- lapply(1:nrow(CompareGene.df), function(x)UpdateGeneDUPE(CompareGene.df,x))  %>% as.data.frame() %>% t
row.names(GeneExp_Temp.df) <- UpGeneName2.df
GeneExp.df <- GeneExp_Temp.df
CompareGene.df <- cbind(CompareGene.df, UpGeneName2.df)
colnames(CompareGene.df) <- c("Ori","UpGeneName","DUPEGene")
row.names(CompareGene.df) <- seq(1:nrow(CompareGene.df))

CompareGene_Sum.df <- data.frame(
  A = sum(CompareGene.df[,1] == CompareGene.df[,2]),
  B = sum(CompareGene.df[,1] != CompareGene.df[,2]),
  C =sum(CompareGene.df[,2] != CompareGene.df[,3])
)


rm(GeneExp_Temp.df,UpGeneName.df,UpGeneName2.df)

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

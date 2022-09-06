##### Update the genename ####
# #### Test ####
# UpdateSymbolList("SEPT1")
# A <- UpdateSymbolList(row.names(GeneExp.df))
# B <- row.names(GeneExp.df)
# # sum(c("a","c")==c("a","b"))
# sum(A==B)
# summary(A==B)
#
## Error: Timeout was reached: [rest.genenames.org] Operation timed out after 10005 milliseconds with 0 bytes received

# ## https://rdrr.io/github/vertesy/Seurat.utils/src/Development/Functions/Seurat.update.gene.symbols.HGNC.R
# HGNC.EnforceUniquet("SEPT1")

# ## Update the genename ##* Take very long time
# UpdateGene <- "No"  # UpdateGene <- c("Yes","No")
# if(UpdateGene == "Yes"){
#   row.names(GeneExp.df) <- UpdateSymbolList(row.names(GeneExp.df))
# }


# #### Test2 ####
# ## Ref: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
# library(limma)
# alias2Symbol("SEPT1", species = "Hs", expand.symbols = FALSE)


## Ref: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
library(limma)

UpdateGene <- function(TestGeneName, Species = "Hs") {
  UpdateGeneName <- alias2Symbol(TestGeneName, species = Species, expand.symbols = FALSE)
  if( length(UpdateGeneName) == 0 ){
    TestGeneName <- TestGeneName
  }else{
    TestGeneName <- UpdateGeneName
  }

  if(length(TestGeneName) > 1){
    TestGeneName <- TestGeneName[1]
  }
  return(TestGeneName)
}

# ## Test UpdateGene function
# TestGene <- "SEPT1"
# TestGene <- "HBII-52-46"
# TestGene <- UpdateGene(TestGene)

TTT <- GeneExp.df
# EX1 <- lapply(1:3, function(x)x^x)
# TTT2 <- lapply(row.names(TTT), function(x)UpdateGene(x))
# # TTT2 <- lapply(1:3, function(x)UpdateGene(x))
# TTT3 <- unlist(TTT2) %>% as.data.frame()
TTT3 <- lapply(row.names(TTT), function(x)UpdateGene(x))  %>% unlist() %>% as.data.frame()

df <- cbind(row.names(TTT),TTT3[,1]) %>% as.data.frame()
sum(df[,1] != df[,2])

## Find Duplicate name
## Ref: http://guangzheng.name/2017/10/07/%E5%A6%82%E4%BD%95%E6%9F%A5%E6%89%BE%E6%95%B0%E6%8D%AE%E6%A1%86%E4%B8%AD%E9%87%8D%E5%A4%8D%E7%9A%84%E6%95%B0%E6%8D%AE/
library(dplyr)
df %>% group_by(V2) %>%
  mutate(index = n()) %>%
  filter(index > 1) %>%
  select(2) %>%
  ungroup() %>%
  unique() %>%
  unlist() -> Dup.set

# x=11
UpdateGeneDUPE <- function(df,x) {
  if( (df[x,2] %in% Dup.set)== TRUE ){
    df[x,1] = df[x,1]
  }else{
    df[x,1] = df[x,2]
  }
  return(df[x,1])
}

df2 <- lapply(1:nrow(df), function(x)UpdateGeneDUPE(df,x))  %>% as.data.frame()
row.names(TTT) <- df2






## Update the genename ##
library(limma)
UpdateGene <- "No"  # UpdateGene <- c("Yes","No")
if(UpdateGene == "Yes"){
  ## Error ##  row.names(GeneExp.df) <- alias2Symbol(row.names(GeneExp.df), species = "Hs", expand.symbols = FALSE)
  row.names(GeneExp.df) <- alias2Symbol(row.names(GeneExp.df), species = "Hs", expand.symbols = FALSE)

}

## https://www.nature.com/articles/s41588-020-0669-3#Sec18
## Genenames.org: the HGNC and VGNC resources in 2021
## https://academic.oup.com/nar/article/49/D1/D939/5957168

# HGNC
# https://www.genenames.org/

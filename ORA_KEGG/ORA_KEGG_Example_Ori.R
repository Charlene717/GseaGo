## Ref https://sa123.cc/e8ntyq917kc26kellk68.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####
  # if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  # BiocManager::install("pathview")
  library(pathview)

##### Load Data  #####    
  data(gse16873.d)
  head(gse16873.d)
# 讀取自己的檔案可以使用類似下面的命令
# gse16873.d <- read.table("filename", sep="t", row.names=1, check.names=F, header=T)
  
##### Current path and new folder setting  ##### 
  Version = paste0(Sys.Date(),"_","KEGG")
  Save.Path = paste0(getwd(),"/",Version)
  dir.create(Save.Path) 
  
  
##### KEGG Analysis  ##### 
  p1 <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
                out.suffix = "gse16873", kegg.native = T,kegg.dir = Save.Path)
  
  p2 <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
                out.suffix = "gse16873.2layer", kegg.native = T, same.layer = F, kegg.dir = Save.Path)
  
  p3 <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
                out.suffix = "gse16873", kegg.native = F, sign.pos="bottomleft", kegg.dir = Save.Path)
  
  p4 <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
                out.suffix = "gse16873.2layer", kegg.native = F, sign.pos="bottomleft", 
                same.layer = F, kegg.dir = Save.Path)
  p5 <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
                out.suffix = "gse16873.split", kegg.native = F, sign.pos="bottomleft",
                split.group = T, kegg.dir = Save.Path)
  

##### KEGG Intedration Analysis  ##### 
  sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000)
  data(cpd.simtypes)
  head(sim.cpd.data)
  p6 <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
                pathway.id ="00640", species = "hsa", out.suffix = "gse16873.cpd",
                keys.align = "y", kegg.native = T, key.pos = "topright", kegg.dir = Save.Path)
  
  p7 <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
                pathway.id ="04110", species = "hsa", out.suffix = "gse16873.cpd",
                keys.align = "y", kegg.native = T, key.pos = "topright", kegg.dir = Save.Path)
  
  p8 <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data, pathway.id ="00640",
                species = "hsa", out.suffix = "gse16873.cpd", keys.align = "y", kegg.native = F,
                key.pos = "bottomright", sign.pos ="topright", cpd.lab.offset =-0.8)

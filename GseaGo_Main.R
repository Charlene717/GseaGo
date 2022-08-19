##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  if(!require("Seurat")) install.packages("Seurat")
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!require("Seurat")) install.packages("SeuratData")
  if(!require("tidyverse")) install.packages("patchwork")
  library(tidyverse)
  library(Seurat)
  library(SeuratData)
  library(patchwork)

##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")

##### Current path and new folder setting* #####
  ProjectName = "ifnb"
  Sampletype = "PBMC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

  ## Import information
  InputFolder = "Input_files_10x"
  InputAnno = "PBMC_Ano.csv"

  InputGSEA = "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"

##### Load dataset* #####
  # install dataset
  InstallData("ifnb")

  # load dataset
  LoadData("ifnb")

  # split the dataset into a list of two seurat objects (stim and CTRL)
  ifnb.list <- SplitObject(ifnb, split.by = "stim")

  # normalize and identify variable features for each dataset independently
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = ifnb.list)

##### Perform integration #####
  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
  # this command creates an 'integrated' data assay
  immune.combined <- IntegrateData(anchorset = immune.anchors)

##### Perform an integrated analysis #####
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(immune.combined) <- "integrated"

  # Run the standard workflow for visualization and clustering
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindClusters(immune.combined, resolution = 0.5)

  # Visualization
  p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
  p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
  p1 + p2

  p1 %>% BeautifyggPlot(.,LegPos = c(1, 0.5))

  DimPlot(immune.combined, reduction = "umap", split.by = "stim")

##### Identify conserved cell type markers #####
  # For performing differential expression after integration, we switch back to the original
  # data
  DefaultAssay(immune.combined) <- "RNA"
  nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
  head(nk.markers)

  FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                                            "CCL2", "PPBP"), min.cutoff = "q9")

  immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
                                  `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
                                  `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono Mk Doublets", `14` = "HSPC")
  immune.combined$celltype <- Idents(immune.combined)

  DimPlot(immune.combined, label = TRUE) %>% BeautifyggPlot(.,LegPos = c(1, 0.5))


  Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("HSPC", "Mono/Mk Doublets",
                                                                        "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
                                                                        "CD4 Naive T", "CD4 Memory T"))
  markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                       "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                       "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
  DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
    RotatedAxis()

##### Find Marker in different Cell type and VolcanoPlot (SPA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")
  immune.combined$celltype.Stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
  Idents(immune.combined) <- "celltype.Stim"

  DefaultAssay(immune.combined) <- "RNA"


  Idents(immune.combined) <- "celltype.Stim"
  #CellType.list <- as.character(unique(immune.combined@meta.data[["celltype"]]))
  dir.create(paste0(Save.Path,"/PBMC_SPA_FindMarkers"))

  CellType.list <- as.character(unique(immune.combined@meta.data[["celltype"]]))
  CCMarker_SPA.lt <- list()
  for(i in c(1:length(CellType.list))){
    try({
      CCMarker_SPA.lt[[i]] <- Find_Markers(immune.combined,
                                           paste0(CellType.list[i],"_STIM"),
                                           paste0(CellType.list[i],"_CTRL"),
                                           CellType.list[i],
                                           Path = Save.Path,
                                           ResultFolder = "PBMC_SPA_FindMarkers")

      # names(CCMarker_SPA.lt)[[i]] <- paste0("CCMarker_SPA.lt.",CellType.list[i])
      names(CCMarker_SPA.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i)

  CCMarker_SPA.lt <- CCMarker_SPA.lt[!unlist(lapply(CCMarker_SPA.lt,is.null))]


  ## Generate pdf and tif file for VolcanoPlot
  dir.create(paste0(Save.Path,"/PBMC_SPA_VolcanoPlot/"))

  pdf(file = paste0(Save.Path,"/PBMC_SPA_VolcanoPlot/PBMC_SPA_VolcanoPlot.pdf"),width = 7, height = 7 )
  for (i in 1:length(CellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_SPA.lt[[i]][["CCMarker.S"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+
              ggtitle(paste0("PBMC_",CellType.list[i]))
      )
    })
  }
  # graphics.off()
  dev.off()
  rm(i)

  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Save.Path,"/PBMC_SPA_VolcanoPlot/PBMC_SPA_VolcanoPlot",CellType.list[i],".tif"), width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_SPA.lt[[i]][["CCMarker.S"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Neg_List"]])+ ggtitle(paste0("PBMC_",CellType.list[i]))
      )

      graphics.off()
    })
  }
  rm(i)

  #### Save RData ####
  save.image(paste0(Save.Path,"/GseaGo.RData"))





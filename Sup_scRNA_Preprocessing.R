##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!require("Seurat")) install.packages("Seurat")
  if(!require("SeuratData")) install.packages("SeuratData")
  if(!require("patchwork")) install.packages("patchwork")

  library(tidyverse)
  library(Seurat)
  library(SeuratData)
  library(patchwork)

##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")
  # source("FUN_Find_Markers.R")
  # source("FUN_VolcanoPlot.R")
  # source("FUN_GSEA_LargeGeneSet.R")
  # source("FUN_GSEA_ggplot.R")

##### Current path and new folder setting* #####
  ProjectName = "ifnb"
  Sampletype = "PBMC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){ dir.create(Save.Path)}

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

  immune.combined <- RenameIdents(immune.combined, `0` = "CD14_Mono", `1` = "CD4_Naive_T", `2` = "CD4_Memory_T",
                                  `3` = "CD16_Mono", `4` = "B", `5` = "CD8_T", `6` = "NK", `7` = "T_activated", `8` = "DC", `9` = "B_Activated",
                                  `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono_Mk_Doublets", `14` = "HSPC")
  immune.combined$celltype <- Idents(immune.combined)

  DimPlot(immune.combined, label = TRUE) %>% BeautifyggPlot(.,LegPos = c(1, 0.5))
  DimPlot(immune.combined, label = TRUE,group.by = "celltype") %>% BeautifyggPlot(.,LegPos = c(1, 0.5))



  Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("HSPC", "Mono_Mk_Doublets",
                                                                        "pDC", "Eryth", "Mk", "DC", "CD14_Mono", "CD16_Mono", "B_Activated", "B", "CD8_T", "NK", "T_activated",
                                                                        "CD4_Naive_T", "CD4_Memory_T"))
  markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                       "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                       "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
  DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
    RotatedAxis()

##---------------------------------------------------------------------------------------------------------------##


# ##### Find Marker in different Cell type and VolcanoPlot (SPA) ########
#   ### Define group by different phenotype ###
#   source("FUN_Find_Markers.R")
#   immune.combined$celltype.Stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
#   Idents(immune.combined) <- "celltype.Stim"
#
#   DefaultAssay(immune.combined) <- "RNA"
#
#
#   Idents(immune.combined) <- "celltype.Stim"
#   #CellType.list <- as.character(unique(immune.combined@meta.data[["celltype"]]))
#   dir.create(paste0(Save.Path,"/PBMC_SPA_FindMarkers"))
#
#   CellType.list <- as.character(unique(immune.combined@meta.data[["celltype"]]))
#   CCMarker_SPA.lt <- list()
#   for(i in c(1:length(CellType.list))){
#     try({
#       CCMarker_SPA.lt[[i]] <- Find_Markers(immune.combined,
#                                            paste0(CellType.list[i],"_STIM"),
#                                            paste0(CellType.list[i],"_CTRL"),
#                                            CellType.list[i],
#                                            Path = Save.Path,
#                                            ResultFolder = "PBMC_SPA_FindMarkers")
#
#       # names(CCMarker_SPA.lt)[[i]] <- paste0("CCMarker_SPA.lt.",CellType.list[i])
#       names(CCMarker_SPA.lt)[[i]] <- paste0(CellType.list[i])
#     })
#   }
#   rm(i)
#
#   CCMarker_SPA.lt <- CCMarker_SPA.lt[!unlist(lapply(CCMarker_SPA.lt,is.null))]
#
#
#   ## Generate pdf and tif file for VolcanoPlot
#   dir.create(paste0(Save.Path,"/PBMC_SPA_VolcanoPlot/"))
#
#   pdf(file = paste0(Save.Path,"/PBMC_SPA_VolcanoPlot/PBMC_SPA_VolcanoPlot.pdf"),width = 7, height = 7 )
#   for (i in 1:length(CellType.list)) {
#     try({
#       print(VolcanoPlot(CCMarker_SPA.lt[[i]][["TarMarker.S"]],
#                         CCMarker_SPA.lt[[i]][["TarMarker.S_Pos_List"]],
#                         CCMarker_SPA.lt[[i]][["TarMarker.S_Neg_List"]], ShowGeneNum = 6)+
#               ggtitle(paste0("PBMC_",CellType.list[i]))
#       )
#     })
#   }
#   # graphics.off()
#   dev.off()
#   rm(i)
#
#   for (i in 1:length(CellType.list)) {
#     try({
#       tiff(file = paste0(Save.Path,"/PBMC_SPA_VolcanoPlot/PBMC_SPA_VolcanoPlot",CellType.list[i],".tif"), width = 17, height = 17, units = "cm", res = 200)
#       print(VolcanoPlot(CCMarker_SPA.lt[[i]][["TarMarker.S"]],
#                         CCMarker_SPA.lt[[i]][["TarMarker.S_Pos_List"]],
#                         CCMarker_SPA.lt[[i]][["TarMarker.S_Neg_List"]])+ ggtitle(paste0("PBMC_",CellType.list[i]))
#       )
#
#       graphics.off()
#     })
#   }
#   rm(i)
#
#
#
#
# ##### 09_1 GSEA Analysis (SPA) #####
#
#   ## Load the GSEA Dataset
#   load("GSEA_Analysis_Geneset.RData")
#   InputGSEA = "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"
#
#
#   ## Geneset from GSEA
#   # Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
#   Pathway.all <- read.delim2(paste0(getwd(),"/",InputGSEA),
#                              col.names = 1:max(count.fields(paste0(getwd(),"/",InputGSEA))),
#                              header = F,sep = "\t")
#
#
#   ## Run GSEA
#   GSEA_Large <- list()
#   GSEA_Large.df <- as.data.frame(matrix(nrow=0,ncol=10))
#   colnames(GSEA_Large.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
#   GSEA_Large.df.TOP <- GSEA_Large.df
#
#   dir.create(paste0(Save.Path,"/PBMC_GSEA"))
#
#
#   pdf(file = paste0(Save.Path, "/PBMC_GSEA/PBMC_GSEA_SPA.pdf"),width = 15, height = 7 )
#
#   for(i in 1:length(CellType.list)){
#
#     gseaDat <- CCMarker_SPA.lt[[paste0(CellType.list[i])]][["TarMarker.All"]]
#     gseaDat <- data.frame(row.names(gseaDat),gseaDat)
#     colnames(gseaDat)[[1]] <- c("Gene")
#     ranks <- gseaDat$avg_log2FC
#     names(ranks) <- gseaDat$Gene
#     # head(ranks)
#     # barplot(sort(ranks, decreasing = T))
#
#     GSEA_Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all,10)
#
#     fgseaRes <- GSEA_Large.Output[["fgseaRes"]]
#     # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
#
#     pathwaysH <- GSEA_Large.Output[["Pathway.all.list"]]
#
#     # plot.new()
#     # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
#
#     topPathways <- GSEA_Large.Output[["topPathways"]]
#
#     library(ggplot2)
#     plot.new()
#     plotGseaTable(pathwaysH[topPathways$pathway],
#                   ranks,
#                   fgseaRes,
#                   gseaParam = 0.5) + title( paste0("SC.",CellType.list[i]), adj = 0, line =3)
#
#     plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("SC.",CellType.list[i],": ",as.character(topPathways[1,1])))
#     #plotEnrichment_Pos1
#     plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("SC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
#     #plotEnrichment_Neg1
#
#     Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
#     names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
#     GSEA_Large[[i]] <- Sum
#     names(GSEA_Large)[[i]] <- paste0(CellType.list[i])
#
#     fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
#     colnames(fgseaRes2)[[1]] <- c("PhenoType")
#     GSEA_Large.df <- rbind(GSEA_Large.df,fgseaRes2 )
#
#     topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
#     colnames(topPathways2)[[1]] <- c("PhenoType")
#     GSEA_Large.df.TOP <- rbind(GSEA_Large.df.TOP, topPathways2)
#
#     rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
#
#   }
#
#   dev.off()
#
#   ## GSEA_Large.Sum.TOP ##
#   GSEA_Large.Sum.TOP <- rbind(GSEA_Large.df.TOP)
#   GSEA_Large.Sum.TOP <- GSEA_Large.Sum.TOP[,!colnames(GSEA_Large.Sum.TOP) %in% c("leadingEdge")]
#   write.table(GSEA_Large.Sum.TOP, file=paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Pathway_LargeTOP_SPA.txt"),sep="\t",
#               row.names=F, quote = FALSE)
#
#   ##### Bubble plot #####
#   library(ggplot2)
#   library(scales)
#   GSEA_Color.lt = list(high = "#ef476f",mid = "white",low = "#0077b6")
#
#   GSEA_Large.Sum.TOP$PhenoType <- factor(GSEA_Large.Sum.TOP$PhenoType,
#                                          levels = Cell_Type_Order.set)
#
#   GSEA_ggplot_SPA.lt <- GSEA_ggplot(GSEA_Large.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
#   GSEA_Large.Sum.TOP.S <- GSEA_ggplot_SPA.lt[["GSEA_TOP.df"]]
#   # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$NES) > 1,]
#   # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$padj) < 0.05,]
#   # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$padj) < 0.25,]
#   # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$pval) < 0.05,]
#
#   pdf(file = paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Bubble_SPA.pdf"),width = 17, height = 12 )
#   GSEA_ggplot_SPA.lt[["BBPlot_Ori"]]
#   GSEA_ggplot_SPA.lt[["BBPlot"]]
#   GSEA_ggplot_SPA.lt[["BBPlot2"]]
#   GSEA_ggplot_SPA.lt[["BBPlotB1"]]
#   GSEA_ggplot_SPA.lt[["BBPlotB1"]]
#   dev.off()
#
#
#   ##### Extract SubType #####
#   ## Mac
#   GSEA_Mac.df <- GSEA_Large.Sum.TOP.S[grep("Mac",GSEA_Large.Sum.TOP.S$PhenoType),]
#
#   BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
#     geom_point() +
#     scale_size_area(max_size = 5)+
#     scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
#                            guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
#     theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#
#   BBPlot_Mac
#
#   BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,OL_Thick = 1.5,
#                                                XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
#
#   BBPlot_MacB1 <- BBPlot_MacB %>%
#     insert_left(GSEA_ggplot_SPA.lt[["Y_Order"]],width = 0.2)
#   BBPlot_MacB1
#
#   pdf(file = paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Bubble_SPA_SubType_Mac.pdf"),width = 17, height = 20 )
#   BBPlot_MacB
#   BBPlot_MacB1
#   dev.off()
#
#
#   rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
#      df1.1,df1,BBPlot,BBPlot_Fib,BBPlot_FibB,BBPlot_T,BBPlot_TB)

  #### Save RData ####
  save.image(paste0(Save.Path,"/GseaGo_",ProjectName,"_",Sampletype,".RData"))





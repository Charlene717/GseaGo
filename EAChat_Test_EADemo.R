## Ref: https://zhuanlan.zhihu.com/p/432284358
## Ref: https://zhuanlan.zhihu.com/p/377356510 #With GOCluster


################################################################################
##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

Set_TarGenename <- "TOP2A"
Save_Path <- paste0("20231213_EAChat_Test_",Set_TarGenename)
if (!dir.exists(Save_Path)){dir.create(Save_Path)}   ## Create new folder

SetExport_Name <- "20231213_3"

##### Load Packages #####
source("FUN_Package_InstLoad.R")
PKG_Basic.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT")
PKG_BiocManager.set <- c("clusterProfiler","enrichplot","pathview","limma")

FUN_Package_InstLoad(Basic.set = PKG_Basic.set, BiocManager.set = PKG_BiocManager.set)

##### Function setting #####
## Call function
source("FUN_DistrPlot_GE.R")
source("FUN_Beautify_ggplot.R")
# source("FUN_Find_Markers.R")
# source("FUN_GSEA_LargeGeneSet.R")
# source("FUN_GSEA_ggplot.R")
source("FUN_ggPlot_vline.R")
source("FUN_GSEA_ANAL.R")
source("FUN_VolcanoPlot.R")

##### Import setting and data loading* #####
## Set Import path
SetInputPath_FOL <- "Input_TCGA"  # Input Folder Name
SetInput_GE <- "Xena_TCGA_LGG_GE"
SetInput_Meta <- "TCGA.LGG.sampleMap_LGG_clinicalMatrix"

## Load Gene expression file
GeneExp.df <- read.table(paste0(SetInputPath_FOL,"/",SetInput_GE), header=T, row.names = 1, sep="\t")
colnames(GeneExp.df) <-  gsub("\\.", "-", colnames(GeneExp.df))
## Load Metadata file
Metadata.df <- read.table(paste0(SetInputPath_FOL,"/",SetInput_Meta), header=T, sep="\t")
row.names(Metadata.df) <- Metadata.df[,1]

## Reorder the Metadata.df
Metadata.df <- left_join(data.frame("sampleID"=colnames(GeneExp.df)),
                         Metadata.df)
row.names(Metadata.df) <- Metadata.df[,1]


##### Conditions setting* #####
Set_Species <- "Homo sapiens"
Set_GroupMode <- "GoupByGeneExp"  # "GoupByPheno", "GoupByGeneExp"

if(Set_GroupMode == "GoupByPheno"){
  Set_GroupCond <-  list(GroupType = "sample_type",
                         GroupPair = c("Primary Tumor","Recurrent Tumor"))
}else if(Set_GroupMode == "GoupByGeneExp"){
  Set_TarGene_name = Set_TarGenename
  Set_TarGene <-  list(TarGeneName = Set_TarGene_name,
                       GEGroupMode = "Mean1SD", # c("Mean","Mean1SD","Mean2SD","Mean3SD","Median","Quartile","Customize"))
                       UpCutoff = 1, LowerCutoff = 1)

  Set_GroupCond <- list(GroupType = Set_TarGene_name, GroupPair = c("High","Low") )   ## DEG by GeneExp group
}else{
  print("Please set Set_GroupMode by GoupByPheno or GoupByGeneExp")
}

##### Grouping by GeneExp #####
if(Set_GroupMode == "GoupByGeneExp"){

  source("FUN_Group_GE.R")
  GeneExp_group.set <- FUN_Group_GE(GeneExp.df, Metadata.df,
                                    TarGeneName = Set_TarGene_name, GroupSet = Set_TarGene,
                                    Save.Path = Save_Path, ExportName = SetExport_Name)
  Metadata.df <- GeneExp_group.set[["AnnoNew.df"]]
  GeneExp_high.set <- GeneExp_group.set[["GeneExp_high.set"]]
  GeneExp_low.set <- GeneExp_group.set[["GeneExp_low.set"]]

}

## Set DEG Analysis
Set_DEGThr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) )
## - [] Add Metric for ranking gene
## - [] Modify DGE

##### Run Differential Expression Gene(DEG) analysis in R #####
#### Run DEG ####
source("RUN_DEG_Analysis.R")
DEG_Extract.df <- DEG_ANAL.lt[["DEG_Extract.df"]]
DEG_Extract_S.df <- DEG_Extract.df[abs(DEG_Extract.df$logFC) > 2 & DEG_Extract.df$FDR < 0.01,]

dif_mat <- DEG_Extract_S.df
# dif_mat <- GeneExp.df[1:1000,]

################################################################################
library(clusterProfiler) # Enrichment analysis R package
library(stringr) # Line break for labels
library(AnnotationDbi)
library(org.Hs.eg.db) # Reference species gene annotation database (hg19)
library(DOSE)
library(ggplot2) # Plotting
library(ggrepel) # Label-related

# Convert gene symbols to Entrez IDs to prevent analysis errors
id_list <- mapIds(org.Hs.eg.db, rownames(dif_mat), "ENTREZID", "SYMBOL")
# Remove cases where the conversion was not successful, i.e., id = NA
id_list <- na.omit(id_list)
# GO enrichment analysis
go <- enrichGO(gene = id_list, # Entrez ID list
               OrgDb = org.Hs.eg.db, # Specify species database
               keyType = "ENTREZID", # Specify the given name type
               ont = "ALL", # Optional: BP (Biological Process) / CC (Cellular Component) / MF (Molecular Function) / ALL (Specify all)
               pAdjustMethod = "BH", # P-value adjustment method, can also be FDR
               pvalueCutoff = 0.05, qvalueCutoff = 0.2, # p/q value thresholds
               readable = TRUE # Convert ID to symbol
)
go.res <- data.frame(go) # Convert GO results to a data frame for subsequent analysis (optional, depends on personal preference)
# write.csv(go.res, "Table_GO_result.csv", quote = FALSE) # Output GO enrichment analysis results
# Create a barplot for GO enrichment analysis. By default, select the top ten terms sorted by q-value for plotting.
goBP <- subset(go.res, subset = (ONTOLOGY == "BP"))[1:15,]
goCC <- subset(go.res, subset = (ONTOLOGY == "CC"))[1:15,]
goMF <- subset(go.res, subset = (ONTOLOGY == "MF"))[1:15,]
go.df <- rbind(goBP, goCC, goMF)
# Ensure that the order of GO terms in the plot matches the input
go.df$Description <- factor(go.df$Description, levels = rev(go.df$Description))
go.df <- go.df[!is.na(go.df$ONTOLOGY),]
# Create the plot
go_bar_ori <- ggplot(data = go.df, # Data for plotting
                     aes(x = reorder(Description, Count), y = Count)) +
  geom_bar(stat = "identity", width = 0.9) + # Create a bar chart with specified width
  coord_flip() + theme_bw() + # Flip the x and y axes and remove background color
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) + # Wrap term names when they are too long
  labs(x = "GO terms", y = "Gene Number", title = "Barplot of Enriched GO Terms") + # Set axis labels and title
  theme(axis.title = element_text(size = 13), # Axis title size
        axis.text = element_text(size = 11), # Axis label size
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # Title settings
        legend.title = element_text(size = 13), # Legend title size
        legend.text = element_text(size = 11), # Legend label size
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) # Plot margins
go_bar_ori


# Create the plot
go_bar <- ggplot(data = go.df, # Data for plotting
                 aes(x = Description, y = Count, fill = ONTOLOGY)) + # X-axis and color fill classification
  geom_bar(stat = "identity", width = 0.9) + # Create a bar chart with specified width
  coord_flip() + theme_bw() + # Flip the x and y axes and remove background color
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) + # Wrap term names when they are too long
  labs(x = "GO terms", y = "Gene Number", title = "Barplot of Enriched GO Terms") + # Set axis labels and title
  theme(axis.title = element_text(size = 13), # Axis title size
        axis.text = element_text(size = 11), # Axis label size
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # Title settings
        legend.title = element_text(size = 13), # Legend title size
        legend.text = element_text(size = 11), # Legend label size
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) # Plot margins
go_bar


## Sort
library(dplyr)

go.df <- go.df %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

go_bar <- ggplot(data = go.df, # Data for plotting
                 aes(x = Description, y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.9) +
  coord_flip() + theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(x = "GO terms", y = "Gene Number", title = "Barplot of Enriched GO Terms") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
go_bar


go_bar_ori[["data"]] %>% head()
write.table(go_bar_ori[["data"]], file = paste0(Save_Path,"/GO_Result_",SetExport_Name,".tsv"),
            sep="\t", row.names= F, quote = FALSE)


pdf(file = paste0(Save_Path,"/BarPlot_",SetExport_Name,".pdf"),
    width = 10,  height = 8)
go_bar_ori %>% print()
go_bar %>% print()
dev.off()

save.image(paste0(Save_Path,"/EAChat_",SetExport_Name,".RData"))

# ggsave(go_bar, filename = "GO_Barplot.pdf", width = 9, height = 7)
####################################################################################################

# Load the necessary library if not already loaded
# install.packages("data.table")
library(data.table)

# Define the terms
terms <- c(
  "vesicle transport along actin filament",
  "actin filament-based transport",
  "GPI anchor biosynthetic process",
  "GPI anchor metabolic process",
  "glycolipid biosynthetic process",
  "protein lipidation",
  "cell adhesion mediated by integrin",
  "phosphatidylinositol biosynthetic process",
  "phosphatidylinositol metabolic process",
  "glycerophospholipid biosynthetic process",
  "sodium:potassium-exchanging ATPase complex",
  "messenger ribonucleoprotein complex",
  "complex of collagen trimers",
  "SCF ubiquitin ligase complex",
  "microvillus",
  "myosin complex",
  "integrin complex",
  "protein complex involved in cell adhesion",
  "actin-based cell projection",
  "external side of plasma membrane",
  "chromatin-protein adaptor activity",
  "complement component C3b binding",
  "voltage-gated chloride channel activity",
  "insulin-like growth factor I binding",
  "voltage-gated monoatomic anion channel activity",
  "P-type ion transporter activity",
  "opsonin binding",
  "insulin-like growth factor binding",
  "microfilament motor activity",
  "integrin binding"
)

# Create groups based on common themes or categories
groups <- list(
  "Actin Filament" = c(
    "vesicle transport along actin filament",
    "actin filament-based transport",
    "microfilament motor activity"
  ),
  "GPI Anchor" = c(
    "GPI anchor biosynthetic process",
    "GPI anchor metabolic process",
    "glycolipid biosynthetic process",
    "protein lipidation"
  ),
  "Cell Adhesion" = c(
    "cell adhesion mediated by integrin",
    "protein complex involved in cell adhesion",
    "integrin binding"
  ),
  "Phosphatidylinositol" = c(
    "phosphatidylinositol biosynthetic process",
    "phosphatidylinositol metabolic process",
    "glycerophospholipid biosynthetic process"
  ),
  "Protein Complex" = c(
    "sodium:potassium-exchanging ATPase complex",
    "messenger ribonucleoprotein complex",
    "complex of collagen trimers",
    "SCF ubiquitin ligase complex",
    "myosin complex",
    "integrin complex"
  ),
  "Membrane" = c(
    "microvillus",
    "external side of plasma membrane",
    "voltage-gated chloride channel activity",
    "voltage-gated monoatomic anion channel activity",
    "P-type ion transporter activity"
  ),
  "Binding" = c(
    "chromatin-protein adaptor activity",
    "complement component C3b binding",
    "insulin-like growth factor I binding",
    "opsonin binding",
    "insulin-like growth factor binding"
  )
)

# Create a dataframe with an additional "Group" column
go.df$Group <- unlist(lapply(terms, function(term) {
  for (group in names(groups)) {
    if (term %in% groups[[group]]) {
      return(group)
    }
  }
  return("Other") # Set to NA if term does not belong to any group
}))

# Print the updated dataframe
print(go.df)


## Sort
library(dplyr)

go.df <- go.df %>%
  arrange(Group, Count) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

go_bar_GPT <- ggplot(data = go.df, # Data for plotting
                 aes(x = Description, y = Count, fill = Group)) +
  geom_bar(stat = "identity", width = 0.9) +
  coord_flip() + theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(x = "GO terms", y = "Gene Number", title = "Barplot of Enriched GO Terms") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
go_bar_GPT



###############################################################################
#### Cancer Related ####
# Load the necessary library if not already loaded
# install.packages("data.table")
library(data.table)

# Define the terms
terms <- c(
  "vesicle transport along actin filament",
  "actin filament-based transport",
  "GPI anchor biosynthetic process",
  "GPI anchor metabolic process",
  "glycolipid biosynthetic process",
  "protein lipidation",
  "cell adhesion mediated by integrin",
  "phosphatidylinositol biosynthetic process",
  "phosphatidylinositol metabolic process",
  "glycerophospholipid biosynthetic process",
  "sodium:potassium-exchanging ATPase complex",
  "messenger ribonucleoprotein complex",
  "complex of collagen trimers",
  "SCF ubiquitin ligase complex",
  "microvillus",
  "myosin complex",
  "integrin complex",
  "protein complex involved in cell adhesion",
  "actin-based cell projection",
  "external side of plasma membrane",
  "chromatin-protein adaptor activity",
  "complement component C3b binding",
  "voltage-gated chloride channel activity",
  "insulin-like growth factor I binding",
  "voltage-gated monoatomic anion channel activity",
  "P-type ion transporter activity",
  "opsonin binding",
  "insulin-like growth factor binding",
  "microfilament motor activity",
  "integrin binding"
)

# Create groups based on cancer-related mechanisms
cancer_mechanisms <- list(
  "Cell Adhesion and Migration" = c(
    "cell adhesion mediated by integrin",
    "protein complex involved in cell adhesion",
    "integrin binding"
  ),
  "Signal Transduction" = c(
    "chromatin-protein adaptor activity",
    "complement component C3b binding",
    "insulin-like growth factor I binding",
    "opsonin binding",
    "insulin-like growth factor binding"
  ),
  "Cellular Transport" = c(
    "vesicle transport along actin filament",
    "actin filament-based transport",
    "sodium:potassium-exchanging ATPase complex",
    "microvillus",
    "myosin complex",
    "integrin complex",
    "external side of plasma membrane",
    "voltage-gated chloride channel activity",
    "voltage-gated monoatomic anion channel activity",
    "P-type ion transporter activity"
  ),
  "Protein Modification" = c(
    "GPI anchor biosynthetic process",
    "GPI anchor metabolic process",
    "glycolipid biosynthetic process",
    "protein lipidation"
  ),
  "Cellular Complexes" = c(
    "messenger ribonucleoprotein complex",
    "complex of collagen trimers",
    "SCF ubiquitin ligase complex"
  ),
  "Cytoskeletal Changes" = c(
    "actin-based cell projection",
    "microfilament motor activity"
  )
)

# Create a dataframe with an additional "Cancer Mechanism" column
go.df$Cancer_Mechanism <- unlist(lapply(terms, function(term) {
  for (mechanism in names(cancer_mechanisms)) {
    if (term %in% cancer_mechanisms[[mechanism]]) {
      return(mechanism)
    }
  }
  return("Other") # Set to NA if term does not belong to any cancer mechanism group
}))

# Print the updated dataframe
print(go.df)

## Sort
library(dplyr)

go.df <- go.df %>%
  arrange(Cancer_Mechanism, Count) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

go_bar_GPT_Cancer <- ggplot(data = go.df, # Data for plotting
                     aes(x = Description, y = Count, fill = Cancer_Mechanism)) +
  geom_bar(stat = "identity", width = 0.9) +
  coord_flip() + theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(x = "GO terms", y = "Gene Number", title = "Barplot of Enriched GO Terms") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
go_bar_GPT_Cancer

################################################################################
# Load the necessary library if not already loaded
# install.packages("data.table")
library(data.table)

# Define the terms
terms <- c(
  "vesicle transport along actin filament",
  "actin filament-based transport",
  "GPI anchor biosynthetic process",
  "GPI anchor metabolic process",
  "glycolipid biosynthetic process",
  "protein lipidation",
  "cell adhesion mediated by integrin",
  "phosphatidylinositol biosynthetic process",
  "phosphatidylinositol metabolic process",
  "glycerophospholipid biosynthetic process",
  "sodium:potassium-exchanging ATPase complex",
  "messenger ribonucleoprotein complex",
  "complex of collagen trimers",
  "SCF ubiquitin ligase complex",
  "microvillus",
  "myosin complex",
  "integrin complex",
  "protein complex involved in cell adhesion",
  "actin-based cell projection",
  "external side of plasma membrane",
  "chromatin-protein adaptor activity",
  "complement component C3b binding",
  "voltage-gated chloride channel activity",
  "insulin-like growth factor I binding",
  "voltage-gated monoatomic anion channel activity",
  "P-type ion transporter activity",
  "opsonin binding",
  "insulin-like growth factor binding",
  "microfilament motor activity",
  "integrin binding"
)

# Create groups based on COVID-19 related mechanisms
covid_mechanisms <- list(
  "Viral Entry and Attachment" = c(
    "cell adhesion mediated by integrin",
    "protein complex involved in cell adhesion",
    "integrin binding"
  ),
  "Host-Pathogen Interaction" = c(
    "chromatin-protein adaptor activity",
    "complement component C3b binding",
    "opsonin binding",
    "insulin-like growth factor binding"
  ),
  "Cellular Transport" = c(
    "vesicle transport along actin filament",
    "actin filament-based transport",
    "sodium:potassium-exchanging ATPase complex",
    "microvillus",
    "myosin complex",
    "integrin complex",
    "external side of plasma membrane",
    "voltage-gated chloride channel activity",
    "voltage-gated monoatomic anion channel activity",
    "P-type ion transporter activity"
  ),
  "Protein Modification" = c(
    "GPI anchor biosynthetic process",
    "GPI anchor metabolic process",
    "glycolipid biosynthetic process",
    "protein lipidation"
  ),
  "Cellular Complexes" = c(
    "messenger ribonucleoprotein complex",
    "complex of collagen trimers",
    "SCF ubiquitin ligase complex"
  ),
  "Cytoskeletal Changes" = c(
    "actin-based cell projection",
    "microfilament motor activity"
  )
)

# Create a dataframe with an additional "Covid Mechanism" column
go.df$Covid_Mechanism <- unlist(lapply(terms, function(term) {
  for (mechanism in names(covid_mechanisms)) {
    if (term %in% covid_mechanisms[[mechanism]]) {
      return(mechanism)
    }
  }
  return("Other") # Set to NA if term does not belong to any COVID-19 mechanism group
}))

# Print the updated dataframe
print(go.df)

## Sort
library(dplyr)

go.df <- go.df %>%
  arrange(Covid_Mechanism, Count) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

go_bar_GPT_Covid <- ggplot(data = go.df, # Data for plotting
                            aes(x = Description, y = Count, fill = Covid_Mechanism)) +
  geom_bar(stat = "identity", width = 0.9) +
  coord_flip() + theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(x = "GO terms", y = "Gene Number", title = "Barplot of Enriched GO Terms") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
go_bar_GPT_Covid

################################################################################
#### Preliminary verification ####







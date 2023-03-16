## Volcano Plot
FUN_VolcanoPlot <- function(Marker.df,
                            DiffThr = list("log2FC",-2,2),
                            StatsTestThr = list("PValue",0.05),
                            color = c(High = "#ef476f",Mid = "gray",Low = "#0077b6"),     # color = c(High = "red",Mid = "gray",Low = "blue")
                            ShowGeneNumPos = 7, ShowGeneNumNeg = 7,
                            SizePoint = 3,  SizeAxisTitle = 16, SizeAxisText = 14, SizeLableText = 5,
                            ThkFrameLine = 2 , ThkThrLine = 0.8){


    #### Load Packages  #####
    if(!require("tidyverse")) install.packages("tidyverse")
    library(tidyverse)
    if(!require("cowplot")) install.packages("cowplot")
    library(cowplot)
    if(!require("ggrepel")) install.packages("ggrepel")
    library(ggrepel)
    # if(!require("ggplot2")) install.packages("ggplot2")
    # library(ggplot2)

    #### Add Gene name to column ####
    if("Gene" %in% colnames(Marker.df)){

    }else{
      Marker.df <- data.frame(row.names(Marker.df),Marker.df)
      colnames(Marker.df)[[1]] <- c("Gene")
    }

    #### Avoid 0 value in the Stats Test ####
    Marker.df[,StatsTestThr[[1]]] <- Marker.df[,StatsTestThr[[1]]] + 1.0e-300  # Marker.df$p_val <- Marker.df$p_val+1.0e-300

    #### Set color ####
    Marker.df$color <- ifelse(Marker.df[,StatsTestThr[[1]]] < StatsTestThr[[2]] & Marker.df[,DiffThr[[1]]] >= DiffThr[[3]],'High',
                              ifelse(Marker.df[,StatsTestThr[[1]]] < StatsTestThr[[2]] & Marker.df[,DiffThr[[1]]] <= DiffThr[[2]],'Low','Mid'))

    ## Old version (Symmetrical thresholding in logFC)
    # Marker.df$color <- ifelse(Marker.df$p_val< PValue & abs(Marker.df$avg_log2FC)>= log2FC,ifelse(Marker.df$avg_log2FC > log2FC,'High','Low'),'Mid')


    #### Arrange and filter ####
    Marker.df <- Marker.df %>% arrange(desc(Marker.df[,DiffThr[[1]]]))

    Pos.List <- Marker.df[rowSums(Marker.df[DiffThr[[1]]] > DiffThr[[3]] & Marker.df[StatsTestThr[[1]]] < StatsTestThr[[2]]) > 0, ] %>% rownames() # Pos.List <- Marker.df[rowSums(Marker.df["logFC"] >= 1) > 0, ] %>% rownames()
    Neg.List <- Marker.df[rowSums(Marker.df[DiffThr[[1]]] < DiffThr[[2]] & Marker.df[StatsTestThr[[1]]] < StatsTestThr[[2]]) > 0, ] %>% rownames()

    ## Record the name of the genes being displayed.
    # ShowGene_Pos.List <- row.names(Marker.df)[1:ShowGeneNumPos]
    # ShowGene_Neg.List <- row.names(Marker.df)[(nrow(Marker.df)-ShowGeneNumNeg+1):nrow(Marker.df)]


    if (length(Pos.List) >= ShowGeneNumPos && length(Neg.List) >= ShowGeneNumNeg) {
        Marker.df$genelabels <- factor(Marker.df$Gene,
                                        levels = c(Pos.List[1:ShowGeneNumPos],
                                                   Neg.List[(length(Neg.List)-(ShowGeneNumNeg-1)):length(Neg.List)]))
      }else if(length(Pos.List) >= ShowGeneNumPos && length(Neg.List) < ShowGeneNumNeg){
        Marker.df$genelabels <- factor(Marker.df$Gene,
                                        levels = c(Pos.List[1:ShowGeneNumPos],
                                                   Neg.List[(length(Neg.List)-(length(Neg.List)-1)):length(Neg.List)]))
      }else if(length(Pos.List) < ShowGeneNumPos && length(Neg.List) >= ShowGeneNumNeg){
        Marker.df$genelabels <- factor(Marker.df$Gene,
                                        levels = c(Pos.List[1:length(Pos.List)],
                                                   Neg.List[(length(Neg.List)-(ShowGeneNumNeg-1)):length(Neg.List)]))
      }else {
        Marker.df$genelabels <- factor(Marker.df$Gene,
                                        levels = c(Pos.List[1:length(Pos.List)],
                                                   Neg.List[(length(Neg.List)-(length(Neg.List)-1)):length(Neg.List)]))
      }

    ## Redefine levels:
    # Marker.df$genelabels <- factor(Marker.df$Gene, levels = c(Pos.List,Neg.List))

    #### Set Thr line ####
    Xintercept = c(DiffThr[[2]], DiffThr[[3]])  # Xintercept = c(-log2FC, log2FC)
    Yintercept = -log10(StatsTestThr[[2]]) # Yintercept = -log10(PValue)

    #### Volcano Plot ####
    library(ggrepel)
    VolcanoPlot <- ggplot(Marker.df, aes(Marker.df[,DiffThr[[1]]], -log10(Marker.df[,StatsTestThr[[1]]]), label = genelabels, col = color)) +
      theme_bw() + # Set to white background
      scale_color_manual(values = color) +
      geom_point(size = SizePoint) +
      geom_hline(yintercept = Yintercept, lty=8,col="black",lwd=ThkThrLine) +
      geom_vline(xintercept = Xintercept, lty=8,col="black",lwd=ThkThrLine) +
      theme(legend.position = "none",
            panel.grid=element_blank(),
            # text = element_text(size = 15),
            axis.title = element_text(size = SizeAxisTitle),
            axis.text = element_text(size = SizeAxisText)
            ) +
      labs(x=DiffThr[[1]], y= paste0("-log10 (",StatsTestThr[[1]],")")) + # labs(x="log2 (fold change)",y="-log10 (p-value)") +
      geom_text_repel(col = "#14213d", na.rm = TRUE,size = SizeLableText, box.padding = unit(0.45, "lines"), hjust = 1)+
      #geom_label(nudge_y = 2, alpha = 0.5)+
      theme(aspect.ratio=1) +
      theme(panel.border = element_rect(fill=NA,color="black", size= ThkFrameLine, linetype="solid"))  ## Ref: https://www.cnblogs.com/liujiaxin2018/p/14257944.html

    VolcanoPlot



return(VolcanoPlot_2)
}

#### To-Do List
## -[V] Clean up PKG Set
## -[V] Modify word size setting
## -[V] Clean up the code
## -[V] Annotation

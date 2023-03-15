FUN_VolcanoPlot <- function(Marker.df,
                            DiffThr = list("log2FC",-2,2),
                            StatsTestThr = list("PValue",0.05),
                            color = c(High = "#ef476f",Mid = "gray",Low = "#0077b6"),
                            ShowGeneNumPos = 7, ShowGeneNumNeg = 7){

    Xintercept = c(DiffThr[[2]], DiffThr[[3]])  # Xintercept = c(-log2FC, log2FC)
    Yintercept = -log10(StatsTestThr[[2]]) # Yintercept = -log10(PValue)
    #
    library(ggplot2)
    library(cowplot)

    Marker.df <- Marker.df %>% arrange(desc(Marker.df[,DiffThr[[1]]]))

    Pos.List <- Marker.df[rowSums(Marker.df[DiffThr[[1]]] > DiffThr[[3]] & Marker.df[StatsTestThr[[1]]] < StatsTestThr[[2]]) > 0, ] %>% rownames() # Pos.List <- Marker.df[rowSums(Marker.df["logFC"] >= 1) > 0, ] %>% rownames()
    Neg.List <- Marker.df[rowSums(Marker.df[DiffThr[[1]]] < DiffThr[[2]] & Marker.df[StatsTestThr[[1]]] < StatsTestThr[[2]]) > 0, ] %>% rownames()

    # ShowGene_Pos.List <- row.names(Marker.df)[1:ShowGeneNumPos]
    # ShowGene_Neg.List <- row.names(Marker.df)[(nrow(Marker.df)-ShowGeneNumNeg+1):nrow(Marker.df)]


    ##-------------- Volcano Plot --------------##

    if("Gene" %in% colnames(Marker.df)){

    }else{
      Marker.df <- data.frame(row.names(Marker.df),Marker.df)
      colnames(Marker.df)[[1]] <- c("Gene")
    }

    Marker.df[,StatsTestThr[[1]]] <- Marker.df[,StatsTestThr[[1]]] + 1.0e-300  # Marker.df$p_val <- Marker.df$p_val+1.0e-300

    Marker.df$color <- ifelse(Marker.df[,StatsTestThr[[1]]] < StatsTestThr[[2]] & Marker.df[,DiffThr[[1]]] >= DiffThr[[3]],'High',
                              ifelse(Marker.df[,StatsTestThr[[1]]] < StatsTestThr[[2]] & Marker.df[,DiffThr[[1]]] <= DiffThr[[2]],'Low','Mid'))

    # Marker.df$color <- ifelse(Marker.df[,StatsTestThr[[1]]] < StatsTestThr[[2]] & abs(Marker.df[,DiffThr[[1]]])>= DiffThr[[3]],ifelse(Marker.df[,DiffThr[[1]]] > DiffThr[[3]],'High','Low'),'Mid')
    # Marker.df$color <- ifelse(Marker.df$p_val< PValue & abs(Marker.df$avg_log2FC)>= log2FC,ifelse(Marker.df$avg_log2FC > log2FC,'High','Low'),'Mid')
    # color <- c(High = "red",Mid = "gray",Low = "blue")
    # color <- c(High = "#ef476f",Mid = "gray",Low = "#0077b6")

    # redefine levels:
    # Marker.df$genelabels <- factor(Marker.df$Gene, levels = c(Pos.List,Neg.List))

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


    library(ggrepel)
    VolcanoPlot <- ggplot(Marker.df, aes(Marker.df[,DiffThr[[1]]], -log10(Marker.df[,StatsTestThr[[1]]]), label = genelabels, col = color)) +
      geom_point(size = 3) +
      theme_bw() +
      scale_color_manual(values = color) +
      labs(x=DiffThr[[1]], y= paste0("-log10 (",StatsTestThr[[1]],")")) +
      # labs(x="log2 (fold change)",y="-log10 (p-value)") +
      geom_hline(yintercept = Yintercept, lty=8,col="black",lwd=0.8) +
      geom_vline(xintercept = Xintercept, lty=8,col="black",lwd=0.8) +
      theme(legend.position = "none",
            panel.grid=element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            text = element_text(size = 15)) +
      geom_point() +
      geom_text_repel(col = "#14213d", na.rm = TRUE,size = 5, box.padding = unit(0.45, "lines"), hjust = 1)+
      #geom_label(nudge_y = 2, alpha = 0.5)+
      theme(aspect.ratio=1)

    VolcanoPlot

    VolcanoPlot_2 <- VolcanoPlot + theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
    VolcanoPlot_2


    # https://www.cnblogs.com/liujiaxin2018/p/14257944.html

    # UMAP3 <- FeaturePlot(PBMC.combined, features = Pos.List[1], split.by = SplitBy, max.cutoff = 3,
    #                      cols = c("grey","#de3767", "red"), ncol = 2)
    # UMAP4 <- FeaturePlot(PBMC.combined, features = Neg.List[length(Neg.List)], split.by = SplitBy, max.cutoff = 3,
    #                      cols = c("grey", "blue"), ncol = 2)

return(VolcanoPlot_2)
}

#### To-Do List
## -[] Clean up PKG Set
## -[] Modify word size setting



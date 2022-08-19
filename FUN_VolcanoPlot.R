VolcanoPlot <- function(Marker.df, Pos.List, Neg.List,
                        color = c(red = "#ef476f",gray = "gray",blue = "#0077b6"),
                        log2FC = 1,PValue = 0.05,
                        ShowGeneNum = 5){

    Xintercept = c(-log2FC, log2FC)
    Yintercept = -log10(PValue)  
    #
    library(ggplot2)
    library(cowplot)
    
    ##-------------- Volcano Plot --------------##
    Marker.df2 <- data.frame(row.names(Marker.df),Marker.df)
    colnames(Marker.df2)[[1]] <- c("Gene")
    
    Marker.df2$p_val <- Marker.df2$p_val+1.0e-300
    
    Marker.df2$color <- ifelse(Marker.df2$p_val< PValue & abs(Marker.df2$avg_log2FC)>= log2FC,ifelse(Marker.df2$avg_log2FC > log2FC,'red','blue'),'gray')
    # color <- c(red = "red",gray = "gray",blue = "blue")
    # color <- c(red = "#ef476f",gray = "gray",blue = "#0077b6")
    
    # redefine levels:
    # Marker.df2$genelabels <- factor(Marker.df2$Gene, levels = c(Pos.List,Neg.List))
    
    if (length(Pos.List) >= ShowGeneNum && length(Neg.List) >= ShowGeneNum) {
        Marker.df2$genelabels <- factor(Marker.df2$Gene, 
                                        levels = c(Pos.List[1:ShowGeneNum],
                                                   Neg.List[(length(Neg.List)-(ShowGeneNum-1)):length(Neg.List)]))
      }else if(length(Pos.List) >= ShowGeneNum && length(Neg.List) < ShowGeneNum){
        Marker.df2$genelabels <- factor(Marker.df2$Gene, 
                                        levels = c(Pos.List[1:ShowGeneNum],
                                                   Neg.List[(length(Neg.List)-(length(Neg.List)-1)):length(Neg.List)]))
      }else if(length(Pos.List) < ShowGeneNum && length(Neg.List) >= ShowGeneNum){
        Marker.df2$genelabels <- factor(Marker.df2$Gene, 
                                        levels = c(Pos.List[1:length(Pos.List)],
                                                   Neg.List[(length(Neg.List)-(ShowGeneNum-1)):length(Neg.List)]))
      }else {
        Marker.df2$genelabels <- factor(Marker.df2$Gene, 
                                        levels = c(Pos.List[1:length(Pos.List)],
                                                   Neg.List[(length(Neg.List)-(length(Neg.List)-1)):length(Neg.List)]))
      }
    
    
    library(ggrepel)
    VolcanoPlot <- ggplot(Marker.df2, aes(avg_log2FC, -log10(p_val), label = genelabels, col = color)) +
      geom_point(size = 3) +
      theme_bw() +
      scale_color_manual(values = color) +
      labs(x="log2 (fold change)",y="-log10 (p-value)") +
      geom_hline(yintercept = Yintercept, lty=8,col="black",lwd=0.8) +
      geom_vline(xintercept = Xintercept, lty=8,col="black",lwd=0.8) +
      theme(legend.position = "none",
            panel.grid=element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            text = element_text(size = 15)) + 
      geom_point() +
      geom_text_repel(col = "#14213d", na.rm = TRUE,size = 6, box.padding = unit(0.45, "lines"), hjust = 1)+ 
      #geom_label(nudge_y = 2, alpha = 0.5)+
      theme(aspect.ratio=1)
    
    VolcanoPlot_2 <- VolcanoPlot + theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
    # https://www.cnblogs.com/liujiaxin2018/p/14257944.html 
    
    
    # UMAP3 <- FeaturePlot(PBMC.combined, features = Pos.List[1], split.by = SplitBy, max.cutoff = 3,
    #                      cols = c("grey","#de3767", "red"), ncol = 2)
    # UMAP4 <- FeaturePlot(PBMC.combined, features = Neg.List[length(Neg.List)], split.by = SplitBy, max.cutoff = 3,
    #                      cols = c("grey", "blue"), ncol = 2)

return(VolcanoPlot_2)
}




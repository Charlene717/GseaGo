ggPlot_vline = function( TGeneDen.p,
                         data,
                         Line.clr = list(rect="#ecbdfc", line="#994db3",text="#6a3b7a" ),
                         Line1 = Target_gene_Mean+Target_gene_SD,
                         Line2 = Target_gene_Mean,
                         Line3 = Target_gene_Mean-Target_gene_SD,
                         Text.set = c("Mean+SD","Mean","Mean-SD"),
                         Text.yPos = 0.5,
                         rectP = list(xWidth=0.03, yminP=0.45, ymaxP=0.55,alpha=0.8) 
                        ){
  
  TGeneDen.p +
    # Add line for Mean+SD
    geom_vline(aes(xintercept=c(Line1)),
               colour = Line.clr[["line"]], linetype="dashed")+
    # Add line for Mean
    geom_vline(aes(xintercept=c(Line2)),
               colour = Line.clr[["line"]],linetype="solid") +
    # Add block for Mean-SD
    geom_vline(aes(xintercept=c(Line3)),
               colour = Line.clr[["line"]],linetype="dashed")+
    
    # Add block for Mean+SD
    annotate(geom = "rect", 
             xmin = Line1 - rectP[["xWidth"]]*max(data),
             xmax = Line1 + rectP[["xWidth"]]*max(data),
             ymin = max(density(data$value)$y*rectP[["yminP"]]), 
             ymax = max(density(data$value)$y*rectP[["ymaxP"]]),
             fill = Line.clr[["rect"]], alpha = rectP[["alpha"]])+ 
    # Add block for Mean
    annotate(geom = "rect", 
             xmin = Line2 - rectP[["xWidth"]]*max(data),
             xmax = Line2 + rectP[["xWidth"]]*max(data),
             ymin = max(density(data$value)$y*rectP[["yminP"]]), 
             ymax = max(density(data$value)$y*rectP[["ymaxP"]]),
             fill = Line.clr[["rect"]], alpha = rectP[["alpha"]])+ 
    # Add block for Mean-SD
    annotate(geom = "rect", 
             xmin = Line3 - rectP[["xWidth"]]*max(data),
             xmax = Line3 + rectP[["xWidth"]]*max(data),
             ymin = max(density(data$value)$y*rectP[["yminP"]]), 
             ymax = max(density(data$value)$y*rectP[["ymaxP"]]),
             fill = Line.clr[["rect"]], alpha = rectP[["alpha"]])+
    
    # Add Annotation for Mean+SD
    geom_text(aes(x = Line1, 
                  y = max(density(data$value)$y)*Text.yPos, 
                  label = paste0(Text.set[1],"\n",
                                 round(Line1, digits = 2))),
              colour=Line.clr[["text"]])+
    # Add Annotation for Mean
    geom_text(aes(x = Line2, 
                  y = max(density(data$value)$y)*Text.yPos, 
                  label = paste0(Text.set[2],"\n",
                                 round(Line2, digits = 2))),
              colour=Line.clr[["text"]])+
    # Add Annotation for Mean-SD
    geom_text(aes(x = Line3, 
                  y = max(density(data$value)$y)*Text.yPos, 
                  label = paste0(Text.set[3],"\n",
                                 round(Line3, digits = 2))),
              colour=Line.clr[["text"]]) -> TGeneDen_SD.p

  return(TGeneDen_SD.p)
}

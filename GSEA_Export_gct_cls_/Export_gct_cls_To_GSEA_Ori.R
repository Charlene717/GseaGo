##### Presetting ######
  ## Clear variables
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Version information ######  
  # platform       x86_64-w64-mingw32          
  # arch           x86_64                      
  # os             mingw32                     
  # system         x86_64, mingw32             
  # status                                     
  # major          4                           
  # minor          1.2                         
  # year           2021                        
  # month          11                          
  # day            01                          
  # svn rev        81115                       
  # language       R                           
  # version.string R version 4.1.2 (2021-11-01)
  # nickname       Bird Hippie  
  
##### Load library #####  
  library(data.table)
  library("plyr") 
  library(dplyr) # Basic data manupilation tools
  library(tidyverse)
  library(eoffice) # Export plot to PPT
  
##### Function setting #####
  source("FUN_Beautify_ggplot.R")
  
##### Files setting and import * ##### 
  ## File setting*
  FileName <- "Xena_TCGA_LGG_GE"

  ## Import genetic data file
  GeneExp.df <- read.table(FileName, header=T, row.names = 1, sep="\t")
  
##### Conditions setting* ##### 
  Target_gene_name <- "TOP2A"
  Mode_Group <- list(Mode="Mean",SD=1) # Mode_Group <- list(Mode="Quartile",Q2="Only")
  
##### Current path and new folder setting* ##### 
  Result_Folder_Name <- paste0(Target_gene_name,"_",Sys.Date()) ## Generate output folder automatically
  dir.create(Result_Folder_Name)
  
  
  
##### Extract Target gene and Statistics ####
  # Extract data with Target_gene_name
  Target_gene_Mean <- GeneExp.df[Target_gene_name,] %>%
    as.numeric() %>% 
    mean()
  
  #rowMeans(data.matrix(Target_gene))
  Target_gene_SD <- GeneExp.df[Target_gene_name,] %>%
    as.numeric() %>% 
    sd()
  
  # Quartile
  Target_gene_Q <- GeneExp.df[Target_gene_name,] %>%
    as.numeric() %>% 
    quantile()


##### Group the expression matrix according to the expression level of Target gene ####  
  if(Mode_Group$Mode=="Mean"){
    if(Mode_Group$SD==0){
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD*(Mode_Group$SD)]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] < Target_gene_Mean-Target_gene_SD*(Mode_Group$SD)]
    }else{
    GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD*(Mode_Group$SD)]
    GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] <= Target_gene_Mean-Target_gene_SD*(Mode_Group$SD)]
    }
    #rm(Target_gene_Mean, Target_gene_SD)
    
  }else{
    if(Mode_Group$Q2=="Only"){ # Mode="Quartile"
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Q[3]]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] < Target_gene_Q[3]]
    }else{
      GeneExp_high.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] >= Target_gene_Q[4]]
      GeneExp_low.set <- colnames(GeneExp.df)[GeneExp.df[Target_gene_name,] <= Target_gene_Q[2]]
    }
    #rm(Target_gene_Q)
    
  }
  
##### Build Expression matrix for GSEA #####    
  GeneExp_GSEA.df <- cbind(
    NAME=row.names(GeneExp.df), 
    Description = rep("na", nrow(GeneExp.df)),
    GeneExp.df[,c(GeneExp_high.set, GeneExp_low.set)]
  )
  
  GSEA_SampleCol.df <- data.frame(t(colnames(GeneExp_GSEA.df)), stringsAsFactors=FALSE)
  colnames(GSEA_SampleCol.df) <- GSEA_SampleCol.df
  
  GeneExp_GSEA.df <- rbind(GSEA_SampleCol.df,GeneExp_GSEA.df)
  
  GeneExp_GSEA.df <- data.frame(
    "NAME" = c("#1.2",nrow(GeneExp.df)),
    "Description" = c('',length(c(GeneExp_high.set, GeneExp_low.set)))
  ) %>% 
    rbind.fill(GeneExp_GSEA.df)
  rm(GSEA_SampleCol.df)
  
##### Build Group Files #####
  ## Set the group array 
  Pheno_sum.df <- c(ncol(GeneExp_GSEA.df)-2,2,1) %>% t() %>% data.frame() %>% 
    rbind.fill(c(paste0("#",Target_gene_name,"_high"),paste0(Target_gene_name,"_Low")) %>% t() %>% data.frame(stringsAsFactors=FALSE)) %>% 
    rbind.fill(c(rep(0,length(GeneExp_high.set)),rep(1,length(GeneExp_low.set))) %>% t() %>% data.frame())
  
  
##### Export Result ##### 
  if(Mode_Group$Mode=="Mean"){  
    write.table(
      GeneExp_GSEA.df,
      file=paste0(Result_Folder_Name,"/",FileName,"_",
                  Mode_Group$Mode,Mode_Group$SD,"SD_",
                  Target_gene_name,"_collapsed.gct"),
      quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t'
    )
    write.table(
      Pheno_sum.df,
      file=paste0(Result_Folder_Name,"/",FileName,"_",
                  Mode_Group$Mode,Mode_Group$SD,"SD_",
                  Target_gene_name,".cls"),
      quote = FALSE,row.names = FALSE, na = "",col.names = FALSE
    )
  }else{
    write.table(
      GeneExp_GSEA.df,
      file=paste0(Result_Folder_Name,"/",FileName,"_",
                  Mode_Group$Mode,"Q2",Mode_Group$Q2,"_",
                  Target_gene_name,"_collapsed.gct"),
      quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t'
    )
    write.table(
      Pheno_sum.df,
      file=paste0(Result_Folder_Name,"/",FileName,"_",
                  Mode_Group$Mode,"Q2",Mode_Group$Q2,"_",
                  Target_gene_name,".cls"),
      quote = FALSE,row.names = FALSE, na = "",col.names = FALSE
    )
    
  }

##### Visualization #####  
  ## https://www.jianshu.com/p/9e5b7ffcf80f
  
  data <- reshape2::melt(GeneExp.df[Target_gene_name,]%>%
                           as.numeric())
  TGeneDen.p <- ggplot(data,aes(value,fill=value, color=value)) + 
    xlab("Expression level") + 
    geom_density(alpha = 0.6, fill = "lightgray") + 
    geom_rug() + theme_bw()
  
  ## Set the color
  Mean_SD.clr <- list(rect="#ecbdfc", line="#994db3",text="#6a3b7a" )
  Mean_Q.clr <- list(rect="#abede1", line="#12705f",text="#12705f" )
  
  TGeneDen.p 
  TGeneDen.p +
    # Add line for Mean+SD
    geom_vline(aes(xintercept=c(Target_gene_Mean+Target_gene_SD)),
               colour = Mean_SD.clr[["line"]], linetype="dashed")+
    # Add line for Mean
    geom_vline(aes(xintercept=c(Target_gene_Mean)),
               colour = Mean_SD.clr[["line"]],linetype="solid") +
    # Add block for Mean-SD
    geom_vline(aes(xintercept=c(Target_gene_Mean-Target_gene_SD)),
               colour = Mean_SD.clr[["line"]],linetype="dashed")+
    
    # Add block for Mean+SD
    annotate(geom = "rect", 
             xmin = Target_gene_Mean+Target_gene_SD-0.03*max(data),
             xmax = Target_gene_Mean+Target_gene_SD+0.03*max(data),
             ymin = max(density(data$value)$y*0.45), 
             ymax = max(density(data$value)$y*0.55),
             fill = Mean_SD.clr[["rect"]], alpha = 0.8)+ 
    # Add block for Mean
    annotate(geom = "rect", 
             xmin = Target_gene_Mean-0.02*max(data),
             xmax = Target_gene_Mean+0.02*max(data),
             ymin = max(density(data$value)$y*0.45), 
             ymax = max(density(data$value)$y*0.55),
             fill = Mean_SD.clr[["rect"]], alpha = 0.8)+ 
    # Add block for Mean-SD
    annotate(geom = "rect", 
             xmin = Target_gene_Mean-Target_gene_SD-0.03*max(data),
             xmax = Target_gene_Mean-Target_gene_SD+0.03*max(data),
             ymin = max(density(data$value)$y*0.45), 
             ymax = max(density(data$value)$y*0.55),
             fill = Mean_SD.clr[["rect"]], alpha = 0.8)+
    
    # Add Annotation for Mean+SD
    geom_text(aes(x = Target_gene_Mean+Target_gene_SD, 
                  y = max(density(data$value)$y)/2, 
                  label = paste0("Mean+SD\n",
                                 round(Target_gene_Mean+Target_gene_SD, digits = 2))),
              colour=Mean_SD.clr[["text"]])+
    # Add Annotation for Mean
    geom_text(aes(x = Target_gene_Mean, 
                  y = max(density(data$value)$y)/2, 
                  label = paste0("Mean\n",
                                 round(Target_gene_Mean, digits = 2))),
              colour=Mean_SD.clr[["text"]])+
    # Add Annotation for Mean-SD
    geom_text(aes(x = Target_gene_Mean-Target_gene_SD, 
                  y = max(density(data$value)$y)/2, 
                  label = paste0("Mean-SD\n",
                                 round(Target_gene_Mean-Target_gene_SD, digits = 2))),
              colour=Mean_SD.clr[["text"]]) -> TGeneDen_SD.p
  TGeneDen_SD.p + labs(title= Target_gene_name,
                       x ="Expression level", y = "Density")
  
  TGeneDen.p +
    # Add line for Q1
    geom_vline(aes(xintercept=c(Target_gene_Q[2])),
               colour="#12705f",linetype="dashed")+
    
    # Add line for Q2
    geom_vline(aes(xintercept=c(Target_gene_Q[3])),
               colour="#12705f",linetype="solid") +
    
    # Add line for Q3
    geom_vline(aes(xintercept=c(Target_gene_Q[4])),
               colour="#12705f",linetype="dashed")+
    
    # Add block for Q1
    annotate(geom = "rect", 
             xmin = Target_gene_Q[2]-0.015*max(data),
             xmax = Target_gene_Q[2]+0.015*max(data),
             ymin = max(density(data$value)$y*0.45), 
             ymax = max(density(data$value)$y*0.55),
             fill = "#abede1", alpha = 0.8)+
    # Add block for Q2
    annotate(geom = "rect", 
             xmin = Target_gene_Q[3]-0.015*max(data),
             xmax = Target_gene_Q[3]+0.015*max(data),
             ymin = max(density(data$value)$y*0.45), 
             ymax = max(density(data$value)$y*0.55),
             fill = "#abede1", alpha = 0.8)+ 
    # Add block for Q3
    annotate(geom = "rect", 
             xmin = Target_gene_Q[4]-0.015*max(data),
             xmax = Target_gene_Q[4]+0.015*max(data),
             ymin = max(density(data$value)$y*0.45), 
             ymax = max(density(data$value)$y*0.55),
             fill = "#abede1", alpha = 0.8)+ 
    
    # Add Annotation for Q1
    geom_text(aes(x = Target_gene_Q[2], 
                  y = max(density(data$value)$y)/2, 
                  label = paste0("Q1\n",
                                 round(Target_gene_Q[2], digits = 2))),
              colour="#12705f")+
    # Add Annotation for Q2
    geom_text(aes(x = Target_gene_Q[3], 
                  y = max(density(data$value)$y)/2, 
                  label = paste0("Q2\n",
                                 round(Target_gene_Q[3], digits = 2))),
              colour="#12705f")+
    
    # Add Annotation for Q3
    geom_text(aes(x = Target_gene_Q[4], 
                  y = max(density(data$value)$y)/2, 
                  label = paste0("Q3\n",
                                 round(Target_gene_Q[4], digits = 2))),
              colour="#12705f") -> TGeneDen_Q.p
  TGeneDen_Q.p + labs(title= Target_gene_name,
                      x ="Expression level", y = "Density")
  
  TGeneDen_SD.p+
    # Add line for Q1
    geom_vline(aes(xintercept=c(Target_gene_Q[2])),
               colour="#12705f",linetype="dashed")+
    
    # Add line for Q2
    geom_vline(aes(xintercept=c(Target_gene_Q[3])),
               colour="#12705f",linetype="solid") +
    
    # Add line for Q3
    geom_vline(aes(xintercept=c(Target_gene_Q[4])),
               colour="#12705f",linetype="dashed")+
    
    # Add block for Q1
    annotate(geom = "rect", 
             xmin = Target_gene_Q[2]-0.015*max(data),
             xmax = Target_gene_Q[2]+0.015*max(data),
             ymin = max(density(data$value)$y*0.3), 
             ymax = max(density(data$value)$y*0.4),
             fill = "#abede1", alpha = 0.8)+
    # Add block for Q2
    annotate(geom = "rect", 
             xmin = Target_gene_Q[3]-0.015*max(data),
             xmax = Target_gene_Q[3]+0.015*max(data),
             ymin = max(density(data$value)$y*0.3), 
             ymax = max(density(data$value)$y*0.4),
             fill = "#abede1", alpha = 0.8)+ 
    # Add block for Q3
    annotate(geom = "rect", 
             xmin = Target_gene_Q[4]-0.015*max(data),
             xmax = Target_gene_Q[4]+0.015*max(data),
             ymin = max(density(data$value)$y*0.3), 
             ymax = max(density(data$value)$y*0.4),
             fill = "#abede1", alpha = 0.8)+ 
    
    # Add Annotation for Q1
    geom_text(aes(x = Target_gene_Q[2], 
                  y = max(density(data$value)$y)*0.35, 
                  label = paste0("Q1\n",
                                 round(Target_gene_Q[2], digits = 2))),
              colour="#12705f")+
    # Add Annotation for Q2
    geom_text(aes(x = Target_gene_Q[3], 
                  y = max(density(data$value)$y)*0.35, 
                  label = paste0("Q2\n",
                                 round(Target_gene_Q[3], digits = 2))),
              colour="#12705f")+
    
    # Add Annotation for Q3
    geom_text(aes(x = Target_gene_Q[4], 
                  y = max(density(data$value)$y)*0.35, 
                  label = paste0("Q3\n",
                                 round(Target_gene_Q[4], digits = 2))),
              colour="#12705f") -> TGeneDen_SD_Q.p
  TGeneDen_SD_Q.p
  
  
  pdf(
    file = paste0(Result_Folder_Name,"/",FileName,"_",Target_gene_name,"_DensityPlot.pdf"),
    width = 10,  height = 8
  )
  
    TGeneDen_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) + 
      labs(title= Target_gene_name,
           x ="Expression level", y = "Density")
    TGeneDen_SD.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) + 
      labs(title= Target_gene_name,
           x ="Expression level", y = "Density")
    TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7) + 
      labs(title= Target_gene_name,
           x ="Expression level", y = "Density")
  
  dev.off()
  
  # Export PPT
  TGeneDen_SD_Q.p  %>% BeautifyggPlot(LegPos = c(0.9, 0.8),AxisTitleSize=1.7,
                                      OL_Thick = 1.5) + 
    labs(title= Target_gene_name,
         x ="Expression level", y = "Density") -> TGeneDen_SD_Q2.p
  
  topptx(TGeneDen_SD_Q2.p,paste0(Result_Folder_Name,"/",FileName,"_",Target_gene_name,"_DensityPlot.pptx"))
  
  ## Finding Peak Values For a Density Distribution
  # http://ianmadd.github.io/pages/PeakDensityDistribution.html
  which.max(density(data$value)$y)
  max(density(data$value)$y)
  
  ## Plot multiple gene
  
    
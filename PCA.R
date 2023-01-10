#####GGPLOT FOR PCA
##Add Color Brewer function for 
ggpca=function(raw_dt,pheno_factor1,pheno_factor2,pheno_title1=NULL,pheno_title2=NULL,logT=F,scale_plot=NULL){
require(ggsci)
if(logT==T){
  exp_raw <- log2(raw_dt)
  tit1="PCA plot of the log-transformed raw expression data"
}else{
  exp_raw=raw_dt
  tit1="PCA plot of the expression data"

}
   
  PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
  
  percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  
  
  
  dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                       Factor1 = pheno_factor1,
                       Factor2 = pheno_factor2)
  
  if(missing(pheno_title1) & missing(pheno_title2)){
    pheno_title1="Factor1"
    pheno_title2="Factor2"
  }
  
  fcol=length(unique(dataGG$Factor2))
  if(fcol>10){
    fcol=10
  }
  mypal = pal_npg("nrc", alpha = 0.8)(fcol)
  
  
 #c("darkorange2", "dodgerblue4")
  plot1=ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(shape = Factor1, colour = Factor2)) +
    ggtitle(tit1) +
    xlab(paste0("PC1, Variance Explained: ", percentVar[1], "%")) +
    ylab(paste0("PC2, Variance Explained: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_shape_manual(pheno_title1,values = c(1:length(unique(dataGG$Factor1)))) + 
    scale_color_manual(pheno_title2,values = mypal)
    if(hasArg(scale_plot)){
    plot1 + coord_fixed(ratio = sd_ratio) 
  }else{
    plot1
  }                 
}

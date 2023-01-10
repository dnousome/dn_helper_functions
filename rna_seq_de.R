
##Run RSEM DESeq2Run the FXN loop
defxn=function(counts,covar=covar,variable,test=NULL,readscutoff=5,num="Yes",denom="No"){
  dds = DESeqDataSetFromMatrix(counts, covar, formula(paste("~",variable)))
  keep = rowSums(counts >=10 ) >= readscutoff
  dds=dds[keep,]
  
  if(missing(test)){
    dds = DESeq(dds)
  }else if(test=="LRT"){
    dds = DESeq(dds,test = "LRT",reduced = ~1)
  }
  
  results(dds,contrast=c(variable,num,denom)) %>% data.frame() %>% 
    rownames_to_column("Gene") %>% separate(Gene,into=c("gencode","symbol"),extra = 'merge',sep="_")  %>%
    arrange(padj)
}


##Run LIMMA-VOOM
voomfxn=function(counts,covar,variable,matform,cutoff=1){
  
  d0 <- DGEList(counts)
  dge <- calcNormFactors(d0,method="TMM")
  
  drop <- apply(cpm(dge), 1, max) < cutoff
  d <- dge[!drop,] 
  mm <- model.matrix(~CBene1,covar)
  
  #colnames(mm)=make.names(colnames(mm))
  v = voom(d, mm,plot = F)
  fit = lmFit(v, mm)
  
  
  fit.cont <- contrasts.fit(fit, coef=2)
  fit2 = eBayes(fit.cont)
  
  
  
  
  topTable(fit2, sort.by = "P", n = Inf) %>% 
    select(everything(),pvalue=P.Value,padj=adj.P.Val,log2FoldChange=logFC) %>%
    rownames_to_column("Gene") %>% 
    separate(Gene,into=c("gencode","symbol"),extra = 'merge',sep="_")  %>%
    arrange(pvalue) 
}


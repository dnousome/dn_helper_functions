make_oncoplot3=function(maf.filtered, cohort_freq_thresh = 0.1, auto_adjust_threshold=T,
         use_clinvar_anno=F, title_text="",
         genes_to_plot=NULL, include_all=F,
         show_sample_names=NULL,
         clin_data=NULL,clin_data_colors=NULL,
         custom_column_order=NULL, full_output=F,
         savename=NULL,plot_height=800,plot_width=700,pointsize=12) {
  
  ### Read in MAF file
  # maf.filtered <- read.maf(maf_file)
  if (! is.null(genes_to_plot)) {
    if (class(genes_to_plot)=="character") {
      if (length(genes_to_plot)==1) {
        ## Then it's either a file name or a single gene
        if (file.exists(genes_to_plot)) {
          ### Need to parse file type and read accordingly; assuming tsv for now
          gene_data <- read.table(genes_to_plot,sep="\t", header=T)
          if (sum(c("Hugo_Symbol","Reason") %in% colnames(gene_data)) != 2) {
            stop("Can't find Hugo Symbol or Reason in custom gene input.")
          }
          genes_for_oncoplot <- gene_data
        } else {
          stop(paste0("Can't find file: ",genes_to_plot))
        }
      } else {
        genes_for_oncoplot <- data.frame(Hugo_Symbol=genes_to_plot,Reason="Selected Genes")
      }
    } else if (class(genes_to_plot)=="data.frame") {
      genes_for_oncoplot <- genes_to_plot
      if (! "Reason" %in% colnames(genes_for_oncoplot)) {
        genes_for_oncoplot$Reason <- "Selected Genes"
      }
    } else {
      stop(paste0("Don't know what to do with 'genes_to_plot' of class: ",class(genes_to_plot)))
    }
    
    genes_for_oncoplot <- genes_for_oncoplot[,c("Hugo_Symbol","Reason")]
    genes_for_oncoplot <- data.frame(apply(genes_for_oncoplot,2,as.character), stringsAsFactors = F)
    genes_for_oncoplot <- genes_for_oncoplot[genes_for_oncoplot$Hugo_Symbol %in% maf.filtered@gene.summary$Hugo_Symbol, ]
  } else {
    genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), Reason=c())
  }
  
  frac_mut <- data.frame(Hugo_Symbol=maf.filtered@gene.summary$Hugo_Symbol,
                         frac_mut=(maf.filtered@gene.summary$MutatedSamples/as.numeric(maf.filtered@summary$summary[3])),
                         stringsAsFactors = F)
  if (! is.null(cohort_freq_thresh)) {
    ### Structure info about the fraction of the cohort that has each gene mutated
    # browser()
    # target_frac = sort(frac_mut$frac_mut, decreasing = T)[min(50,nrow(frac_mut))]
    # cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
    ### Select genes based on the frequency threshold
    freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut > cohort_freq_thresh]
    if (length(freq_genes) == 0) {
      stop("No genes to plot; change the frequency threshold to include more genes.")
    }
    if (length(freq_genes) > 200) {
      # target_frac = round(sort(frac_mut$frac_mut, decreasing = T)[min(50,nrow(frac_mut))],2)
      ngene_max=25
      # auto_adjust_threshold=T
      target_frac = sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))]
      if (auto_adjust_threshold) {
        cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
      }
      warning(paste0("Too many genes for oncoplot. Setting the cohort mutated fraction to > ", target_frac))
      freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut >= cohort_freq_thresh]
      freq_genes <- freq_genes[1:ngene_max]
      # return(NA)
    }
    gene_list <- list(freq_genes)
    reasons <- paste0("Cohort Freq > ",cohort_freq_thresh)
    
    ### Collect genes to plot
    
    for (i in 1:length(gene_list)) {
      if (is.na(gene_list[[i]][1])) {
        next
      }
      genes_for_oncoplot <- rbind(genes_for_oncoplot,
                                  data.frame(Hugo_Symbol=gene_list[[i]],
                                             Reason=reasons[i]))
    }
    genes_for_oncoplot <- cbind(genes_for_oncoplot,
                                frac=frac_mut$frac_mut[match(genes_for_oncoplot$Hugo_Symbol, frac_mut$Hugo_Symbol)])
    
    genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$Reason, -genes_for_oncoplot$frac),]
  }
  # browser()
  ### Split the oncoplot based on the reason for picking the gene
  ###   Here, we're only picked based on the frequency
  ###   But this framework is useful for plotting genes picked using various criteria
  split_idx=factor(genes_for_oncoplot$Reason)
  split_colors <- rainbow(length(levels(split_idx)))
  names(split_colors) <- as.character(genes_for_oncoplot$Reason[!duplicated(genes_for_oncoplot$Reason)])
  split_colors <- list(Reason=split_colors)
  
  # source("scripts/helper_functions.oncoplot.R")
  ### Make matrix to plot, and order it correctly
  if (use_clinvar_anno) {
    oncomat <- createOncoMatrix_CLINSIG(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = include_all)$oncoMatrix
  } else {
    oncomat <- createOncoMatrix(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = include_all)$oncoMatrix
  }
  oncomat <- oncomat[match(genes_for_oncoplot$Hugo_Symbol,rownames(oncomat)), , drop=F]
  onco_genes <- rownames(oncomat)
  
  # browser()
  if (include_all) {
    ### createOncoMatrix drops empty samples, so this adds them back in
    # all_wes_samples <- as.character(sample_info.exome$Tumor_Sample_Barcode[!is.na(sample_info.exome$Tumor_Sample_Barcode)])
    all_wes_samples <- levels(maf.filtered@variants.per.sample$Tumor_Sample_Barcode)
    extra_samples <- setdiff(all_wes_samples, colnames(oncomat) )
    print(paste0("Adding back ", length(extra_samples), " samples with no reported mutations..."))
    empty_data <- matrix(data = "", nrow=nrow(oncomat), ncol=length(extra_samples), dimnames=list(rownames(oncomat), extra_samples))
    oncomat <- cbind(oncomat, empty_data)
  }
  
  if (!is.null(custom_column_order)) {
    custom_order <- match(custom_column_order, colnames(oncomat), nomatch=0)
    oncomat <- oncomat[,custom_order]
  }
  oncomat.plot <- oncomat
  
  ### Set the height of the plot based on number of genes
  onco_height=NULL
  if (is.null(onco_height)) {
    onco_height=max(round(0.2*nrow(oncomat.plot),0),6)
  }
  
  ### Make the mutation type names prettier by removing the underscore
  # my_mut_col <- mutation_colors
  # names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
  oncomat.plot <- gsub("_"," ",oncomat.plot)
  # browser()
  if (ncol(oncomat.plot) < 1) {
    stop("No samples to plot.")
  }
  ### Column labels get cluttered if too many samples
  if (!is.logical(show_sample_names)) {
    show_sample_names=T
    if (ncol(oncomat.plot) > 20) {
      show_sample_names=F
    }
  }
  
  # if (show_burden) {
  variant_type_data <- data.frame(maf.filtered@variant.classification.summary)
  unmut_samples <- setdiff(levels(variant_type_data$Tumor_Sample_Barcode),as.character(variant_type_data$Tumor_Sample_Barcode))
  
  if (length(unmut_samples) > 0) {
    unmut_data <- data.frame(unmut_samples, matrix(0, nrow=length(unmut_samples), ncol=ncol(variant_type_data)-1))
    colnames(unmut_data) <- colnames(variant_type_data)
    variant_type_data <- rbind(variant_type_data,unmut_data)
  }
  
  rownames(variant_type_data) <- variant_type_data$Tumor_Sample_Barcode
  colnames(variant_type_data) <- gsub("_"," ",colnames(variant_type_data))
  variant_type_data <- variant_type_data[,c(-1,-ncol(variant_type_data)), drop=F]
  variant_type_data <- variant_type_data[match(colnames(oncomat.plot), rownames(variant_type_data)),
                                         rev(order(colSums(variant_type_data))), drop=F]
  # browser()
  var_anno_colors <- mutation_colors[match(colnames(variant_type_data), names(mutation_colors))]
  ha = HeatmapAnnotation(Burden = anno_barplot(variant_type_data, gp = gpar(fill = var_anno_colors), border = F))
  # }
  
  pct_anno <- paste0(prettyNum(frac_mut$frac_mut[match(onco_genes, frac_mut$Hugo_Symbol)]*100,digits=1),"%")
  left_ha = rowAnnotation("Cohort Pct"=anno_text(pct_anno,gp = gpar(cex=0.7)), show_annotation_name=F)
  
  myanno=NULL
  if (!is.null(clin_data)) {
    # browser()
    anno_data <- data.frame(clin_data[match(colnames(oncomat.plot), clin_data$Tumor_Sample_Barcode, nomatch=0),],stringsAsFactors = F)
    extra_ids <-setdiff(colnames(oncomat.plot), rownames(anno_data))
    emptydat <- data.frame(matrix(nrow=length(extra_ids), ncol=ncol(anno_data)), stringsAsFactors = F)
    colnames(emptydat) <- colnames(anno_data)
    emptydat[,"Tumor_Sample_Barcode"] <- extra_ids
    anno_data <- rbind(anno_data,
                       emptydat)
    
    row.names(anno_data) <- anno_data$Tumor_Sample_Barcode
    anno_data <- anno_data[,!colnames(anno_data) %in% "Tumor_Sample_Barcode", drop=F]
    if (ncol(anno_data) > 0) {
      myanno <- HeatmapAnnotation(df=anno_data,col = clin_data_colors)
    }
  }
  
  ### Make the oncoplot
  # browser()
  
  onco_base_default <- oncoPrint(oncomat.plot, alter_fun = alter_fun, col=mutation_colors, row_order=1:nrow(oncomat.plot),
                                 name="oncoplot",
                                 row_title=title_text,
                                 show_pct = F,
                                 row_split=split_idx,
                                 bottom_annotation = myanno,
                                 left_annotation = rowAnnotation(Reason = split_idx, col=split_colors, annotation_width = unit(0.3, "mm"), show_annotation_name=F),
                                 right_annotation = left_ha,
                                 top_annotation = ha,
                                 column_order=custom_column_order,
                                 show_column_names = show_sample_names)#,
  # column_names_rot = 30,
  # column_gap = unit(0.0001,"npc"),
  # width = unit(0.75, "npc"))
  # print(dim(oncomat))
  ### Save the oncoplot
  if (full_output) {
    return_val <- list(onco_base_default, oncomat.plot)
  } else {
    return_val <- onco_base_default
  }
  
  if ( ! is.null(savename) ) {
    # save_name <- paste0(out_dir,"/oncoplot.",cohort_freq_thresh,".pdf")
    #onco_height=max(round(0.15*nrow(oncomat.plot),0),6)
    #onco_width=max(c(onco_height*0.75,max(round(0.4*ncol(oncomat.plot),0))))
    
    png(file = savename,height=plot_height,width=plot_width,units="px",pointsize=pointsize)
    draw(onco_base_default)
    dev.off()
  }
  
  invisible(return_val)
}

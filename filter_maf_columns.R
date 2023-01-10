###Copy the MAFDASH output similar to make sure that we have the same fields
generate_maf_fields_annovar=function(input_maf_dt,germline=F){
  
  if(germline==T){
    ##NORMAL NOT ANNOTATED WITH ANNOVAR YET
    input_maf_dt$tumor_genotype=apply(input_maf_dt[, c("Tumor_Seq_Allele1", 
                                                       "Tumor_Seq_Allele2")], 1, paste, collapse = "/")
    input_maf_dt$tumor_freq <- as.numeric(as.character(input_maf_dt$t_alt_count))/as.numeric(as.character(input_maf_dt$t_depth))
    
    
    cols_for_table <- c(`Hugo Symbol` = "Hugo_Symbol", `Sample ID` = "Tumor_Sample_Barcode", 
                        `Variant Classification` = "Variant_Classification", 
                        `Variant Type` = "Variant_Type", 
                        Consequence = "Consequence", 
                        Chromosome = "Chromosome", `Start Position` = "Start_Position", 
                        `End Position` = "End_Position", #Strand = "Strand", 
                        `Reference Allele` = "Reference_Allele", 
                        `Genotype` = 'tumor_genotype', 
                        `Known Effects ClinVar` = "CLIN_SIG", `Transcript Change` = "HGVSc", 
                        `Protein Change` = "HGVSp_Short", 
                        #`Normal Depth` = "n_depth", 
                        #`Normal Ref Depth` = "n_ref_count", `Normal Alt Depth` = "n_alt_count", 
                        `Ref Depth` = "t_depth", `Normal Ref Depth` = "t_ref_count", 
                        `Alt Depth` = "t_alt_count", `Alt Frequency` = "tumor_freq", 
                        #`Existing Annotation` = "Existing_variation", 
                        `gnomAD Frequency` = "AF", 
                        #`ExAC Frequency` = "ExAC_ALL", 
                        `1000Genomes Frequency` = "AF", 
                        `Effect Prediction - SIFT` = "SIFT", `Effect Prediction - PolyPhen` = "PolyPhen")
    
  }else{
    input_maf_dt$tumor_genotype=apply(input_maf_dt[, c("Tumor_Seq_Allele1", 
                                           "Tumor_Seq_Allele2")], 1, paste, collapse = "/")
    input_maf_dt$normal_genotype <- apply(input_maf_dt[, c("Match_Norm_Seq_Allele1", 
                                               "Match_Norm_Seq_Allele2")], 1, paste, collapse = "/")
    input_maf_dt$tumor_freq <- as.numeric(as.character(input_maf_dt$t_alt_count))/as.numeric(as.character(input_maf_dt$t_depth))
    
    
    cols_for_table <- c(`Hugo Symbol` = "Hugo_Symbol", `Sample ID` = "Tumor_Sample_Barcode", 
                        `Variant Classification` = "Variant_Classification", 
                        `Variant Type` = "Variant_Type", 
                        Consequence = "Consequence", 
                        Chromosome = "Chromosome", `Start Position` = "Start_Position", 
                        `End Position` = "End_Position", Strand = "Strand", `Reference Allele` = "Reference_Allele", 
                        `Tumor Genotype` = "tumor_genotype", `Normal Genotype` = "normal_genotype", 
                        `Known Effects ClinVar` = "CLIN_SIG", `Transcript Change` = "txChange", 
                        `Protein Change` = "HGVSp_Short", `Normal Depth` = "n_depth", 
                        `Normal Ref Depth` = "n_ref_count", `Normal Alt Depth` = "n_alt_count", 
                        `Tumor Depth` = "t_depth", `Tumor Ref Depth` = "t_ref_count", 
                        `Tumor Alt Depth` = "t_alt_count", `Tumor Alt Frequency` = "tumor_freq", 
                        #`Existing Annotation` = "Existing_variation", 
                        `gnomAD Frequency` = "gnomAD_AF", 
                        #`ExAC Frequency` = "ExAC_ALL", 
                        `1000Genomes Frequency` = "AF", 
                        `Effect Prediction - SIFT` = "SIFT", `Effect Prediction - PolyPhen` = "PolyPhen")
  }
  ##crate gene alterd infrac
  
  n_samples <- length(unique(input_maf_dt$Tumor_Sample_Barcode))
  gene_frac <- input_maf_dt %>% 
    dplyr::group_by(Hugo_Symbol) %>% 
    dplyr::mutate(`Gene Altered in Cohort frac`=length(unique(Tumor_Sample_Barcode))/n_samples) %>%
    select(`Gene Altered in Cohort frac`) %>%
    distinct()
  var_frac <- input_maf_dt %>% 
    dplyr::group_by(Hugo_Symbol, HGVSp_Short) %>%
    dplyr::mutate(`Variant in Cohort frac`=length(unique(Tumor_Sample_Barcode))/n_samples)%>%
    select(`Variant in Cohort frac`) %>%
    distinct()
  
  input_maf_dt %>% dplyr::select(all_of(cols_for_table)) %>%
    left_join(.,gene_frac,by=c(`Hugo Symbol`='Hugo_Symbol')) %>%
    left_join(.,var_frac,by=c(`Hugo Symbol`="Hugo_Symbol",`Protein Change`="HGVSp_Short")) %>%
    select(`Hugo Symbol`,`Sample ID`, `Gene Altered in Cohort frac`,`Variant in Cohort frac`,everything())
}


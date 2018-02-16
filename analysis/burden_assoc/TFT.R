# TFT function for testing associations

TFT = function(data){
  
  p = NA; OR = NA
  
  fisher_elements = as.numeric(data)
  
  if (fisher_elements[1] > 0 && fisher_elements[3] > 0 && fisher_elements[2] >= 0 && fisher_elements[4] >= 0){
    test.table = matrix(as.numeric(fisher_elements), nrow=2)
    f.test = fisher.test(test.table, alternative = "greater")
    OR = f.test$estimate
    p = f.test$p.value
  }
  result_row = c(fisher_elements,p,OR)
  
  return(result_row)
  
}

run_TFT = function(data, AF_thres = 0.01){
  data=data[data$Freq > AF_thres,]
  num_genes = length(unique(data$Gene))
  # burden test: TFT: http://slideplayer.com/slide/8660600/
  stats = matrix(,nrow=num_genes,ncol=7)
  
  for (i in 1:num_genes){
    gene = unique(data$Gene)[i]
    #data_g = data[data$gene_symbol==gene,]
    data_g = c(53105,data$ExAC_Count[data$Gene == gene],data$Sample_size[1],data$Count[data$Gene == gene])
    gene_stat = TFT(data_g)
    stats[i,] = c(as.character(gene), gene_stat)
  }
  colnames(stats) = c("Gene","ExAC_nonTCGA_AN","ExAC_nonTCGA_AC","cohort_AN","cohort_AC","P","OR")
  stats = data.frame(stats,stringsAsFactors = F)
  stats[,2:7] = sapply(stats[,2:7],as.numeric,2)
  
  #stats$P = as.numeric(stats$P)
  #stats$FDR = p.adjust(stats$P, method="BH")
  #stats = stats[order(stats$P),]
  return(stats)
}

run_TFT_against_others = function(data_c, PCA_count_byCancer, all_cancer_stat_m_suggest, AF_thres = 0.01){
  data_c=data_c[data_c$Freq > AF_thres,]
  cancer = data_c$Cancer[1]
  num_genes = length(unique(data_c$Gene))
  # burden test: TFT: http://slideplayer.com/slide/8660600/
  stats = matrix(,nrow=num_genes,ncol=7)
  
  for (i in 1:num_genes){
    gene = unique(data_c$Gene)[i]
    
    sig_cancers = all_cancer_stat_m_suggest$Cancer[all_cancer_stat_m_suggest$Gene==gene]
    data_other_c = PCA_count_byCancer[!(PCA_count_byCancer$Cancer %in% c(sig_cancers,as.character(cancer))),]
    
    data_other_c_g = data_other_c[data_other_c$Gene == gene,]
    #data_c_g = data_c[data_c$gene_symbol==gene,]
    data_c_g = c(sum(data_other_c_g$Sample_size),sum(data_other_c_g$Count),data_c$Sample_size[i],data_c$Count[data_c$Gene == gene])
    gene_stat = TFT(data_c_g)
    stats[i,] = c(as.character(gene), gene_stat)
  }
  colnames(stats) = c("Gene","other_cancers_AN","other_cancers_AC","cohort_AN","cohort_AC","P","OR")
  stats = data.frame(stats,stringsAsFactors = F)
  stats[,2:7] = sapply(stats[,2:7],as.numeric,2)
  
  #stats$P = as.numeric(stats$P)
  #stats$FDR = p.adjust(stats$P, method="BH")
  #stats = stats[order(stats$P),]
  return(stats)
}

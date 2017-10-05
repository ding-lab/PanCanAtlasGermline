##### burden test functions #####
TFT = function(data){
  ExAC_nonTCGA_AN = max(data$ExAC_nonTCGA_AN, na.rm=T)
  TCGA_AN = max(data$TCGA_AN, na.rm=T)
  ExAC_nonTCGA_AC = sum(data$ExAC_nonTCGA_AC, na.rm=T)
  TCGA_AC = sum(data$TCGA_AC, na.rm=T)
  
  p = NA; OR = NA
  
  fisher_elements = as.numeric(c(ExAC_nonTCGA_AN,ExAC_nonTCGA_AC,TCGA_AN,TCGA_AC))
  
  if (fisher_elements[1] > 0 && fisher_elements[3] > 0 && fisher_elements[2] >= 0 && fisher_elements[4] >= 0){
    test.table = matrix(as.numeric(fisher_elements), nrow=2)
    f.test = fisher.test(test.table)
    OR = f.test$estimate
    p = f.test$p.value
  }
  result_row = c(fisher_elements,p,OR)
  
  return(result_row)
}

run_TFT = function(data, AF_thres = 0.01){
  # some clean-ups
  data = data[data$ExAC_AC/data$ExAC_AN < AF_thres,]
  data = data[data$ExAC_AN > 114600,] # require enough samples with observed data; 114600 = first quantile
  num_genes = length(unique(data$gene_symbol))
  # burden test: TFT: http://slideplayer.com/slide/8660600/
  stats = matrix(,nrow=num_genes,ncol=7)
  
  for (i in 1:num_genes){
    gene = unique(data$gene_symbol)[i]
    data_g = data[data$gene_symbol==gene,]
    gene_stat = TFT(data_g)
    stats[i,] = c(gene, gene_stat)
  }
  colnames(stats) = c("gene","nonTCGA_AN","total_nonTCGA_AC","TCGA_AN","total_TCGA_AC","P","OR")
  stats = data.frame(stats,stringsAsFactors = F)
  stats[,2:7] = sapply(stats[,2:7],as.numeric,2)
  
  #stats$P = as.numeric(stats$P)
  stats$FDR = p.adjust(stats$P, method="BH")
  stats = stats[order(stats$P),]
  return(stats)
}

plot_burden_result = function(stats){
  p = ggplot(data=stats, aes(x=total_nonTCGA_AC,y=total_TCGA_AC, color=OR))
  p = p + geom_point(aes(size=-log(FDR)),alpha=0.5)
  p = p + geom_text(aes(label=ifelse(FDR<0.05,gene, NA)))
  p = p + labs(x = "Total nonTCGA variant counts", y = "Total TCGA variant counts") + theme_bw()
  p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14),
                axis.text.y = element_text(colour="black", size=14))
  p = p + geom_abline(slope=15202/106210, alpha=0.8)
  return(p)
}

run_plot_burden = function(data,AF_thres=0.01){
  data_name = deparse(substitute(data))
  data_stats = run_TFT(data,AF_thres = AF_thres)
  # write result
  tn = paste("out/burden",data_name, "AF", AF_thres, "TFT_stats.tsv", sep="_")
  write.table(data_stats, file=tn, quote=F, sep = '\t', col.names=NA)
  # plot result
  plot_burden_result(data_stats)
  fn = paste("out/burden",data_name, "AF", AF_thres, "point.pdf", sep="_")
  ggsave(file=fn, h=6,w=6, useDingbats=FALSE)
}
##### somaticDriver_germline_fisher.R #####
# Kuan-lin Huang @ WashU 2016 June
# updated 2018
# Conduct fisher's exact test to find mutual occurence/mutual exclusivity of germline/somatic variants
# somatic driver only

### dependencies ###
#setwd("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/somatic_germline")
source("../global_aes_out.R")
source("../dependency_files.R")
source("../load_somatic.R")

### function ###
sg_fisher_test = function(all_samples, germline_carriers, somatic_carriers){
  p = NA; OR = NA
  
  fisher_elements = c(sum(!all_samples %in% c(germline_carriers,somatic_carriers)),sum((all_samples %in% germline_carriers) & !(all_samples %in% somatic_carriers)),
                      sum(!(all_samples %in% germline_carriers) & (all_samples %in% somatic_carriers)),sum((all_samples %in% germline_carriers) & (all_samples %in% somatic_carriers)))
  if (fisher_elements[1] > 0 && fisher_elements[3] > 0 && fisher_elements[2] >= 0 && fisher_elements[4] >= 0){
    test.table = matrix(as.numeric(fisher_elements), nrow=2)
    f.test = fisher.test(test.table)
    OR = f.test$estimate
    p = f.test$p.value
    
    count00 = test.table[1,1]
    count10 = test.table[2,1]
    count01 = test.table[1,2]
    count11 = test.table[2,2]
    return(list("p"=p, "OR"=OR, "count00" = count00, "count00" = count00
                , "count10" = count10, "count01" = count01, "count11" = count11))
  }
}

### get input date and files ###

### germline ###
germ_list = names(table(pathVarP$HUGO_Symbol)[table(pathVarP$HUGO_Symbol)>9]) # limit analyses to germline genes with equal to or more than 10 variants

### germline somatic interaction ###
# samples being tested
test_samples = intersect(unique(somatic$bcr_patient_barcode), clin$bcr_patient_barcode)
cat("Running germline-somatic interaction in ",length(test_samples),"samples.\n")

# run through function
out_table=character(0);
for (g_gene in germ_list){
  for (s_gene in somaticDriver299){
    germline_carriers = unique(pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol==g_gene])
    somatic_carriers = unique(somatic_likelyfunctional_driver$bcr_patient_barcode[somatic_likelyfunctional_driver$Hugo_Symbol==s_gene])
    if(length(germline_carriers) > 9 & length(somatic_carriers) > 19){
      t_result = sg_fisher_test(all_samples = test_samples, germline_carriers, somatic_carriers)
      p = t_result$p
      OR = t_result$OR
      count00 = t_result$count00
      count10 = t_result$count10
      count01 = t_result$count01
      count11 = t_result$count11
    
    
      out_row = c(g_gene, s_gene, count00, count10, count01, count11, OR, p)
      out_table = rbind(out_table,out_row)
    }
  }
}
rownames(out_table)=NULL
colnames(out_table) = c("GermlineGene", "SomaticGene", "germline0somatic0", "germline1somatic0", "germline0somatic1", "germline1somatic1","OR", "P")
out_table = data.frame(out_table)
out_table[,"P"] = as.numeric(as.character(out_table[,"P"]))

FDR = p.adjust(out_table[,"P"],method="fdr")
out_table = cbind(out_table,FDR)
out_table = out_table[order(out_table[,"P"]),]

fn = paste("out/germline_somatic_driver_fisher.tsv")
write.table(out_table,file=fn,quote=FALSE,row.names=FALSE,sep="\t")
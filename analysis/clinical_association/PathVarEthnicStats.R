##### plotPathVarEthnic.R #####
# Kuan-lin Huang @ WashU 201711
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/clinical_association"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")
source("plotPathVarEthnic.R")

out_table=NULL

for (cancer in unique(count.F_25$Cancer)){
  count_c = count.F_25[count.F_25$Cancer==cancer,]
  ethnicities = unique(count_c$Ethnicity)
  caucasian_all = count_c[!is.na(count_c$Ethnicity) & count_c$Ethnicity=="Europian","CohortCount"][1]
  
  if ("Asian" %in% ethnicities){
    asian_all = count_c[!is.na(count_c$Ethnicity) & count_c$Ethnicity=="Asian","CohortCount"][1]
    for (gene in unique(count_c$Gene)){
      count_c_g = count_c[count_c$Gene==gene,]
      caucasian_carrier = count_c_g[!is.na(count_c_g$Ethnicity) & count_c_g$Ethnicity=="Europian","CarrierCount"]
      asian_carrier = count_c_g[!is.na(count_c_g$Ethnicity) & count_c_g$Ethnicity=="Asian","CarrierCount"]
      if (!length(caucasian_carrier)){caucasian_carrier=0}
      if (!length(asian_carrier)){asian_carrier=0}
      dist_nums = c(caucasian_all-caucasian_carrier,caucasian_carrier,asian_all-asian_carrier,asian_carrier)
      test.t = matrix(dist_nums,nrow=2)
      f.test = fisher.test(test.t)
      OR = f.test$estimate
      p = f.test$p.value
      
      row_stat = cbind(cancer, gene, "Europian", "Asian", caucasian_all-caucasian_carrier,caucasian_carrier,asian_all-asian_carrier,asian_carrier, OR, p)
      out_table = rbind(out_table, row_stat)
    }
  }
  if ("African American" %in% ethnicities){
    asian_all = count_c[!is.na(count_c$Ethnicity) & count_c$Ethnicity=="African American","CohortCount"][1]
    for (gene in unique(count_c$Gene)){
      count_c_g = count_c[count_c$Gene==gene,]
      caucasian_carrier = count_c_g[!is.na(count_c_g$Ethnicity) & count_c_g$Ethnicity=="Europian","CarrierCount"]
      asian_carrier = count_c_g[!is.na(count_c_g$Ethnicity) & count_c_g$Ethnicity=="African American","CarrierCount"]
      if (!length(caucasian_carrier)){caucasian_carrier=0}
      if (!length(asian_carrier)){asian_carrier=0}
      dist_nums = c(caucasian_all-caucasian_carrier,caucasian_carrier,asian_all-asian_carrier,asian_carrier)
      test.t = matrix(dist_nums,nrow=2)
      f.test = fisher.test(test.t)
      OR = f.test$estimate
      p = f.test$p.value
      
      row_stat = cbind(cancer, gene, "Europian", "African American", caucasian_all-caucasian_carrier,caucasian_carrier,asian_all-asian_carrier,asian_carrier, OR, p)
      out_table = rbind(out_table, row_stat)
    }
  }
  if ("American" %in% ethnicities){
    asian_all = count_c[!is.na(count_c$Ethnicity) & count_c$Ethnicity=="American","CohortCount"][1]
    for (gene in unique(count_c$Gene)){
      count_c_g = count_c[count_c$Gene==gene,]
      caucasian_carrier = count_c_g[!is.na(count_c_g$Ethnicity) & count_c_g$Ethnicity=="Europian","CarrierCount"]
      asian_carrier = count_c_g[!is.na(count_c_g$Ethnicity) & count_c_g$Ethnicity=="]American","CarrierCount"]
      if (!length(caucasian_carrier)){caucasian_carrier=0}
      if (!length(asian_carrier)){asian_carrier=0}
      dist_nums = c(caucasian_all-caucasian_carrier,caucasian_carrier,asian_all-asian_carrier,asian_carrier)
      test.t = matrix(dist_nums,nrow=2)
      f.test = fisher.test(test.t)
      OR = f.test$estimate
      p = f.test$p.value
      
      row_stat = cbind(cancer, gene, "Europian", "American", caucasian_all-caucasian_carrier,caucasian_carrier,asian_all-asian_carrier,asian_carrier, OR, p)
      out_table = rbind(out_table, row_stat)
    }
  }
}
out_table = as.data.frame(out_table)
colnames(out_table) = c("cancer","gene","EthnicityA","EthnicityB","EthnicityA_noncarrier","EthnicityA_carrier"
                        ,"EthnicityB_noncarrier","EthnicityB_carrier","OR","P")
out_table = out_table[as.numeric(as.character(out_table$EthnicityB_carrier))>1,]
out_table$FDR = p.adjust(out_table[,"P"], method="fdr") # MAW new, calculates FDR based on the method from,
out_table=out_table[order(out_table$P, decreasing=FALSE),]
tn = "out/pathVar_2ethni_TFT_assoc.txt"
write.table(out_table, quote=F, sep="\t", file = tn, row.names = F)

##### compile_compare_samples.R #####
# Kuan-lin Huang @ WashU 2017 Oct.
# make a clinical supplementary table for pan-germline manuscript

setwd("/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/sample_listing")
source("../global_aes_out.R")

clin_f = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/PanCan_ClinicalData_V4_20170428.txt"
clin_full = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)
clin = clin_full[,c("bcr_patient_barcode", "type","age_at_initial_pathologic_diagnosis","gender","race")]
colnames(clin) = c("sample","cancer","age_at_onset","gender","ethnicity")

clin$ethnicity[clin$ethnicity %in% c("","[Not Available]","[Not Evaluated]","[Unknown]")]=NA
sample_cancer_clin$gender[sample_cancer_clin$gender==""]=NA

s_c_list_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/sampleQC/pca_table.20171019.wclin.tsv"
sample_cancer = read.table(header=T, quote = "", sep="\t", file = s_c_list_f, stringsAsFactors=FALSE)
sample_cancer = sample_cancer[,c("bcr_patient_barcode", "cancer")]
colnames(sample_cancer) = c("sample","cancer")

sample_cancer_clin = merge(sample_cancer,clin, by = c("sample","cancer"), all.x =T)
mean(sample_cancer_clin$age_at_onset, na.rm=T)

cancer_count = data.frame(table(data = sample_cancer_clin$cancer ))
cancer_ethni_count = as.data.frame(table(sample_cancer_clin$cancer,sample_cancer_clin$ethnicity))
cancer_ethni_count_d = dcast(cancer_ethni_count, formula = Var1 ~ Var2, value.var = "Freq")
cancer_gender_count = as.data.frame(table(sample_cancer_clin$cancer,sample_cancer_clin$gender))
cancer_gender_count_d = dcast(cancer_gender_count, formula = Var1 ~ Var2, value.var = "Freq")
cancer_gender_count_d$FemalePercent=cancer_gender_count_d$FEMALE/(cancer_gender_count_d$FEMALE+ cancer_gender_count_d$MALE)
cancer_aao = data.frame(aggregate(data = sample_cancer_clin, age_at_onset~cancer, FUN = "mean"))
cancer_aao_sd = data.frame(aggregate(data = sample_cancer_clin, age_at_onset~cancer, FUN = "sd"))
colnames(cancer_count) = c("Cancer","Sample size")
colnames(cancer_ethni_count_d)[1] = c("Cancer")
colnames(cancer_aao) = c("Cancer", "Average_AAO")
colnames(cancer_aao_sd) = c("Cancer", "AAO_SD")
colnames(cancer_gender_count_d)[1]="Cancer"
cancer_gender_count_d_sum = cancer_gender_count_d[,c("Cancer","FemalePercent")]

cancer_count_gender = merge(cancer_count,cancer_gender_count_d_sum,by="Cancer")
cancer_count_wethni = merge(cancer_count_gender,cancer_ethni_count_d, by="Cancer") 
cancer_count_wethni_waao = merge(cancer_count_wethni, cancer_aao, by = "Cancer")
cancer_count_wethni_waao_sd = merge(cancer_count_wethni_waao, cancer_aao_sd, by = "Cancer")
cancer_count_wethni_waao_sd$AAO = paste(round(cancer_count_wethni_waao_sd$Average_AAO,1), "+/-", round(cancer_count_wethni_waao_sd$AAO_SD,1))
colnames(cancer_count_wethni_waao_sd) = tolower(colnames(cancer_count_wethni_waao_sd))
colnames(cancer_count_wethni_waao_sd) = paste(toupper(substring(colnames(cancer_count_wethni_waao_sd), 1,1)), substring(colnames(cancer_count_wethni_waao_sd), 2),sep="")
all_sum = sum(as.numeric(cancer_count_wethni_waao_sd[,2]))
all_sum_ethni = sapply(cancer_count_wethni_waao_sd[,c(4:8)],sum)
gender_ratio = sum(sample_cancer_clin$gender=="FEMALE",na.rm=T)/(sum(sample_cancer_clin$gender=="FEMALE",na.rm=T)+sum(sample_cancer_clin$gender=="MALE",na.rm=T))
all_aao = round(mean(sample_cancer_clin$age_at_onset, na.rm=T),1)
all_aao_sd = round(sd(sample_cancer_clin$age_at_onset, na.rm=T),1)
all_row = c("All", all_sum, gender_ratio, all_sum_ethni,all_aao, all_aao_sd, paste(all_aao, "+/-", all_aao_sd))
cancer_count_wethni_waao_sd$Cancer = as.character(cancer_count_wethni_waao_sd$Cancer)
cancer_count_wethni_waao_sd = rbind(cancer_count_wethni_waao_sd,all_row)

cancer_count_wethni_waao_sd_p = cancer_count_wethni_waao_sd[,!(colnames(cancer_count_wethni_waao_sd) %in% c("Average_aao","Aao_sd"))]
#colnames(cancer_count_wethni_waao_sd_p) = c("Cancer", "Sample size", "American indian", "Asian", "African american", "Pacific islander", "White", "Age at onset")

tn = "out/cancer_count_wethni_waao.txt"
write.table(cancer_count_wethni_waao_sd_p, quote=F, sep="\t", file = tn, row.names = F)

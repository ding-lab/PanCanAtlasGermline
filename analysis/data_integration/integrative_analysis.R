##### integrative_analysis.R #####
# Kuan-lin Huang @ WashU 2016-2017
# basic stats on germline variants in the paper

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/functional_clinical_integration"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")
library(plyr)
library(UpSetR)
library(ggrepel)

# pathogenic only
pathVar_TSG = pathVar[pathVar$Gene_Classification=="TSG",]
pathVar_onco = pathVar[pathVar$Gene_Classification == "Oncogene",]

pathVar_onco$exp_q.cut = cut(pathVar_onco$expression_quantile, breaks=seq(0,1,0.25))
df = as.data.frame(with(pathVar_onco, table(exp_q.cut,Final_Classification2)))
# next: compute percentages per group
df = ddply(df, .(Final_Classification2), transform, p = Freq/sum(Freq))
tn = paste(td, "pan8000_MAF0.05_germline_onco_exp.txt",sep="_")
write.table(df, quote=F, sep="\t", file = tn, row.names = F)

cat("Pathogenic defined as any variants that are high-confidence pathogenic, pathogenic, and likely pathogenic\n")

cat("Pathogenic oncogene variants\n")
pathVar_onco = pathVar[pathVar$Gene_Classification=="Oncogene",]
cat("Percentage top 25% expression\n")
sum(pathVar_onco$expression_quantile > 0.75, na.rm=T)/sum(!is.na(pathVar_onco$expression_quantile))
cat("Percentage Significant LOH\n")
sum(as.character(pathVar_onco$LOH_Sig)=="Significant")/length(pathVar_onco$LOH_Sig)
cat("Number of unique oncogenes:",length(unique(pathVar_onco$gene_name)), "\n")
table(pathVar_onco$gene_name)[order(table(pathVar_onco$gene_name),decreasing = T)]

cat("Pathogenic TSG variants\n")
pathVar_tsg = pathVar[pathVar$Gene_Classification=="TSG",]
cat("Percentage bottom 25% expression\n")
sum(pathVar_tsg$expression_quantile < 0.25, na.rm=T)/sum(!is.na(pathVar_tsg$expression_quantile))
cat("Percentage Significant LOH\n")
sum(as.character(pathVar_tsg$LOH_Sig)=="Significant")/length(pathVar_tsg$LOH_Sig)
cat("Number of unique TSGs:",length(unique(pathVar_tsg$gene_name)), "\n")

cat("Expression quantile difference\n")
ks.test(pathVar_onco$expression_quantile,pathVar_tsg$expression_quantile)

# overall fraction showing association with cancer
num_total = nrow(pathVar)
num_assoc = nrow(pathVar[!is.na(pathVar$ExAC_assoc_P) & pathVar$ExAC_assoc_P < 0.05,])
cat("Num pathogenic association:",num_assoc, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing association:",num_assoc/num_total,"\n")
cat("\n")

cat("Oncogenes:\n")
num_total = nrow(pathVar_onco)
num_assoc = nrow(pathVar_onco[!is.na(pathVar_onco$ExAC_assoc_P) & pathVar_onco$ExAC_assoc_P < 0.05,])
cat("Num pathogenic association:",num_assoc, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing association:",num_assoc/num_total,"\n")
cat("\n")

cat("Tumor Suppressor genes:\n")
num_total = nrow(pathVar_tsg)
num_assoc = nrow(pathVar_tsg[!is.na(pathVar_tsg$ExAC_assoc_P) & pathVar_tsg$ExAC_assoc_P < 0.05,])
cat("Num pathogenic association:",num_assoc, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing association:",num_assoc/num_total,"\n")
cat("\n")


# overall fraction showing AI
num_total = nrow(pathVar)
num_AI = sum(pathVar$LOH_Sig =="Significant",na.rm=T)
cat("Num pathogenic AI:",num_AI, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing AI:",num_AI/num_total,"\n")
cat("\n")

cat("Oncogenes:\n")
num_total = nrow(pathVar_onco)
num_AI = sum(pathVar_onco$LOH_Sig =="Significant",na.rm=T)
cat("Num pathogenic AI:",num_AI, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing AI:",num_AI/num_total,"\n")
cat("\n")

cat("Tumor Suppressor genes:\n")
num_total = nrow(pathVar_tsg)
num_AI = sum(pathVar_tsg$LOH_Sig =="Significant",na.rm=T)
cat("Num pathogenic AI:",num_AI, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing AI:",num_AI/num_total,"\n")
cat("\n")

# overall fraction showing compound het
 # check if i need to lump in likely pathogenic here
num_total = nrow(pathVar)
num_CH = 34 #27 pathogenic variants
cat("Num pathogenic compound het:",num_CH, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing compound het:",num_CH/num_total,"\n")

cat("\n")

# overall fraction showing co-localizing
num_total = nrow(pathVar)
num_coloc = nrow(pathVar[pathVar$colocalized_somatic_mutation_count >=3 & pathVar$normal_VAF > 20,])
cat("Num pathogenic colocalizing:",num_coloc, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing co-localizing:",num_coloc/num_total,"\n")
cat("\n")
cat("Pathogenic variants gene distribution:")
table(pathVar[pathVar$colocalized_somatic_mutation_count >=3 & pathVar$normal_VAF > 20,]$HUGO_Symbol)
# cat("Likely pathogenic variants gene distribution:")
# table(pathVar[pathVar$colocalized_somatic_mutation_count >=3 & pathVar$normal_VAF > 20 & pathVar$Final_Classification2=="Likely Pathogenic",]$HUGO_Symbol)
cat("\n")

cat("Oncogenes:\n")
num_total = nrow(pathVar_onco)
num_coloc = nrow(pathVar_onco[pathVar_onco$colocalized_somatic_mutation_count >=3,])
cat("Num pathogenic colocalizing:",num_coloc, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing colocalizing:",num_coloc/num_total,"\n")
cat("\n")

cat("TSGs:\n")
num_total = nrow(pathVar_tsg)
num_coloc = nrow(pathVar_tsg[pathVar_tsg$colocalized_somatic_mutation_count >=3,])
cat("Num pathogenic colocalizing:",num_coloc, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing colocalizing:",num_coloc/num_total,"\n")
cat("\n")

# overall fraction showing clustering
num_total = nrow(pathVar)
num_clus = nrow(pathVar[pathVar$clustered,])
cat("Num pathogenic clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing clustering:",num_clus/num_total,"\n")
cat("\n")
cat("Pathogenic variants gene distribution:")
table(pathVar[pathVar$clustered,]$HUGO_Symbol)
# cat("Likely pathogenic variants gene distribution:")
# table(pathVar[pathVar$clustered & pathVar$normal_VAF > 20 & pathVar$Final_Classification2=="Likely Pathogenic",]$HUGO_Symbol)
cat("\n")

cat("Oncogenes:\n")
num_total = nrow(pathVar_onco)
num_clus = nrow(pathVar_onco[pathVar_onco$clustered,])
cat("Num pathogenic clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing clustering:",num_clus/num_total,"\n")
cat("\n")

cat("TSGs:\n")
num_total = nrow(pathVar_tsg)
num_clus = nrow(pathVar_tsg[pathVar_tsg$clustered,])
cat("Num pathogenic clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing clustering:",num_clus/num_total,"\n")
cat("\n")


# overall fraction showing co-localizing or clustering
num_total = nrow(pathVar)
num_clus = nrow(pathVar[(pathVar$clustered | (pathVar$colocalized_somatic_mutation_count >=3)),])
cat("Num pathogenic co-localizing or clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing co-localizing or clustering:",num_clus/num_total,"\n")
cat("\n")

cat("Oncogenes:\n")
num_total = nrow(pathVar_onco)
num_clus = nrow(pathVar_onco[(pathVar_onco$clustered | (pathVar_onco$colocalized_somatic_mutation_count >=3)),])
cat("Num pathogenic co-localizing or clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing co-localizing or clustering:",num_clus/num_total,"\n")
cat("\n")

cat("TSGs:\n")
num_total = nrow(pathVar_tsg)
num_clus = nrow(pathVar_tsg[(pathVar_tsg$clustered | (pathVar_tsg$colocalized_somatic_mutation_count >=3)),])
cat("Num pathogenic co-localizing or clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing co-localizing or clustering:",num_clus/num_total,"\n")
cat("\n")
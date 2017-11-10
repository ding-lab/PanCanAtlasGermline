##### pathVar_integrative_analysis.R #####
# Kuan-lin Huang @ WashU 201711
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/data_integration"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

# library(plyr)
# library(UpSetR)
# library(ggrepel)

source("../dependency_files.R")

# define subsets of data
pathVar_tsg = pathVar[!is.na(pathVar$Gene_Classification) & pathVar$Gene_Classification=="TSG",]
pathVar_onco = pathVar[!is.na(pathVar$Gene_Classification) & pathVar$Gene_Classification == "Oncogene",]

cat("Pathogenic defined as any variants that are high-confidence pathogenic, pathogenic, and likely pathogenic\n")

cat("Pathogenic oncogene variants\n")
cat("Percentage top 25% expression\n")
sum(pathVar_onco$expressionQuantile > 0.75, na.rm=T)/sum(!is.na(pathVar_onco$expressionQuantile))
cat("Percentage significant LOH\n")
sum(as.character(pathVar_onco$LOH_Sig)=="Significant")/length(pathVar_onco$LOH_Sig)
cat("Number of unique oncogenes:",length(unique(pathVar_onco$HUGO_Symbol)), "\n")
table(pathVar_onco$HUGO_Symbol)[order(table(pathVar_onco$HUGO_Symbol),decreasing = T)]

cat("Pathogenic TSG variants\n")
cat("Percentage bottom 25% expression\n")
sum(pathVar_tsg$expressionQuantile < 0.25, na.rm=T)/sum(!is.na(pathVar_tsg$expressionQuantile))
cat("Percentage significant LOH\n")
sum(as.character(pathVar_tsg$LOH_Sig)=="Significant")/length(pathVar_tsg$LOH_Sig)
cat("Number of unique TSGs:",length(unique(pathVar_tsg$HUGO_Symbol)), "\n")
table(pathVar_tsg$HUGO_Symbol)[order(table(pathVar_tsg$HUGO_Symbol),decreasing = T)]

cat("Expression quantile difference:\n")
ks.test(pathVar_onco$expressionQuantile,pathVar_tsg$expressionQuantile)

### frequency comparison plots###
p = ggplot(pathVar)
p = p + facet_grid(.~CharGer_Classification)#, scale = "free", space = "free", drop=T)
p = p + geom_density(aes(x = Allele_Frequency_Plot, fill = CharGer_Classification))
p = p + theme_bw() + theme_nogrid()
p = p + guides(fill=FALSE) + xlim(0,0.05) + labs(x="MAF in ExAC (%)")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90),
              axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
fn = 'out/pathVar_MAF_ExAC_by_classification.pdf'
ggsave(file=fn, useDingbats=FALSE)

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
num_AI = nrow(pathVar[pathVar$LOH_Sig =="Significant",])
cat("Num pathogenic AI:",num_AI, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing AI:",num_AI/num_total,"\n")
cat("\n")

cat("Oncogenes:\n")
num_total = nrow(pathVar_onco)
num_AI = nrow(pathVar_onco[pathVar_onco$LOH_Sig =="Significant",])
cat("Num pathogenic AI:",num_AI, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing AI:",num_AI/num_total,"\n")
cat("\n")

cat("Tumor Suppressor genes:\n")
num_total = nrow(pathVar_tsg)
num_AI = nrow(pathVar_tsg[pathVar_tsg$LOH_Sig =="Significant",])
cat("Num pathogenic AI:",num_AI, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing AI:",num_AI/num_total,"\n")
cat("\n")

# tn = "pan8000_pathogenic_oncogene.txt"
# write.table(pathVar_onco, quote=F, sep="\t", file = tn, row.names = F)

# overall fraction showing compound het
 # check if i need to lump in likely pathogenic here
num_total = nrow(pathVar)
num_CH = 56 #XX pathogenic variants
cat("Num pathogenic compound het:",num_CH, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing compound het:",num_CH/num_total,"\n")

cat("\n")

# overall fraction showing co-localizing
num_total = nrow(pathVar)
num_coloc = nrow(pathVar[pathVar$colocalized_somatic_mutation_count >=3,])
cat("Num pathogenic colocalizing:",num_coloc, "Num all pathogenic:", num_total,"\n")
cat("Percentage of all pathogenic variant showing co-localizing:",num_coloc/num_total,"\n")
cat("\n")
cat("Pathogenic variants gene distribution:")
table(pathVar[pathVar$colocalized_somatic_mutation_count >=3,]$HUGO_Symbol)
# cat("Likely pathogenic variants gene distribution:")
# table(pathVar[pathVar$colocalized_somatic_mutation_count >=3 & pathVar$CharGer_Classification=="Likely Pathogenic",]$HUGO_Symbol)
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

# # overall fraction showing clustering
# num_total = nrow(pathVar)
# num_clus = nrow(pathVar[pathVar$clustered,])
# cat("Num pathogenic clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
# cat("Percentage of all pathogenic variant showing clustering:",num_clus/num_total,"\n")
# cat("\n")
# cat("Pathogenic variants gene distribution:")
# table(pathVar[pathVar$clustered,]$HUGO_Symbol)
# # cat("Likely pathogenic variants gene distribution:")
# # table(pathVar[pathVar$clustered & pathVar$CharGer_Classification=="Likely Pathogenic",]$HUGO_Symbol)
# cat("\n")
# 
# cat("Oncogenes:\n")
# num_total = nrow(pathVar_onco)
# num_clus = nrow(pathVar_onco[pathVar_onco$clustered,])
# cat("Num pathogenic clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
# cat("Percentage of all pathogenic variant showing clustering:",num_clus/num_total,"\n")
# cat("\n")
# 
# cat("TSGs:\n")
# num_total = nrow(pathVar_tsg)
# num_clus = nrow(pathVar_tsg[pathVar_tsg$clustered,])
# cat("Num pathogenic clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
# cat("Percentage of all pathogenic variant showing clustering:",num_clus/num_total,"\n")
# cat("\n")
# 
# 
# # overall fraction showing co-localizing or clustering
# num_total = nrow(pathVar)
# num_clus = nrow(pathVar[(pathVar$clustered | (pathVar$colocalized_somatic_mutation_count >=3)),])
# cat("Num pathogenic co-localizing or clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
# cat("Percentage of all pathogenic variant showing co-localizing or clustering:",num_clus/num_total,"\n")
# cat("\n")
# 
# cat("Oncogenes:\n")
# num_total = nrow(pathVar_onco)
# num_clus = nrow(pathVar_onco[(pathVar_onco$clustered | (pathVar_onco$colocalized_somatic_mutation_count >=3)),])
# cat("Num pathogenic co-localizing or clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
# cat("Percentage of all pathogenic variant showing co-localizing or clustering:",num_clus/num_total,"\n")
# cat("\n")
# 
# cat("TSGs:\n")
# num_total = nrow(pathVar_tsg)
# num_clus = nrow(pathVar_tsg[(pathVar_tsg$clustered | (pathVar_tsg$colocalized_somatic_mutation_count >=3)),])
# cat("Num pathogenic co-localizing or clustering:",num_clus, "Num all pathogenic:", num_total,"\n")
# cat("Percentage of all pathogenic variant showing co-localizing or clustering:",num_clus/num_total,"\n")
# cat("\n")


# # split out truncations & missenses as statistical tests were conducted separately
# truncations = pathVar[pathVar$binary_type=="Truncation",]
# missenses = pathVar[pathVar$binary_type=="Missense",]
# 
# ##### Truncations #####
# ### stats ###
# get_low_exp_stat = function(term){
#   quantile_thres = 0.25
#   num_total = sum(truncations$CharGer_Classification==term)
#   num_low_exp = nrow(truncations[truncations$expression_quantile <= quantile_thres & truncations$CharGer_Classification==term,])
#   num_AI = nrow(truncations[truncations$Sig =="Significant" & truncations$CharGer_Classification==term,])
#   low_exp_fraction = num_low_exp/num_total
#   AI_fraction = num_AI/num_total
#   cat("Truncation Clinical_classification",term, "Number_of_variants",num_total,"\n")
#   cat("Truncation Expression quantile threshold",quantile_thres,"\n")
#   cat("Truncation Clinical_classification",term, "Number_of_low_exp_variants",num_low_exp,"\n")
#   cat("Truncation Clinical_classification",term, "Fraction_of_low_exp_variants",low_exp_fraction,"\n")
#   cat("Truncation Clinical_classification",term, "Number_of_AI_variants",num_AI,"\n")
#   cat("Truncation Clinical_classification",term, "Fraction_of_AI_variants",AI_fraction,"\n")
# }
# #get_low_exp_stat("High-confidence Pathogenic")
# get_low_exp_stat("Pathogenic")
# get_low_exp_stat("Likely Pathogenic")
# get_low_exp_stat("Uncertain Significance")
# 
# ### get stats for high-confidence variants specifically
# term = "High-confidence Pathogenic"
# quantile_thres = 0.2
# num_total = sum(truncations$Final_Classification==term)
# num_low_exp = nrow(truncations[truncations$expression_quantile <= quantile_thres & truncations$Final_Classification==term,])
# num_AI = nrow(truncations[truncations$Sig =="Significant" & truncations$Final_Classification==term,])
# low_exp_fraction = num_low_exp/num_total
# AI_fraction = num_AI/num_total
# cat("Truncation Clinical_classification",term, "Number_of_variants",num_total,"\n")
# cat("Truncation Expression quantile threshold",quantile_thres,"\n")
# cat("Truncation Clinical_classification",term, "Number_of_low_exp_variants",num_low_exp,"\n")
# cat("Truncation Clinical_classification",term, "Fraction_of_low_exp_variants",low_exp_fraction,"\n")
# cat("Truncation Clinical_classification",term, "Number_of_AI_variants",num_AI,"\n")
# cat("Truncation Clinical_classification",term, "Fraction_of_AI_variants",AI_fraction,"\n")
# 
# 
# # ### statistical test ###
# # hpath = truncations[truncations$Final_Classification=="High-confidence Pathogenic",]
# # path = truncations[truncations$CharGer_Classification=="Pathogenic",]
# # lpath = truncations[truncations$CharGer_Classification=="Likely Pathogenic",]
# # uncer = truncations[truncations$CharGer_Classification=="Uncertain Significance",]
# # cat("Truncation: K-S test for difference in expression effect across clinical classification")
# # ks.test(hpath$expression_quantile, uncer$expression_quantile)
# # ks.test(path$expression_quantile, uncer$expression_quantile)
# # ks.test(lpath$expression_quantile, uncer$expression_quantile)
# # #ks.test(path$expression_quantile, lpath$expression_quantile)
# # 
# # cat("Truncation: K-S test for difference in AI across clinical classification")
# # ks.test(hpath$tumor_by_normal_VAF, uncer$tumor_by_normal_VAF)
# # ks.test(path$tumor_by_normal_VAF, uncer$tumor_by_normal_VAF)
# # ks.test(lpath$tumor_by_normal_VAF, uncer$tumor_by_normal_VAF)
# # #ks.test(path$tumor_by_normal_VAF, lpath$tumor_by_normal_VAF)
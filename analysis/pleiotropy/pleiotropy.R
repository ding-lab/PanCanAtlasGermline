##### examine_gene_list.R #####
# Kuan-lin Huang @ WashU 201802
# examine the curated gene list

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/pleiotropy"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

pathVarP_otherSymptoms = pathVarP[!pathVarP$cancer_term_trait & pathVarP$binary_type=="Missense" & pathVarP$ClinVar_Pathogenicity == "Pathogenic",]
pathVarP_otherSymptoms_sele = pathVarP_otherSymptoms[,c(2,15,23,27,108,109)]
pathVarP_otherSymptoms_sele_uni = pathVarP_otherSymptoms_sele[!duplicated(pathVarP_otherSymptoms_sele$HGVSp),]

dim(pathVarP_otherSymptoms_sele_uni)

fn = "out/pleitropy_vars_inCPG.tsv"
write.table(pathVarP_otherSymptoms_sele_uni, file=fn, quote=F, sep="\t", col.names=T, row.names=F)

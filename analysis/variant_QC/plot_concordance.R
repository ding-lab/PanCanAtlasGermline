##### plot_pseq_stats.R #####
# Kuan-lin Huang @ WashU 2017 Aug./Sep.

setwd("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/variant_QC/")
source("../global_aes_out.R")
system("mkdir out")

##### individual level stats #####
fileName = "out/run_calc_vcf_concordance.chr22only.sampleID.GT.txt"
concordance = read.table(header=F, sep=" ", file=fileName)
colnames(concordance) = c("sample","validated","unvalidated")
concordance$validated[concordance$sample=="All"]/(concordance$validated[concordance$sample=="All"] + concordance$unvalidated[concordance$sample=="All"])
concordance =concordance[-c(9236:9240),]


clin_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/all.clin.merged.picked.txt"
clin = read.table(header=T, quote = "", sep="\t", row.names =NULL, file = clin_f, stringsAsFactors=FALSE)
clin = clin[,c("sample","cancer","race")]
colnames(clin) = c("sample","cancer","ethnicity")
clin = clin[!(clin$cancer %in% c("GBMLGG","COAD","KIPAN")),]

concordance_clin = merge(concordance,clin,by="sample",all.x=T)

p = ggplot(data=concordance_clin,aes(x = validated, y=unvalidated, color = cancer))
p = p + geom_point(alpha=0.1,stroke=0) 
p = p + theme_bw() + theme(legend.position="bottom")
p = p + geom_label(aes(label=ifelse(unvalidated > validated | unvalidated > 50, as.character(sample),NA)),size=1)
p
fn = paste('out/chr22_validated_vs_unvalidated.pdf',sep=".")
ggsave(file=fn, useDingbats=FALSE,limitsize=FALSE)

concordance_clin[concordance_clin$unvalidated > concordance_clin$validated | concordance_clin$unvalidated > 50,]

# ##### compare to pan8000 #####
# miss = read.table(header=F, sep="\t", "out/PCA_pan8000_missense_concordance.txt")
# trun = read.table(header=F, sep="\t", "out/PCA_pan8000_truncation_concordance.txt")
# PCA_sample = as.vector(t(read.table(header=F, sep="\t", "out/PCA_samples.txt")))
# colnames(miss) = c("Var","Sample","Validated")
# colnames(trun) = c("Var","Sample","Validated")
# 
# PCA_sample_ID = substr(PCA_sample,1,12)
# miss$inPCA = miss$Sample %in% PCA_sample_ID
# trun$inPCA = trun$Sample %in% PCA_sample_ID
# 
# table(miss$Validated,miss$inPCA)
# table(trun$Validated,trun$inPCA)

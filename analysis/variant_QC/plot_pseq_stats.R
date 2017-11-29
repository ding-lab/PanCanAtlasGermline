##### plot_pseq_stats.R #####
# Kuan-lin Huang @ WashU 2017 Aug./ updated Nov.

setwd("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/variant_QC/")
source("../global_aes_out.R")
source("../dependency_files.R")

##### individual level stats #####

fileName = "PCA.r1.TCGAbarcode.merge.exon.vcf.istats.tsv"
istat = read.table(header=TRUE, sep="\t", file=fileName, fill=T)

##### plot NVAR by ethnicity and cancer type #####
istat$bcr_patient_barcode = substr(istat$ID,1,12)
istat_clin = merge(istat,clin,by="bcr_patient_barcode")
# # filter for low call samples
#low_call_samples = istat_clin[istat_clin$NVAR <= 15000,]

# istat_clin[istat_clin$bcr_patient_barcode %in% swap_samples,c("NVAR","type","consensus_call")]

# nrow(low_call_samples)
# table(low_call_samples$type)
# tn = "out/low_call_samples.txt"
# write.table(low_call_samples, quote=F, sep="\t", file = tn, row.names = F)

istat_clin = istat_clin[istat_clin$NVAR > 15000,]

cat("Total number of variants in Exons: ", sum(istat_clin$NVAR),"\n")
cat("Average number of variants in Exons: ", mean(istat_clin$NVAR),"\n")

for (eth in unique(istat_clin$consensus_call)){
  cat("Average number of variants in Exons (", eth,"): ", mean(istat_clin$NVAR[!is.na(istat_clin$consensus_call) & istat_clin$consensus_call==eth]),"\n")
}

p = ggplot(data=istat_clin,aes(x = type, y=NVAR))
#p = p + geom_violin(alpha=0.1, fill=NA) 
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
             geom = "crossbar", width = 0.8)
p = p + geom_dotplot(aes(fill= consensus_call),binaxis="y", stackdir = "center",colour=NA,binwidth=100)
#p = p + geom_jitter(aes(color= consensus_call),height=0, alpha = 0.1, size=1) 
p = p  + theme_bw() #+ theme(legend.position="none")
p = p + labs(x = "Cancer", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p
fn = 'out/PCA10467.pseq.istat.NVAR.dotplot.pdf'
ggsave(file=fn, w = 10, useDingbats=FALSE,limitsize=FALSE)

# p = ggplot(data=istat_clin,aes(x = type, y=DP))
# #p = p + geom_violin(alpha=0.1, fill=NA) 
# p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#                      geom = "crossbar", width = 0.8)
# p = p + geom_dotplot(aes(fill= consensus_call),binaxis="y", stackdir = "center",colour=NA,binwidth=1000)
# #p = p + geom_jitter(aes(color= consensus_call),height=0, alpha = 0.1, size=1) 
# p = p  + theme_bw() #+ theme(legend.position="none")
# p = p + labs(x = "Cancer", y = "DP")
# p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
# p
# fn = 'out/PCA10467.pseq.istat.DP.dotplot.pdf'
# ggsave(file=fn, w = 10, useDingbats=FALSE,limitsize=FALSE)
# 
# p = ggplot(data=istat_clin,aes(x = type, y=TITV))
# #p = p + geom_violin(alpha=0.1, fill=NA) 
# p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#                      geom = "crossbar", width = 0.8)
# p = p + geom_dotplot(aes(fill= consensus_call),binaxis="y", stackdir = "center",colour=NA,binwidth=1000)
# #p = p + geom_jitter(aes(color= consensus_call),height=0, alpha = 0.1, size=1) 
# p = p  + theme_bw() #+ theme(legend.position="none")
# p = p + labs(x = "Cancer", y = "TiTv")
# p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
# p
# fn = 'out/PCA10467.pseq.istat.TITV.dotplot.pdf'
# ggsave(file=fn, w = 10, useDingbats=FALSE,limitsize=FALSE)
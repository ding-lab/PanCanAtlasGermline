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
#istat_clin = istat_clin[istat_clin$NVAR > 15000,]

cat("Total number of variants in Exons: ", sum(istat_clin$NVAR),"\n")
cat("Average number of variants in Exons: ", mean(istat_clin$NVAR),"\n")

for (eth in unique(istat_clin$consensus_call)){
  cat("Average number of variants in Exons (", eth,"): ", mean(istat_clin$NVAR[!is.na(istat_clin$consensus_call) & istat_clin$consensus_call==eth]),"\n")
}

p = ggplot(data=istat_clin,aes(x = type, y=NVAR))
#p = p + geom_violin(alpha=0.1, fill=NA) 
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
             geom = "crossbar", width = 0.8)
#p = p + geom_dotplot(aes(fill= consensus_call),binaxis="y", stackdir = "center",colour=NA,binwidth=100)
p = p + geom_jitter(aes(color= consensus_call),height=0, alpha = 0.1, size=1) 
p = p  + theme_bw() #+ theme(legend.position="none")
p = p + labs(x = "Cancer", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10467.pseq.istat.dotplot.pdf'
ggsave(file=fn, w = 10, useDingbats=FALSE,limitsize=FALSE)


##### plot NVAR by ethnicity and cancer type #####
fileName = "PCA.r1.TCGAbarcode.merge.4.1.istats.tsv"
istatAll = read.table(header=TRUE, sep="\t", file=fileName, fill=T)
istatAll$bcr_patient_barcode = substr(istatAll$ID,1,12)
istatAll_clin = merge(istatAll,clin,by="bcr_patient_barcode")

cat("Total number of variants in all WES: ", sum(istatAll_clin$NVAR),"\n")
cat("Average number of variants in all WES: ", mean(istatAll_clin$NVAR),"\n")

for (eth in unique(istatAll_clin$consensus_call)){
  cat("Average number of variants in all WES (", eth,"): ", mean(istatAll_clin$NVAR[!is.na(istatAll_clin$consensus_call) & istatAll_clin$consensus_call==eth]),"\n")
}

p = ggplot(data=istatAll_clin,aes(x = type, y=log10(NVAR)))
#p = p + geom_violin(alpha=0.1, fill=NA) 
#p = p + geom_dotplot(aes(fill= consensus_call),binaxis="y", stackdir = "center",colour=NA,binwidth=100)
p = p + geom_jitter(aes(color= consensus_call),height=0, alpha = 0.1, size=0.5) 
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8, size=0.02)
p = p  + theme_bw() #+ theme(legend.position="none")
p = p + labs(x = "Cancer", y = "log10(Count of total variants)")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10467.pseq.istatAll.NVAR.log10.pdf'
ggsave(file=fn, w = 10, useDingbats=FALSE,limitsize=FALSE)

##### coverage #####
coverageF = "coverage_by_sample.wclin_10467.20171201.txt"
coverage = read.table(header=TRUE, sep="\t", file=coverageF, fill=T)
colnames(coverage)[1] = "bcr_patient_barcode"
coverage_clin = merge(coverage,clin,by="bcr_patient_barcode")

p = ggplot(data=coverage_clin,aes(x = type, y=EstdDepthOfCoverage ))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8)
p = p + geom_jitter(aes(color= consensus_call),height=0, alpha = 0.1, size=1) 
p = p  + theme_bw() 
p = p + labs(x = "Cancer", y = "Coverage of cancer predisposition genes")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10467.coverage.dotplot.pdf'
ggsave(file=fn, w = 10, useDingbats=FALSE,limitsize=FALSE)

##### coverage x NVAR #####
coverage_istat = merge(coverage,istat,by="bcr_patient_barcode")
coverage_istat_clin = merge(coverage_istat,clin, by="bcr_patient_barcode")
p = ggplot(data=coverage_istat_clin,aes(x = EstdDepthOfCoverage, y=NVAR ))
p = p + geom_point(aes(color= type),alpha = 0.1, size=1,stroke=0) #
p = p  + theme_bw() + getPCACancerColor()
p = p + labs(x = "Coverage of cancer predisposition genes", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
#p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10467.coverage.vs.NVAR.pdf'
ggsave(file=fn, useDingbats=FALSE,limitsize=FALSE)

##### # pathogenic variants #####
pathCount = data.frame(table(pathVarP$bcr_patient_barcode))
colnames(pathCount) = c("bcr_patient_barcode","pathVarCount")
pathCount_clin = merge(pathCount,clin, by="bcr_patient_barcode",all.y=T)
pathCount_clin$pathVarCount[is.na(pathCount_clin$pathVarCount)]=0
p = ggplot(data=pathCount_clin,aes(x = type, y=pathVarCount ))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8)
p = p + geom_jitter(aes(color= consensus_call),height=0, alpha = 0.1, size=1) 
p = p  + theme_bw() 
p = p + labs(x = "Cancer", y = "Coverage of cancer predisposition genes")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10467.pathVarCount.pdf'
ggsave(file=fn, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=pathCount_clin,aes(x = pathVarCount ))
p = p + facet_grid(.~type)
p = p + geom_histogram(aes(fill= consensus_call),binwidth = 1) 
p = p  + theme_bw() 
p = p + labs(x = "Pathogenic variant counts", y = "Number of Samples")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) + xlim(-0.5,2.5)
p
fn = 'out/PCA10467.pathVarCount.pdf'
ggsave(file=fn, w=20,useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=pathCount_clin,aes(x = pathVarCount ))
p = p + facet_grid(.~type)
p = p + geom_density(aes(fill= consensus_call),alpha=0.3) 
p = p  + theme_bw() 
p = p + labs(x = "Pathogenic variant counts", y = "Density of sample countss")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) + xlim(-0.5,2.5)
p
fn = 'out/PCA10467.pathVarCount.density.pdf'
ggsave(file=fn, w=20,useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=pathCount_clin,aes( x = type, y = pathVarCount ))
#p = p + facet_grid(.~type)
p = p + geom_violin() 
p = p + geom_jitter(aes(color=consensus_call), width =0.2, height = 0, alpha=0.1) 
p = p  + theme_bw() 
p = p + labs(x = "Pathogenic variant counts", y = "Number of Samples")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + ylim(0,2.5)
p
fn = 'out/PCA10467.pathVarCount_violin.pdf'
ggsave(file=fn, w=10,useDingbats=FALSE,limitsize=FALSE)

##### 
pathCount_clin_istat = merge(pathCount_clin, istat,by="bcr_patient_barcode")
p = ggplot(data=pathCount_clin_istat[pathCount_clin_istat$pathVarCount<3,],aes( x = factor(pathVarCount), y = NVAR ))
p = p + facet_grid(.~type)
p = p + geom_violin() 
p = p + geom_jitter(aes(color=consensus_call), width =0.2, height = 0, alpha=0.2,stroke=0,shape=19,size=0.5) 
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8,size=0.02, alpha=0.1)
p = p  + theme_bw() 
p = p + labs(x = "Pathogenic Variant Counts", y = "Total Exonic Variant Counts")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14))
p = p + theme(legend.position = "bottom")
p = p + ylim(y=0,45000)
p
fn = 'out/PCA10467.pathVarCount_NVAR_violin.pdf'
ggsave(file=fn, w=20,useDingbats=FALSE,limitsize=FALSE)

# see if there is a correlation between number of exonic variants and pathogenic variants
cor.test(pathCount_clin_istat$NVAR,pathCount_clin_istat$pathVarCount)
fit = lm(data = pathCount_clin_istat, pathVarCount ~ NVAR)
summary(fit)
fit = lm(data = pathCount_clin_istat, pathVarCount ~ NVAR + type)
summary(fit)
for (cancer in unique(pathCount_clin_istat$type)){
  m = pathCount_clin_istat[pathCount_clin_istat$type==cancer,]
  cat("#####",cancer,"#####\n")
  print(cor.test(m$NVAR,m$pathVarCount))
}
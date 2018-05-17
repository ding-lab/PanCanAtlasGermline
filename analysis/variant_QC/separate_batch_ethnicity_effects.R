##### separate_batch_ethnicity_effects.R #####
# Kuan-lin Huang @ WashU 2017 Aug./ updated Nov.

#setwd("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/variant_QC/")
source("../global_aes_out.R")
source("../dependency_files.R")

# implementation of LMG to find independent contribution of regressors
# library("relaimpo") # couldn't get this to work, seems to require all quantitative variables
library("hier.part")

##### individual level stats #####

fileName = "PCA.r1.TCGAbarcode.merge.exon.vcf.istats.tsv"
istat = read.table(header=TRUE, sep="\t", file=fileName, fill=T)
istat$TSS = gsub("TCGA-(..)-.*","\\1",istat$ID)
istat$Center = gsub(".*(..)","\\1",istat$ID)
istat$Analyte = gsub("TCGA-..-....-...-..(.)-.*","\\1",istat$ID)
# table(istat$Analyte)
# table(istat$TSS)
# table(istat$Center) #https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes
istat$Center[istat$Center=="01"] = "BI" # both BI
istat$Center[istat$Center=="08"] = "BI" # both BI
istat$Center[istat$Center=="09"] = "WUSM"
istat$Center[istat$Center=="10"] = "BCM"
istat$Center[istat$Center=="32"] = "SANGER"

##### check overall variant count as a function of ancestry and/or batch #####

##### merge all data #####
istat$bcr_patient_barcode = substr(istat$ID,1,12)
istat_clin = merge(istat,clin,by="bcr_patient_barcode")

coverageF = "coverage_by_sample.wclin_10467.20171201.txt"
coverage = read.table(header=TRUE, sep="\t", file=coverageF, fill=T)
colnames(coverage)[1] = "bcr_patient_barcode"
istat_clin_coverage = merge(coverage,istat_clin,by="bcr_patient_barcode")

pathCount = data.frame(table(pathVarP$bcr_patient_barcode))
colnames(pathCount) = c("bcr_patient_barcode","pathVarCount")
istat_clin_coverage_pathCount = merge(pathCount,istat_clin_coverage, by="bcr_patient_barcode",all.y=T)
colnames(istat_clin_coverage_pathCount)[which(colnames(istat_clin_coverage_pathCount)=="consensus_call")] = "Ancestry"
colnames(istat_clin_coverage_pathCount)[which(colnames(istat_clin_coverage_pathCount)=="type")] = "Cancer"
istat_clin_coverage_pathCount$pathVarCount[is.na(istat_clin_coverage_pathCount$pathVarCount)]=0
istat_clin_coverage_pathCount_brief = istat_clin_coverage_pathCount[,c(1:20,55:58)]
fn = "out/10389_data_for_batch_effect_analysis.tsv"
write.table(istat_clin_coverage_pathCount_brief, file=fn, quote=F, sep="\t", col.names=T, row.names=F)

# reporting stats
cat("Total number of variants in Exons: ", sum(istat_clin_coverage_pathCount$NVAR),"\n")
cat("Average number of variants in Exons: ", mean(istat_clin_coverage_pathCount$NVAR),"\n")
for (eth in unique(istat_clin_coverage_pathCount$Ancestry)){
  cat("Average number of variants in Exons (", eth,"): ", mean(istat_clin_coverage_pathCount$NVAR[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry==eth]),"\n")
}

# implement LMG (1980) to identify the independent contribution from regressions
istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("NVAR","Center","Analyte","Cancer","Ancestry")]),c("NVAR","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$NVAR, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")

istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("pathVarCount","Center","Analyte","Cancer","Ancestry")]),c("pathVarCount","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$pathVarCount, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")

# check for overall exomic variants
fit = lm(data = istat_clin_coverage_pathCount, NVAR ~ Ancestry + Cancer + Center + Analyte) # the order of importance found using LMG
summary(fit)
results = data.frame(anova(fit))
colnames(results) = gsub("\\.","",colnames(results))
results = rbind(results, colSums(results))
results$PercSumSq = results$SumSq / 100*results$SumSq[length(results$SumSq)]
fn = "out/NVAR_model.tsv"
write.table(results, file=fn, quote=F, sep="\t", col.names=T, row.names=T)

# check for pathogenic variants
fit = lm(data = istat_clin_coverage_pathCount, pathVarCount ~ Cancer + Center + Analyte + Ancestry) # the order of importance found using LMG
summary(fit)
results = data.frame(anova(fit))
colnames(results) = gsub("\\.","",colnames(results))
results = rbind(results, colSums(results))
results$PercSumSq = results$SumSq / 100*results$SumSq[length(results$SumSq)]
fn = "out/pathVarCount_model.tsv"
write.table(results, file=fn, quote=F, sep="\t", col.names=T, row.names=T)

# titv 
istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("TITV","Center","Analyte","Cancer","Ancestry")]),c("TITV","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$TITV, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")

fit = lm(data = istat_clin_coverage_pathCount, TITV ~ Cancer + Ancestry + Center + Analyte)
summary(fit)
results = data.frame(anova(fit))
colnames(results) = gsub("\\.","",colnames(results))
results = rbind(results, colSums(results))
results$PercSumSq = results$SumSq / 100*results$SumSq[length(results$SumSq)]
fn = "out/titv_model.tsv"
write.table(results, file=fn, quote=F, sep="\t", col.names=T, row.names=T)


##### implement a stricter regression test on pathogenic variant test #####
istat_clin_coverage_pathCount_eur = istat_clin_coverage_pathCount[istat_clin_coverage_pathCount$Ancestry=="eur",]
pathVarP_eur = pathVarP[pathVarP$AIM_ethnicity=="eur",]
# load overall potential association cancer results 
tn = "../burden_assoc/out/ExAC_vs_cancer_pathogenic_variants_burden.tsv"
all_cancer_stat_m = read.table(header=TRUE, sep="\t", file=tn)
all_cancer_stat_m_suggest = all_cancer_stat_m[all_cancer_stat_m$FDR < 0.15,]

tn = "../burden_assoc/out/cancer_vs_otherCancer_pathogenic_variants_burden.tsv"
all_cancer_stat_m_by_cancer = read.table(header=TRUE, sep="\t", file=tn)

test_genes = names(table(pathVarP_eur$HUGO_Symbol)[table(pathVarP_eur$HUGO_Symbol)>1])
all_cancer_stat_byCancer = vector("list")
i = "a"
for (gene in test_genes){
  pathVarP_eur_g = pathVarP_eur[!is.na(pathVarP_eur$HUGO_Symbol) & pathVarP_eur$HUGO_Symbol == gene,]
  cancer_to_exclude = all_cancer_stat_m_suggest$Cancer[all_cancer_stat_m_suggest$Gene == gene]
  for (cancer in unique(pathVarP_eur_g$cancer)){
    pathVarP_eur_g$carrier = 1
    pathVarP_eur_g_for_merge = pathVarP_eur_g[,c("bcr_patient_barcode","carrier")]
    istat_clin_coverage_pathCount_eur_overall = merge(istat_clin_coverage_pathCount_eur, pathVarP_eur_g_for_merge, by="bcr_patient_barcode",all.x=T)
    istat_clin_coverage_pathCount_eur_overall$carrier[is.na(istat_clin_coverage_pathCount_eur_overall$carrier)] = 0
    istat_clin_coverage_pathCount_eur_overall$inCancer = 0
    istat_clin_coverage_pathCount_eur_overall$inCancer[istat_clin_coverage_pathCount_eur_overall$Cancer == cancer] = 1
    test_set = istat_clin_coverage_pathCount_eur_overall[istat_clin_coverage_pathCount_eur_overall$inCancer==1 | 
                                                         !(istat_clin_coverage_pathCount_eur_overall$Cancer %in% cancer_to_exclude),]
    table(test_set$inCancer)
    # Model: carrier ~ cancer + analyte + center 
    fit = glm(data = test_set, carrier ~ inCancer + Analyte + Center)
    stats = as.data.frame(t(summary(fit)$coef[2,]))
    stats$Cancer = cancer
    stats$Gene = gene
    all_cancer_stat_byCancer[[i]]=stats
    i = paste(i,"a",sep="")
  }
}

all_cancer_stat_byCancer_m = do.call(rbind,all_cancer_stat_byCancer)
row.names(all_cancer_stat_byCancer_m) = NULL
colnames(all_cancer_stat_byCancer_m) = c("logic_Estimate","logic_StdError","logic_tvalue","logic_P","Cancer","Gene")
all_cancer_stat_byCancer_m$logic_FDR = p.adjust(all_cancer_stat_byCancer_m$logic_P, method="BH")
all_cancer_stat_byCancer_m = all_cancer_stat_byCancer_m[order(all_cancer_stat_byCancer_m$logic_P),]

all_cancer_stat_m_by_cancer_merged = merge(all_cancer_stat_byCancer_m,all_cancer_stat_m_by_cancer,by=c("Cancer","Gene"),all=T)
all_cancer_stat_m_by_cancer_merged = all_cancer_stat_m_by_cancer_merged[order(all_cancer_stat_m_by_cancer_merged$P),]

tn = "out/cancer_vs_otherCancer_pathogenic_variants_logic_and_burden.tsv"
write.table(all_cancer_stat_m_by_cancer_merged, quote=F, sep="\t", file = tn, row.names = F)

##### plotting by cancer Cancer and the specific variable of interest #####

### all variants in exome ###
p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry != "mix",],aes(x = Ancestry, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Ancestry))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Ancestry", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.NVAR.by.ancestry.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Center),],aes(x = Center, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.NVAR.by.center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Analyte),],aes(x = Analyte, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Analyte))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Analyte", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.NVAR.by.Analyte.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

# TITV
p = ggplot(data=istat_clin_coverage_pathCount,aes(x = Center, y=TITV))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Ti/Tv")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.titv.by.center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

# the six cancers
six_cancers = c("SKCM","STAD","THCA","BRCA","UCEC","LIHC")
istat_clin_coverage_pathCount_sixc = istat_clin_coverage_pathCount[istat_clin_coverage_pathCount$Cancer %in% six_cancers,]
p = ggplot(data=istat_clin_coverage_pathCount_sixc[!is.na(istat_clin_coverage_pathCount_sixc$Ancestry) & istat_clin_coverage_pathCount_sixc$Ancestry != "mix",],aes(x = Ancestry, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Ancestry))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Ancestry", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.NVAR.by.ancestry.sixC.pdf'
ggsave(file=fn, h=3,w = 6, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_sixc[!is.na(istat_clin_coverage_pathCount_sixc$Center),],aes(x = Center, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.NVAR.by.center.sixC.pdf'
ggsave(file=fn, h=3,w = 6, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_sixc[!is.na(istat_clin_coverage_pathCount_sixc$Analyte),],aes(x = Analyte, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Analyte))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Analyte", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.NVAR.by.Analyte.sixC.pdf'
ggsave(file=fn, h=3,w = 6, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_sixc,aes(x = Center, y=TITV))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Ti/Tv")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.titv.by.center.sixC.pdf'
ggsave(file=fn, h=3,w = 6, useDingbats=FALSE,limitsize=FALSE)

# ### pathogenic variants ###
# p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry != "mix",],aes(x = Ancestry, y=pathVarCount))
# p = p + facet_grid(.~Cancer, drop = T)
# p = p + geom_violin(alpha=0.8, aes(fill = Ancestry))
# p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#                      geom = "crossbar", width = 0.8)
# p = p + theme_bw() + theme(legend.position="none")
# p = p + labs(x = "Ancestry", y = "Count of predisposing variants")
# p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
# p = p + theme(legend.position = "bottom")
# p = p + expand_limits(y=0)
# p
# fn = 'out/PCA10389.pathVarCount.by.ancestry.pdf'
# ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)
# 
# p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Center),],aes(x = Center, y=pathVarCount))
# p = p + facet_grid(.~Cancer, drop = T)
# p = p + geom_violin(alpha=0.8, aes(fill = Center))
# p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#                      geom = "crossbar", width = 0.8)
# p = p + theme_bw() + theme(legend.position="none")
# p = p + labs(x = "Center", y = "Count of predisposing variants")
# p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
# p = p + theme(legend.position = "bottom")
# p = p + expand_limits(y=0)
# p
# fn = 'out/PCA10389.pathVarCount.by.center.pdf'
# ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)
# 
# p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Analyte),],aes(x = Analyte, y=pathVarCount))
# p = p + facet_grid(.~Cancer, drop = T)
# p = p + geom_violin(alpha=0.8, aes(fill = Analyte))
# p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#                      geom = "crossbar", width = 0.8)
# p = p + theme_bw() + theme(legend.position="none")
# p = p + labs(x = "Analyte", y = "Count of predisposing variants")
# p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
# p = p + theme(legend.position = "bottom")
# p = p + expand_limits(y=0)
# p
# fn = 'out/PCA10389.pathVarCount.by.Analyte.pdf'
# ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)


# ####
# p = ggplot(data=istat_clin_coverage_pathCount[istat_clin_coverage_pathCount$pathVarCount<3,],aes( x = factor(pathVarCount), y = NVAR ))
# p = p + facet_grid(.~Cancer)
# p = p + geom_violin()
# p = p + geom_jitter(aes(color=Ancestry), width =0.2, height = 0, alpha=0.2,stroke=0,shape=19,size=0.5)
# p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8,size=0.02, alpha=0.8)
# p = p  + theme_bw()
# p = p + labs(x = "Pathogenic Variant Counts", y = "Total Exonic Variant Counts")
# p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14))
# p = p + theme(legend.position = "bottom")
# p = p + ylim(y=0,45000)
# p
# fn = 'out/PCA10389.pathVarCount_NVAR_violin.pdf'
# ggsave(file=fn, h=3.5,w=15,useDingbats=FALSE,limitsize=FALSE)


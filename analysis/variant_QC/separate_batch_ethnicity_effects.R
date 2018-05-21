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

coverage152F = "coverage_by_gene.revert_id_swap.filtered_by_PCAgermline_paper_analysisIds2.tsv.gz"
coverage152 = read.table(header=TRUE, sep="\t", file=gzfile(coverage152F), fill=T)

colnames(coverage152)[2] = "bcr_patient_barcode"
colnames(coverage152)[6:length(colnames(coverage152))] = paste("coverage",colnames(coverage152)[6:length(colnames(coverage152))],sep="_")
coverage152$coverage152g =  rowSums(coverage152[,6:length(colnames(coverage152))])/152
coverage152$coverage152g_over10 =  rowSums(coverage152[,6:length(colnames(coverage152))]>10)/152
istat_clin_coverage_pathCount = merge(coverage152,istat_clin_coverage_pathCount,by="bcr_patient_barcode")

length(istat_clin_coverage_pathCount$coverage152g[istat_clin_coverage_pathCount$coverage152g>=10])

## 152 gene stats ###
fileName = "PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.152g.41.vcf.gz.pseq.istats.tsv"
istat_152 = read.table(header=TRUE, sep="\t", file=fileName, fill=T)
colnames(istat_152) = paste("Predispose152Genes",colnames(istat_152),sep="_")
istat_152$bcr_patient_barcode = substr(istat_152$Predispose152Genes_ID,1,12)
istat_clin_coverage_pathCount = merge(istat_clin_coverage_pathCount,istat_152, by="bcr_patient_barcode",all.x=T)

## 99 gene stats ###
fileName = "10389.99g.pseq.istats.tsv"
istat_99 = read.table(header=TRUE, sep="\t", file=fileName, fill=T)
colnames(istat_99) = paste("Predispose99Genes",colnames(istat_99),sep="_")
istat_99$bcr_patient_barcode = substr(istat_99$Predispose99Genes_ID,1,12)
istat_clin_coverage_pathCount = merge(istat_clin_coverage_pathCount,istat_99, by="bcr_patient_barcode",all.x=T)

## 21 gene stats ###
fileName = "10389.21g.pseq.istats.tsv"
istat_21 = read.table(header=TRUE, sep="\t", file=fileName, fill=T)
colnames(istat_21) = paste("Predispose21Genes",colnames(istat_21),sep="_")
istat_21$bcr_patient_barcode = substr(istat_21$Predispose21Genes_ID,1,12)
istat_clin_coverage_pathCount = merge(istat_clin_coverage_pathCount,istat_21, by="bcr_patient_barcode",all.x=T)

colnames(istat_clin_coverage_pathCount)[which(colnames(istat_clin_coverage_pathCount)=="consensus_call")] = "Ancestry"
colnames(istat_clin_coverage_pathCount)[which(colnames(istat_clin_coverage_pathCount)=="type")] = "Cancer"
istat_clin_coverage_pathCount$pathVarCount[is.na(istat_clin_coverage_pathCount$pathVarCount)]=0
# istat_clin_coverage_pathCount_brief = istat_clin_coverage_pathCount[,c(1:20,55:58)]
# fn = "out/10389_data_for_batch_effect_analysis.tsv"
# write.table(istat_clin_coverage_pathCount_brief, file=fn, quote=F, sep="\t", col.names=T, row.names=F)

# reporting stats
cat("Total number of variants in Exons: ", sum(istat_clin_coverage_pathCount$NVAR),"\n")
cat("Average number of variants in Exons: ", mean(istat_clin_coverage_pathCount$NVAR),"\n")
for (eth in unique(istat_clin_coverage_pathCount$Ancestry)){
  cat("Average number of variants in Exons (", eth,"): ", mean(istat_clin_coverage_pathCount$NVAR[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry==eth]),"\n")
}
mean(istat_clin_coverage_pathCount_eur$TITV,na.rm=T)

### implement LMG (1980) to identify the independent contribution from regressions
# coverage
cat("LMG: 152 gene coverage\n")
pdf("out/152coverage_effect_lmg_iEffects.pdf")
istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("coverage152g","Center","Analyte","Cancer","Ancestry")]),c("coverage152g","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$coverage152g, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")
dev.off()

cat("LMG: 152 gene coverage >= 10\n")
pdf("out/152coverage_gt10_effect_lmg_iEffects.pdf")
istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("coverage152g_over10","Center","Analyte","Cancer","Ancestry")]),c("coverage152g_over10","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$coverage152g_over10, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")
dev.off()

# all variants
cat("LMG: CDS\n")
pdf("out/CDS_variantCounts_lmg_iEffects.pdf")
istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("NVAR","Center","Analyte","Cancer","Ancestry")]),c("NVAR","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$NVAR, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")
dev.off()

# variants in 152 genes
cat("LMG: 152 genes\n")
pdf("out/152g_variantCounts_lmg_iEffects.pdf")
istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("Predispose152Genes_NVAR","Center","Analyte","Cancer","Ancestry")]),c("Predispose152Genes_NVAR","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$Predispose152Genes_NVAR, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")
dev.off()

# variants in 99 genes
cat("LMG: 99 genes\n")
pdf("out/99g_variantCounts_lmg_iEffects.pdf")
istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("Predispose99Genes_NVAR","Center","Analyte","Cancer","Ancestry")]),c("Predispose99Genes_NVAR","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$Predispose99Genes_NVAR, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")
dev.off()

# variants in 21 genes
cat("LMG: 21 genes\n")
pdf("out/21g_variantCounts_lmg_iEffects.pdf")
istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("Predispose21Genes_NVAR","Center","Analyte","Cancer","Ancestry")]),c("Predispose21Genes_NVAR","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$Predispose21Genes_NVAR, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")
dev.off()

# pathogenic variants
cat("LMG: predisposing variants\n")
pdf("out/Presiposing_variantCounts_lmg_iEffects.pdf")
istat_clin_coverage_pathCount_test = istat_clin_coverage_pathCount[complete.cases(istat_clin_coverage_pathCount[,c("pathVarCount","Center","Analyte","Cancer","Ancestry")]),c("pathVarCount","Center","Analyte","Cancer","Ancestry")]
hier.part(istat_clin_coverage_pathCount_test$pathVarCount, istat_clin_coverage_pathCount_test[,-1], family = "gaussian", gof = "Rsqu")
dev.off()

# check for overall exomic variants
fit = lm(data = istat_clin_coverage_pathCount, NVAR ~ Ancestry + Cancer + Center + Analyte) # the order of importance found using LMG
summary(fit)
results = data.frame(anova(fit))
colnames(results) = gsub("\\.","",colnames(results))
results = rbind(results, colSums(results))
results$PercSumSq = results$SumSq / 100*results$SumSq[length(results$SumSq)]
fn = "out/NVAR_model.tsv"
write.table(results, file=fn, quote=F, sep="\t", col.names=T, row.names=T)

# check for variants in 99 genes
fit = lm(data = istat_clin_coverage_pathCount, Predispose99Genes_NVAR ~ Ancestry + Cancer + Center + Analyte) # the order of importance found using LMG
summary(fit)
results = data.frame(anova(fit))

# check for pathogenic variants
fit = lm(data = istat_clin_coverage_pathCount, pathVarCount ~ Cancer + Center + Analyte + Ancestry) # the order of importance found using LMG
summary(fit)
results = data.frame(anova(fit))
colnames(results) = gsub("\\.","",colnames(results))
results = rbind(results, colSums(results))
results$PercSumSq = results$SumSq / 100*results$SumSq[length(results$SumSq)]
fn = "out/pathVarCount_model.tsv"
write.table(results, file=fn, quote=F, sep="\t", col.names=T, row.names=T)

# one-way anova:
cat("One-way anova tests \n")
fit = lm(data = istat_clin_coverage_pathCount, pathVarCount ~ Cancer) # the order of importance found using LMG
anova(fit)
fit = lm(data = istat_clin_coverage_pathCount, pathVarCount ~ Center) # the order of importance found using LMG
anova(fit)
fit = lm(data = istat_clin_coverage_pathCount, pathVarCount ~ Analyte) # the order of importance found using LMG
anova(fit)
fit = lm(data = istat_clin_coverage_pathCount, pathVarCount ~ Ancestry) # the order of importance found using LMG
anova(fit)

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

##### plotting by cancer and the specific variable of interest #####
### coverage in 152 genes ###
p = ggplot(data=istat_clin_coverage_pathCount,aes(x = Center, y=coverage152g))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Coverage")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) + geom_hline(alpha=0.5, yintercept = 10)
p
fn = 'out/PCA10389.coverage.by.Center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

### all variants in exome ###
p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry != "mix",],aes(x = Ancestry, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Ancestry))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
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
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
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
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
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
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Ti/Tv")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.titv.by.center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

# the six cancers
six_cancers = c("SKCM","STAD","THCA","BRCA","UCEC","LIHC")
istat_clin_coverage_pathCount_sixc = istat_clin_coverage_pathCount_eur[istat_clin_coverage_pathCount_eur$Cancer %in% six_cancers,]
istat_clin_coverage_pathCount_sixc$Cancer = factor(istat_clin_coverage_pathCount_sixc$Cancer,levels = c("SKCM","STAD","THCA","BRCA","UCEC","LIHC"))
p = ggplot(data=istat_clin_coverage_pathCount_sixc[!is.na(istat_clin_coverage_pathCount_sixc$Ancestry) & istat_clin_coverage_pathCount_sixc$Ancestry != "mix",],aes(x = Ancestry, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Ancestry))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Ancestry", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.NVAR.by.ancestry.sixC.pdf'
ggsave(file=fn, h=3,w = 6, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_sixc[!is.na(istat_clin_coverage_pathCount_sixc$Center),],aes(x = Cancer, y=NVAR))
#p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.NVAR.by.center.sixC.pdf'
ggsave(file=fn, h=3,w = 4, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_sixc[!is.na(istat_clin_coverage_pathCount_sixc$Center),],aes(x = Cancer, y=Predispose152Genes_NVAR))
#p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of variants in 152 genes")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.Predispose152Genes_NVAR.by.center.sixC.pdf'
ggsave(file=fn, h=3,w = 4, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_sixc[!is.na(istat_clin_coverage_pathCount_sixc$Analyte),],aes(x = Analyte, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Analyte))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Analyte", y = "Count of exonic variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.NVAR.by.Analyte.sixC.pdf'
ggsave(file=fn, h=3,w = 6, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_sixc,aes(x = Cancer, y=TITV))
#p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="bottom")
p = p + labs(x = "Center", y = "Ti/Tv")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.titv.by.center.sixC.pdf'
ggsave(file=fn, h=3,w = 4, useDingbats=FALSE,limitsize=FALSE)

### stats

# all still significant
fit = lm(data =istat_clin_coverage_pathCount_sixc, TITV ~ Center )
anova(fit)

fit = lm(data =istat_clin_coverage_pathCount_sixc, NVAR ~ Center )
anova(fit)

fit = lm(data =istat_clin_coverage_pathCount_sixc, Predispose152Genes_NVAR ~ Center )
anova(fit)

fit = lm(data =istat_clin_coverage_pathCount_sixc, pathVarCount ~ Center )
anova(fit)


### all together
p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Center),],aes(x = Center, y=NVAR))
p = p + geom_violin(aes(fill = Center),color=NA,alpha=0.2)
#p = p + geom_jitter(aes(color = Ancestry),width =0.2, height = 0, alpha=0.1,stroke=0,shape=19,size=0.5)
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of CDS variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.NVAR.by.center.all.pdf'
ggsave(file=fn, h=3,w = 3, useDingbats=FALSE,limitsize=FALSE)

### 152 gene variants ###
p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry != "mix",],aes(x = Ancestry, y=Predispose152Genes_NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Ancestry))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Ancestry", y = "Count of variants in 152 genes")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.152geneVarCount.by.ancestry.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount,aes(x = Center, y=Predispose152Genes_NVAR))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of variants in 152 genes")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.152geneVarCount.by.center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

### Caucasians only ###
# all exome
p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry == "eur",],aes(x = Center, y=NVAR))
p = p + facet_grid(.~Cancer, drop = T, space = "free_x" , scale= "free_x")
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of variants in CDS")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.eurExomeVarCount.by.center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

# 152 genes
p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry == "eur",],aes(x = Center, y=Predispose152Genes_NVAR))
p = p + facet_grid(.~Cancer, drop = T, space = "free_x" , scale= "free_x")
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of variants in 152 genes")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.eur152geneVarCount.by.center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry == "eur",],aes(x = Center, y=Predispose21Genes_NVAR))
p = p + facet_grid(.~Cancer, drop = T, space = "free_x" , scale= "free_x")
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of variants in 21 genes")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.eur21geneVarCount.by.center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

# predisposing variants
p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry == "eur",],aes(x = Center, y=pathVarCount))
p = p + facet_grid(.~Cancer, drop = T, space = "free_x" , scale= "free_x")
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of predisposing variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.eurPathVarCount.by.center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$Center),],aes(x = Center, y=NVAR))
p = p + geom_violin(aes(fill = Center),color=NA,alpha=0.2)
p = p + geom_jitter(aes(color = Center),width =0.2, height = 0, alpha=0.2,stroke=0,shape=19,size=0.5)
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of CDS variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.eurNVAR.by.center.all.pdf'
ggsave(file=fn, h=3,w = 3, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$Center),],aes(x = Center, y=Predispose152Genes_NVAR))
p = p + geom_violin(aes(fill = Center),color=NA,alpha=0.2)
p = p + geom_jitter(aes(color = Center),width =0.2, height = 0, alpha=0.2,stroke=0,shape=19,size=0.5)
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of variants in 152 predisposition genes")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.eurPredispose152Genes_NVAR.by.center.all.pdf'
ggsave(file=fn, h=3,w = 3, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$Center),],aes(x = Center, y=Predispose21Genes_NVAR))
p = p + geom_violin(aes(fill = Center),color=NA,alpha=0.2)
p = p + geom_jitter(aes(color = Center),width =0.2, height = 0, alpha=0.2,stroke=0,shape=19,size=0.5)
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of variants in 21 associated genes")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.eurPredispose21Genes_NVAR.by.center.all.pdf'
ggsave(file=fn, h=3,w = 3, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$Center),],aes(x = Center, y=pathVarCount))
p = p + geom_violin(aes(fill = Center),color=NA,alpha=0.2)
p = p + geom_jitter(aes(color = Center),width =0.2, height = 0.1, alpha=0.2,stroke=0,shape=19,size=0.5)
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                      geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of predisposing variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.eurPathVarCount.by.center.all.pdf'
ggsave(file=fn, h=3,w = 3, useDingbats=FALSE,limitsize=FALSE)

sd(istat_clin_coverage_pathCount_eur$NVAR,na.rm=T)/mean(istat_clin_coverage_pathCount_eur$NVAR,na.rm=T)
sd(istat_clin_coverage_pathCount_eur$Predispose152Genes_NVAR,na.rm=T)/mean(istat_clin_coverage_pathCount_eur$Predispose152Genes_NVAR,na.rm=T)
sd(istat_clin_coverage_pathCount_eur$Predispose21Genes_NVAR,na.rm=T)/mean(istat_clin_coverage_pathCount_eur$Predispose21Genes_NVAR,na.rm=T)

### density plots
p = ggplot(data=istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$Center),],aes(fill = Center, x=NVAR))
p = p + geom_density(aes(fill = Center),color=NA,alpha=0.6)
p = p + theme_nogrid() + theme(legend.position="bottom")
p = p + labs( x = "Count of CDS variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(x=0)
p
fn = 'out/PCA10389.eurNVAR.by.center.all.density.pdf'
ggsave(file=fn, h=3,w = 3, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$Center),],aes(fill = Center, x=Predispose152Genes_NVAR))
p = p + geom_density(aes(fill = Center),color=NA,alpha=0.6)
p = p + theme_nogrid() + theme(legend.position="bottom")
p = p + labs( x = "Count of variants in 152 predisposition genes")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(x=0)
p
fn = 'out/PCA10389.eurPredispose152Genes_NVAR.by.center.all.density.pdf'
ggsave(file=fn, h=3,w = 3, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$Center),],aes(fill = Center, x=Predispose21Genes_NVAR))
p = p + geom_density(aes(fill = Center),color=NA,alpha=0.6)
p = p + theme_nogrid() + theme(legend.position="bottom")
p = p + labs( x = "Count of variants in 21 associated genes")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(x=0)
p
fn = 'out/PCA10389.eurPredispose21Genes_NVAR.by.center.all.density.pdf'
ggsave(file=fn, h=3,w = 3, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$Center),],aes(fill = Center, x=pathVarCount))
p = p + geom_density(aes(fill = Center),color=NA,alpha=0.6)
p = p + theme_nogrid() + theme(legend.position="bottom")
p = p + labs( x = "Count of predisposing variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + expand_limits(x=0)
p
fn = 'out/PCA10389.eurPathVarCount.by.center.all.density.pdf'
ggsave(file=fn, h=3,w = 3, useDingbats=FALSE,limitsize=FALSE)

### percentage of variance across centers in European samples ###
NVAR_diff = vector("list")
G152_diff = vector("list")
Path_diff = vector("list")
for (cancer in unique(istat_clin_coverage_pathCount_eur$Cancer)){
  istat_clin_coverage_pathCount_eur_c = istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$Cancer) & istat_clin_coverage_pathCount_eur$Cancer==cancer,]
  if (length(unique(istat_clin_coverage_pathCount_eur_c$Center))>1){
    center1 = unique(istat_clin_coverage_pathCount_eur_c$Center)[1]
    center2 = unique(istat_clin_coverage_pathCount_eur_c$Center)[2]
    avg1= mean(istat_clin_coverage_pathCount_eur_c$NVAR[istat_clin_coverage_pathCount_eur_c$Center==center1])
    avg2= mean(istat_clin_coverage_pathCount_eur_c$NVAR[istat_clin_coverage_pathCount_eur_c$Center==center2])
    NVAR_diff[[cancer]] = 100*abs(avg1-avg2)/max(avg1,avg2)

    avg1= mean(istat_clin_coverage_pathCount_eur_c$Predispose152Genes_NVAR[istat_clin_coverage_pathCount_eur_c$Center==center1])
    avg2= mean(istat_clin_coverage_pathCount_eur_c$Predispose152Genes_NVAR[istat_clin_coverage_pathCount_eur_c$Center==center2])
    G152_diff[[cancer]] = 100*abs(avg1-avg2)/max(avg1,avg2)

    avg1= mean(istat_clin_coverage_pathCount_eur_c$pathVarCount[istat_clin_coverage_pathCount_eur_c$Center==center1])
    avg2= mean(istat_clin_coverage_pathCount_eur_c$pathVarCount[istat_clin_coverage_pathCount_eur_c$Center==center2])
    Path_diff[[cancer]] = 100*abs(avg1-avg2)/max(avg1,avg2)
  }
}
NVAR_diff_m = cbind(NVAR_diff)
G152_diff_m = cbind(G152_diff)
Path_diff_m = cbind(Path_diff)
all_diff = cbind(NVAR_diff_m,G152_diff_m,Path_diff_m)

### pathogenic variants ###
p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Ancestry) & istat_clin_coverage_pathCount$Ancestry != "mix",],aes(x = Ancestry, y=pathVarCount))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Ancestry))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Ancestry", y = "Count of predisposing variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.pathVarCount.by.ancestry.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Center),],aes(x = Center, y=pathVarCount))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Center))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Center", y = "Count of predisposing variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.pathVarCount.by.center.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_clin_coverage_pathCount[!is.na(istat_clin_coverage_pathCount$Analyte),],aes(x = Analyte, y=pathVarCount))
p = p + facet_grid(.~Cancer, drop = T)
p = p + geom_violin(alpha=0.8, aes(fill = Analyte))
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Analyte", y = "Count of predisposing variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) 
p
fn = 'out/PCA10389.pathVarCount.by.Analyte.pdf'
ggsave(file=fn, h=3,w = 15, useDingbats=FALSE,limitsize=FALSE)

##### the % carrier in each sub-categories #####
### assign less visual weight (maybe size) to smaller sample sizes ###
# Ancestry
all_counts = data.frame(table(istat_clin_coverage_pathCount$Ancestry,istat_clin_coverage_pathCount$Cancer))
counts = data.frame(table(istat_clin_coverage_pathCount$Ancestry,istat_clin_coverage_pathCount$Cancer,istat_clin_coverage_pathCount$pathVarCount!=0))
merged_counts = merge(all_counts,counts, by = c("Var1","Var2"))
colnames(merged_counts) = c("Ancestry","Cancer","Count","Carrier","CarrierCount")
merged_counts_carrier = merged_counts[merged_counts$Carrier==TRUE,]
merged_counts_carrier$PercentageCarrier = merged_counts_carrier$CarrierCount/merged_counts_carrier$Count

p = ggplot(data=merged_counts_carrier,aes(x = Cancer, y=PercentageCarrier))
p = p + facet_grid(.~Cancer, drop = T, space = "free_x", scale = "free_x")
p = p + geom_point(alpha=0.8, aes(color = Ancestry,size = Count))
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Cancer", y = "% Carrier of Predisposing Variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) + scale_size_area(max_size = 6) + ylim(0,0.25)
p
fn = 'out/PCA10389.pathVarCount.by.Ancestry.dotplot.pdf'
ggsave(file=fn, h=3,w = 10, useDingbats=FALSE,limitsize=FALSE)

# Analyte
all_counts = data.frame(table(istat_clin_coverage_pathCount$Analyte,istat_clin_coverage_pathCount$Cancer))
counts = data.frame(table(istat_clin_coverage_pathCount$Analyte,istat_clin_coverage_pathCount$Cancer,istat_clin_coverage_pathCount$pathVarCount!=0))
merged_counts = merge(all_counts,counts, by = c("Var1","Var2"))
colnames(merged_counts) = c("Analyte","Cancer","Count","Carrier","CarrierCount")
merged_counts_carrier = merged_counts[merged_counts$Carrier==TRUE,]
merged_counts_carrier$PercentageCarrier = merged_counts_carrier$CarrierCount/merged_counts_carrier$Count

p = ggplot(data=merged_counts_carrier,aes(x = Cancer, y=PercentageCarrier))
p = p + facet_grid(.~Cancer, drop = T, space = "free_x", scale = "free_x")
p = p + geom_point(alpha=0.8, aes(color = Analyte,size = Count))
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Cancer", y = "% Carrier of Predisposing Variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) + scale_size_area(max_size = 6)
p
fn = 'out/PCA10389.pathVarCount.by.Analyte.dotplot.pdf'
ggsave(file=fn, h=3,w = 10, useDingbats=FALSE,limitsize=FALSE)

# Center
all_counts = data.frame(table(istat_clin_coverage_pathCount$Center,istat_clin_coverage_pathCount$Cancer))
counts = data.frame(table(istat_clin_coverage_pathCount$Center,istat_clin_coverage_pathCount$Cancer,istat_clin_coverage_pathCount$pathVarCount!=0))
merged_counts = merge(all_counts,counts, by = c("Var1","Var2"))
colnames(merged_counts) = c("Center","Cancer","Count","Carrier","CarrierCount")
merged_counts_carrier = merged_counts[merged_counts$Carrier==TRUE,]
merged_counts_carrier$PercentageCarrier = merged_counts_carrier$CarrierCount/merged_counts_carrier$Count

p = ggplot(data=merged_counts_carrier,aes(x = Cancer, y=PercentageCarrier))
p = p + facet_grid(.~Cancer, drop = T, space = "free_x", scale = "free_x")
p = p + geom_point(alpha=0.8, aes(color = Center,size = Count))
p = p + theme_nogrid() + theme(legend.position="none")
p = p + labs(x = "Cancer", y = "% Carrier of Predisposing Variants")
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "bottom")
p = p + expand_limits(y=0) + scale_size_area(max_size = 6)
p
fn = 'out/PCA10389.pathVarCount.by.Center.dotplot.pdf'
ggsave(file=fn, h=3,w = 10, useDingbats=FALSE,limitsize=FALSE)

##### NVAR vs. pathVar###
p = ggplot(data=istat_clin_coverage_pathCount[istat_clin_coverage_pathCount$pathVarCount<3,],aes( x = factor(pathVarCount), y = NVAR ))
p = p + facet_grid(.~Cancer)
p = p + geom_violin()
p = p + geom_jitter(aes(color=Ancestry), width =0.2, height = 0, alpha=0.2,stroke=0,shape=19,size=0.5)
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.8,size=0.02, alpha=0.8)
p = p  + theme_nogrid()
p = p + labs(x = "Pathogenic Variant Counts", y = "Total Exonic Variant Counts")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14))
p = p + theme(legend.position = "bottom")
p = p + ylim(y=0,45000)
p
fn = 'out/PCA10389.pathVarCount_NVAR_violin.pdf'
ggsave(file=fn, h=3.5,w=15,useDingbats=FALSE,limitsize=FALSE)

# eur only
p = ggplot(data=istat_clin_coverage_pathCount_eur[!is.na(istat_clin_coverage_pathCount_eur$pathVarCount),],aes( x = factor(pathVarCount==0), y = NVAR ))
p = p + facet_grid(.~Cancer)
p = p + geom_violin(aes(fill = factor(pathVarCount==0)),color=NA,alpha=0.2)
p = p + geom_jitter(width =0.2, height = 0, alpha=0.2,stroke=0,shape=19,size=0.5)
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.8,size=0.02, alpha=0.8)
p = p  + theme_nogrid()
p = p + labs(x = "Pathogenic Variant Carrier", y = "Total Exonic Variant Counts")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14))
p = p + theme(legend.position = "none")
p = p + expand_limits(y=0)
p
fn = 'out/PCA10389.eur_pathVarCount_NVAR_violin.pdf'
ggsave(file=fn, h=3.5,w=15,useDingbats=FALSE,limitsize=FALSE)


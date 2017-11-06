##### label_onco_var.R #####
# Kuan-lin Huang @ WashU 2017 Oct
# find non-cancer pathogenic variant in the Pan8000 cohort

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/pathogenic_variants"
setwd(bdir)
source("../global_aes_out.R")


fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/process_files/germline/local/ChargedSample/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.lowAF.sele.tsv"
variants = read.table(sep="\t",header=T,file=fn, stringsAsFactors=FALSE, quote = "",fill=TRUE)

gene_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/reference_files/20160713_Rahman_KJ_KH_152_gene_table_list.txt"
predisposition_genes = as.vector(t(read.table(sep="\t",header=F,file=gene_fn, stringsAsFactors=FALSE, quote = "")))

samples_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/sampleQC/pca_table.20171019.wclin.tsv"
samples = read.table(sep="\t",header=T,file=samples_fn, stringsAsFactors=FALSE, quote = "")
# samples_brief = samples[,1:5]
# colnames(samples_brief)[2] = "Sample"

bam_location_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/sampleQC/pca_table_with_tumorBams.20171019.patched_28oct.tsv"
bam_location = read.table(sep="\t",header=T,file=bam_location_fn, stringsAsFactors=FALSE, quote = "", fill=T)
bam_location_s = bam_location[bam_location$sbarcode %in% samples$bcr_patient_barcode,]
bam_location_s_brief = bam_location_s[,c(1:5,8:10)]
colnames(bam_location_s_brief)[4] = "Sample"
colnames(bam_location_s_brief)[1:2] = c("bcr_patient_barcode","cancer")

# subset to PCA germline project sample
#variants$bcr_patient_barcode = substring(variants$Sample,1,12)
cat("Original count of variants: ",nrow(variants),"\n")
variants_pca = merge(variants,bam_location_s_brief,by="Sample",all=F)
cat("PCA sample sets of 10,467, # of variants left: ",nrow(variants_pca),"\n")

##### Add readcounts and readcount filter #####
# snp_rc_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/process_files/germline/local/ChargedSample/single_sub.rc_anno.tsv"
# snp_rc = read.table(sep="\t", header=T, file=snp_rc_fn, stringsAsFactors=FALSE, quote = "",fill=T)
# del_rc_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/process_files/germline/local/ChargedSample/del.rc_anno.tsv"
# del_rc = read.table(sep="\t", header=T, file=del_rc_fn , stringsAsFactors=FALSE, quote = "",fill=T)
# all_rc = rbind(snp_rc,del_rc)
all_rc_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/process_files/germline/local/ChargedSample/all.rc_anno"
all_rc = read.table(sep="\t", header=T, file=all_rc_fn , stringsAsFactors=FALSE, quote = "",fill=T)
all_rc[all_rc==-1]=NA # swap out -1 for NA to avoid confusion
all_rc$altBase = gsub("\\-","del",all_rc$altBase)
all_rc$altBase = gsub("\\+","ins",all_rc$altBase)# excel friendly format
variants_pca = merge(variants_pca,all_rc,all.x = T, by=c("Sample","HGVSg"))
variants_pca$normalVAF = variants_pca$normalAltCnt/variants_pca$normalDepth
variants_pca$tumorVAF = variants_pca$tumorAltCnt/variants_pca$tumorDepth

p = ggplot(data=variants_pca)
p = p + geom_point(aes(x=normalVAF,y=tumorVAF), alpha = 0.02,stroke=0)
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
p
fn = "out/normal_tumor_VAF_prefilter.pdf"
ggsave(file=fn, h=5,w=5,useDingbats=FALSE)

p = ggplot(data=variants_pca)
p = p + geom_point(aes(x=normalAltCnt,y=tumorAltCnt), alpha = 0.02,stroke=0)
p = p + theme_bw() + xlim(0,50) + ylim(0,50)
p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
p
fn = "out/normal_tumor_altcnt_prefilter.pdf"
ggsave(file=fn, h=5,w=5,useDingbats=FALSE)

variants_pca = variants_pca[is.na(variants_pca$normalAltCnt) | (variants_pca$normalAltCnt > 4),]
cat("Applying normal Alt readcount >= 5, # of variants left: ",nrow(variants_pca),"\n")
variants_pca = variants_pca[is.na(variants_pca$tumorAltCnt) | (variants_pca$tumorAltCnt > 4),]
cat("Applying tumor Alt readcount >= 5, # of variants left: ",nrow(variants_pca),"\n")

variants_pca = variants_pca[is.na(variants_pca$normalVAF) | (variants_pca$normalVAF >= 0.2),]
cat("Applying normal Alt readcount >= 0.2, # of variants left: ",nrow(variants_pca),"\n")
variants_pca = variants_pca[is.na(variants_pca$tumorVAF) | (variants_pca$tumorVAF >= 0.2),]
cat("Applying tumor Alt readcount >= 0.2, # of variants left: ",nrow(variants_pca),"\n")

p = ggplot(data=variants_pca)
p = p + geom_point(aes(x=normalAltCnt,y=tumorAltCnt), alpha = 0.02,stroke=0)
p = p + theme_bw() + xlim(0,50) + ylim(0,50)
p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
p
fn = "out/normal_tumor_altcnt_postfilter.pdf"
ggsave(file=fn, h=5,w=5,useDingbats=FALSE)

p = ggplot(data=variants_pca)
p = p + geom_point(aes(x=normalVAF,y=tumorVAF), alpha = 0.02,stroke=0)
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
p
fn = "out/normal_tumor_VAF_postfilter.pdf"
ggsave(file=fn, h=5,w=5,useDingbats=FALSE)

##### classify whether a variant is cancer-relevant #####
cancer_terms = c("tumor","cancer","neoplasia")

variants_pca$predisposition_gene = F
variants_pca$predisposition_gene[variants_pca$HUGO_Symbol %in% predisposition_genes] = TRUE
variants_pca$cancer_term_trait = FALSE
for (term in cancer_terms){
  variants_pca$cancer_term_trait[grep(term,tolower(variants_pca$ClinVar_Traits))] = TRUE
}
variants_pca$cancer_term_trait[grep("oma$",tolower(variants_pca$ClinVar_Traits))] = TRUE

table(variants_pca$predisposition_gene,variants_pca$cancer_term_trait)
variants_pca$cancer_related = F
variants_pca$cancer_related[variants_pca$predisposition_gene | variants_pca$cancer_term_trait] = T

table(variants_pca$ClinVar_Traits[variants_pca$cancer_term_trait])[table(variants_pca$ClinVar_Traits[variants_pca$cancer_term_trait])>3]

# variant frequency annotation
var_freq = data.frame(table(variants_pca$HGVSg))
colnames(var_freq) = c("HGVSg","Cohort_AC")
variants_pca_frq = merge(variants_pca,var_freq,by="HGVSg")
variants_pca_frq$Cohort_AF = variants_pca_frq$Cohort_AC/nrow(samples)
# dim(variants_pca_frq[variants_pca_frq$Cohort_AF > 0.01,])

tn = "charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.cleaned.expanded.AFcorrected.lowAF.sele.labeled.tsv"
write.table(variants_pca_frq, quote=F, sep="\t", file = tn, row.names = F)

# table(variants_pca_frq$Cohort_AF < 0.0005,variants_pca_frq$cancer_related) 
# sum(variants_pca_frq$Cohort_AF < 0.0005)
# table(variants_pca_frq$Cohort_AF < 0.01,variants_pca_frq$cancer_related)

# all truncations and missense with less cohort frequency
# table(variants_pca_frq$Consequence)
# table(variants_pca_frq$Variant_Classification)

# ##### check against previous variants #####
# pan8000_fn = "/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/functional_clinical_integration/tables/2016-07-18/2016-07-18_KH_pan8000_MAF0.05_germline_class_func_clin_somatic_PLP_only.txt"
# pan8000 = read.table(sep="\t",header=T,file=pan8000_fn, stringsAsFactors=FALSE, quote = "")
# pan8000$pan8000_sample_var = paste(pan8000$sample,pan8000$Start)  
# 
# pca_var = paste(variants_pca_frq$bcr_patient_barcode,variants_pca_frq$Start)
# length(pan8000$pan8000_sample_var)
# length(pan8000$pca_var)
# sum(pan8000$pan8000_sample_var %in% pca_var)
# sum(pca_var %in% pan8000$pan8000_sample_var)
# not_found = pan8000[!(pan8000$pan8000_sample_var %in% pca_var),]
# not_found = not_found[order(not_found$start),]
# not_found_trun =not_found[not_found$Variant_Classification=="frameshift_variant",]

##### rare frequency filter #####
variants_pca_frq_rare_1000G = variants_pca_frq[is.na(variants_pca_frq$GMAF) | variants_pca_frq$GMAF < 0.0005,]
variants_pca_frq_rare_1000G_ExAC = variants_pca_frq_rare_1000G[is.na(variants_pca_frq_rare_1000G$ExAC_adj_AF) | variants_pca_frq_rare_1000G$ExAC_adj_AF < 0.0005,]
# have to manually check these variants
unique(variants_pca_frq_rare_1000G_ExAC$HGVSg[variants_pca_frq_rare_1000G_ExAC$Cohort_AF > 0.001 & variants_pca_frq_rare_1000G_ExAC$Variant_Classification != "missense_variant" & variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF==0])
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual = NA
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="11:g.22646872_22646873delAG" ] = 	0.00008343
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="12:g.121432117_121432118insC"] = 0.001854
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="13:g.32914438_32914438delT"] =	0.0002651
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="15:g.91303325"] = 0.4311
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="16:g.2096132_2096133insG" ] =		0.0001745
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="16:g.89805421_89805427insAACAG"] =	0.002389
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="17:g.41209079_41209080insG" ] =		0.0001565
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="17:g.63532584_63532585insC" ] =	0.0005561
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="2:g.48030639_48030642insCC"] =		0.000008310
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="2:g.48033791"] =		0.001293
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="6:g.32007887G>T"] =		0.01042
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="12:g.12044534"] =	1 # low coverage
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="18:g.42456664_42456665insTCTC"] =	1 # low coverage
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="2:g.48030639_48030642insCC"] =	0.001753
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="2:g.241808314C>G"] =	0.1545
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="2:g.190670539_190670543insAAA"] =	0.2434
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="3:g.142274739_142274742insTT"] =	0.003416
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="9:g.36882049_36882050insG"] =	0.0001290
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="9:g.37020768_37020769insC" ] =		0.0001402
variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual[variants_pca_frq_rare_1000G_ExAC$HGVSg=="9:g.98209616_98209617insG" ] =		0.0001554

variants_pca_frq_rare_1000G_ExAC_pass = variants_pca_frq_rare_1000G_ExAC[is.na(variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual) | variants_pca_frq_rare_1000G_ExAC$ExAC_adj_AF_Manual < 0.0005,]
cat("Number of variants passing 1000G/ExAC 0.05% cut-off:",nrow(variants_pca_frq_rare_1000G_ExAC_pass),"\n")
## don't do cohort freq for now
# variants_pca_frq_rare_1000G_ExAC_pass_001 = variants_pca_frq_rare_1000G_ExAC_pass[variants_pca_frq_rare_1000G_ExAC_pass$Cohort_AF< 0.01,]
# tail(table(variants_pca_frq_rare_1000G_ExAC_pass_001$HGVSg)[order(table(variants_pca_frq_rare_1000G_ExAC_pass_001$HGVSg))])[2:5]

variants_pca_frq_rare_1000G_ExAC_pass_cancer = variants_pca_frq_rare_1000G_ExAC_pass[variants_pca_frq_rare_1000G_ExAC_pass$cancer_related,]
cat("Number of these variants that are cancer-relevant:",nrow(variants_pca_frq_rare_1000G_ExAC_pass_cancer),"\n")

# p = ggplot(data=variants_pca_frq_rare_1000G_ExAC_pass_cancer)
# p = p + geom_histogram(aes(x=Cohort_AF, fill=Variant_Classification),binwidth = 0.0005)
# #p = p + geom_text(aes(x=Cancer,y=Count, label=ifelse(Technology=="WXS",Count,NA)),color="black", vjust = -0.5, size=2)
# p = p + theme_bw()# + xlim(0,0.1) #+ ylim(0,0.1)
# p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
# p
# fn = "out/variants_pca_frq_rare_1000G_ExAC_pass_cancer_cohortAF.pdf"
# ggsave(file=fn, h=5,w=6,useDingbats=FALSE)


##### for manual review ##### 

## previous manual review results ##
prev_truncation_f = "/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/variant_files/201509_pan8000_0.05percentMAF/truncation/pan8000.truncation.var.1percentMAF.remove-low-coverage_0.05percentMAF_624_genes_rc.txt"
prev_pass_truncation_f = "/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/variant_files/201509_pan8000_0.05percentMAF/truncation/pan8000.truncation.var.1percentMAF.remove-low-coverage_0.05percentMAF_624_genes_rc_oldManualRev_manualRev_cleaned_germline_noNonstop.txt"
prev_truncation = read.table(sep="\t",header=F,file=prev_truncation_f, stringsAsFactors=FALSE, quote = "",fill=TRUE)
prev_pass_truncation = read.table(sep="\t",header=T,file=prev_pass_truncation_f, stringsAsFactors=FALSE, quote = "",fill=TRUE)
all = paste(prev_truncation$V22,prev_truncation$V2,sep=":")
pass = paste(prev_pass_truncation$sample,prev_pass_truncation$start,sep=":")
fail = all[!(all %in% pass)]

variants_pca_frq_rare_1000G_ExAC_pass_cancer$previous_manual_review = NA
variants_pca_frq_rare_1000G_ExAC_pass_cancer$previous_manual_review[paste(variants_pca_frq_rare_1000G_ExAC_pass_cancer$bcr_patient_barcode,variants_pca_frq_rare_1000G_ExAC_pass_cancer$Start,sep=":") %in% pass] = "Pass"
variants_pca_frq_rare_1000G_ExAC_pass_cancer$previous_manual_review[paste(variants_pca_frq_rare_1000G_ExAC_pass_cancer$bcr_patient_barcode,variants_pca_frq_rare_1000G_ExAC_pass_cancer$Start,sep=":") %in% fail] = "Fail"


tn = "charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.cleaned.expanded.AFcorrected.lowAF.sele.labeled.rare.cancer.tsv"
write.table(variants_pca_frq_rare_1000G_ExAC_pass_cancer, quote=F, sep="\t", file = tn, row.names = F)

##### combine with manual review result #####
RJ_manual_review_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/pathogenic_variants/Germline_PCA_RJ/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.cleaned.expanded.AFcorrected.lowAF.sele.labeled.rare.cancer.forReview.RJ.txt"
RJ_manual_review = read.table(sep="\t",header=T,file=RJ_manual_review_f, stringsAsFactors=FALSE, quote = "",fill=TRUE)
RJ_manual_review_brief =  RJ_manual_review[,c("Sample","HGVSg","reyka_manual_review")]

final_var = merge(variants_pca_frq_rare_1000G_ExAC_pass_cancer,RJ_manual_review_brief, by = c("Sample","HGVSg"))
final_var_pass = final_var[final_var$reyka_manual_review %in% c("pass" , "pass - nearby mutation" ),]
cat("Reyka's manual review result count for cancer variants:","\n")
table(final_var_pass$reyka_manual_review)

tn = "charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.cleaned.expanded.AFcorrected.lowAF.sele.labeled.rare.cancer.pass.tsv"
write.table(final_var_pass, quote=F, sep="\t", file = tn, row.names = F)
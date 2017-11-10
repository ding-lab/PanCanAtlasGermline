##### combine_data.R #####
# Kuan-lin Huang @ WashU 2017 Oct
# combine files to get all attributes for germline variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/data_integration"
setwd(bdir)
library(reshape2)
source("../global_aes_out.R")

##### variant files #####
var_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/pathogenic_variants/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.cleaned.expanded.AFcorrected.lowAF.sele.labeled.rare.cancer.pass.tsv"
var = read.table(sep="\t",header=T,quote="",file=var_fn, stringsAsFactors=FALSE)
var$HGVSp_short = gsub(".*:","",var$HGVSp)
var$HGVSp_short[var$HGVSp_short=="p."] = NA

cat("Starting out with matrix size dimension:","\n")
dim(var)
# process to get some attributes of interest

##### merge association results #####
assoc_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/association_test/assoc_results/PCA.pathVar.ExAC.r1.sites.vep.biallelic.combine.fisher.anno.v2.tsv"
assoc_f = read.table(sep="\t",header=T,file=assoc_fn, fill=T,stringsAsFactors=FALSE)
assoc_f = assoc_f[assoc_f$ExAC_AC<61,]
assoc_f = assoc_f[order(assoc_f$P, decreasing = F),]
assoc_f = assoc_f[!duplicated(assoc_f$HGVSg),]
assoc_f_brief = assoc_f[,c("HGVSg","Var","OR","P")]

# var_split = strsplit(assoc_f_brief$Var,':')
# assoc_f_brief$pos = colsplit(string=assoc_f_brief$Var, pattern=":", names=c("CHR", "Start","ID","REF","ALT"))[,2]
# assoc_f_brief$ALT_clean = colsplit(string=assoc_f_brief$Var, pattern=":", names=c("CHR", "Start","ID","REF","ALT"))[,5]
# 
colnames(assoc_f_brief) = paste("ExAC_assoc",colnames(assoc_f_brief),sep="_")
colnames(assoc_f_brief)[1] = "HGVSg"
# colnames(assoc_f_brief)[5] = "Start"
# colnames(assoc_f_brief)[5] ="ALT_clean"
# var_assoc = merge(var,assoc_f_brief,by=c("Start","ALT_clean"),all.x=T,all.y=F)

var_assoc = merge(var,assoc_f_brief,by="HGVSg",all.x=T,all.y=F)
cat("After merging association:","\n")
dim(var_assoc)

##### merge clustering results #####
#hotspot_fn = ""

##### merge expression #####
exp_score_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/expression_effect/2016-06-21_KH_pancan_exp_log2RSEM_all_uniq.tsv.gz"
exp_score = read.table(sep="\t",header=F,file=gzfile(exp_score_f), stringsAsFactors=FALSE)
exp_quantile_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/expression_effect/2016-06-21_KH_pancan_exp_quantile_all_uniq.tsv.gz"
exp_quantile = read.table(sep="\t",header=F,file=gzfile(exp_quantile_f), stringsAsFactors=FALSE)

colnames(exp_score)=c("HUGO_Symbol","bcr_patient_barcode","log2RSEM","cancer")
colnames(exp_quantile)=c("HUGO_Symbol","bcr_patient_barcode","expressionQuantile","cancer")

var_assoc_exp = merge(var_assoc,exp_score,by=c("HUGO_Symbol","bcr_patient_barcode","cancer"),all.x=T,all.y=F)
var_assoc_exp_q = merge(var_assoc_exp,exp_quantile,by=c("HUGO_Symbol","bcr_patient_barcode","cancer"),all.x=T,all.y=F)
cat("After merging expression:","\n")
dim(var_assoc_exp_q)

### merge the respective transcript length, and thus their relative position in the transcript ###
var_assoc_exp_q$transcript_name = gsub(":.*","",var_assoc_exp_q$HGVSc)
transcript_l_f = "/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/reference_files/gene_transcript_length"
transcript_l = read.table(header=F, quote = "", sep="\t", fill =T, file = transcript_l_f , stringsAsFactors=FALSE)
colnames(transcript_l) = c("HUGO_Symbol","transcript_name","transcript length")
var_assoc_exp_q_t = merge(var_assoc_exp_q, transcript_l, by = c("HUGO_Symbol","transcript_name"), all.x=T,all.y=F)
cat("After merging transcript length:","\n")
dim(var_assoc_exp_q_t)

##### merge LOH #####
## read in AI files and pre-process
ai_missense_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/LOH/LOH.missense.out.tsv"
ai_missense = read.table(sep="\t",header=T,file=ai_missense_fn, stringsAsFactors=FALSE)
ai_truncation_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/LOH/LOH.truncation.out.tsv"
ai_truncation = read.table(sep="\t",header=T,file=ai_truncation_fn, stringsAsFactors=FALSE)
ai_truncation$binary_type = "Truncation"
ai_missense$binary_type = "Missense"
ai_f = rbind(ai_missense,ai_truncation)
ai_f$TumorByNormalVAF = ai_f$TumorVAF/ai_f$NormalVAF
ai_f_brief = ai_f[,c(3,14:21)]
colnames(ai_f_brief)[1] = "Start"
colnames(ai_f_brief)[3:7] = paste("LOH",colnames(ai_f_brief)[3:7],sep="_")
# merge AI files
var_assoc_exp_q_t_ai = merge(var_assoc_exp_q_t,ai_f_brief,by=c("Sample","Start"),all.x=T,all.y=F)
cat("After merging LOH:","\n")
dim(var_assoc_exp_q_t_ai)

##### merge clinical file #####
clin_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/PanCan_ClinicalData_V4_wAIM.txt"
clin = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)
#colnames(clin) = c("sample","age_at_onset","ethnicity","cancer")
clin = clin[,-c(38:39)]
colnames(clin)[colnames(clin) %in% colnames(var_assoc_exp_q_t_ai)]
colnames(clin)[2] = "cancer"
colnames(clin)[38] = "AIM_ethnicity"

germ_clin = merge(var_assoc_exp_q_t_ai,clin, by =c("bcr_patient_barcode","cancer"), all.x=T,all.y=F)
cat("After merging clinical data:","\n")
dim(germ_clin)

##### merge PCGP information #####
truncation_f_PCGP = paste("/Users/khuang/Box Sync/PhD/germline/pan8000_germline_clinical/CharGer_analysis/pan8000/05_07_2016/05_07_2016_pan8000_MAF0.05_truncation_charger_PCGP.txt",sep="")
missense_f_PCGP = paste("/Users/khuang/Box Sync/PhD/germline/pan8000_germline_clinical/CharGer_analysis/pan8000/05_07_2016/05_07_2016_pan8000_MAF0.05_missense_charger_PCGP.txt",sep="")

trun_PCGP = read.table(header=T, quote = "", sep="\t", fill =T, file = truncation_f_PCGP)
miss_PCGP = read.table(header=T, quote = "", sep="\t", fill =T, file = missense_f_PCGP)

germ_clin$PCGP = FALSE
germ_clin[germ_clin$HGVSg %in% trun_PCGP$HGVSg,]$PCGP = TRUE 
germ_clin[germ_clin$HGVSg %in% miss_PCGP$HGVSg,]$PCGP = TRUE
cat("After merging PCGP:","\n")
dim(germ_clin)
##### merge somatic information #####
somatic_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/somatic/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz"
somatic = read.table(header=T, quote = "", sep="\t", file = gzfile(somatic_f), stringsAsFactors=FALSE)
colnames(somatic) = c("HUGO_Symbol","Somatic_Variant_Classification","TumorSample","Somatic_HGVSp") #Hugo_Symbol     Variant_Classification  Tumor_Sample_Barcode    HGVSp_Short
somatic$bcr_patient_barcode = gsub("(^TCGA-[A-Z0-9][A-Z0-9]-[A-Z0-9][A-Z0-9][A-Z0-9][A-Z0-9])-.*","\\1",somatic$TumorSample)
somatic = somatic[somatic$Somatic_Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                                                                "Missense_Mutation","Nonsense_Mutation","Splice_Site"),]
somatic_class_agg = aggregate(somatic[c('Somatic_Variant_Classification','Somatic_HGVSp')], by=somatic[c('TumorSample',"HUGO_Symbol","bcr_patient_barcode")], paste, collapse = ",")
germ_clin_so = merge(germ_clin,somatic_class_agg, by =c("bcr_patient_barcode","HUGO_Symbol"), all.x=T, all.y=F)
cat("After merging somatic mutation:","\n")
dim(germ_clin_so)

# find number of co-localizing somatic mutations
germ_clin_so$var = paste(germ_clin_so$HUGO_Symbol,germ_clin_so$HGVSp_short)
somatic_colocalized_var = somatic[paste(somatic$HUGO_Symbol,somatic$Somatic_HGVSp) %in% germ_clin_so$var,]
somatic_occurence = data.frame(table(paste(somatic_colocalized_var$HUGO_Symbol,somatic_colocalized_var$Somatic_HGVSp)))
colnames(somatic_occurence) = c("var","colocalized_somatic_mutation_count")
germ_clin_so = merge(germ_clin_so,somatic_occurence, by ="var", all.x=T, all.y=F)
germ_clin_so$colocalized_somatic_mutation_count[is.na(germ_clin_so$colocalized_somatic_mutation_count)] = 0
germ_clin_so = germ_clin_so[,-which(colnames(germ_clin_so)=="var")]

##### write output files #####
cat("Ending with matrix size dimension:","\n")
dim(germ_clin_so)

fn = "out/PCA_pathVar_integrated.tsv"
write.table(germ_clin_so, quote=F, sep="\t", file = fn, row.names = F)

## additional output files to focus on specific cases ##

# cases with more than 1 variants
duplicated_sample = germ_clin_so$bcr_patient_barcode[duplicated(germ_clin_so$bcr_patient_barcode)]
germ_clin_so_ds = germ_clin_so[germ_clin_so$bcr_patient_barcode %in% duplicated_sample,]
fn = "out/PCA_pathVar_integrated_multi-var_carrier.tsv"
write.table(germ_clin_so_ds, quote=F, sep="\t", file = fn, row.names = F)

# variants found in more than 1 samples
duplicated_HGVSg = germ_clin_so$HGVSg[duplicated(germ_clin_so$HGVSg)]
germ_clin_so_dv = germ_clin_so[germ_clin_so$HGVSg %in% duplicated_HGVSg,]
fn = "out/PCA_pathVar_integrated_multi-carrier_var.tsv"
write.table(germ_clin_so_dv, quote=F, sep="\t", file = fn, row.names = F)

# potential compound hets
germ_clin_so_ComHet = germ_clin_so[!is.na(germ_clin_so$Somatic_Variant_Classification),]
fn = "out/PCA_pathVar_integrated_comHet.tsv"
write.table(germ_clin_so_ComHet, quote=F, sep="\t", file = fn, row.names = F)
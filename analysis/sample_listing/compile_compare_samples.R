##### compile_compare_samples.R #####
# Kuan-lin Huang @ WashU 2017 Oct.
# compile sample lists and compare PCA germline samples to MC3 sample set

setwd("/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/sample_listing")
source("../global_aes_out.R")
system("mkdir out")
library(UpSetR)

# read and preprocess files
ISB_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/ISB-manifest/CGHub_legacyGDC_DNA_bams.csv.gz"
ISB = read.table(header=T, sep = ',',file=ISB_fn)
ISB_normal = ISB[substr(ISB$sample_barcode,14,14)==1,]
ISB_normal$inferred_ctype = gsub(".*TCGA.([A-Z]+)/.*","\\1",toupper(ISB_normal$bam_gcs_url))

ISB_normal_WXS = ISB_normal[ISB_normal$experimental_strategy == "WXS",]
ISB_normal_WXS_uniq = ISB_normal_WXS[!duplicated(ISB_normal_WXS$case_barcode),]

ISB_normal_WGS = ISB_normal[ISB_normal$experimental_strategy == "WGS",]
ISB_normal_WGS_uniq = ISB_normal_WGS[!duplicated(ISB_normal_WGS$case_barcode),]

WXS_data = data.frame(table(ISB_normal_WXS_uniq$inferred_ctype))
WGS_data = data.frame(table(ISB_normal_WGS_uniq$inferred_ctype))
WXS_data$Technology = "WXS"
WGS_data$Technology = "WGS"
ISB_data = rbind(WXS_data,WGS_data)#merge(WXS_data,WGS_data,by="Var1")
colnames(ISB_data)[1:2] = c("Cancer","Count")
ISB_data$Cancer = as.character(ISB_data$Cancer)
#ISB_data$yPos = ISB_data$Count

p = ggplot(data=ISB_data)
p = p + geom_bar(aes(x=Cancer, y=Count, fill=Technology),stat = "identity")
p = p + geom_text(aes(x=Cancer,y=Count, label=ifelse(Technology=="WXS",Count,NA)),color="black", vjust = -0.5, size=2)
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
p
fn = "out/cancer_datatype_distribution.pdf"
ggsave(file=fn, h=5,w=6,useDingbats=FALSE)

completed_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/sampleQC/GCS_listing.02jun2016.tcga.all.DNA.WXS.normals.tsv.noDupBam.fastqcMaxFail_2_exclKmerCont.minGoodCvg20x.analysisID_w_barcodes"
completed = read.table(header=F, sep = '\t',file=completed_fn )
colnames(completed) = c("uuid","bcr_sample_barcode")
completed$bcr_patient_barcode = substring(completed$bcr_sample_barcode,1,12)

# clinical data quick look:
# 10957 /Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/clinical_PANCAN_patient_with_followup.tsv
# 11160 /Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/PanCan_ClinicalData_V4_20170428.txt
# 14505 /Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/all.clin.merged.picked.txt

PCA_clin_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/PanCan_ClinicalData_V4_20170428.txt"
PCA_clin = read.table(header=T, sep = '\t',file=PCA_clin_fn,fill=T )
PCA_clin = PCA_clin[,1:6]
MC3_clin_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/clinical_PANCAN_patient_with_followup.tsv"
MC3_clin = read.table(header=T, sep = '\t',file=MC3_clin_fn,quote="" )
MC3_clin = MC3_clin[,1:6]

##### compare our WGS samples to MC3 and PCA samples #####
all_cases_data = data.frame(unique(c(as.character(ISB_normal_WXS_uniq$case_barcode),as.character(PCA_clin$bcr_patient_barcode),as.character(MC3_clin$bcr_patient_barcode),as.character(completed$bcr_patient_barcode))))
colnames(all_cases_data) = "bcr_patient_barcode"
all_cases_data$inISB = all_cases_data$bcr_patient_barcode %in% as.character(ISB_normal_WXS_uniq$case_barcode)
all_cases_data$inPCA = all_cases_data$bcr_patient_barcode %in% as.character(PCA_clin$bcr_patient_barcode)
all_cases_data$inMC3 = all_cases_data$bcr_patient_barcode %in% as.character(MC3_clin$bcr_patient_barcode)
all_cases_data$completed9401 = all_cases_data$bcr_patient_barcode %in% as.character(completed$bcr_patient_barcode)
all_cases_data[!all_cases_data] = 0
all_cases_data[all_cases_data] = 1
fn = "out/201710_sample_upset.pdf"
pdf(fn, useDingbats = F)
upset(all_cases_data, sets = c("inISB", "inPCA", "inMC3","completed9401"),order.by = "freq")
dev.off()

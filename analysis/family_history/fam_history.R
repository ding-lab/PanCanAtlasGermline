##### fam_history.R #####
# Kuan-lin Huang @ WashU 2017 updated Nov.

setwd("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/family_history/")
source("../global_aes_out.R")
source("../dependency_files.R")

##### individual level stats #####

fileName = "clinical_PANCAN_patient_cancerhistory.111817.tsv"
famHist = read.table(header=TRUE, sep="\t", file=fileName, fill=T, quote="",stringsAsFactors = F)
str(famHist)

fileNameB = "clinical_PANCAN_patient_cancerhistory.2col.111817.tsv"
famHistB = read.table(header=TRUE, sep="\t", file=fileNameB, fill=T, quote="",stringsAsFactors = F)
str(famHistB)

fam_hist_samples = famHistB$bcr_patient_barcode[famHistB[,2]=="Yes"]

famHistPos = famHist[famHist$bcr_patient_barcode %in% fam_hist_samples,]

for (col in colnames(famHistPos)[2:14]){
  print(col)
  print(table(famHistPos[,col]))
}

fam_var_merge = merge(pathVarP, famHistPos, by="bcr_patient_barcode")

for (col in colnames(fam_var_merge)[182:194]){
  print(col)
  print(table(fam_var_merge[,col]))
}

fam_var_merge$fam_tag = "Other"
fam_var_merge$fam_tag[grep("reast",fam_var_merge$relative_family_cancer_hx_text)] = "Breast cancer"
fam_var_merge$fam_tag[grep("rostate",fam_var_merge$relative_family_cancer_hx_text)] = "Prostate cancer"
fam_var_merge$fam_tag[fam_var_merge$family_history_of_stomach_cancer=="YES"] = "Stomach cancer"
fam_var_merge$fam_tag[fam_var_merge$number_of_first_degree_relatives_with_cancer_diagnosis==1] = "First degree relatives"
table(fam_var_merge$fam_tag)
table(fam_var_merge$HUGO_Symbol)
table(fam_var_merge$HUGO_Symbol,fam_var_merge$fam_tag)
table(fam_var_merge$cancer)

fam_var_merge_noNA = fam_var_merge[!is.na(fam_var_merge$HGVSp_short),]
p = ggplot(fam_var_merge_noNA,aes(x = HGVSp_short, fill = cancer))
p = p + facet_grid(.~HUGO_Symbol, scale="free", space="free")
p = p + geom_bar() + theme_bw() + theme_nogrid() + getPCACancerFill()
p = p + labs(x = "Variant", y="Number of carriers with family history")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + theme(legend.position = "top")
p
fn = 'out/family_var_count.pdf'
ggsave(file=fn, w=15, h =6 ,useDingbats=FALSE)

dup_var = fam_var_merge[duplicated(fam_var_merge$HGVSp_short),c("HUGO_Symbol","HGVSp_short","cancer")]
cat("Duplicated variants","\n")
dup_var 
dup_genes = fam_var_merge$HUGO_Symbol[duplicated(fam_var_merge$HUGO_Symbol)]
fam_var_merge_g = fam_var_merge[fam_var_merge$HUGO_Symbol %in% dup_genes,]

p = ggplot(fam_var_merge_g,aes(y = HUGO_Symbol, x = cancer, color=fam_tag))
#p = p + facet_grid(Gene_Classification~., scale="free", space="free")
p = p + geom_jitter(alpha=0.5,size=1.5, height = 0.2,width = 0.22,shape=16, stroke=0)
p = p + geom_text_repel(aes(label=as.character(HGVSp_short),size=3,angle=0,vjust=1.5))
p = p + theme_bw()
#p = p + geom_hline(yintercept = -log10(0.05),alpha=0.3)
#p = p + xlim(0,6)
#p = p + getVarColorScale()
p = p + labs( x="TCGA case cancer type", y = "Gene") + scale_colour_discrete(name = "Case family history")
p = p + theme(axis.text.x = element_text(colour="black", size=14, angle=90, vjust = 0.5),axis.text.y = element_text(colour="black", size=14))
#p = p + coord_equal()
p
fn = "out/fam_history_var.pdf"
ggsave(fn, h=6, w = 8,useDingbat=F)

p = ggplot(fam_var_merge_g,aes(y = HUGO_Symbol, x = cancer, color=fam_tag))
#p = p + facet_grid(Gene_Classification~., scale="free", space="free")
p = p + geom_jitter(alpha=0.4,size=3, height = 0.1,width = 0.1,shape=16, stroke=0)
p = p + geom_text(aes(label=ifelse(as.character(fam_var_merge_g$HGVSp_short) %in% dup_var$HGVSp_short,HGVSp_short,NA),size=3,angle=0,vjust=1.5))
p = p + theme_bw()
#p = p + geom_hline(yintercept = -log10(0.05),alpha=0.3)
#p = p + xlim(0,6)
#p = p + getVarColorScale()
p = p + labs( x="TCGA case cancer type", y = "Gene") + scale_colour_discrete(name = "Case family history")
p = p + theme(axis.text.x = element_text(colour="black", size=14, angle=90, vjust = 0.5),axis.text.y = element_text(colour="black", size=14))
#p = p + coord_equal()
p
fn = "out/fam_history_var_dupVarOnly.pdf"
ggsave(fn, h=6, w = 8,useDingbat=F)

p = ggplot(fam_var_merge_g,aes(x = HUGO_Symbol, fill = fam_tag))
p = p + facet_grid(.~cancer, scale="free", space="free")
p = p + geom_bar() + theme_bw() + theme_nogrid() #+ getPCACancerFill()
p = p + labs(x = "Gene", y="Number of carriers with family history") + scale_fill_discrete(name = "Case family history")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + theme(legend.position = "top")
p
fn = "out/fam_history_var_famtag.pdf"
ggsave(fn, h=4, w = 8,useDingbat=F)

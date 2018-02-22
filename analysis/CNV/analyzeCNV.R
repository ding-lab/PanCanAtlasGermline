##### analyzeCNV.R #####
# Kuan-lin Huang @ WashU 2018 Feb
# conduct association of CNV with Expression

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/CNV"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")
### MAIN ###

###  using array data ###
# read input files
cnv_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/CNV/data/int_xhmm10389_array_by_gene_expq_cleaned.txt"
cnv = read.table(header=T, quote = "", sep=" ", fill =T, file = cnv_f, stringsAsFactors=FALSE)

burden_f ="/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/burden_assoc/out/cancer_vs_otherCancer_pathogenic_variants_burden.tsv"
burden = read.table(header=T, quote = "", sep="\t", fill =T, file = burden_f, stringsAsFactors=FALSE)
burden_sig = burden[burden$FDR<0.15,]

cnv$associated = FALSE
cnv$associated[paste(cnv$cancer,cnv$gene) %in% paste(burden_sig$Cancer,burden_sig$Gene)] = T
# cnv$associated[cnv$gene %in% burden_sig$Gene]= T

p = ggplot(cnv[cnv$associated,],aes(x=depth_xhmm, y =qt,color=cancer))
p = p + geom_point(alpha=0.3,shape=19)#, scale = "free", space = "free", drop=T)
p = p + geom_text_repel(aes(label=paste(gene)))
p + theme_bw() + geom_vline(xintercept = 0, alpha=0.5) + getPCACancerColor()+ ylim(0,1)
fn = 'out/overlapped_CNV_qt_assocGene.pdf'
ggsave(file=fn, useDingbats=FALSE)

###  using array data ###
# read input files
array_cnv_f = "data/array.cnv.filtered.rare.seg.refseq.expq.annotated.rahman.txt"
array_cnv = read.table(header=T, quote = "", sep="\t", fill =T, file = array_cnv_f, stringsAsFactors=FALSE)
colnames(array_cnv)[2]= "bcr_patient_barcode"
array_cnv$CNV = "None"
array_cnv$CNV[array_cnv$seg_mean < -0.1] = "DEL"
array_cnv$CNV[array_cnv$seg_mean > 0.1] = "DUP"
array_cnv = array_cnv[!duplicated(paste(array_cnv$bcr_patient_barcode,array_cnv$gene)),]

array_cnv$associated = FALSE
#array_cnv$associated[array_cnv$gene %in% burden_sig$Gene]= "Off target" #off target doesn't work that well 
array_cnv$associated[paste(array_cnv$cancer,array_cnv$gene) %in% paste(burden_sig$Cancer,burden_sig$Gene)] = TRUE

p = ggplot(array_cnv[array_cnv$associated !=F,],aes(x=seg_mean, y =qt,color=cancer))
p = p + geom_point(alpha=0.1,shape=19)#, scale = "free", space = "free", drop=T)
p = p + geom_text_repel(aes(label=paste(gene)))
p + theme_bw() + geom_vline(xintercept = 0, alpha=0.5) + getPCACancerColor()+ ylim(0,1)
fn = 'out/SNParray_CNV_qt_assocGene.pdf'
ggsave(file=fn, useDingbats=FALSE)

p = ggplot(array_cnv[array_cnv$CNV == "Deletion",],aes(x=associated, y =qt, fill = associated))
#p = p + geom_jitter(aes(color = associated), alpha=0.1,shape=19, height=0)#, scale = "free", space = "free", drop=T)
p = p + geom_violin(alpha=0.2)
p = p + geom_text_repel(aes(label=ifelse(associated,paste(cancer,gene),NA)))
p = p + ggtitle("Copy number deletions (SNP array)")
p = p + theme_bw() + geom_vline(xintercept = 0, alpha=0.5)+ ylim(0,1) + labs(x="Associated cancer-gene pair", y="Expression quantile")
p + theme(legend.position = "none")
fn = 'out/SNParray_CNV_qt_deletion.pdf'
ggsave(file=fn, useDingbats=FALSE)

### xhmm ###
xhmm_cnv_f = "data/xhmm.auto.rare.q60.rm3sd.expq.annotated.rahman.txt"
xhmm_cnv = read.table(header=T, quote = "", sep="\t", fill =T, file = xhmm_cnv_f, stringsAsFactors=FALSE)
colnames(xhmm_cnv)[1]= "bcr_patient_barcode"
xhmm_cnv = xhmm_cnv[!duplicated(paste(xhmm_cnv$bcr_patient_barcode,xhmm_cnv$gene)),]

xhmm_cnv$associated = FALSE
#array_cnv$associated[array_cnv$gene %in% burden_sig$Gene]= "Off target" #off target doesn't work that well 
xhmm_cnv$associated[paste(xhmm_cnv$cancer,xhmm_cnv$gene) %in% paste(burden_sig$Cancer,burden_sig$Gene)] = TRUE

p = ggplot(xhmm_cnv[xhmm_cnv$associated !=F,],aes(x=MEAN_RD, y =qt,color=cancer))
p = p + geom_point(alpha=0.1,shape=19)#, scale = "free", space = "free", drop=T)
p = p + geom_text_repel(aes(label=paste(cancer,gene)))
p + theme_bw() + geom_vline(xintercept = 0, alpha=0.5)+ getPCACancerColor() + ylim(0,1)
fn = 'out/xhmm_CNV_qt_assocGene.pdf'
ggsave(file=fn, useDingbats=FALSE)

##### combined counts #####
colnames(array_cnv)[colnames(array_cnv)=="seg_mean"] = "CNV_value"
colnames(xhmm_cnv)[colnames(xhmm_cnv)=="MEAN_RD"] = "CNV_value"
overlappled_colnames = colnames(array_cnv)[colnames(array_cnv) %in% colnames(xhmm_cnv)]
array_cnv_sele = array_cnv[array_cnv$associated !=F & array_cnv$CNV == "DEL" & array_cnv$gene %in% all_TSGs,overlappled_colnames]
array_cnv_sele$detection = "SNP array"
xhmm_cnv_sele = xhmm_cnv[xhmm_cnv$associated !=F & xhmm_cnv$CNV == "DEL" & xhmm_cnv$gene %in% all_TSGs,overlappled_colnames]
xhmm_cnv_sele$detection = "Whole-exome sequencing"
cnv_sele = rbind(array_cnv_sele,xhmm_cnv_sele)

tn = "out/cnv_ontarget_TSG_dele.txt"
write.table(cnv_sele, quote=F, sep="\t", file = tn, row.names = F)

cnv_sele_count = data.frame(table(cnv_sele$gene,cnv_sele$cancer,cnv_sele$detection))
colnames(cnv_sele_count) = c('gene',"cancer","detection","count")
cnv_sele_count= cnv_sele_count[cnv_sele_count$count!=0,]
p = ggplot(cnv_sele_count,aes(x=gene, y = cancer))
p = p + facet_grid(.~detection, space="free",scale="free")
#p = p + geom_jitter(aes(color = associated), alpha=0.1,shape=19, height=0)#, scale = "free", space = "free", drop=T)
p = p + geom_point(aes(size=count)) #+ scale_size_area()
#p = p + ggtitle("Copy number deletions")
p = p + theme_bw() + labs(x="Associated cancer-gene pair", y="Expression quantile")+ getPCACancerFill()
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))
p #+ theme(legend.position = "none")
fn = 'out/both_CNV_qt_deletion.pdf'
ggsave(file=fn,h=3,w=5, useDingbats=FALSE)

p = ggplot(cnv_sele,aes(x=CNV_value, y =qt, color=cancer))
p = p + facet_grid(.~detection,scale="free")
#p = p + geom_jitter(aes(color = associated), alpha=0.1,shape=19, height=0)#, scale = "free", space = "free", drop=T)
p = p + geom_jitter(height=0,alpha=0.5)
p = p + geom_text_repel(aes(label=paste(cancer,gene))) + getPCACancerColor()
p = p + theme_bw() + labs(x="CNV value", y="Expression quantile") + expand_limits(x=0)
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))
p = p + geom_vline(xintercept = 0, alpha=0.5)
p = p + geom_hline(yintercept = 0.25, alpha=0.5)
p #+ theme(legend.position = "none")
fn = 'out/both_CNV_qt_deletion.pdf'
ggsave(file=fn,h=5,w=7, useDingbats=FALSE)



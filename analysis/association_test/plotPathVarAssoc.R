##### plotPathVarAssoc.R #####
# Kuan-lin Huang @ WashU 201711
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/association_test"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

# attach NFE only association
assoc_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/association_test/assoc_results/ExAC.r1.sites.vep.biallelic.combine.fisher.anno.NFE.152gene.tsv"
assoc_f = read.table(sep="\t",header=T,file=assoc_fn, fill=T,stringsAsFactors=FALSE)
assoc_f = assoc_f[assoc_f$ExAC_AC<61,]
assoc_f = assoc_f[order(assoc_f$P, decreasing = F),]
assoc_f = assoc_f[!duplicated(assoc_f$Var),]
assoc_f_brief = assoc_f

colnames(assoc_f_brief) = paste("ExAC_NFEassoc",colnames(assoc_f_brief),sep="_")
colnames(assoc_f_brief)[12] = "HGVSp.1"

pathVar_assoc = merge(pathVar,assoc_f_brief,by="HGVSp.1",all.x=T,all.y=F)
pathVar_assoc$ExAC_NFEassoc_OR[!is.na(pathVar_assoc$ExAC_NFEassoc_TCGA_AC) & pathVar_assoc$Cohort_AC < pathVar_assoc$ExAC_NFEassoc_TCGA_AC] = NA
pathVar_assoc$ExAC_NFEassoc_P[!is.na(pathVar_assoc$ExAC_NFEassoc_TCGA_AC) & pathVar_assoc$Cohort_AC < pathVar_assoc$ExAC_NFEassoc_TCGA_AC] = NA

pathVar_assoc = pathVar_assoc[!duplicated(pathVar_assoc$HGVSp.1),]
pathVar_assoc_sig = pathVar_assoc[!is.na(pathVar_assoc$ExAC_NFEassoc_P) & pathVar_assoc$ExAC_NFEassoc_P < 0.05,]
pathVar_assoc_sig = pathVar_assoc_sig[order(pathVar_assoc_sig$ExAC_NFEassoc_P),c(1,3,6:15,117,187:197)]

tn = "out/pathVarSingleAssoc.txt"
write.table(pathVar_assoc_sig[pathVar_assoc_sig$Overall_Classification != "Prioritized VUS",], quote=F, sep="\t", file = tn, row.names = F)

tn = "out/prioritizedVUSSingleAssoc.txt"
write.table(pathVar_assoc_sig[pathVar_assoc_sig$Overall_Classification == "Prioritized VUS",], quote=F, sep="\t", file = tn, row.names = F)

pathVarP_assoc = pathVar_assoc[pathVar_assoc$Overall_Classification %in% c("Likely Pathogenic","Pathogenic"),]
gene = pathVarP_assoc$HUGO_Symbol[!is.na(pathVarP_assoc$ExAC_NFEassoc_P) & pathVarP_assoc$ExAC_NFEassoc_P < 0.05]
pathVar_g_p_g = pathVarP_assoc[pathVarP_assoc$HUGO_Symbol %in% gene,]
#pathVar_g_p$minusLogP[!is.na(pathVar_g_p$minusLogP) & pathVar_g_p$minusLogP >6 ] = 6

featGenesAssoc = unique(c(pathVarP_assoc$HUGO_Symbol[!is.na(pathVarP_assoc$ExAC_NFEassoc_P) & pathVarP_assoc$ExAC_NFEassoc_P < 0.01],"BUB1B","PALB2","VHL","MET"))
featGenesAssoc = featGenesAssoc[featGenesAssoc!="BRCA2"] #somehow that sneaked in

p = ggplot(pathVar_g_p_g[pathVar_g_p_g$HUGO_Symbol %in% featGenesAssoc,],aes(x = HUGO_Symbol, y = -log10(ExAC_NFEassoc_P), color=Overall_Classification))
#p = p + facet_grid(Gene_Classification~., scale="free", space="free")
p = p + geom_jitter(alpha=0.5,size=1.5, height = 0,width = 0.25,shape=16, stroke=0)
p = p + geom_text_repel(aes(label=ifelse( (-log10(ExAC_NFEassoc_P) > 1.30103 & HUGO_Symbol %in% c("BUB1B","PALB2","VHL","MET")) | (-log10(ExAC_NFEassoc_P) > -log10(0.01)), 
                                    gsub("p.","",as.character(HGVSp_short)),NA)),size=3,angle=0)
p = p + theme_bw()
p = p + geom_hline(yintercept = -log10(0.05),alpha=0.3)
#p = p + xlim(0,6)
p = p + getVarColorScale()
p = p + labs( x="Gene", y = "-log10(P)") #+ guide_legend(title ="dosage")
p = p + theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5))
p = p + coord_equal()
p
fn = "out/germline_NFE_assoc_result_path.pdf"
ggsave(fn,height = 5,width=5, useDingbat=F)

pathVar_g_p_g$assocFDR = p.adjust(pathVar_g_p_g$ExAC_NFEassoc_P, method = "BH")
sum(pathVar_g_p_g$assocFDR < 0.05, na.rm=T)

# tn = "out/pathVarSingleAssoc.txt"
# write.table(pathVar_g_p_g, quote=F, sep="\t", file = tn, row.names = F)

p = ggplot(pathVar_g_p_g[pathVar_g_p_g$HUGO_Symbol %in% featGenes,],aes(x = HUGO_Symbol, y = -log10(assocFDR), color=Overall_Classification))
#p = ggplot(pathVar_g_p_p,aes(y = HUGO_Symbol, x = minusLogP, color=Final_Classification))
#p = p + facet_grid(Gene_Classification~., scale="free", space="free")
p = p + geom_point(alpha=0.5,size=1.5)
p = p + geom_text_repel(aes(label=ifelse(-log10(assocFDR) > -log10(0.05), as.character(HGVSp_short),NA)),size=3)
p = p + theme_bw()
p = p + getVarColorScale()
#p = p + xlim(0,30)
p = p + labs( x="Gene", y= "-log10(FDR)") #+ guide_legend(title ="dosage")
p = p + theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5))
p = p + coord_equal()
p
fn = "out/germline_NFE_assoc_result_pathFDR.pdf"
ggsave(fn,height = 5,width=10, useDingbat=F)

p = ggplot(pathVar_assoc[pathVar_assoc$HUGO_Symbol %in% featGenes,],aes(x = HUGO_Symbol, y = -log10(ExAC_NFEassoc_P), color=Overall_Classification))
#p = p + facet_grid(.~Gene_Classification, scale="free", space="free")
p = p + geom_jitter(alpha=0.2,size=0.5,height=0)
p = p + geom_text(aes(label=ifelse( -log10(ExAC_NFEassoc_P) > -log10(0.05) & Overall_Classification=="Prioritized VUS", 
                                    as.character(HGVSp_short),NA)),size=2,angle=0,vjust=-0.3)
p = p + theme_bw()
p = p + geom_hline(yintercept = -log10(0.05),alpha=0.3)
p = p + ylim(0,6)
p = p + getVarColorScale2()
p = p + labs( x="Gene", y = "-log(P)") #+ guide_legend(title ="dosage")
p = p + theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5))
p
fn = "out/germline_NFE_assoc_result_vus_included.pdf"
ggsave(fn,height = 5,width=6, useDingbat=F)


# pathVar_g_p_sig = pathVar_g_p[!is.na(pathVar_g_p$Pvalue) & pathVar_g_p$Pvalue < 0.05,]
# tn = paste(td, "sig_germline_in_select_genes.txt",sep="_")
# write.table(pathVar_g_p_sig, quote=F, sep="\t", file = tn, row.names = F)
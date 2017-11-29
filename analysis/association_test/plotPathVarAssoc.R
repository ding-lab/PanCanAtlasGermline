##### plotPathVarAssoc.R #####
# Kuan-lin Huang @ WashU 201711
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/association_test"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")
library(ggrepel)

# take a look
pathVarP$minusLogP = -log10(pathVarP$ExAC_assoc_P)
pathVar_g = pathVarP[pathVarP$HUGO_Symbol %in% featGenes,]
pathVar_g_p = pathVar_g[!duplicated(pathVar_g$HGVSg),]
pathVar_g_p$minusLogP[!is.na(pathVar_g_p$minusLogP) & pathVar_g_p$minusLogP >6 ] = 6
# p = ggplot(pathVar_g,aes(x = HUGO_Symbol, y = minusLogP, color=Overall_Classification))
# p = p + facet_grid(.~Gene_Classification, scale="free", space="free")
# p = p + geom_point(alpha=0.2,size=0.5)
# p = p + geom_text(aes(label=ifelse(minusLogP > 2.30103, as.character(HGVSp_brief),NA)),size=2,angle=0,vjust=-0.8)
# #p = p + geom_text(aes(label=ifelse(minusLogP > 2.30103, as.character(amino_acid_change),NA)),size=2)
# p = p + theme_bw()
# #p = p + xlim(0,30)
# p = p + geom_hline(yintercept = -log10(0.05),alpha=0.3)
# p = p + labs( x="Gene", y = "-log(P)") #+ guide_legend(title ="dosage")
# p = p + theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5))
# p
# fn = paste(pd,"germline_assoc_result.pdf",sep="_")
# ggsave(fn,height = 5,width=7, useDingbat=F)

#pathVar_g_p$HUGO_Symbol = factor(pathVar_g_p$HUGO_Symbol, levels = rev(gene_order))

pathVarP_assoc = pathVarP[!is.na(pathVarP$ExAC_assoc_P) & pathVarP$ExAC_assoc_P < 0.05,]
gene = unique(pathVarP_assoc$HUGO_Symbol)
pathVar_g_p_g =  pathVar_g_p[pathVar_g_p$HUGO_Symbol %in% gene,]

p = ggplot(pathVar_g_p_g,aes(x = HUGO_Symbol, y = minusLogP, color=Overall_Classification))
#p = p + facet_grid(Gene_Classification~., scale="free", space="free")
p = p + geom_jitter(alpha=0.5,size=1.5, height = 0,width = 0.25,shape=16, stroke=0)
p = p + geom_text_repel(aes(label=ifelse( (minusLogP > 1.30103 & HUGO_Symbol %in% c("BUB1B","PALB2","VHL","MET")) | (minusLogP > -log10(0.01)), 
                                    as.character(HGVSp_short),NA)),size=3,angle=0,vjust=1.5)
p = p + theme_bw()
p = p + geom_hline(yintercept = -log10(0.05),alpha=0.3)
#p = p + xlim(0,6)
p = p + getVarColorScale()
p = p + labs( x="Gene", y = "-log10(P)") #+ guide_legend(title ="dosage")
p = p + theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5))
p = p + coord_equal()
p
fn = "out/germline_assoc_result_path.pdf"
ggsave(fn,height = 5,width=5, useDingbat=F)

pathVar_g_p_g$assocFDR = p.adjust(pathVar_g_p_g$ExAC_assoc_P, method = "BH")
sum(pathVar_g_p_g$assocFDR < 0.05, na.rm=T)

tn = "out/pathVarSingleAssoc.txt"
write.table(pathVar_g_p_g, quote=F, sep="\t", file = tn, row.names = F)

p = ggplot(pathVar_g_p_g,aes(x = HUGO_Symbol, y = -log10(assocFDR), color=Overall_Classification))
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
fn = "out/germline_assoc_result_pathFDR.pdf"
ggsave(fn,height = 5,width=5, useDingbat=F)

p = ggplot(pathVar_g_p,aes(x = HUGO_Symbol, y = minusLogP, color=Overall_Classification))
#p = p + facet_grid(.~Gene_Classification, scale="free", space="free")
p = p + geom_point(alpha=0.2,size=0.5)
p = p + geom_text(aes(label=ifelse( (minusLogP > 1.30103 & Gene_Classification=="Oncogene") | (minusLogP > -log10(0.01) & Gene_Classification=="TSG"), 
                                    as.character(HGVSp_short),NA)),size=2,angle=0,vjust=-0.3)
p = p + theme_bw()
p = p + geom_hline(yintercept = -log10(0.05),alpha=0.3)
p = p + ylim(0,6)
p = p + getVarColorScale2()
p = p + labs( x="Gene", y = "-log(P)") #+ guide_legend(title ="dosage")
p = p + theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5))
p
fn = "out/germline_assoc_result_vus_included.pdf"
ggsave(fn,height = 5,width=6, useDingbat=F)


# pathVar_g_p_sig = pathVar_g_p[!is.na(pathVar_g_p$Pvalue) & pathVar_g_p$Pvalue < 0.05,]
# tn = paste(td, "sig_germline_in_select_genes.txt",sep="_")
# write.table(pathVar_g_p_sig, quote=F, sep="\t", file = tn, row.names = F)
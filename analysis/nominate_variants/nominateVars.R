##### nominateVars.R #####
# Kuan-lin Huang @ WashU 2018
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/nominate_variants"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

ptm = read.table("../PTMs/danbo/variant_ptm_overlap.final.substitutionsOnly.txt", quote="",fill=T,sep="\t", header=T)

### investigation into VUS ###
pathVar_vus = pathVar[pathVar$Overall_Classification=="Prioritized VUS",]
# pathVar_vus = pathVar[(((pathVar$colocalized_somatic_mutation_count >=3) | #pathVar$clustered | 
#                                     pathVar$LOH_Sig =="Uncertain Significance" | 
#                                     (!is.na(pathVar$ExAC_assoc_P) & (pathVar$ExAC_assoc_P < 0.05) ) |
#                                     (!is.na(pathVar$expressionQuantile) & !is.na(pathVar$Gene_Classification) & (pathVar$Gene_Classification=="Oncogene") & (pathVar$expressionQuantile >= 0.75)) | 
#                                     (!is.na(pathVar$expressionQuantile) & !is.na(pathVar$Gene_Classification) & (pathVar$Gene_Classification=="TSG") & (pathVar$expressionQuantile <= 0.25)))) &
#                                   (as.character(pathVar$Overall_Classification) == "Uncertain Significance"),]

pathVar_vus$colocBool=0
pathVar_vus$LOHBool=0
pathVar_vus$LowExpBool=0
pathVar_vus$assocBool = 0
pathVar_vus$HighExpBool = 0 
pathVar_vus$ptmBool = 0 
pathVar_vus$colocBool[(pathVar_vus$colocalized_somatic_mutation_count >=3) ] = 1
pathVar_vus$LOHBool[pathVar_vus$LOH_classification %in% c("Significant")] = 1
pathVar_vus$LowExpBool[pathVar_vus$Gene_Classification=="TSG" & pathVar_vus$expressionQuantile <= 0.25] = 1
pathVar_vus$HighExpBool[pathVar_vus$Gene_Classification=="Oncogene" & (pathVar_vus$expressionQuantile >= 0.75)] = 1
pathVar_vus$ptmBool[pathVar_vus$HGVSp.1 %in% ptm$HGVSp] =1
pathVar_vus$assocBool[pathVar_vus$ExAC_assoc_P < 0.05] = 1

pathVar_vus$countBool = rowSums(pathVar_vus[,grep("Bool",colnames(pathVar_vus))])
# pathVar_vus$class = "None"
# pathVar_vus$class[pathVar_vus$Gene_Classification=="TSG" & pathVar_vus$expressionQuantile <= 0.25] = "Low_Expression"
# pathVar_vus$class[pathVar_vus$Gene_Classification=="Oncogene" & (pathVar_vus$expressionQuantile >= 0.75)] = "High_Expression"
# #pathVar_vus[(pathVar_vus$clustered & pathVar_vus$colocalized_somatic_mutation_count < 3),]$class = "Co-cluster"
# pathVar_vus[(pathVar_vus$colocalized_somatic_mutation_count >=3),]$class = "Co-localize"
# pathVar_vus$class[!is.na(pathVar_vus$LOH_Sig) & pathVar_vus$LOH_Sig=="Significant"] = "LOH"
# pathVar_vus$class[pathVar_vus$countBool > 2] = "Multiple"

fn = "out/PCA_germline_vus_upset.pdf"
pdf(fn,useDingbats=FALSE)
upset(pathVar_vus, sets = c("colocBool", "LOHBool", "LowExpBool","assocBool","HighExpBool","ptmBool"),
      mb.ratio = c(0.8, 0.2))#, order.by = "freq")
dev.off()

# fn = "out/PCA_germline_vus_upset_oncogene.pdf"
# pdf(fn,useDingbats=FALSE)
# upset(pathVar_vus[pathVar_vus$Gene_Classification=="Oncogene",], sets = c("clusterBool", "colocBool", "LOHBool","HighExpBool","assocBool"),
#       mb.ratio = c(0.8, 0.2), order.by = "freq")
# dev.off()
# 
# fn = "out/PCA_germline_vus_upset_tsg.pdf"
# pdf(fn,useDingbats=FALSE)
# upset(pathVar_vus[pathVar_vus$Gene_Classification=="TSG",], sets = c("clusterBool", "colocBool", "LOHBool", "LowExpBool","assocBool"),
#       mb.ratio = c(0.8, 0.2), order.by = "freq")
# dev.off()

pathVar_vus_int = pathVar_vus[(pathVar_vus$countBool > 1),]

tn = "out/PCA_2_evidence_nominated_prioritizedVUS.txt"
write.table(pathVar_vus_int, quote=F, sep="\t", file = tn, row.names = F)
# 
# pathVar_vus_int_tier1 = pathVar_vus_int[pathVar_vus_int$ExAC_assoc_P < 0.01,]
# tn = "PCA_2_evidence_nominated_VUS_all_tier1.txt"
# write.table(pathVar_vus_int_tier1, quote=F, sep="\t", file = tn, row.names = F)
# 
# pathVar_vus_int_onco = pathVar_vus_int[!is.na(pathVar_vus_int$Gene_Classification) & pathVar_vus_int$Gene_Classification=="Oncogene",]
# tn = "PCA_2_evidence_nominated_VUS_oncogene.txt"
# write.table(pathVar_vus_int_onco, quote=F, sep="\t", file = tn, row.names = F)

# for plotting
pathVar_vus_int$LOH_FDR[is.na(pathVar_vus_int$LOH_FDR)] = 1
pathVar_vus_int$ExAC_assoc_P[is.na(pathVar_vus_int$ExAC_assoc_P)] = 1

# #p = ggplot(pathVar_vus_int,aes(x = tumorByNormalVAF, y = expressionQuantile, color=class))
# #p = ggplot(pathVar_vus_int,aes(x = tumorByNormalVAF, y = -log10(FDR), color=class))
# p = ggplot(pathVar_vus_int,aes(x = -log10(FDR), y = HUGO_Symbol, color=countBool))
# p = p + facet_grid(.~Overall_Classification,drop=T)
# p = p + geom_point()
# p = p + theme_bw() + labs(x="LOH analysis: -log10(FDR)",y="gene")
# p = p + geom_text(aes(label=ifelse(countBool > 2,paste(HUGO_Symbol, HGVSp_short),NA)),vjust=0.3,size=3)
# p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(colour="black", size=10), axis.text.x = element_text(colour="black", size=10, angle=90, vjust=0.5))
# #p = p + coord_trans(x="log2")
# p = p + geom_vline(xintercept=-log10(0.05), alpha=0.5)
# p = p + xlim(0,14)
# p
# fn = "out/PCA_germline_vus_candidates_LOH.pdf"
# ggsave(file=fn, width = 8, height=7,useDingbats=FALSE)

#top_genes = names(table(pathVar_vus_int$HUGO_Symbol)[order(table(pathVar_vus_int$HUGO_Symbol),decreasing = T)][1:23])
#pathVar_vus_int_g = pathVar_vus_int[pathVar_vus_int$HUGO_Symbol %in% top_genes,]
pathVar_vus_int_OT = pathVar_vus_int[!is.na(pathVar_vus_int$Gene_Classification),]
gene_ordered = table(pathVar_vus_int_OT$HUGO_Symbol)[order(table(pathVar_vus_int_OT$HUGO_Symbol),decreasing = T)]
top_genes = names(gene_ordered)
pathVar_vus_int_OT_g = pathVar_vus_int_OT[pathVar_vus_int_OT$HUGO_Symbol %in% top_genes,]
pathVar_vus_int_OT_g$HUGO_Symbol = factor(pathVar_vus_int_OT_g$HUGO_Symbol,levels=top_genes)

p = ggplot(pathVar_vus_int_OT_g,aes(x = HUGO_Symbol))
p = p + facet_grid(.~Gene_Classification,scale="free",space="free")
#p = p + facet_grid(Gene_Classification~.,scale="free",space="free")
p = p + geom_bar() #+ ylim(0,16)
p = p + theme_bw() #+ scale_y_continuous(trans='log2')
# p = p + geom_text_repel(aes(label=ifelse((HighExpBool==1 & spatialBool & Gene_Classification=="Oncogene") | (countBool >3 & Gene_Classification=="TSG"),paste(HUGO_Symbol, HGVSp_short),NA),size=2))
# p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(colour="black", size=10), axis.text.x = element_text(colour="black", size=10, angle=90, vjust=0.5))
# p = p + labs(x="Expression Quantile", y = "Tumor VAF / Normal VAF")
# p = p + ylim(-0.8,2.8)
# p = p + geom_hline(yintercept=1, alpha=0.5)
p = p + theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5))
p = p  + labs(x="Gene") #+ coord_flip()
p
fn = "out/nominated_vus_count_oncogene.TSG.pdf"
ggsave(file=fn, width = 7, height=5,useDingbats=FALSE)

p = ggplot(pathVar_vus_int_OT,aes(x = TumorByNormalVAF, y = expressionQuantile))
#p = ggplot(pathVar_vus_int_g,aes(x = -log10(FDR), y = HUGO_Symbol, color=class))
#p = ggplot(pathVar_vus_int_g,aes(x = -log10(FDR), y = HUGO_Symbol, color=class))
#p = p + facet_grid(Gene_Classification~Overall_Classification,drop=T, scale="free", space="free")
p = p + facet_grid(.~Gene_Classification,drop=T, scale="free", space="free")
p = p + geom_point(alpha=0.3,size=0.5)
p = p + theme_bw() #+ labs(x="LOH analysis: -log10(FDR)",y="Association analysis: -log10(P)")
#p = p + geom_text(aes(label=ifelse(((countBool) >2),paste(HUGO_Symbol, HGVSp_short),NA)),vjust=0.3,size=2,angle=-45)
#p = p + geom_text_repel(aes(label=ifelse((((countBool) >2 & Gene_Classification=="Oncogene") | ((countBool) >3 & Gene_Classification=="TSG")),paste(HUGO_Symbol, HGVSp_short),NA)),size=2)
p = p + geom_text_repel(aes(label=ifelse((Gene_Classification=="Oncogene") | 
                                           (countBool >2 & Gene_Classification=="TSG"),paste(HUGO_Symbol, HGVSp_short),NA)),size=3)#,vjust=-0.6)
p = p + theme(axis.title = element_text(size=12), axis.text.y = element_text(colour="black", size=10), axis.text.x = element_text(colour="black", size=10, angle=90, vjust=0.5))
p = p + labs(y="Expression Quantile", x = "Tumor VAF / Normal VAF")
p = p + xlim(-0.3,2.3)
p = p + geom_vline(xintercept=1, alpha=0.5)
p
fn = "out/PCA_germline_vus_candidates_oncogene.TSG.pdf"
ggsave(file=fn, width = 7, height=5,useDingbats=FALSE)

cat("All VUS studied: ", sum(pathVar$Overall_Classification == "Uncertain Significance"),"\n")
cat("VUS showing potential enrichment in cancer:","\n")
sum(pathVar_vus$assocBool==1,na.rm=T)
cat("Nominated alleles based on an additional evidence:","\n")
dim(pathVar_vus_int)[1]
cat("In ", length(unique(pathVar_vus_int$HUGO_Symbol)), " genes\n")
cat("By gene counts:\n")
table(pathVar_vus_int$HUGO_Symbol)[order(table(pathVar_vus_int$HUGO_Symbol),decreasing = T)][1:30]
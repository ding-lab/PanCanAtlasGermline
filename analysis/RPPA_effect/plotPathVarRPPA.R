##### plotPathVarExpression.R #####
# Kuan-lin Huang @ WashU 201711
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/RPPA_effect"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

RPPA = read.table("out/pancan_RPPA_quantile_all.tsv",header=T, stringsAsFactors = F, quote = "", sep = "\t")

RPPA_g = RPPA[RPPA$genes %in% pathVarP$HUGO_Symbol,]
RPPA_g_m = RPPA_g[,c("marker","genes","sample_l","bcr_patient_barcode","expression","quantile")]
RPPA_g_m$marker = gsub(".*\\|","",RPPA_g_m$marker)
colnames(RPPA_g_m)[2] = "HUGO_Symbol" 
pathVarP_RPPA = merge(pathVarP,RPPA_g_m,by=c("HUGO_Symbol","bcr_patient_barcode"))
pathVarP_RPPA_fg = pathVarP_RPPA[pathVarP_RPPA$HUGO_Symbol %in% featGenes,]
pathVarP_RPPA_fg = pathVarP_RPPA_fg[!is.na(pathVarP_RPPA_fg$binary_type),]

p = ggplot(pathVarP_RPPA_fg,aes(x=marker,y=quantile, fill=binary_type))
p = p + facet_grid(.~Gene_Classification, scale = "free", space = "free", drop=T)
p = p + geom_dotplot(dotsize=1.2,binwidth=.01, binaxis= "y",colour=NA,stackdir ="centerwhole")
p = p + geom_text(aes(label=ifelse(Gene_Classification=="Oncogene" & quantile>0.75, gsub("p.","",HGVSp_short),NA)),size=2.5)
p = p + theme_bw() 
p = p + ylab("RPPA Expression Quantile") + xlab("Protein with Germline Variant")
p = p + scale_y_continuous(breaks = seq(0,1, by= 0.25))
#p = p + getVarColorScale()
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))
p
fn = "out/pathVarRPPAExpression_byGene.pdf"
ggsave(file=fn, w=6, h=6,useDingbats=FALSE)

# ### somatic information ###
# somatic_f = "/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/somatic_germline/pancan.merged.v0.2.4.filtered.maf.gene_vclass_HGVSp_sample.txt"
# somatic = read.table(header=T, quote = "", sep="\t", file = somatic_f, stringsAsFactors=FALSE)
# colnames(somatic) = c("HUGO_Symbol","Somatic_Variant_Classification","sample","Somatic_HGVSp")
# somatic$sample = gsub("(^TCGA-[A-Z0-9][A-Z0-9]-[A-Z0-9][A-Z0-9][A-Z0-9][A-Z0-9])-.*","\\1",somatic$sample)
# # one entry per sample
# somatic_class_agg = aggregate(somatic[c('Somatic_Variant_Classification','Somatic_HGVSp')], by=somatic[c('sample',"HUGO_Symbol")], paste, collapse = ",")
# germ_so = merge(germ_clin_abb,somatic_class_agg, by =c("sample","HUGO_Symbol"), all=T)
# 
# ### clustering information###
# cluster_var_f = "/Users/khuang/Box Sync/PhD/germline/pan8000_germline_clinical/germline_hotspot/20161010_germline_ARD_ASD_run/pan8000_somatic_germline_combined.maf.3D_Proximity.pairwise.singleprotein.collapsed.l0.p0.05.r10.clusters_wcounts.tsv"
# cluster_var = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = cluster_var_f)
# 
# cluster_germsoma_f = "/Users/khuang/Box Sync/PhD/germline/pan8000_germline_clinical/germline_hotspot/20161010_germline_ARD_ASD_run/pan8000_somatic_germline_combined.maf.3D_Proximity.pairwise.singleprotein.collapsed.l0.p0.05.r10.clusters.summary_wcounts_cc10.3_wgermsoma.tsv"
# cluster_germsoma = read.table(header=T, quote = "", sep="\t", fill =T, file = cluster_germsoma_f, stringsAsFactors=FALSE)
# 
# cluster_var_in_hybrid = cluster_var[cluster_var$Cluster %in% cluster_germsoma$Cluster_ID,]
# germ_clin_var = germ_clin[paste(germ_clin$chromosome_name, germ_clin$start, germ_clin$stop) %in% 
#                             paste(cluster_var_in_hybrid$Chromosome, cluster_var_in_hybrid$Start, cluster_var_in_hybrid$Stop),]



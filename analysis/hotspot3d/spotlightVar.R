##### spotlightVar.R #####
# Kuan-lin Huang @ WashU 201711
# plot special variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/hotspot3d"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

# check if there is statistical enrichment of overlaps
# $ awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' all_CDS_and_ncRNA_24Chroms_Contigs_1BasedStart_2bpFlanks_ForMusic_merged
# 49586385
exonSize = 49586385
nPath = nrow(pathVarP)

### somatic mutation
# $ gzcat /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/somatic/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz | cut -f1,4 | sort | uniq -c | awk '$1 >2 && $3 != "."' | wc -l
# 68537
numSomaticOverlap = nrow(pathVarP[pathVarP$colocalized_somatic_mutation_count > 2,])
somaticMutRate = 68537/exonSize  
poisson.test(numSomaticOverlap, T = nPath, r = somaticMutRate, conf.level = 0.95, alternative = "greater")

### PCGP germ var
# /charged.2015_stJude_germline_nejm_S4_AD_varOnly.vep.tsv 
# 551 /Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/pathogenic_variants/PCGP/charged.2015_stJude_germline_nejm_S4_AD_varOnly.vep.tsv
# vpn-10-1-24-5:analysis khuang$ wc -l /Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/pathogenic_variants/PCGP/charged.2015_stJude_germline_nejm_S4_AR_varOnly.vep.tsv 
# 239 /Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/pathogenic_variants/PCGP/charged.2015_stJude_germline_nejm_S4_AR_varOnly.vep.tsv

numPCGPOverlap = nrow(pathVarP[pathVarP$PCGP,])
pcgpMutRate = (551 + 239)/exonSize
poisson.test(numPCGPOverlap, T = nPath, r = pcgpMutRate, conf.level = 0.95, alternative = "greater")

# write file
pathVarP_hot = pathVarP[pathVarP$colocalized_somatic_mutation_count > 2 | pathVarP$PCGP,]
write.table(file = "out/colocalize_var.tsv", pathVarP_hot,quote=F, sep = '\t',row.names = F)

# plot
pathVarPOT_hot = pathVarPOT[pathVarPOT$colocalized_somatic_mutation_count > 2 | pathVarPOT$PCGP,]
table(pathVarPOT_hot$HUGO_Symbol)
pathVarPOT_hot$somatic_count_plot = pathVarPOT_hot$colocalized_somatic_mutation_count
pathVarPOT_hot$HGVSp_short_plot = gsub("p.","",pathVarPOT_hot$HGVSp_short)
#pathVarPOT_hot$somatic_count_plot[pathVarPOT_hot$somatic_count_plot> 100 ]  = 100

p = ggplot(pathVarPOT_hot,aes(y=HUGO_Symbol, x =somatic_count_plot, color = PCGP))
#p = p + facet_grid(PCGP~Gene_Classification,drop=T,scale="free",space="free")
p = p + facet_grid(Gene_Classification~ .,drop=T,scale="free_y",space="free_y")
p = p + geom_point(stroke=0) + theme_bw()  #+ guides(color=FALSE)
#p = p + geom_abline(intercept = 0, slope=1, alpha=0.2) #+ geom_density2d(alpha=0.5)
p = p + geom_text_repel(aes(label=ifelse(duplicated(HGVSp_short),NA,HGVSp_short_plot)))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + scale_x_log10()
p = p + expand_limits(x = 0)
#p = p + coord_equal() + getLOHColorScale()
p = p + labs(x = "Co-localizing somatic mutation count", y = "Gene")
p
fn = "out/pathVarP_spotlight.pdf"
ggsave(file=fn, width=7, h =5, useDingbats=FALSE)

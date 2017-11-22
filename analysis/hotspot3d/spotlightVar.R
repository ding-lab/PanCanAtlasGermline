##### spotlightVar.R #####
# Kuan-lin Huang @ WashU 201711
# plot special variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/hotspot3d"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

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

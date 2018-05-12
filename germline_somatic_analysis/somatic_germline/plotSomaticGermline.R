##### plotSomaticGermline.R #####
# Kuan-lin Huang 2018

#setwd("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/somatic_germline/")

source("../global_aes_out.R")
source("../dependency_files.R")

tn = "out/germline_somatic_driver_fisher.tsv"
tt = read.table(sep="\t",header=T,file=tn, stringsAsFactors=FALSE)

germlineG = unique(tt$GermlineGene[tt$P<0.05])
somaticG = unique(tt$SomaticGene[tt$P<0.01])
ttG = tt[tt$GermlineGene %in% germlineG & tt$SomaticGene %in% somaticG,]

# pre-plotting
ttG$minusLogP = -log10(ttG$P) 
ttG$plotP = ttG$minusLogP
ttG$plotP[ttG$OR < 1] = -ttG$plotP[ttG$OR < 1] # opposite effect size (mutual exclusivity)
ttG$plotP[ttG$plotP > 5] = 5
ttG$plotP[ttG$plotP < -5] = -5

### plotting ###

#getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
p = ggplot(data=ttG)
p = p + geom_tile(data=ttG,aes(y=SomaticGene, x=GermlineGene, fill= plotP), linetype="blank") + scale_fill_gradientn(name= "-log10(P)", colours=getPalette(100), na.value=NA, limit=c(-5,5))
#p = p + geom_text(data=ttG,aes(x=cancer, y=signature, label = coefficient), color="black", size=3)
#p = p + geom_tile(data=ttG[ttG$FDR < 0.15,],aes(y=cancer, x=signature), color="grey",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
p = p + geom_tile(data=ttG[ttG$FDR < 0.05,],aes(y=SomaticGene, x=GermlineGene, fill= plotP), color="black",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p = p + labs(title="Germline-somatic Interaction: Pan-cancer",x="Germline variant carrier",y = "Somatic mutation carrier")
p
fn = 'out/pan10389_germlineAssocWithSomaticHeatmap.pdf'
ggsave(fn,useDingbat=F)
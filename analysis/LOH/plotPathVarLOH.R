##### plotPathVarLOH.R #####
# Kuan-lin Huang @ WashU 201711
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/LOH"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

getLOHColorScale = function() {
  colors = c("#e31a1c", "#b2df8a", "#1f78b4") #positive is dark grey       
  color.names = c("Significant","Suggestive","None")
  names(colors) = color.names
  color.scale = scale_color_manual(name="Loss of Heterozygosity", values=colors)
  return(color.scale)
}

getLOHFillScale = function() {
  colors = c("#e31a1c", "#b2df8a", "#1f78b4") #positive is dark grey       
  color.names = c("Significant","Suggestive","None")
  names(colors) = color.names
  color.scale = scale_fill_manual(name="Loss of Heterozygosity", values=colors)
  return(color.scale)
}

#pathVarOT$LOH_Sig = factor(pathVarOT$LOH_Sig, levels=c("None","Suggestive","Significant"))

p = ggplot(pathVarOT,aes(x=normalVAF, y =tumorVAF, color = LOH_Sig))
p = p + facet_grid(.~Gene_Classification,drop=T)
p = p + geom_point(alpha=0.3, stroke=0) + theme_bw() + theme_nogrid() #+ guides(color=FALSE)
p = p + geom_abline(intercept = 0, slope=1, alpha=0.2) #+ geom_density2d(alpha=0.5)
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + expand_limits(x = 0, y = 0)
p = p + coord_equal() + getLOHColorScale()
p = p + labs(x = "Normal VAF", y = "Tumor VAF")
p
fn = "out/pathVarLOH.pdf"
ggsave(file=fn, width=6, useDingbats=FALSE)

# find LOH gene percentage
geneLOH = data.frame(table(pathVarOT$HUGO_Symbol,pathVarOT$LOH_Sig))
colnames(geneLOH) = c("Gene","LOH_category","Count")

geneLOHsig = geneLOH[geneLOH$LOH_category=="Significant",]
sig_gene = geneLOHsig[geneLOHsig$Count>2,]$Gene
sig_gene_order = geneLOHsig$Gene[order(geneLOHsig$Count,decreasing = T)]
geneLOH_g = geneLOH[geneLOH$Gene %in% sig_gene,]

geneLOH_g$Gene = factor(geneLOH_g$Gene, levels=sig_gene_order)
geneLOH_g$LOH_category = factor(geneLOH_g$LOH_category, levels=c("None","Suggestive","Significant"))

p = ggplot(geneLOH_g,aes(x = Gene, y = Count, fill = LOH_category))
p = p + geom_bar(stat = "identity") + theme_bw() + theme_nogrid() 
p = p + labs(x = "Gene", y="Count of variants")
p = p + getLOHFillScale()
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "top")
p
fn = 'out/LOH_var_count_by_gene.pdf'
ggsave(file=fn, height = 6, width = 6, useDingbats=FALSE)

pathVarOT$TumorByNormalVAFPlot = pathVarOT$TumorByNormalVAF
pathVarOT$TumorByNormalVAFPlot[pathVarOT$TumorByNormalVAFPlot> 2 ] = 2
pathVarOT$TumorByNormalVAFPlot[pathVarOT$TumorByNormalVAFPlot< 0.5 ] = 0.5
p = ggplot(pathVarOT[pathVarOT$HUGO_Symbol %in% sig_gene,],aes(x = HUGO_Symbol, y = TumorByNormalVAFPlot, color = LOH_Sig, fill = HUGO_Symbol))
p = p + geom_point(position = position_jitter(w = 0.2, h = 0), alpha = 0.3)
#p = p + geom_violin(alpha = 0.3) 
p = p + labs(x = "Gene", y="Tumor VAF / Normal VAF") + theme_bw()
#p = p + getVarFillScale()
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "none")
p
fn = 'out/LOH_var_VAFratio_by_gene.pdf'
ggsave(file=fn, height = 6, width = 6, useDingbats=FALSE)

p = ggplot(pathVarOT[pathVarOT$HUGO_Symbol %in% sig_gene,],aes(x=HUGO_Symbol,y = TumorByNormalVAFPlot, color = LOH_Sig, fill = LOH_Sig))
#p = p + facet_grid(.~Gene_Classification, scale = "free", space = "free", drop=T)
p = p + geom_dotplot(dotsize=0.6,binwidth=.015, binaxis= "y",stackdir ="centerwhole",alpha=0.5)
p = p + theme_bw() 
p = p + labs(x = "Gene", y="Tumor VAF / Normal VAF")
p = p + scale_y_continuous(breaks = seq(0.5,2.0, by= 0.5))
#p = p + getVarColorScale()
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "none")
p
fn = 'out/LOH_var_VAFratio_by_gene_dotplot.pdf'
ggsave(file=fn, height = 6, width = 6, useDingbats=FALSE)
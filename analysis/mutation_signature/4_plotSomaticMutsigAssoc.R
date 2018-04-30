##### plotPathVarMutsigAssoc.R #####
# Kuan-lin Huang 2018

source("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/global_aes_out.R")
source("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/dependency_files.R")

g_tn = "out/pathVarMutsigAssoc.txt"
g_tt = read.table(sep="\t",header=T,file=g_tn, stringsAsFactors=FALSE)

tn = "out/somaticMutMutsigAssoc.txt"
tt = read.table(sep="\t",header=T,file=tn, stringsAsFactors=FALSE)

### plotting ###
tt$signature = factor(tt$signature)
tt$signature = factor(tt$signature,levels = c("Signature-1","Signature-2","Signature-3","Signature-4","Signature-5","Signature-6"
                                              ,"Signature-7","Signature-8","Signature-9","Signature-10","Signature-11","Signature-12"
                                              ,"Signature-13","Signature-14","Signature-15","Signature-16","Signature-17","Signature-18"
                                              ,"Signature-19","Signature-20","Signature-21","Signature-22","Signature-23","Signature-24"
                                              ,"Signature-25","Signature-26","Signature-27","Signature-28","Signature-29","Signature-30"))
tt$association = "None"
tt$association[tt$FDR<0.15] = "Suggestive"
tt$association[tt$FDR<0.05] = "Significant"
tt$gene = as.character(tt$gene)
tt$FDR_plot = -log(tt$FDR)
tt$FDR_plot[tt$FDR_plot > 5 ] = 5
#uniqG = unique(tt$gene[tt$FDR<0.05])
uniqG = unique(g_tt$gene[g_tt$FDR<0.05]) #  plot just the germline genes for now
ttG = tt[tt$gene %in% uniqG,]

getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
p = ggplot(data=ttG)
p = p + facet_grid(gene~.,space="free",scale="free")
p = p + geom_tile(data=ttG,aes(y=cancer, x=signature, fill= coefficient), linetype="blank") + scale_fill_gradientn(name= "Coefficient", colours=getPalette(100), na.value=NA, limit=c(0,NA))
#p = p + geom_text(data=ttG,aes(x=cancer, y=signature, label = coefficient), color="black", size=3)
#p = p + geom_tile(data=ttG[ttG$FDR < 0.15,],aes(y=cancer, x=signature), color="grey",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
p = p + geom_tile(data=ttG[ttG$FDR < 0.05,],aes(y=cancer, x=signature), color="black",fill=NA, size=1.5) #+ scale_color_gradientn(name= "Sig", colours=sigColors)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p + labs(x="Signature",y = "Cancer")

fn = 'out/SomaticWithmutSignatureHeatmap.pdf'
ggsave(fn,h=25,useDingbat=F)

# plot by gene
p = ggplot(data=ttG,aes(x=coefficient,y=cancer,color = cancer))
p = p + facet_grid(gene~.,space="free",scale="free")
# p = p + geom_point(aes(y=-log10(wilcoxFDR),x= coefficient,color = cancer),alpha=0.5)
p = p + geom_point(aes(shape=association),alpha=0.3,size=2)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.0005,signature,NA)))
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "-log10(FDR)")
p = p + geom_vline(xintercept = 0, alpha=0.5) #+ xlim(-3.1,3.1)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p + labs(x = "coefficient",y="cancer")
fn = 'out/SomaticWithmutSignatureByGene.pdf'
ggsave(fn,h=28,useDingbat=F)

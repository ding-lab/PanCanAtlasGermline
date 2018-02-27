##### plotPathVarExpressAssoc.R #####
# Kuan-lin Huang @ WashU 2016 August updated 2017
# conduct association of pathVarPline variants with AAO

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/expression_effect"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

tn = "out/pathVarExpressAssoc.txt"
tt = read.table(sep="\t",header=T,file=tn, stringsAsFactors=FALSE)

### plotting ###
tt$gene = as.character(tt$gene)
tt$FDR_plot = tt$FDR

# # using GLM test result
# tt$FDR_plot[tt$FDR_plot<10^(-6)]= 0.95*10^(-6)
# p = ggplot(data=tt)
# p = p + geom_point(aes(y=-log10(FDR_plot),x= coefficient,color = cancer),alpha=0.5)
# #p = p + geom_text_repel(aes(y=-log10(FDR),x= cohort_AF,label=ifelse(FDR<0.05, Gene,NA),color = Cancer))#,alpha=1.3)
# #p = p + geom_point(aes(y=cohort_AF,x=Cancer,size=-log10(FDR),color = Cancer))
# p = p + geom_text_repel(aes(y=-log10(FDR_plot),x= coefficient,color = cancer,label=ifelse(FDR<0.05, gene,NA)))
# p = p + getPCACancerColor()
# p = p + labs(x="Coefficient",y= "-log10(FDR)")
# p = p + geom_vline(xintercept = 0, alpha=0.5)
# p = p  + theme_bw() +
#   theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# fn = 'out/geneExpressAssocVolcanoGLM.pdf'
# ggsave(fn,w = 5, h = 5, useDingbat=F)

# # using the Wilcox test result
# p = ggplot(data=tt)
# p = p + geom_point(aes(y=-log10(wilcoxFDR),x= coefficient,color = cancer),alpha=0.5)
# #p = p + geom_text_repel(aes(y=-log10(FDR),x= cohort_AF,label=ifelse(FDR<0.05, Gene,NA),color = Cancer))#,alpha=1.3)
# #p = p + geom_point(aes(y=cohort_AF,x=Cancer,size=-log10(FDR),color = Cancer))
# p = p + geom_text_repel(aes(y=-log10(wilcoxFDR),x= coefficient,color = cancer,label=ifelse(FDR<0.05, gene,NA)))
# p = p + getPCACancerColor()
# p = p + labs(x="Coefficient",y= "-log10(FDR)")
# p = p + geom_vline(xintercept = 0, alpha=0.5)
# p = p  + theme_bw() +
#   theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# fn = 'out/geneExpressAssocVolcanoWCOX.pdf'
# ggsave(fn,w = 5, h = 5, useDingbat=F)

tt$association = "None"
tt$association[tt$FDR<0.15] = "Suggestive"
tt$association[tt$FDR<0.05] = "Significant"
#tt$association = factor(tt$association,level=c("None","Suggestive","Significant"))

# using the Wilcox test result: plot by gene
p = ggplot(data=tt,aes(x=coefficient,y=cancer,color = cancer))
# p = p + geom_point(aes(y=-log10(wilcoxFDR),x= coefficient,color = cancer),alpha=0.5)
p = p + geom_point(aes(shape=association),alpha=0.3,size=2)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.05,gene,NA)))
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "-log10(FDR)")
p = p + geom_vline(xintercept = 0, alpha=0.5) + xlim(-3.1,3.1)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/geneExpressAssocByGene.pdf'
ggsave(fn,w = 5, h = 5, useDingbat=F)

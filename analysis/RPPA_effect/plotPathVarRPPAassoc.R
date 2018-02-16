##### plotPathVarExpressAssoc.R #####
# Kuan-lin Huang @ WashU 2016 August updated 2017
# conduct association of pathVarPline variants with AAO

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/RPPA_effect"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

getPCACancerColor = function() {
  # according to google spreadsheet: https://docs.google.com/spreadsheets/d/1Nb9mMkonAhZR1_2OI9nv4ylCei0LZjAf-2vYTRQcXKw/edit#gid=1704872109
  colors = c(
    "#C1A72F",
    "#FAD2D9",
    "#ED2891",
    "#F6B667",
    "#104A7F",
    "#9EDDF9",
    "#3953A4",
    "#007EB5",
    "#B2509E",
    "#97D1A9",
    "#ED1C24",
    "#F8AFB3",
    "#EA7075",
    "#754C29",
    "#D49DC7",
    "#CACCDB",
    "#D3C3E0",
    "#A084BD",
    "#542C88",
    "#D97D25",
    "#6E7BA2",
    "#E8C51D",
    "#7E1918",
    "#DAF1FC",
    "#00A99D",
    "#BBD642",
    "#00AEEF",
    "#BE1E2D",
    "#F9ED32",
    "#CEAC8F",
    "#FBE3C7",
    "#F89420",
    "#009444")    
  color.names = c("ACC",
                  "BLCA",
                  "BRCA",
                  "CESC",
                  "CHOL",
                  "COAD",
                  "DLBC",
                  "ESCA",
                  "GBM",
                  "HNSC",
                  "KICH",
                  "KIRC",
                  "KIRP",
                  "LAML",
                  "LGG",
                  "LIHC",
                  "LUAD",
                  "LUSC",
                  "MESO",
                  "OV",
                  "PAAD",
                  "PCPG",
                  "PRAD",
                  "READ",
                  "SARC",
                  "SKCM",
                  "STAD",
                  "TGCT",
                  "THCA",
                  "THYM",
                  "UCEC",
                  "UCS",
                  "UVM")
  names(colors) = color.names
  color.scale = scale_color_manual(name="Cancer", values=colors)
  return(color.scale)
}

tn = "out/pathVarRPPAAssoc.txt"
tt = read.table(sep="\t",header=T,file=tn, stringsAsFactors=FALSE)

### plotting ###
tt$gene = as.character(tt$gene)
tt$marker = as.character(tt$marker)
tt$FDR_plot = tt$FDR

# # using GLM test result
# tt$FDR_plot[tt$FDR_plot<10^(-6)]= 0.95*10^(-6)
# p = ggplot(data=tt)
# p = p + geom_point(aes(y=-log10(FDR_plot),x= coefficient,color = cancer),alpha=0.5)
# #p = p + geom_text_repel(aes(y=-log10(FDR),x= cohort_AF,label=ifelse(FDR<0.05, Gene,NA),color = Cancer))#,alpha=1.3)
# #p = p + geom_point(aes(y=cohort_AF,x=Cancer,size=-log10(FDR),color = Cancer))
# p = p + geom_text_repel(aes(y=-log10(FDR_plot),x= coefficient,color = cancer,label=ifelse(FDR<0.05, marker,NA)))
# p = p + getPCACancerColor()
# p = p + labs(x="Coefficient",y= "-log10(FDR)")
# p = p + geom_vline(xintercept = 0, alpha=0.5)
# p = p  + theme_bw() +
#   theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# fn = 'out/RPPAAssocVolcanoGLM.pdf'
# ggsave(fn,w = 5, h = 5, useDingbat=F)

# using the Wilcox test result
p = ggplot(data=tt)
p = p + geom_point(aes(y=-log10(wilcoxFDR),x= coefficient,color = cancer),alpha=0.5)
#p = p + geom_text_repel(aes(y=-log10(FDR),x= cohort_AF,label=ifelse(FDR<0.05, Gene,NA),color = Cancer))#,alpha=1.3)
#p = p + geom_point(aes(y=cohort_AF,x=Cancer,size=-log10(FDR),color = Cancer))
p = p + geom_text_repel(aes(y=-log10(wilcoxFDR),x= coefficient,color = cancer,label=ifelse(FDR<0.05, marker,NA)))
p = p + getPCACancerColor()
p = p + labs(x="Coefficient",y= "-log10(FDR)")
p = p + geom_vline(xintercept = 0, alpha=0.5) + xlim(-1.6,1.6)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/geneExpressAssocVolcanoWCOX.pdf'
ggsave(fn,w = 5, h = 5, useDingbat=F)

# using the Wilcox test result: plot by gene
p = ggplot(data=tt,aes(x=coefficient,y=cancer,color = cancer))
# p = p + geom_point(aes(y=-log10(wilcoxFDR),x= coefficient,color = cancer),alpha=0.5)
p = p + geom_point(aes(size=-log10(FDR)))
p = p + geom_text_repel(aes(label=ifelse(FDR<0.05,marker,NA)))
p = p + getPCACancerColor()
p = p + labs(x="Coefficient",y= "-log10(FDR)")
p = p + geom_vline(xintercept = 0, alpha=0.5) + xlim(-1.6,1.6)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/geneExpressAssocVolcanoWCOX_bygene.pdf'
ggsave(fn,w = 5, h = 5, useDingbat=F)
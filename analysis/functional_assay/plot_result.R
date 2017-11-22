##### plot_result.R #####
# Kuan-lin Huang @ WashU 2017 July
# plot RET functional experiment results
library(ggplot2)
setwd("/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/functional_assay/")

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  fn = args[1]
}
#fn = "20170717_CO/All_Results_RET_Gel1_75_4_071717.txt"
outFile = paste("output/",gsub(".*/","",fn),".pdf",sep="")
outFile2 = paste("output/",gsub(".*/","",fn),"2.pdf",sep="")
  
results = read.table(fn, sep="\t", header=T)

colnames(results) = gsub("\\.","_",colnames(results))
#colnames(results)
results$Ligand = gsub(".*_","",results$Sample)
results$Mut = gsub("_.*","",results$Sample)

results$Sample = factor(make.names(results$Sample,unique = T), levels=make.names(results$Sample,unique = T))

p = ggplot(results,aes(x=Sample, y=MAPK_RET_GAPDH_WT, fill=Ligand))
p = p + facet_grid(.~Ligand, drop=T, space="free",scale="free")
p = p + geom_bar(stat="identity") + theme_bw() #+ theme_nogrid()
p = p + labs(x = "Sample", y="MAPK/RET/GAPDH (normalized to WT)")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position="bottom")
p = p + geom_hline(yintercept = 1, alpha=0.3)
#p = p + geom_text(aes(label=Number_of_Phosphosites), vjust=-0.25)
p
#fn = "output/All_Results_RET_Gel1_75_4_071717.pdf"
ggsave(file=outFile, useDingbats=FALSE)

p = ggplot(results,aes(x=Sample, y=MAPK_GAPDH_WT, fill=Ligand))
p = p + facet_grid(.~Ligand, drop=T, space="free",scale="free")
p = p + geom_bar(stat="identity") + theme_bw() #+ theme_nogrid()
p = p + labs(x = "Sample", y="MAPK/RET/GAPDH (normalized to WT)")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position="bottom")
p = p + geom_hline(yintercept = 1, alpha=0.3)
#p = p + geom_text(aes(label=Number_of_Phosphosites), vjust=-0.25)
p

ggsave(file=outFile2, useDingbats=FALSE)
##### plot_pseq_stats.R #####
# Kuan-lin Huang @ WashU 2017 Aug.

setwd("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/variant_QC/")
source("../global_aes_out.R")
system("mkdir out")

##### individual level stats #####
fileNames = Sys.glob("filtered_VCF_pseq_stats/*istats*9401*")

istat_list = vector("list")
i=1
for (fileName in fileNames) {
  fileName1 = strsplit(fileName, split="/")[[1]][2]
  fileName2 = gsub(".anno.whitelist.vcf.gz.pseq.istats.tsv.filtered.tsv.9401.tsv","",fileName1)
  chrName = gsub("PCA.merge.tcgaBarcode.","",fileName2)
  #cat(paste(cancer,"\n"))
  #exp_table = read.table(row.names=1,header=TRUE, sep="\t", file=fileName)
  istat = read.table(header=TRUE, sep="\t", file=fileName, fill=T)
  istat_filt = istat[1:9401,]
  istat_filt$chr = chrName
  istat_list[[i]] = istat_filt
  i = i + 1
}

istat_all = do.call(rbind,istat_list)
istat_all$chr[istat_all$chr=="PCA.merge.ExAConly.tcgaBarcode"] = "ExAC"
istat_all$chr = gsub("chr","",istat_all$chr)
#istat_all_m = melt(istat_all, id.vars = c("ID","chr"))
istat_all$chr_plot = factor(istat_all$chr, levels = c(sort(unique(as.numeric(istat_all$chr))),"X","Y","ExAC"))

# for (col in colnames(istat_all)[2:12]){
#   #p = ggplot(data=istat_all,aes_string(x_string = "chr", y_string=col), fill= "chr") # can't get aes_string to work
#   p = ggplot(data=istat_all,aes(x = chr, y=NALT), fill= chr)
#   #p = p + facet_grid(variable~., drop=T, space="free", scale="free")
#   p = p + geom_boxplot() 
#   p = p  + theme_bw() + theme(legend.position="none")
#   #theme(axis.title = element_text(size=18), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5), axis.text.y = element_text(colour="black", size=14),axis.ticks = element_blank())#element_text(colour="black", size=14))
#   p
#   fn = paste('out/pseq.istat.9401',col,'boxplot.pdf',sep=".")
#   ggsave(file=fn, useDingbats=FALSE,limitsize=FALSE)
# }

p = ggplot(data=istat_all,aes(x = chr_plot, y=NVAR), fill= chr_plot)
p = p + geom_boxplot(alpha=0.1) 
p = p  + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Chromosome")
p
fn = paste('out/pseq.istat.9401.NVAR.boxplot.pdf',sep=".")
ggsave(file=fn, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=istat_all,aes(x = chr_plot, y=TITV), fill= chr_plot)
p = p + geom_boxplot(alpha=0.1) 
p = p  + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Chromosome")
p
fn = paste('out/pseq.istat.9401.TITV.boxplot.pdf',sep=".")
ggsave(file=fn, useDingbats=FALSE,limitsize=FALSE)

##### plot NVAR by ethnicity and cancer type #####
istat_all_NVAR_ID = aggregate(NVAR ~ ID, data=istat_all, sum)
# use previous clinical table first
clin_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/all.clin.merged.picked.txt"
clin = read.table(header=T, quote = "", sep="\t", row.names =NULL, file = clin_f, stringsAsFactors=FALSE)
clin = clin[,c("sample","cancer","race")]
colnames(clin) = c("sample","cancer","ethnicity")
clin = clin[!(clin$cancer %in% c("GBMLGG","COADREAD","KIPAN")),]
istat_all_NVAR_ID$sample = substring(istat_all_NVAR_ID$ID,1,12)
istat_all_NVAR_ID_clin = merge(istat_all_NVAR_ID,clin,by="sample")

p = ggplot(data=istat_all_NVAR_ID_clin,aes(x = cancer, y=NVAR))
#p = p + geom_violin(alpha=0.1, fill=NA) 
# p = p + geom_dotplot(binaxis="y", stackdir = "center",binwidth=100) 
p = p + geom_jitter(aes(color= ethnicity),height=0, alpha = 0.1, size=1) 
p = p  + theme_bw() #+ theme(legend.position="none")
p = p + labs(x = "Cancer", y = "Count of rare variants in ROI")
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p
fn = paste('out/pseq.istat.9401.NVAR.dotplot.pdf',sep=".")
ggsave(file=fn, w = 10, useDingbats=FALSE,limitsize=FALSE)

##### chromosome level stats #####

v_fileNames = Sys.glob("filtered_VCF_pseq_stats/*vstats*")

vstat_list = vector("list")
i=1
for (v_fileName in v_fileNames) {
  v_fileName1 = strsplit(v_fileName, split="/")[[1]][2]
  v_fileName2 = gsub(".anno.whitelist.vcf.gz.pseq.vstats.tsv.filtered.tsv","",v_fileName1)
  chrName = gsub("PCA.merge.tcgaBarcode.","",v_fileName2)
  #cat(paste(cancer,"\n"))
  #exp_table = read.table(row.names=1,header=TRUE, sep="\t", v_file=v_fileName)
  vstat = read.table(header=F, sep="\t", file=v_fileName, fill=T)
  vstat$V3 = chrName
  
  vstat_list[[i]] = vstat
  i = i + 1
}

vstat_all = do.call(rbind,vstat_list)
colnames(vstat_all) = c("attribute","value","chr")
vstat_all$chr[vstat_all$chr=="PCA.merge.ExAConly.tcgaBarcode"] = "ExAC"
vstat_all$chr = gsub("chr","",vstat_all$chr)
#vstat_all_m = melt(vstat_all, id.vars = c("ID","chr"))
vstat_all_c = dcast(vstat_all, chr~attribute)
vstat_all_c$chr_plot = factor(vstat_all_c$chr, levels = c(sort(unique(as.numeric(vstat_all_c$chr))),"X","Y","ExAC"))

p = ggplot(data=vstat_all_c,aes(x = chr_plot, y=NVAR), fill= chr_plot)
p = p + geom_bar(stat="identity") 
p = p  + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Chromosome")
p
fn = paste('out/pseq.vstat.NVAR.barplot.pdf',sep=".")
ggsave(file=fn, useDingbats=FALSE,limitsize=FALSE)

p = ggplot(data=vstat_all_c,aes(x = chr_plot, y=TITV), fill= chr_plot)
p = p + geom_bar(stat="identity") 
p = p  + theme_bw() + theme(legend.position="none")
p = p + labs(x = "Chromosome")
p
fn = paste('out/pseq.vstat.TITV.barplot.pdf',sep=".")
ggsave(file=fn, useDingbats=FALSE,limitsize=FALSE)
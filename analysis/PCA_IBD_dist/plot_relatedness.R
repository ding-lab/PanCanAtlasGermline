##### plot_PCA.R #####
# Kuan-lin Huang @ WashU 2017 March
# plot PCA output from plink
# for TCGA samples, label based on reported ethnicity

setwd("/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/PCA_IBD_dist")
source("../global_aes_out.R")

# plink relation file
# pca_f = "plink_out/all.normal.merge.vcf.pca.eigenvec"
# relationship based on Yang J, Lee SH, Goddard ME, Visscher PM (2011) GCTA: A Tool for Genome-wide Complex Trait Analysis 
rel_f = "plink_out/all.normal.merge.indep_50_5_2.vcf.rel.rel"
rel = read.table(header=F, quote = "", sep="\t", file = rel_f)

# merge with sample
sample_f = "plink_out/all.normal.merge.indep_50_5_2.vcf.rel.rel.id"
sample = read.table(header=F, quote = "", sep="\t", row.names =NULL, file = sample_f, stringsAsFactors=FALSE)
samples = sample[,1]
row.names(rel) = samples
colnames(rel) = samples
rel$sample = samples

# make a small version of the table and take a quick look
rel_m = melt(rel[1:10,],id.var="sample")
rel_m[rel_m$sample ==rel_m$variable,]

# plotting
p = ggplot(data=rel_m)
p = p + geom_bar(aes(x=value),stat="bin",bins=100)  
#p = p + scale_colour_gradientn(na.value="grey", colours=getPalette(100))#, limits=c(-4.2,4.2))
p = p + theme_bw() + scale_y_log10() #+ expand_limits(y=1)#+ guides(fill=FALSE)
p = p + geom_vline(xintercept = 0,alpha=.7) + geom_vline(xintercept = 1,alpha=.7)
p
fn = paste(pd, "PanCanAtlas_rel_10samples_hist.pdf", sep="_")
ggsave(file=fn, useDingbats=FALSE)


### IBD
dist_f = "plink_out/all.normal.merge.indep_50_5_2.vcf.dist.dist"
dist = read.table(header=F, quote = "", sep="\t", file = dist_f)

# merge with sample
sample_f = "plink_out/all.normal.merge.indep_50_5_2.vcf.dist.dist.id"
sample = read.table(header=F, quote = "", sep="\t", row.names =NULL, file = sample_f, stringsAsFactors=FALSE)
samples = sample[,1]
row.names(dist) = samples
colnames(dist) = samples
dist$sample = samples

# make a small version of the table and take a quick look
dist_m = melt(dist[1:10,],id.var="sample")
dist_m[dist_m$sample ==dist_m$variable,]

# plotting
p = ggplot(data=dist_m)
p = p + geom_bar(aes(x=value),stat="bin",bins=100)  
#p = p + scale_colour_gradientn(na.value="grey", colours=getPalette(100))#, limits=c(-4.2,4.2))
p = p + theme_bw() + scale_y_log10() #+ expand_limits(y=1)#+ guides(fill=FALSE)
p = p + geom_vline(xintercept = 0,alpha=.7) + geom_vline(xintercept = 1,alpha=.7)
p
fn = paste(pd, "PanCanAtlas_dist_10samples_hist.pdf", sep="_")
ggsave(file=fn, useDingbats=FALSE)
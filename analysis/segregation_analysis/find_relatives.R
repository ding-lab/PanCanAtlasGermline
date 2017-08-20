##### find_relatives.R #####
# Kuan-lin Huang @ WashU 2017 July

setwd("/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/segregation_analysis")
source("../global_aes_out.R")
system("mkdir out")
# ethnicity assignment
ethni_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/PCA_IBD_dist/out/2017-04-24/2017-04-24_GDAN_AIM_PCA_ethnicity_assigned_WashU.tsv"
ethni = read.table(header=T, sep = '\t',file=ethni_fn)
ethni_short = ethni[,c("Case","cancer","washu_assigned_ethnicity")]


#genome file for all the pi HAT
rel_f = "plink_out/all.normal.merge.indep_50_5_2.vcf.ibd.genome.PI_HAT0.05.tsv"
rel = read.table(header=T, quote = "", sep="\t", row.names =NULL, file = rel_f, stringsAsFactors=FALSE)
rel_short = rel[rel$Z1 + rel$Z2 > 0.2,]
rel_short$Sample1 = gsub("(.{12}).*","\\1",as.character(rel_short$FID1))
rel_short$Sample2 = gsub("(.{12}).*","\\1",as.character(rel_short$FID2))

colnames(ethni_short) = c("Sample1","cancer1","assigned_ethnicity1")
rel_short_m1 = merge(rel_short,ethni_short, by="Sample1",all.x=T)
colnames(ethni_short) = c("Sample2","cancer2","assigned_ethnicity2")
rel_short_m2 = merge(rel_short_m1,ethni_short, by="Sample2",all.x=T)
rel_short_m2$same_sample = (rel_short_m2$Sample1 == rel_short_m2$Sample2)

rel_short_m2_same = rel_short_m2[rel_short_m2$assigned_ethnicity2==rel_short_m2$assigned_ethnicity1,]
rel_short_m2_same_withethni = rel_short_m2_same[rel_short_m2_same$assigned_ethnicity1!="unknown",]
rel_short_m2_same_withethni$assigned_ethnicity1 = as.character(rel_short_m2_same_withethni$assigned_ethnicity1)

p = ggplot(data=rel_short_m2_same_withethni,aes(x = Z1, y=Z2, color=same_sample))
p = p + facet_grid(.~assigned_ethnicity1, drop=T, scales="free",space="free")
p = p + geom_point(alpha=0.2)  
p = p + theme_bw()  #+ expand_limits(y=1)#+ guides(fill=FALSE)
p = p + xlim(0,1) + ylim(0,1)
#p = p + geom_vline(xintercept = 0,alpha=.7) + geom_vline(xintercept = 1,alpha=.7)
p
fn = paste(pd, "PanCanAtlas_rel_z1.z2_withinEthni.pdf", sep="_")
ggsave(file=fn, height=5, width=12, useDingbats=FALSE)

rel_short_m2_same_withethni_ofinterest = rel_short_m2_same_withethni[!rel_short_m2_same_withethni$same_sample,]
rel_short_m2_same_withethni_ofinterest = rel_short_m2_same_withethni_ofinterest[(rel_short_m2_same_withethni_ofinterest$Z2 > 0.125 & rel_short_m2_same_withethni_ofinterest$Z1 > 0.25) |
                                                                                  (rel_short_m2_same_withethni_ofinterest$Z1 > 0.20) ,]

tn = paste("out/TCGA_z1_z2_relatives.tsv", sep="_")
write.table(rel_short_m2_same_withethni_ofinterest,quote=F, sep = '\t', row.names = FALSE,file=tn)

rel_short_m2_same_withethni_ofinterest2 = rel_short_m2_same_withethni_ofinterest[(rel_short_m2_same_withethni_ofinterest$Z2 > 0.125 & rel_short_m2_same_withethni_ofinterest$Z1 > 0.25) |
                                                                                   (rel_short_m2_same_withethni_ofinterest$Z1 > 0.625) ,]
#rel_short_m2_same_withethni_ofinterest2[,c("Sample1","Sample2","cancer1","cancer2")]
tn = paste("out/TCGA_z1_z2_relatives_strict.tsv", sep="_")
write.table(rel_short_m2_same_withethni_ofinterest2,quote=F, sep = '\t', row.names = FALSE,file=tn)

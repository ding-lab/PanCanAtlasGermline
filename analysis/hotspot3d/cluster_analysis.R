##### cluster_analysis.R #####
# Kuan-lin Huang @ WashU 201711
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/hotspot3d"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

s_fn = "mc3.v0.2.8.PUBLIC.code.filtered.PCAgenes.hotspot.maf"
somatic = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=s_fn)

fn = "PCA_somatic_germline_combined_missense.maf.3D_Proximity.pairwise.recurrence.l0.r10.clusters"
cluster = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=fn)

colnames(cluster)[2:3] = c("HUGO_Symbol","HGVSp_short")
cluster$type = NA
cluster$type[paste(cluster$HUGO_Symbol,cluster$HGVSp_short) %in% paste(somatic$Hugo_Symbol,somatic$HGVSp_Short)] = "Somatic"
cluster$type[is.na(cluster$type) & (paste(cluster$HUGO_Symbol,cluster$HGVSp_short) %in% paste(pathVarP$HUGO_Symbol,pathVarP$HGVSp_short))] = "Germline"
cluster$type[cluster$type== "Somatic" & (paste(cluster$HUGO_Symbol,cluster$HGVSp_short) %in% paste(pathVar$HUGO_Symbol,pathVar$HGVSp_short))] = "Colocalized"

cluster_w_germ = cluster$Cluster[cluster$type %in% c("Germline","Colocalized")]

cluster_germ = cluster[cluster$Cluster %in% cluster_w_germ,]
cat("Number of sites co-clustered: ","\n")
table(cluster_germ$type)

cat("Number of unique clusters: ",length(unique(cluster_germ$Cluster)),"\n")

table(cluster_germ$HUGO_Symbol[!duplicated(cluster_germ$Cluster)])

tn = "PCA_somatic_germline_combined_missense.maf.3D_Proximity.pairwise.recurrence.l0.r10.clusters_annotated.tsv"
write.table(cluster_germ, quote=F, sep="\t", file = tn, row.names = F)

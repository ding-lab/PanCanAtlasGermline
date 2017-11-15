##### RPPA_effect.R #####
# Kuan-lin Huang @ WashU 2016 May , updated 2017 Nov.
# analyze cohort level RPPA data and convert to different matrices in a sample-gene format

bdir = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/RPPA_effect"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")
ecdf_fun = function(x,perc) ecdf(x)(perc)

### preprocess RPPA file ###
fileNames = Sys.glob("/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/oncogene_signaling/data/RPPA/*.rppa.txt")
all_tables = vector("list")

for (fileName in fileNames) {
  cancer2 = strsplit(fileName, split="/")[[1]][length(strsplit(fileName, split="/")[[1]])]
  cancer = gsub("\\..*","",cancer2)

  exp_table = read.table(header=TRUE, sep="\t", file=fileName)
  exp_table_q = exp_table
  for (i in 1:nrow(exp_table_q)){
    min_RPPA = min(exp_table_q[i,-1],na.rm=T)
    if (min_RPPA<0){ exp_table_q[i,-1] = exp_table_q[i,-1] - min_RPPA}
    exp_table_q[i,-1] = ecdf_fun(unlist(exp_table_q[i,-1]),unlist(exp_table_q[i,-1]))
  }
  exp_table.m = melt(exp_table, id.var="Composite.Element.REF")
  exp_table_q.m = melt(exp_table_q, id.var="Composite.Element.REF")
  colnames(exp_table.m) = c("marker","sample","expression")
  colnames(exp_table_q.m) = c("marker","sample","quantile")
  exp_table_c = merge(exp_table.m,exp_table_q.m,by=c("marker","sample"))
  exp_table_c$cancer = cancer
  all_tables[[cancer]] = exp_table_c
}
RPPA = do.call(rbind,all_tables)
RPPA$sample_l = substr(RPPA$sample, start=0, stop=16)
RPPA$sample_l = gsub("\\.","-",RPPA$sample)
RPPA$bcr_patient_barcode = substr(RPPA$sample_l, start=0, stop=12) 

RPPA$status = "tumor"
status_n = substr(RPPA$sample_l, start=14, stop=14)
RPPA$status[status_n==1] = "normal"
RPPA = RPPA[RPPA$status != "normal",]# exclude normal from BRCA for now

RPPA$genes = gsub("\\|.*","",RPPA$marker)
RPPA$genes[RPPA$genes=="MAPK1 MAPK3"] = "MAPK3"
RPPA$genes[RPPA$genes=="PIK3R1 PIK3R2"] = "PIK3R1"
RPPA$genes[RPPA$genes=="PIK3R1/2"] = "PIK3R1"
RPPA$genes[RPPA$genes=="PIK3CA "] = "PIK3CA"
RPPA$genes[RPPA$genes=="PDK1"] = "PDPK1"

fn = "out/pancan_RPPA_quantile_all.tsv"
write.table(RPPA, file=fn, quote=F, sep="\t", col.names=T, row.names=F)
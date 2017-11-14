##### expression_effect.R #####
# Kuan-lin Huang @ WashU 2016 May , updated 2017 Nov.
# analyze cohort level RNA-Seq data and convert to different matrices in a sample-gene format

bdir = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/expression_effect"
setwd(bdir)
source("../../global_aes_out.R")
source("../../dependency_files.R")


## function ##
unfactorize = function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.numeric(as.character(df[[i]]))
  return(df)
}

ecdf_fun = function(x,perc) ecdf(x)(perc)

expression_effect = function(m){ 
  cat("##### EXPRESSION ANALYSIS #####\n")
  minNum = 5
  m = as.matrix(m)
  num = nrow(m)
  m2 = as.matrix(m[rowSums(!is.na(m)) >= minNum, ])
  num_NA= nrow(m2)
  cat(paste("Original number of markers:", num, "; NA filtered:", num_NA, "\n", sep=" "))
  
  # initiate tables
  outlier = matrix(,nrow=dim(m2)[1],ncol=dim(m2)[2])
  row.names(outlier) = row.names(m2)
  colnames(outlier) = colnames(m2)
  exp_score = outlier
  exp_quantile = outlier
  
  # gene-wise expression score and quantile score
  for (i in 1:nrow(m2)){
    #IQR = quantile(m2[i,], probs=0.75, na.rm=T) - quantile(m2[i,], probs=0.25, na.rm=T) 
    exp_score[i,] = m2[i,]#(m2[i,] - quantile(m2[i,], probs=0.50, na.rm=T))/IQR
    exp_quantile[i,] = ecdf_fun(m2[i,],m2[i,])
  }
  
  return(list("exp_score"=exp_score, "exp_quantile"=exp_quantile))
}

glist_f = read.table(header=FALSE, stringsAsFactors = F, file = "/Users/khuang/Box Sync/PhD/germline/pan8000_germline_clinical/reference_files/CancerGeneListV5-2014-04-18.add-Rahman.add-Fanconi-Gene.txt")
glist = as.vector(t(glist_f))

fileNames = Sys.glob("/Users/khuang/Box\ Sync/PhD/collaborations/premed_2015/data/All_gene_RNASeq/raw_output/*RSEM_hugo.txt")
# # get rid of COADREAD
# CR = "/Users/khuang/Box Sync/PhD/collaborations/premed_2015/data/All_gene_RNASeq/raw_output/COADREAD_RSEM_hugo.txt"
# fileNames = fileNames[-which(fileNames == CR)]

exp_score_tables = vector("list")
exp_quantile_tables = vector("list")
exp_tables = vector("list")

for (fileName in fileNames) {
  cancer2 = strsplit(fileName, split="/")[[1]][11]
  cancer = gsub("_.*","",cancer2)
  #cat(paste(cancer,"\n"))
  #exp_table = read.table(row.names=1,header=TRUE, sep="\t", file=fileName)
  exp_table = read.table(header=TRUE, sep="\t", file=fileName)
  exp_table = exp_table[exp_table$Hybridization.REF %in% glist,]
  
  row.names(exp_table) = make.names(exp_table[,1],unique=T)
  exp_table = exp_table[,-c(1,2)]
  
  # get tumor only
  normal = substr(colnames(exp_table), 14, 14)
  exp_table = exp_table[,normal=="0"] 
  
  exp_table = unfactorize(exp_table)
  exp_table = log2(exp_table+1)
  new_colname = vector()
  for (name in colnames(exp_table)){
    splitted_sname = strsplit(name, split = "\\.")[[1]]
    new_name = paste(splitted_sname[1],splitted_sname[2],splitted_sname[3],sep=".")
    new_colname = c(new_colname,new_name)
  }
  colnames(exp_table) = new_colname

  if (dim(exp_table)[1] == 0 || dim(exp_table)[2] == 0){next;}
  
  exp_results = expression_effect(exp_table)
  
  exp_score_table_m = melt(exp_results$exp_score)
  exp_score_table_m$cancer = cancer
  exp_score_tables[[cancer]] = exp_score_table_m
  
  exp_quantile_table_m = melt(exp_results$exp_quantile)
  exp_quantile_table_m$cancer = cancer
  exp_quantile_tables[[cancer]] = exp_quantile_table_m
  
  exp_tables[[cancer]]
}

exp_score_tables_all = do.call(rbind,exp_score_tables)
colnames(exp_score_tables_all) = c("gene_name","sample","log2RSEM","cancer")
exp_quantile_tables_all = do.call(rbind,exp_quantile_tables)
colnames(exp_quantile_tables_all) = c("gene_name","sample","expression_quantile","cancer")
exp_score_tables_all$sample = gsub("\\.","-",exp_score_tables_all$sample)
exp_quantile_tables_all$sample = gsub("\\.","-",exp_quantile_tables_all$sample)

fn = "out/pancan_exp_log2RSEM_all.tsv"
write.table(exp_score_tables_all, file=fn, quote=F, sep="\t", col.names=T, row.names=F)
fn = "out/pancan_exp_quantile_all.tsv"
write.table(exp_quantile_tables_all, file=fn, quote=F, sep="\t", col.names=T, row.names=F)

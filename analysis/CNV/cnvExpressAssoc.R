##### cnvExpressAssoc.R #####
# Kuan-lin Huang @ WashU 2018 Feb
# conduct association of CNV with Expression

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/CNV"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

### data


#FUNCTION myglm#
myglm=function(z,trait,variant,covar=NA,ytype) {
  if (is.na(covar) | is.null(covar) | nchar(covar)==0  ) { 
    model=formula(paste(trait,"~",variant)) 
  } else {
    model=formula(paste(trait,"~",variant,"+",covar))
  }
  if (ytype=="B") fit=glm(formula=model,data=z,family=binomial(link = "logit"))
  if (ytype=="Q") fit=glm(formula=model,data=z,family=gaussian(link = "identity"))
  fit
}
#END myglm

#FUNCTION plot_violin#
### plot quantitative y in relationship to x
### faceted by covariate covi
plot_violin = function(y,yi,xi,covi){ 
  # change height based on the number of markers
  dat = y[,c(yi, xi, covi)]
  dat[,1] = as.numeric(dat[,1])
  #dat.m = melt(mut_exp, id.var="Status")
  # plot violin plots faceted by marker genes
  n_facet = length(unique(dat[,covi]))
  p = ggplot(data=dat)
  p = p + facet_grid(as.formula(paste(". ~", covi)))
  p = p + geom_violin(aes_string(x=xi, y=yi, fill=xi),alpha=0.5) + guides(fill=FALSE, color =FALSE) 
  p = p + geom_jitter(aes_string(x=xi, y=yi, color=xi), alpha = 0.2) #+ geom_point(aes(x=Status, y=value)) 
  p = p + guides(fill=FALSE, color =FALSE) 
  p = p + labs(x = paste(xi,"pathVarPline Variant"), y = yi) + theme_bw()
  p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8, angle = 90))
  p
  fn = paste(pd, yi, "by", xi, "cross", covi, "violin.pdf", sep="_")
  ggsave(file=fn, width = n_facet, useDingbats=FALSE)
}

run_glm = function(data=NULL, covi="") {
  
  ### data input #####
  data_g = data#[, c("age_at_initial_pathologic_diagnosis","type",gene)]
  
  ######### analysis ##########
  row_stat=NULL
  
  # analysis_type   clinical_data_trait_name    variant/gene_name   covariates  memo
  ytype="Q";yi="log2RSEM";xi="seg_mean";covi=covi;#covi=md[i,4];memo=md[i,5]
  # DEBUG
  #cat(paste("    Processing: yi =", yi, " xi =", xi, " covi =", covi, "\n") )
  
  if (ytype=="Q") {
    test = "F"
  } else if (ytype=="B") {
    test = "Chisq"
  } else {
    stop("Unknown model ytype ", ytype) 
  }
  
  glm = try(myglm(data_g,yi,xi,covi,ytype))  # MAW new
  if(class(glm)[1] == "try-error") {
    cat(paste("    Error caught, continuing.  yi =", yi, " xi =", xi, " covi =", covi, "\n") )
    next
  }
  try(anova(glm,test=test))->fit
  
#   if (xi %in% names(coefficients(glm))){
#     coeff = coefficients(glm)[[xi]]
  if (length(names(coefficients(glm)))>1){
    coeff = coefficients(glm)[[2]]
  } else {coeff = NA}
  
  if (class(fit)[1]!="try-error")
  {
    fit=as.matrix(fit)
    if (xi %in% rownames(fit)) (row_stat = cbind(yi,ytype,xi,as.data.frame(t(fit[xi,])),coeff,covi))
  }
  return(row_stat)
  
}

### MAIN ###

###  using array data ###
# read input files
array_cnv_f = "data/array.cnv.filtered.rare.seg.refseq.expq.annotated.rahman.txt"
array_cnv = read.table(header=T, quote = "", sep="\t", fill =T, file = array_cnv_f, stringsAsFactors=FALSE)
colnames(array_cnv)[2]= "bcr_patient_barcode"
array_cnv$seg_mean[array_cnv$seg_mean < 0] = -1
array_cnv$seg_mean[array_cnv$seg_mean > 0] = 1
array_cnv = array_cnv[!duplicated(paste(array_cnv$bcr_patient_barcode,array_cnv$gene)),]
# gene_sample = data.frame(table(pathVarP$HUGO_Symbol,pathVarP$bcr_patient_barcode))
# colnames(gene_sample) = c("Gene","bcr_patient_barcode","Freq")

exp_score_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/expression_effect/2016-06-21_KH_pancan_exp_log2RSEM_all_uniq.tsv.gz"
exp_score = read.table(sep="\t",header=F,file=gzfile(exp_score_f), stringsAsFactors=FALSE)
colnames(exp_score)=c("Gene","bcr_patient_barcode","log2RSEM","cancer")
exp_score_g = exp_score[exp_score$Gene %in% array_cnv$gene,]
exp_score_g$log2RSEM = as.numeric(exp_score_g$log2RSEM)
##### individual cancer type analysis #####
cancers = unique(array_cnv$cancer)
genes = unique(array_cnv$gene)
# limit runs to cancers with at least 5 likely patho/pathogenic variants 
tt = NULL
for (gene in genes){
  gene_sample_g = array_cnv[array_cnv$gene==gene,]
  exp_score_g_g = exp_score_g[exp_score_g$Gene==gene,]
  var_exp_g = merge(exp_score_g_g,gene_sample_g,by=c("bcr_patient_barcode","cancer"),all.x=T)
  var_exp_g$seg_mean[is.na(var_exp_g$seg_mean)] = 0

  for (cancer in cancers){
    var_exp_g_c = var_exp_g[var_exp_g$cancer %in% cancer,]
    gene_path_count = sum(var_exp_g_c$seg_mean!=0)
    if (gene_path_count > 2){
      # run GLM
      w = wilcox.test(var_exp_g_c$log2RSEM[var_exp_g_c$seg_mean==0],var_exp_g_c$log2RSEM[var_exp_g_c$seg_mean!=0])
      wP = w$p.value
      wWstat = w$statistic
      cancer_gene_stat = run_glm(var_exp_g_c)
      # compile results
      full_cancer_gene_stat = cbind(cancer,gene,gene_path_count,wP,wWstat,cancer_gene_stat)
      tt = rbind(tt, full_cancer_gene_stat)
    }
  }
}

#"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
colnames(tt) = c("cancer","gene","gene_CNV_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F_statistic","p-value","coefficient","covariants");
tt$FDR = p.adjust(tt[,"p-value"], method="fdr") # MAW new, calculates FDR based on the method from,
# Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289–300.
tt$wilcoxFDR = p.adjust(tt[,"wilcoxP"], method="fdr")

tt=tt[order(tt$FDR, decreasing=FALSE),]
tn = "out/arrayCNVExpressAssoc.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

### using XHMM data ###
xhmm_cnv_f = "data/xhmm.auto.rare.q60.rm3sd.expq.annotated.rahman.txt"
xhmm_cnv = read.table(header=T, quote = "", sep="\t", fill =T, file = xhmm_cnv_f, stringsAsFactors=FALSE)
colnames(xhmm_cnv)[1]= "bcr_patient_barcode"
# xhmm_cnv$seg_mean = xhmm_cnv$MEAN_RD
xhmm_cnv$seg_mean = 0
xhmm_cnv$seg_mean[xhmm_cnv$CNV == "DUP"] = 1
xhmm_cnv$seg_mean[xhmm_cnv$CNV == "DEL"] = -1
xhmm_cnv = xhmm_cnv[!duplicated(paste(xhmm_cnv$bcr_patient_barcode,xhmm_cnv$gene)),]

##### individual cancer type analysis #####
cancers = unique(xhmm_cnv$cancer)
genes = unique(xhmm_cnv$gene)
# limit runs to cancers with at least 5 likely patho/pathogenic variants 
tt = NULL
for (gene in genes){
  gene_sample_g = xhmm_cnv[xhmm_cnv$gene==gene,]
  exp_score_g_g = exp_score_g[exp_score_g$Gene==gene,]
  var_exp_g = merge(exp_score_g_g,gene_sample_g,by=c("bcr_patient_barcode","cancer"),all.x=T)
  var_exp_g$seg_mean[is.na(var_exp_g$seg_mean)] = 0
  
  for (cancer in cancers){
    var_exp_g_c = var_exp_g[var_exp_g$cancer %in% cancer,]
    gene_path_count = sum(var_exp_g_c$seg_mean!=0)
    if (gene_path_count > 2){
      # run GLM
      w = wilcox.test(var_exp_g_c$log2RSEM[var_exp_g_c$seg_mean==0],var_exp_g_c$log2RSEM[var_exp_g_c$seg_mean!=0])
      wP = w$p.value
      wWstat = w$statistic
      cancer_gene_stat = run_glm(var_exp_g_c)
      # compile results
      full_cancer_gene_stat = cbind(cancer,gene,gene_path_count,wP,wWstat,cancer_gene_stat)
      tt = rbind(tt, full_cancer_gene_stat)
    }
  }
}

#"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
colnames(tt) = c("cancer","gene","gene_CNV_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F_statistic","p-value","coefficient","covariants");
tt$FDR = p.adjust(tt[,"p-value"], method="fdr") # MAW new, calculates FDR based on the method from,
# Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289–300.
tt$wilcoxFDR = p.adjust(tt[,"wilcoxP"], method="fdr")

tt=tt[order(tt$FDR, decreasing=FALSE),]
tn = "out/xhmmCNVExpressAssoc.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)
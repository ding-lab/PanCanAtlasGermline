##### pathVarP_AAO_assoc.R #####
# Kuan-lin Huang @ WashU 2016 August updated 2017
# conduct association of pathVarPline variants with AAO

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/clinical_association"
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

run_glm = function(data=NULL,gene=NULL, covi="") {
  
  ### data input #####
  data_g = data[, c("age_at_initial_pathologic_diagnosis","type",gene)]
  
  ######### analysis ##########
  row_stat=NULL
  
  # analysis_type   clinical_data_trait_name    variant/gene_name   covariates  memo
  ytype="Q";yi="age_at_initial_pathologic_diagnosis";xi=gene;covi=covi;#covi=md[i,4];memo=md[i,5]
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

# read input files
gene_sample = data.frame(table(pathVarP$HUGO_Symbol,pathVarP$bcr_patient_barcode))
colnames(gene_sample) = c("Gene","bcr_patient_barcode","Freq")
gene_sample_d = dcast(gene_sample, bcr_patient_barcode ~ Gene)

pathVarP_clin = merge(clin,gene_sample_d,by="bcr_patient_barcode",all.x=T)
pathVarP_clin[is.na(pathVarP_clin)]=0
pathVarP_clin = pathVarP_clin[,-c(6:39)]
##### individual cancer type analysis #####
cancers = unique(pathVarP_clin$type)
cols_g = colnames(pathVarP_clin)
genes = cols_g[!(cols_g) %in% c("bcr_patient_barcode","type","age_at_initial_pathologic_diagnosis","gender","race","consensus_call")]

# limit runs to cancers with at least 5 likely patho/pathogenic variants 
tt = NULL
for (cancer in cancers){
  pathVarP_clin_matrix_c = pathVarP_clin[pathVarP_clin$type %in% cancer,]
  for (gene in genes){
    gene_path_count = sum(as.numeric(as.character(pathVarP_clin_matrix_c[,gene])), na.rm=T)
    cancer_count = dim(pathVarP_clin_matrix_c)[1]
    freq = gene_path_count/cancer_count
    if (gene_path_count > 3 & freq >= 0.01){
      # run GLM
      cancer_gene_stat = run_glm(pathVarP_clin_matrix_c, gene)
      # compile results
      full_cancer_gene_stat = cbind(cancer,gene_path_count,freq,cancer_gene_stat)
      tt = rbind(tt, full_cancer_gene_stat)
    }
  }
}

if (!is.null(tt)) {
  #"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
  colnames(tt) = c("cancer","gene_path_count","percentage_carrier","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                                   "F_statistic","p-value","coefficient","covariants");
  #if (ytype=="B") colnames(tt) = c("y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                                   #"p-value","coefficient","covariants","memo");
  tt$FDR = p.adjust(tt[,"p-value"], method="fdr") # MAW new, calculates FDR based on the method from,
  # Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289–300.
  tt=tt[order(tt$FDR, decreasing=FALSE),]
  tn = "out/pathVarAAOassoc.txt"
  write.table(tt, quote=F, sep="\t", file = tn, row.names = F)
}
### plotting ###
genes = unique(as.character(tt[tt$FDR<0.05,]$Gene))
tt_s = tt[tt$FDR<0.05,]
pathVarP_clin_matrix_g = pathVarP_clin[,c(genes,"age_at_initial_pathologic_diagnosis","type")]
for (gene in genes){
  pathVarP_clin_matrix_g[,gene] = as.numeric(as.character(pathVarP_clin_matrix_g[,gene]))*pathVarP_clin_matrix_g[,"age_at_initial_pathologic_diagnosis"]
  pathVarP_clin_matrix_g[,gene][pathVarP_clin_matrix_g[,gene]==0]=NA
}
pathVarP_clin_matrix_g_m = melt(pathVarP_clin_matrix_g, id.var="type")
colnames(pathVarP_clin_matrix_g_m) = c("cancer","carrier","age_at_initial_pathologic_diagnosis")
pathVarP_clin_matrix_g_m$carrier = as.character(pathVarP_clin_matrix_g_m$carrier)
pathVarP_clin_matrix_g_m$carrier[pathVarP_clin_matrix_g_m$carrier=="age_at_initial_pathologic_diagnosis"] = "all"
pathVarP_clin_matrix_g_m$age_at_initial_pathologic_diagnosis[pathVarP_clin_matrix_g_m$age_at_initial_pathologic_diagnosis==0] = NA

for (cancer in cancers){
  #cancer="BRCA"
  tt_sc_gene = as.character(tt_s[tt_s$cancer %in% cancer,]$Gene)
  if (length(tt_sc_gene)>0){
    pathVarP_clin_matrix_gc = pathVarP_clin_matrix_g_m[(pathVarP_clin_matrix_g_m$cancer %in% cancer) & (pathVarP_clin_matrix_g_m$carrier %in% c(tt_sc_gene,"all")), ]
    p = ggplot(data=pathVarP_clin_matrix_gc)
    p = p + facet_grid(. ~cancer)
    p = p + geom_violin(aes(x=carrier, y=age_at_initial_pathologic_diagnosis),alpha=0.5)
    p = p + geom_jitter(aes(x=carrier, y=age_at_initial_pathologic_diagnosis, color=ifelse(carrier!="all",cancer,NA)), alpha = 0.8) #+ geom_point(aes(x=Status, y=value)) 
    p = p + guides(fill=FALSE, color =FALSE) 
    p = p + scale_colour_manual(values = col_vector)
    p = p + labs(x = "", y = "Age at onset") + theme_bw()
    p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
                  axis.text.y = element_text(colour="black", size=14))
    p
    fn = paste("out/PCA",tt_sc_gene, cancer,"AAO_violin.pdf", sep="_")
    ggsave(file=fn, width = 3,height=5,useDingbats=FALSE)
  }
}

##### pan-cancer ##### 
cols_g = colnames(pathVarP_clin)
genes = cols_g[!(cols_g) %in% c("bcr_patient_barcode","type","age_at_initial_pathologic_diagnosis","gender","race","consensus_call")]
tt = NULL
for (gene in genes){
  gene_path_count = sum(as.numeric(as.character(pathVarP_clin[,gene])), na.rm=T)
  if (gene_path_count >= 13){
    cancer_gene_stat = run_glm(pathVarP_clin, gene, covi="type")
    full_cancer_gene_stat = cbind("Pan-cancer",gene_path_count,freq,cancer_gene_stat)
    tt = rbind(tt, full_cancer_gene_stat)
  }
}
if (!is.null(tt)) {
  #"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
  colnames(tt) = c("cancer","gene_path_count","percentage_carrier","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                                   "F_statistic","p-value","coefficient","covariants");
  #if (ytype=="B") colnames(tt) = c("y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
               #                    "p-value","coefficient","covariants","memo");
  tt$FDR = p.adjust(tt[,"p-value"], method="fdr") # MAW new, calculates FDR based on the method from,
  # Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289–300.
  tt=tt[order(tt$FDR, decreasing=FALSE),]
  tn = "out/PanCan_pathVarAAOassoc.txt"
  write.table(tt, quote=F, sep="\t", file = tn, row.names = F)
}

### plotting ###
genes = as.character(tt[tt$FDR<0.05,]$x)
pathVarP_clin_matrix_g = pathVarP_clin[,c(genes,"age_at_initial_pathologic_diagnosis","type")]
for (gene in genes){
  pathVarP_clin_matrix_g[,gene] = as.numeric(as.character(pathVarP_clin_matrix_g[,gene]))*pathVarP_clin_matrix_g[,"age_at_initial_pathologic_diagnosis"]
  pathVarP_clin_matrix_g[,gene][pathVarP_clin_matrix_g[,gene]==0]=NA
}
pathVarP_clin_matrix_g_m = melt(pathVarP_clin_matrix_g, id.var="type")
colnames(pathVarP_clin_matrix_g_m) = c("cancer","carrier","age_at_initial_pathologic_diagnosis")
pathVarP_clin_matrix_g_m$carrier = as.character(pathVarP_clin_matrix_g_m$carrier)
pathVarP_clin_matrix_g_m$carrier[pathVarP_clin_matrix_g_m$carrier=="age_at_initial_pathologic_diagnosis"] = "ALL"
pathVarP_clin_matrix_g_m$age_at_initial_pathologic_diagnosis[pathVarP_clin_matrix_g_m$age_at_initial_pathologic_diagnosis==0] = NA
pathVarP_clin_matrix_g_m$tag = "Pan-cancer"

p = ggplot(data=pathVarP_clin_matrix_g_m)
p = p + facet_grid(. ~ tag)
p = p + geom_violin(aes(x=carrier, y=age_at_initial_pathologic_diagnosis),alpha=0.5)
p = p + geom_jitter(aes(x=carrier, y=age_at_initial_pathologic_diagnosis, color=ifelse(carrier!="ALL",cancer,NA)), alpha = 0.8) #+ geom_point(aes(x=Status, y=value)) 
p = p + guides(fill=FALSE) 
p = p + scale_colour_manual(values = col_vector)
p = p + labs(x = "", y = "Age at onset") + theme_bw()
p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
              axis.text.y = element_text(colour="black", size=14))
#p = p + theme(legend_position="bottom")
p
fn = "out/PCA_pancan_AAO_violin.pdf"
ggsave(file=fn, height=5,useDingbats=FALSE)

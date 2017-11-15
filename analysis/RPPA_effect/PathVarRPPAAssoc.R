##### pathVarP_AAO_assoc.R #####
# Kuan-lin Huang @ WashU 2016 August updated 2017
# conduct association of pathVarPline variants with AAO

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/RPPA_effect"
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
  ytype="Q";yi="expression";xi="Freq";covi=covi;#covi=md[i,4];memo=md[i,5]
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

RPPA = read.table("out/pancan_RPPA_quantile_all.tsv",header=T, stringsAsFactors = F, quote = "", sep = "\t")

RPPA_g = RPPA[RPPA$genes %in% pathVarP$HUGO_Symbol,]
RPPA_g_m = RPPA_g[,c("marker","genes","bcr_patient_barcode","expression","cancer")]
RPPA_g_m$marker = gsub(".*\\|","",RPPA_g_m$marker)

colnames(RPPA_g_m)=c("marker","Gene","bcr_patient_barcode","expression","cancer")
RPPA_g_m = RPPA_g_m[RPPA_g_m$Gene %in% pathVarP$HUGO_Symbol,]
RPPA_g_m$expression = as.numeric(RPPA_g_m$expression)
##### individual cancer type analysis #####
cancers = unique(pathVarP$cancer)
markers = unique(RPPA_g_m$marker)
# limit runs to cancers with at least 3 likely patho/pathogenic variants 
tt = NULL
for (marker in markers){
  RPPA_g_m_g = RPPA_g_m[RPPA_g_m$marker==marker,]
  gene = RPPA_g_m_g$Gene[1]
  gene_sample_g = gene_sample[gene_sample$Gene==gene,]
  var_exp_g = merge(RPPA_g_m_g,gene_sample_g,by="bcr_patient_barcode",all.x=T)
  var_exp_g$Freq[is.na(var_exp_g$Freq)] = 0
  var_exp_g$Freq[var_exp_g$Freq != 0 ] = 1
  for (cancer in cancers){
    var_exp_g_c = var_exp_g[var_exp_g$cancer %in% cancer,]
    gene_path_count = sum(var_exp_g_c$Freq)
    if (gene_path_count > 2){
      # run GLM
      w = wilcox.test(var_exp_g_c$expression[var_exp_g_c$Freq==0],var_exp_g_c$expression[var_exp_g_c$Freq!=0])
      wP = w$p.value
      wWstat = w$statistic
      cancer_gene_stat = run_glm(var_exp_g_c)
      # compile results
      full_cancer_gene_stat = cbind(cancer,marker,gene,gene_path_count,wP,wWstat,cancer_gene_stat)
      tt = rbind(tt, full_cancer_gene_stat)
    }
  }
}

#"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
colnames(tt) = c("cancer","marker","gene","gene_path_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F_statistic","p-value","coefficient","covariants");
tt$FDR = p.adjust(tt[,"p-value"], method="fdr") # MAW new, calculates FDR based on the method from,
# Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289–300.
tt$wilcoxFDR = p.adjust(tt[,"wilcoxP"], method="fdr")

tt=tt[order(tt$FDR, decreasing=FALSE),]
tn = "out/pathVarRPPAAssoc.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

### plotting ###
tt$gene = as.character(tt$gene)
tt$marker = as.character(tt$marker)
tt$FDR_plot = tt$FDR
tt$FDR_plot[tt$FDR_plot<10^(-6)]= 0.95*10^(-6)
p = ggplot(data=tt)
p = p + geom_point(aes(y=-log10(FDR_plot),x= coefficient,color = cancer),alpha=0.5)
#p = p + geom_text_repel(aes(y=-log10(FDR),x= cohort_AF,label=ifelse(FDR<0.05, Gene,NA),color = Cancer))#,alpha=1.3)
#p = p + geom_point(aes(y=cohort_AF,x=Cancer,size=-log10(FDR),color = Cancer))
p = p + geom_text_repel(aes(y=-log10(FDR_plot),x= coefficient,color = cancer,label=ifelse(FDR<0.05, marker,NA)))
p = p + getPCACancerColor()
p = p + labs(x="Coefficient",y= "-log10(FDR)")
p = p + geom_vline(xintercept = 0, alpha=0.5)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/RPPAAssocVolcanoGLM.pdf'
ggsave(fn,w = 5, h = 5, useDingbat=F)

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
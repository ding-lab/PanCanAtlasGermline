##### somaticVsMutationSignature.R #####
# Kuan-lin Huang 2018

source("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/global_aes_out.R")
source("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/dependency_files.R")

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
  ytype="Q";yi="Contribution";xi="Freq";covi=covi;#covi=md[i,4];memo=md[i,5]
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
somaticDriver299_f = "/Users/khuang/Box Sync/HuangLab/reference/Driver_BaileyCell2018/299driverGene.txt"
somaticDriver299 = as.vector(t(read.table(header=F, quote = "", sep="\t", file = somaticDriver299_f, stringsAsFactors=FALSE)))

somatic_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/somatic/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz"
somatic = read.table(header=T, quote = "", sep="\t", file = gzfile(somatic_f), stringsAsFactors=FALSE)
somatic$bcr_patient_barcode = gsub("(^TCGA-[A-Z0-9][A-Z0-9]-[A-Z0-9][A-Z0-9][A-Z0-9][A-Z0-9])-.*","\\1",somatic$Tumor_Sample_Barcode)
somatic_mut_count = data.frame(table(somatic$bcr_patient_barcode))
colnames(somatic_mut_count) = c("bcr_patient_barcode","MutationCount")

table(somatic$Variant_Classification)
likelyFunctionalTypes = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation",
                          "Nonsense_Mutation","Splice_Site","Translation_Start_Site")
somatic_likelyfunctional = somatic[somatic$Hugo_Symbol %in% c(pathVar$HUGO_Symbol,somaticDriver299) & somatic$Variant_Classification %in% likelyFunctionalTypes,]

# read input files
gene_sample = data.frame(table(somatic_likelyfunctional$Hugo_Symbol,somatic_likelyfunctional$bcr_patient_barcode))
colnames(gene_sample) = c("Gene","bcr_patient_barcode","Freq")
gene_sample_m = merge(gene_sample,somatic_mut_count, by="bcr_patient_barcode")

mut_signature_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/signature_profile_sample.txt.gz"
mut_signature = read.table(sep="\t",header=T,file=gzfile(mut_signature_f), stringsAsFactors=FALSE)
mut_signature_tcga = mut_signature[mut_signature$Country=="United States",]
mut_signature_tcga$cancer = gsub("-US","",mut_signature_tcga$project_code)
mut_signature_tcga$Signature = gsub("\\.","-",mut_signature_tcga$Signature)
mut_signature_tcga$bcr_patient_barcode = gsub("\\.","-",mut_signature_tcga$Tumor_Sample_Barcode)

mut_signature_tcga_brief = mut_signature_tcga[,c("cancer" , "Signature", "Contribution", "bcr_patient_barcode")]

##### individual cancer type analysis #####
cancers = unique(pathVarP$cancer)
genes = unique(pathVarP$HUGO_Symbol)
# limit runs to cancers with at least 5 likely patho/pathogenic variants 
tt = NULL
for (gene in genes){
  gene_sample_g = gene_sample_m[gene_sample_m$Gene==gene,]
  var_exp_g = merge(mut_signature_tcga_brief,gene_sample_g,by="bcr_patient_barcode",all.x=T)
  var_exp_g$Freq[is.na(var_exp_g$Freq)] = 0
  var_exp_g$Freq[var_exp_g$Freq != 0 ] = 1
  for (cancer in cancers){
    var_exp_g_c = var_exp_g[var_exp_g$cancer %in% cancer,]
    gene_path_count = sum(var_exp_g_c$Freq[var_exp_g_c$Signature=="Signature-1"])
    if (gene_path_count > 2){
      for (signature in unique(var_exp_g_c$Signature)){
        var_exp_g_c_s = var_exp_g_c[var_exp_g_c$Signature == signature,]
        # run GLM
        w = wilcox.test(var_exp_g_c_s$Contribution[var_exp_g_c_s$Freq==0],var_exp_g_c_s$Contribution[var_exp_g_c_s$Freq!=0])
        wP = w$p.value
        wWstat = w$statistic
        cancer_gene_stat = run_glm(var_exp_g_c_s,covi="MutationCount")
        # compile results
        full_cancer_gene_stat = cbind(cancer,gene,signature,gene_path_count,wP,wWstat,cancer_gene_stat)
        tt = rbind(tt, full_cancer_gene_stat)
      }
    }
  }
}

#"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
colnames(tt) = c("cancer","gene","signature","gene_path_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F_statistic","p-value","coefficient","covariants");
tt$FDR = p.adjust(tt[,"p-value"], method="fdr") # MAW new, calculates FDR based on the method from,
# Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289–300.
tt$wilcoxFDR = p.adjust(tt[,"wilcoxP"], method="fdr")

tt=tt[order(tt$FDR, decreasing=FALSE),]
tn = "out/somaticMutMutsigAssoc.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

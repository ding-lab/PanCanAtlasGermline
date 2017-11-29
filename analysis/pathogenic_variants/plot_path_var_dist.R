##### plot_path_var_dist.R #####
# Kuan-lin Huang @ WashU 2017 Oct
# plot the distribution of pathogenic variants

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/pathogenic_variants"
setwd(bdir)
source("../dependency_files.R")
source("../global_aes_out.R")

## functions 
plot_barplot = function(matrix, x_string, fill_string=NULL, fileName="data.pdf"){
  
  fn = paste("out", fileName,sep ="/")
  
  if (is.null(fill_string)){
    p = ggplot(matrix,aes_string(x = x_string))
  } else{
    p = ggplot(matrix,aes_string(x = x_string, fill = fill_string))
  }
  
  p = p + geom_bar() + theme_bw() + theme_nogrid()
  p = p + labs(x = x_string, y="counts")
  p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
  p
  ggsave(file=fn, useDingbats=FALSE)
}

plot_heatmap = function(data){
  getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
  p = ggplot(data=data)
  p = p + facet_grid(.~classification,drop=T,scales = "free", space = "free")
  p = p + geom_tile(aes(x=cancer, y=gene, fill=count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
  p = p + geom_text(aes(x=cancer, y=gene, label = count, stringsAsFactors=FALSE), color="black", size=3)
  p = p  + theme_bw() + theme_nogrid() +
    theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=14),axis.ticks = element_blank())#element_text(colour="black", size=14))
  return(p)
}

# other files
cancer_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/sample_listing/out/cancer_count_wethni_waao.txt" 
geneList_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/reference_files/20160713_Rahman_KJ_KH_152_gene_table_list.txt"

##### READ inputs #####
variants = pathVarP

cancer_clin = read.table(header=T,sep="\t",file = cancer_f)
cancer_count = cancer_clin[,c("Cancer", "Sample.size")]
colnames(cancer_count)[2] = "sample_size"

geneList_t = read.table(header=F, file = geneList_f)
geneList = as.vector(t(geneList_t))

##### preliminary plotting #####

# plot_barplot(variants,x_string = "cancer", fill_string = "Overall_Classification", 
#              fileName="PCA_PathVar_Classification_by_cancer_bar.pdf")


# Frequency plotting #
# get non-redundant sample entry first
variants = variants[order(variants$Overall_Classification,decreasing = T),]
variants_uni = variants[!duplicated(variants$Sample),]

# plot the bar plot of the counts
all_count = as.data.frame(table(variants_uni$cancer,variants_uni$Overall_Classification))
colnames(all_count) = c("Cancer", "Overall_Classification", "count")
all_count_m = merge(cancer_count, all_count, by="Cancer")
all_count_m$freq = all_count_m$count/all_count_m$sample_size
all_count_m$Overall_Classification = factor(all_count_m$Overall_Classification, levels = c("Pathogenic","Likely Pathogenic"))
all_count_m = all_count_m[order(all_count_m$Overall_Classification),]

cat("Number of pathogenic variant carriers: ",sum(all_count_m$count[all_count_m$Overall_Classification=="Pathogenic"]),"\n")
cat("Additional number of likely pathogenic variant carriers: ",sum(all_count_m$count[all_count_m$Overall_Classification=="Likely Pathogenic"]),"\n")

cat("Percentage of pathogenic variant carriers: ",sum(all_count_m$count[all_count_m$Overall_Classification=="Pathogenic"])/sum(all_count_m$sample_size[all_count_m$Overall_Classification=="Pathogenic"]),"\n")
cat("Additional percentage of likely pathogenic variant carriers: ",sum(all_count_m$count[all_count_m$Overall_Classification=="Likely Pathogenic"])/sum(all_count_m$sample_size[all_count_m$Overall_Classification=="Likely Pathogenic"]),"\n")
tn = "PCA.path.carrier_by_cancer.freq.tsv"
write.table(all_count_m, quote=F, sep="\t", file = tn, row.names = F)

all_count_m_sum = dcast(all_count_m, Cancer ~ Overall_Classification, value.var="freq")
cancer_order = all_count_m_sum$Cancer[order(all_count_m_sum$Pathogenic, decreasing=TRUE)]

all_count_m$Cancer = factor(all_count_m$Cancer, levels=cancer_order)
all_count_m= all_count_m[all_count_m$freq !=0,]
all_count_m$Overall_Classification = factor(all_count_m$Overall_Classification, levels=c("Likely Pathogenic","Pathogenic"))

p = ggplot(all_count_m,aes(x = Cancer, y = freq*100, fill = Overall_Classification))
p = p + geom_bar(stat = "identity") + theme_bw() + theme_nogrid() 
p = p + labs(x = "Cancer Type", y="Percentage Carriers (%)")
p = p + getVarFillScale()
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "top")
p
fn = 'out/PCA_carrierFreq_by_cancer_bar.pdf'
ggsave(file=fn, height = 6, width = 6.5, useDingbats=FALSE)

##### compare counts to Lu et al. #####
# lu et al. discovery strategy
lu_f = "/Users/khuang/Box\ Sync/PhD/germline/pan4000/ncomms10086-s3_lu_et_al_838var.txt"
lu = read.table(header=T,file=lu_f, sep="\t")

lu_34_burden_gene_f = "/Users/khuang/Box\ Sync/PhD/germline/pan4000/burden_test/34_FDR0.15genes_oncoRemoved.list"
lu_34_burden_gene = as.vector(t(read.table(header=F,lu_34_burden_gene_f)))

lu_sample_f = "/Users/khuang/Box\ Sync/PhD/germline/pan4000/mis/figure5c_1.data.tsv"
lu_sample = as.vector(read.table(header=T,lu_sample_f)[,1])

# variants_p_34 = variants_p[(variants_p$sample %in% lu_sample) & (variants_p$gene_name %in% lu_34_burden_gene) & (variants_p$binary_type=="Truncation"),]
variants_in_lu = variants[(variants$bcr_patient_barcode %in% lu_sample) & (variants$Start %in% lu$start),]
cat("Total numbers of pathogenic variants:", nrow(variants),"\n")
cat("Total numbers of pathogenic variants reported in Lu et al:", nrow(variants_in_lu),"\n")

lu_count = as.data.frame(table(variants_in_lu$cancer))
colnames(lu_count) = c("cancer","previously_reported")
huang_count = as.data.frame(table(variants$cancer))
colnames(huang_count) = c("cancer","total")
both_counts = merge(lu_count,huang_count,by="cancer",all=T)
both_counts$previously_reported[is.na(both_counts$previously_reported)]=0
both_counts$novel = both_counts$total - both_counts$previously_reported
both_counts_m = melt(both_counts,id.vars =c("cancer","total"))
both_counts_m$cancer = factor(both_counts_m$cancer, levels=cancer_order)
both_counts_m$variable = factor(both_counts_m$variable, levels=c("novel","previously_reported"))

p = ggplot(both_counts_m,aes(x = cancer, y = value, fill = variable))
p = p + geom_bar(stat = "identity") + theme_bw() + theme_nogrid() 
p = p + labs(x = "Cancer Type", y="Pathogenic Variant Count")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "top")
p
fn = 'out/previous_vs_novel_var_count.pdf'
ggsave(file=fn, height = 6, width = 6, useDingbats=FALSE)

sample_size_order = c(all_count_m$sample_size[all_count_m$Overall_Classification=="Pathogenic"][order(all_count_m$freq[all_count_m$Overall_Classification=="Pathogenic"],decreasing = T)],"45","80")
all_count_m$sample_size = factor(all_count_m$sample_size, levels=sample_size_order)
##### v2: vertical plotting #####
p = ggplot(all_count_m,aes(x = Cancer, y = freq*100, fill = Overall_Classification))
p = p + geom_bar(stat = "identity") + theme_bw() + theme_nogrid() 
p = p + labs(x = "Cancer Type", y="Percentage Carriers (%)")
p = p + getVarFillScale()
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "top")
p = p + coord_flip() 
p
fn = 'out/PCA_carrierFreq_by_cancer_bar_vertical.pdf'
ggsave(file=fn, height = 6, width = 6.5, useDingbats=FALSE)

p = ggplot(both_counts_m,aes(x = cancer, y = value, fill = variable))
p = p + geom_bar(stat = "identity") + theme_bw() + theme_nogrid() 
p = p + labs(x = "Cancer Type", y="Pathogenic Variant Count")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "top")
p = p + coord_flip() 
p
fn = 'out/previous_vs_novel_var_count_verticle.pdf'
ggsave(file=fn, height = 6, width = 6, useDingbats=FALSE)

##### heatmap: count of pathogenic variants by gene/cancer #####
combined_sum_f_added = data.frame(table(variants$HUGO_Symbol, variants$cancer))
colnames(combined_sum_f_added) = c("Gene","Cancer","Count")
combined_sum_f_added = combined_sum_f_added[combined_sum_f_added$Count != 0,]
combined_sum_f_added$Gene_Classification = NA
combined_sum_f_added$Gene_Classification[combined_sum_f_added$Gene %in% all_oncogenes] = "Oncogene"
combined_sum_f_added$Gene_Classification[combined_sum_f_added$Gene %in% all_TSGs] = "TSG"

top = table(variants$HUGO_Symbol)[order(table(variants$HUGO_Symbol),decreasing = T)]
top5 = top[top>4]
top_oncogenes = names(top[top>4])[names(top[top>4]) %in% all_oncogenes]
top_TSGs = names(top[top>13])[names(top[top>13]) %in% all_TSGs]
# write(c(top_oncogenes,top_TSGs), "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/reference_files/PCA_feature_gene_list.txt", sep="\n")
# use the genes defined by burden analysis
gene_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/reference_files/PCA_feature_gene_list.txt"
glist = as.vector(t(read.table(header=FALSE, stringsAsFactors = F, file = gene_fn)))

combined_sum_f_added_g = combined_sum_f_added[combined_sum_f_added$Gene %in% glist,]
#combined_sum_f_added_g$Cancer = factor(combined_sum_f_added_g$Cancer, levels=cancer_order) # site from plot_classification_summary.R
combined_sum_f_added_g$Gene = factor(combined_sum_f_added_g$Gene, levels = rev(gene_order))

getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
combined_sum_f_added_g$Count_plot = combined_sum_f_added_g$Count
combined_sum_f_added_g$Count_plot[combined_sum_f_added_g$Count_plot>10]= 10

p = ggplot(data=combined_sum_f_added_g)
#p = p + facet_grid(Gene_Classification~.,drop=T,scales = "free_y", space = "free_y")
p = p + geom_tile(aes(x=Cancer, y=Gene, fill=Count_plot), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
p = p + geom_text(aes(x=Cancer, y=Gene, label = Count), color="black", size=3)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/PathVar_counts_ggroup_heatmap.pdf'
ggsave(file=fn, height=5, width=10, useDingbats=FALSE)


##### examine_gene_list.R #####
# Kuan-lin Huang @ WashU 201802
# examine the curated gene list

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/gene_list"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

library(readxl)

### annotate 152 gene table with oncogene/TSG
CPG_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/doc/reference/CPG152_table.xlsx"
CPG = data.frame(readxl::read_excel(CPG_fn))
colnames(CPG)[1]="Gene"
CPG$Gene_Classification = "Not classified"
CPG$Gene_Classification[CPG$Gene %in% all_TSGs] = "Tumor Suppressor Gene"
CPG$Gene_Classification[CPG$Gene %in% all_oncogenes] = "Oncogene"
# write.table(CPG[c("Gene","Gene_Classification")], file="/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/doc/reference/CPG152_table_gClass.tsv", quote=F, sep="\t", col.names=T, row.names=F)

### examine the number of pathogenic variants found in each category of curated genes
cat("Source of 152 CPGs","\n")
table(CPG$Source)

gene_var_count = data.frame(table(pathVarP$HUGO_Symbol,pathVarP$Overall_Classification))
colnames(gene_var_count) = c("Gene","Classification","Count")
# write.table(CPG_count[c("Gene","Classification", "Count")], file="/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/doc/reference/CPG152_table_gVarCount.tsv", quote=F, sep="\t", col.names=T, row.names=F)

gene_var_count2= data.frame(table(pathVarP$HUGO_Symbol))
colnames(gene_var_count2) = c("Gene","Count")
CPG_count = merge(CPG,gene_var_count2, by="Gene", all.x=T)
CPG_count$Count[is.na(CPG_count$Count)] = 0
CPG_count$Source[CPG_count$Source=="Cancer Gene Census Germline download 1/5/2016 (http://cancer.sanger.ac.uk/census/ )"] = "Cancer Gene Census Germline"
CPG_count$Source[CPG_count$Source=="Reference (see PMID)"] = "Curated from literature"
CPG_count$Source[CPG_count$Source=="personal communication; related to DICER1"] = "Personal communication"

p = ggplot(CPG_count,aes(x=Source, y = Count ))
#p = p + facet_grid(.~Classification)
p = p + geom_jitter(height =0, width = 0.3, alpha=0.5,aes(fill=Source, color=Source))
p = p + geom_violin(alpha=0.5, stroke=0,aes(fill=Source, color=Source))
p = p + geom_label_repel(aes(label=ifelse(Source != "Rahman 114 CPG" & Count > 5, paste(Gene,Count,sep="-"),NA)))
p = p  + theme_bw() + labs(y="Count of likely pathogenic and pathogenic variant", x = "Source of gene") +
  theme(legend.position = "None", axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p 

fn = 'out/source_var_count.pdf'
ggsave(fn,h=7,w=4,useDingbat=F)
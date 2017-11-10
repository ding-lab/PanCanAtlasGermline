##### plotPathVarEthnic.R #####
# Kuan-lin Huang @ WashU 201711
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/clinical_association"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")


gene_ethni_count = data.frame(table(pathVar$HUGO_Symbol,pathVar$AIM_ethnicity,pathVar$cancer))
colnames(gene_ethni_count) = c("Gene","Ethnicity","Cancer","CarrierCount")

clin_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/PanCan_ClinicalData_V4_wAIM.txt"
clin = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)

s_c_list_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/sampleQC/pca_table.20171019.wclin.tsv"
sample_cancer = read.table(header=T, quote = "", sep="\t", file = s_c_list_f, stringsAsFactors=FALSE)
sample_cancer = sample_cancer[,c("bcr_patient_barcode", "cancer")]

clin_sample_cancer = merge(sample_cancer,clin,by="bcr_patient_barcode",all.x=T)

cancer_ethni_count = data.frame(table(clin_sample_cancer$consensus_call,clin_sample_cancer$cancer))
colnames(cancer_ethni_count) = c("Ethnicity","Cancer","CohortCount")

gene_cancer_ethni_count = merge(gene_ethni_count,cancer_ethni_count,by=c("Ethnicity","Cancer"))

count.F_25 = gene_cancer_ethni_count[gene_cancer_ethni_count$CohortCount >= 25,]

### work on ethnicity label
count.F_25$Ethnicity = as.character(count.F_25$Ethnicity)
count.F_25$Ethnicity[count.F_25$Ethnicity == "afr"] = "African American" 
count.F_25$Ethnicity[count.F_25$Ethnicity == "asian"] = "Asian" 
count.F_25$Ethnicity[count.F_25$Ethnicity == "eur"] = "Europian"
count.F_25$Ethnicity[count.F_25$Ethnicity == "amr"] = "American"
count.F_25$Cohort_num = paste(count.F_25$Ethnicity, "n =", count.F_25$CohortCount)
count.F_25$Frequency = 100*count.F_25$CarrierCount/count.F_25$CohortCount

##### plotting #####
count.F_25_g = count.F_25[count.F_25$Gene %in% featGenes,]
p = ggplot(count.F_25_g,aes(x=Cohort_num, y=Frequency, fill=Gene))
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_bar(stat="identity") + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="% of samples with pathogenic variants")
p = p + scale_fill_manual(values =  col_vector)
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,angle=90, vjust=0.5, hjust=0), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
fn = "out/PathVarCarrierByEthni.pdf"
ggsave(file=fn, height=6, width=15, useDingbats=FALSE)

### bubble plot ###

p = ggplot(count.F_25_g,aes(x=Cohort_num, y=Gene, colour=Ethnicity, size =Frequency))# make this the original ethni
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_point() + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="% carriers")
p = p + scale_colour_manual(values =  c(col_vector[1:(length(unique(count.F_25_g$Gene))+1)]))#brewer.pal(length(unique(count.F_25_g$Gene)),"Set3"))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
fn = "out/PathVarCarrierByEthniBubble.pdf"
ggsave(file=fn, height=7, width=15, useDingbats=FALSE)

top_gene = c()
for (gene in unique(count.F_25$Gene)){
  count = sum(count.F_25$CarrierCount[count.F_25$Gene == gene])
  if (count > 20){top_gene = c(top_gene, gene)}
}

count.F_25$Gene = as.character(count.F_25$Gene)

count.F_25_go = count.F_25
count.F_25_go[!(count.F_25_go$Gene %in% top_gene),]$Gene = "other"
for (cancer in unique(count.F_25_go$Cancer)){
  for (ethni in unique(count.F_25_go[count.F_25_go$Cancer == cancer,]$Ethnicity)){
    all_other = count.F_25_go[count.F_25_go$Cancer == cancer & count.F_25_go$Ethnicity == ethni &
                    count.F_25_go$Gene == "other",]
    s_count = as.numeric(all_other$CohortCount[1])
    num_carrier = sum(as.numeric(all_other$Carrier))
    Cohort_num = paste(ethni, "n =", s_count)
    new_row = c(ethni, cancer, "AllOther", num_carrier, s_count, Cohort_num, num_carrier*100/s_count)
    count.F_25_go = rbind(count.F_25_go, new_row)
  }
}
count.F_25_go = count.F_25_go[count.F_25_go$Gene != "other",]
count.F_25_go$Frequency = as.numeric(count.F_25_go$Frequency)

count.F_25_go = count.F_25_go[count.F_25_go$CarrierCount > 0,]

p = ggplot(count.F_25_go,aes(x=Cohort_num, y=Frequency, fill=Gene))
p = p + facet_grid(.~Cancer, drop=T, space = "free_x",scales = "free_x")#, space = "free", scales = "free")
p = p + geom_bar(stat="identity") + theme_bw() + theme_nogrid()
p = p + labs(x = "", y="% carriers")
#p = p + scale_fill_manual(values =  c("#d9d9d9",col_vector[1:(length(unique(count.F_25$Gene)))]))#brewer.pal(length(unique(count.F_25_g$Gene)),"Set3"))
p = p + scale_fill_manual(values =  c("#d9d9d9","#fc9272","#bcbddc",brewer.pal(length(unique(count.F_25$Gene))-1,"Paired")))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p
fn = 'out/PathVarCarrierByEthniOther.pdf'
ggsave(file=fn, height=6, width=15, useDingbats=FALSE)
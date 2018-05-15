##### somatic_germline_overlap.R #####
# Kuan-lin Huang @ WashU 2018
# Find overlap of genes/variants for somatic/germline variants

### dependencies ###
source("../global_aes_out.R")
source("../dependency_files.R")
source("../load_somatic.R")

# counts of somatic functional mutation by gene
somatic_gene_count = data.frame(table(somatic_likelyfunctional_driver$Hugo_Symbol))
germline_gene_count = data.frame(table(pathVarP$HUGO_Symbol))
colnames(somatic_gene_count) = c("Gene","PredictedFunctionalSomaticMutationCount")
colnames(germline_gene_count) = c("Gene","PathogenicGermlineVariantCount")
gene_count = merge(somatic_gene_count,germline_gene_count,by="Gene",all=T)
gene_count[is.na(gene_count)] = 0
highlight_g = as.character(gene_count$Gene[gene_count$PredictedFunctionalSomaticMutationCount > 400 | gene_count$PathogenicGermlineVariantCount > 20
                                           | (gene_count$PredictedFunctionalSomaticMutationCount > 140 & gene_count$PathogenicGermlineVariantCount > 3)])
gene_count$GeneClass = "Others"
gene_count$GeneClass[gene_count$Gene %in% all_oncogenes] = "Oncogene"
gene_count$GeneClass[gene_count$Gene %in% all_TSGs] = "TSG"

p = ggplot(gene_count,aes(y=PredictedFunctionalSomaticMutationCount, x =PathogenicGermlineVariantCount, color = GeneClass))
p = p + geom_point(stroke=0,alpha = 0.2) + theme_bw()  #+ guides(color=FALSE)
#p = p + geom_abline(intercept = 0, slope=1, alpha=0.2) #+ geom_density2d(alpha=0.5)
p = p + geom_text_repel(aes(label=ifelse(as.character(Gene) %in% highlight_g,as.character(Gene), NA)))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
#p = p + scale_x_log10() + scale_y_log10()
p = p + expand_limits(x = 0,y=0) + ylim(0,1100)
#p = p + coord_equal() + getLOHColorScale()
p
fn = "out/somatic_vs_germline_var_counts_by_gene.pdf"
ggsave(file=fn, width=5, h =5, useDingbats=FALSE)


# by cancer
somatic_cancer_gene_count = data.frame(table(somatic_likelyfunctional_driver$Hugo_Symbol, somatic_likelyfunctional_driver$cancer))
germline_cancer_gene_count = data.frame(table(pathVarP$HUGO_Symbol,pathVarP$cancer))
colnames(somatic_cancer_gene_count) = c("Gene","Cancer","PredictedFunctionalSomaticMutationCount")
colnames(germline_cancer_gene_count) = c("Gene","Cancer","PathogenicGermlineVariantCount")
cancer_gene_count = merge(somatic_cancer_gene_count,germline_cancer_gene_count,by=c("Gene","Cancer"),all=T)
cancer_gene_count[is.na(cancer_gene_count)] = 0
highlight_g = as.character(cancer_gene_count$Gene[cancer_gene_count$PredictedFunctionalSomaticMutationCount > 400 | cancer_gene_count$PathogenicGermlineVariantCount > 20
                                           | (cancer_gene_count$PredictedFunctionalSomaticMutationCount > 140 & cancer_gene_count$PathogenicGermlineVariantCount > 3)])
cancer_gene_count$GeneClass = "Others"
cancer_gene_count$GeneClass[cancer_gene_count$Gene %in% all_oncogenes] = "Oncogene"
cancer_gene_count$GeneClass[cancer_gene_count$Gene %in% all_TSGs] = "TSG"

p = ggplot(cancer_gene_count,aes(y=PredictedFunctionalSomaticMutationCount, x =PathogenicGermlineVariantCount, color = GeneClass))
p = p + facet_wrap(~Cancer)
p = p + geom_point(stroke=0,alpha = 0.2) + theme_bw()  #+ guides(color=FALSE)
#p = p + geom_abline(intercept = 0, slope=1, alpha=0.2) #+ geom_density2d(alpha=0.5)
p = p + geom_text_repel(aes(label=ifelse(as.character(Gene) %in% highlight_g,as.character(Gene), NA)),size=1)
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + scale_x_log10() + scale_y_log10()
p = p + expand_limits(x = 0,y=0) #+ ylim(0,1100)
#p = p + coord_equal() + getLOHColorScale()
p
fn = "out/somatic_vs_germline_var_counts_by_gene_by_cancer.pdf"
ggsave(file=fn, width=10, h =10, useDingbats=FALSE)

##### label_onco_var_ExAC.R #####
# Kuan-lin Huang @ WashU 2017 Oct
# find non-cancer pathogenic variant in the ExAC cohort

bdir = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/burden_assoc"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

fn = "charged.ExAC.r1.sites.vep.biallelic.combine.exon.all.patho.expanded.tsv"
variants = read.table(sep="\t",header=T,file=fn, stringsAsFactors=FALSE, quote = "",fill=TRUE)

gene_fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/reference_files/20160713_Rahman_KJ_KH_152_gene_table_list.txt"
predisposition_genes = as.vector(t(read.table(sep="\t",header=F,file=gene_fn, stringsAsFactors=FALSE, quote = "")))

cat("Original count of variants: ",sum(variants$ExAC_nonTCGA_AC_Adj),"\n")

##### classify whether a variant is cancer-relevant #####
cancer_terms = c("tumor","cancer","neoplasia")

variants$predisposition_gene = F
variants$predisposition_gene[variants$HUGO_Symbol %in% predisposition_genes] = TRUE
variants$cancer_term_trait = FALSE
for (term in cancer_terms){
  variants$cancer_term_trait[grep(term,tolower(variants$ClinVar_Traits))] = TRUE
}
variants$cancer_term_trait[grep("oma$",tolower(variants$ClinVar_Traits))] = TRUE

table(variants$predisposition_gene,variants$cancer_term_trait)
variants$cancer_related = F
variants$cancer_related[variants$predisposition_gene | variants$cancer_term_trait] = T

#table(variants$ClinVar_Traits[variants$cancer_term_trait])[table(variants$ClinVar_Traits[variants$cancer_term_trait])>3]

# variant frequency annotation

#### rare frequency filter #####
variants = variants[variants$ExAC_AF < 0.0005,] # not less

variants_cancer = variants[variants$cancer_related,]
cat("Number of these variants that are cancer-relevant:",sum(variants_cancer$ExAC_nonTCGA_AC_Adj),"\n")

variants_cancer$Overall_Classification = "Uncertain Significance"
variants_cancer$Overall_Classification[variants_cancer$CharGer_Classification=="Pathogenic"] = "Likely Pathogenic"
variants_cancer$Overall_Classification[variants_cancer$ClinVar_Pathogenicity=="Pathogenic"] = "Pathogenic"
variants_cancer$Overall_Classification[grep("PS1",variants_cancer$Positive_Evidence)] = "Pathogenic"
variants_cancerP = variants_cancer[variants_cancer$Overall_Classification %in% c("Likely Pathogenic", "Pathogenic"),]

cat("Number of these variants that are cancer-relevant and pathogenic:",sum(variants_cancerP$ExAC_nonTCGA_AC_Adj),"\n")

# tn = "out/ExAC_pathogenic_variants.tsv"
# write.table(variants_cancerP, quote=F, sep="\t", file = tn, row.names = F)

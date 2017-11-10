##### dependency_files.R #####
# Kuan-lin Huang @ WashU 2017 June
# dependent files for analysis in the PCA Germline project

# gene lists
volg_fn = "/Users/khuang/Box\ Sync/PhD/proteogenomics/reference_files/Volgestin2013Science_125genes_class.txt"
volg_class = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=volg_fn)
colnames(volg_class) = c("gene_name","Gene_Classification","Pathway","Process")
volg_TSGs = volg_class$gene_name[volg_class$Gene_Classification=="TSG"]
volg_oncogenes = volg_class$gene_name[volg_class$Gene_Classification=="Oncogene"]

onco_fn = "/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/reference_files/GSEA_geneLists/oncogenes.txt"
oncogenes = as.vector(t(read.table(header=F, stringsAsFactors = F, file=onco_fn)))

tsg_fn = "/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/reference_files/GSEA_geneLists/tumor_suppressors.txt"
TSGs = as.vector(t(read.table(header=F, stringsAsFactors = F, file=tsg_fn)))

additional_TSGs = c("ATR","BARD1","ERCC1","FANCI","FANCL","FANCM","POLD1","POLE","POLH","RAD50","RAD51","RAD51C","RAD51D","RAD54L")
others = c("ABCB11","ABCC2","AXIN2","CBL","CDKN1B","COL7A1","CYP17A1","CYP1B1","DIS3L2",
           "DKC1","DOCK8","ELANE","FAH","FLCN","GBA","GJB2","HFE","HMBS","LRRK2",
           "MAX","MTAP","MYD88"   
           ,"PRSS1","RHBDF2","RPL5","SDHA","SETBP1","SF3B1","SH2D1A","SLC25A13","SOS1","TMEM127","TRIM37","UROD")

additional_oncogenes = c("AR","STAT3","TERT","MAP2K2","NOTCH3")

all_oncogenes = c(oncogenes,volg_oncogenes,additional_oncogenes)
all_TSGs = c(TSGs,volg_TSGs,additional_TSGs)

gene_fn = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/reference_files/PCA_feature_gene_list.txt"
glist_f = read.table(header=FALSE, stringsAsFactors = F, file = gene_fn)
featGenes = as.vector(t(glist_f))

# get specific gene genera: fanconi; rtk; etc
# fanc: https://humgenomics.biomedcentral.com/articles/10.1186/s40246-015-0054-y
fanc_gene_fn = "/Users/khuang/Box\ Sync/PhD/germline/pan8000_germline_clinical/reference_files/human_fanc_genes.txt"
fanc_gene_f = read.table(header=T, sep="\t", quote="", stringsAsFactors = F, file = fanc_gene_fn)
fanc_genes = gsub("\xa0","",fanc_gene_f[,1])

rtk_fn = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/reference_files/RTKs_list.txt"
rtk_gene_f = read.table(header=FALSE, stringsAsFactors = F, file = rtk_fn)
rtk_genes = as.vector(t(rtk_gene_f))

# ##### all_variants #####

fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/data_integration/out/PCA_pathVar_integrated.tsv"
pathVar = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=fn)

pathVar$ExAC_assoc_P[pathVar$ExAC_adj_AF>=0.001] = NA # not considering anything above ExAC AF of 0.1%; doesn't matter

# volg_fn = "/Users/khuang/Box\ Sync/PhD/proteogenomics/reference_files/Volgestin2013Science_125genes_class.txt"
# volg_class = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=volg_fn)
# colnames(volg_class) = c("gene_name","Gene_Classification","Pathway","Process")
# pathVar = merge(pathVar,volg_class,by="gene_name",all.x=T)
# pathVar$Gene_Classification[is.na(pathVar$Gene_Classification)]="TSG"

# some pre-processing

pathVar$LOH_Sig = "None"
pathVar$LOH_Sig[!is.na(pathVar$LOH_FDR) & (pathVar$LOH_FDR < 0.15)] = "Suggestive"
pathVar$LOH_Sig[!is.na(pathVar$LOH_FDR) & (pathVar$LOH_FDR < 0.05)] = "Significant"
pathVar$LOH_Sig = factor(pathVar$LOH_Sig, levels = c("Significant","Suggestive","None"))

pathVar$Allele_Frequency_Plot = as.numeric(pathVar$Cohort_AF)*100
pathVar$Allele_Frequency_Plot[is.na(pathVar$Allele_Frequency_Plot)] = 0

pathVar$Gene_Classification = NA
pathVar$Gene_Classification[pathVar$HUGO_Symbol %in% all_oncogenes] = "Oncogene"
pathVar$Gene_Classification[pathVar$HUGO_Symbol %in% all_TSGs] = "TSG"

##### subsets #####
pathVarOT = pathVar[!is.na(pathVar$Gene_Classification) & pathVar$Gene_Classification != "None",]

truncations = pathVarOT[pathVarOT$binary_type=="Truncation",]
missenses = pathVarOT[pathVarOT$binary_type=="Missense",]

pathVarFGene = pathVar[pathVar$HUGO_Symbol %in% featGenes,]


### combined classification ###
table(pathVar$CharGer_Classification)
table(pathVar$ACMG_Classification)
table(pathVar$ClinVar_Pathogenicity)

table(pathVar$ClinVar_Pathogenicity,pathVar$CharGer_Classification)
table(pathVar$ACMG_Classification,pathVar$CharGer_Classification)
table(pathVar$ClinVar_Pathogenicity,pathVar$ACMG_Classification)

pathVar$Overall_Classification = paste("CharGer",pathVar$CharGer_Classification)
pathVar$Overall_Classification[pathVar$ClinVar_Pathogenicity=="Pathogenic"] = "Pathogenic"
table(pathVar$Overall_Classification)
pathVar$Overall_Classification[grep("PS1",pathVar$Positive_Evidence)] = "Pathogenic"
table(pathVar$Overall_Classification)
pathVar$CharGer_Classification = pathVar$Overall_Classification
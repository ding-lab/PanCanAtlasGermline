##### analyzeAssocTFT.DDR.R #####
# Kuan-lin Huang @ WashU 2017 Oct. 

setwd("/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/association_test")
source("../global_aes_out.R")
source("TFT_functions.R")

fn = "assoc_results/ExAC.r1.sites.vep.biallelic.combine.fisher.anno.v2.DDRgene.tsv"  
data = read.table( fn, sep="\t",  head=T, fill=T,stringsAsFactors = F)

##### data annotating #####
data$multi_allele = F
#pos = gsub("(^[0-9]+:[0-9]+):.*","\\1",data$Var)
data$pos = sapply(as.character(data$Var),function(x) paste(strsplit(x,":")[[1]][1:2],collapse = ":"))
duplicated_pos = data$pos[duplicated(data$pos)]
data$multi_allele[data$pos %in% duplicated_pos] = T

data$missense = FALSE
data$missense[grep("missense",data$impact)] = TRUE

data$truncation = FALSE
for (type in vep_truncations){
  data$truncation[grep(type,data$impact)] = TRUE
}

data$inframe = FALSE
for (type in vep_inframe){
  data$inframe[grep(type,data$impact)] = TRUE
}

data$variant_type = "other"
data$variant_type[data$missense] = "missense"
data$variant_type[data$inframe] = "inframe"
data$variant_type[data$truncation] = "truncation"
data$HGVSpAbbre = gsub(".*:","",data$HGVSp)

##### some exploration into top candidates #####
data_sig = data[data$P < 0.00001 & data$ExAC_AC/data$ExAC_AN < 0.01 & !data$multi_allele,]

plot_top_counts(data_sig, n=30,x_string="gene_symbol",fill_string="variant_type")
fn = "out/DDR_277gene.p0.00001.dis.top30gene.pdf"
ggsave(file=fn, height=5,w=10, useDingbats=FALSE)

# sele_genes = c("BRCA1","BRCA2","ATM","PALB2","BRIP1","MSH6","FANCI","FANCM")
# data_sig_g = data_sig[data_sig$gene_symbol %in% sele_genes,]
# p = ggplot(data=data_sig_g, aes(x=gene_symbol,y=-log10(P), color=variant_type))
# p = p + geom_point(alpha=0.5)
# p = p + geom_text(aes(label=ifelse(variant_type=="other",as.character(impact), as.character(HGVSpAbbre))))
# p = p + labs(x = "Gene", y = "-log10(P)") + theme_bw()
# p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14),
#               axis.text.y = element_text(colour="black", size=14))
# p
# fn = "out/DDR_seleGene.var.sig.pdf"
# ggsave(file=fn, height=10,w=10,useDingbats=FALSE)

##### burden test #####
DDR_truncation = data[data$variant_type=="truncation" & !data$multi_allele,]
run_plot_burden(DDR_truncation, AF_thres=0.01)
run_plot_burden(DDR_truncation, AF_thres=0.0001)

DDR_truncation_noMulti = data[data$variant_type=="truncation" & !data$multi_allele,]
run_plot_burden(DDR_truncation_noMulti, AF_thres=0.01)
run_plot_burden(DDR_truncation_noMulti, AF_thres=0.0001)

DDR_missense = data[data$variant_type=="missense" & !data$multi_allele,]
run_plot_burden(DDR_missense, AF_thres=0.01)
run_plot_burden(DDR_missense, AF_thres=0.001)
run_plot_burden(DDR_missense, AF_thres=0.0001)

DDR_missense_noMulti = data[data$variant_type=="missense" & !data$multi_allele,]
run_plot_burden(DDR_missense_noMulti, AF_thres=0.01)
run_plot_burden(DDR_missense_noMulti, AF_thres=0.0001)

DDR_FivePrimeUTR_noMulti = data[data$impact=="5_prime_UTR_variant" & !data$multi_allele,]
run_plot_burden(DDR_FivePrimeUTR_noMulti, AF_thres=0.01)
run_plot_burden(DDR_FivePrimeUTR_noMulti, AF_thres=0.0001)

DDR_ThreePrimeUTR_noMulti = data[data$impact=="3_prime_UTR_variant" & !data$multi_allele,]
run_plot_burden(DDR_ThreePrimeUTR_noMulti, AF_thres=0.01)
run_plot_burden(DDR_ThreePrimeUTR_noMulti, AF_thres=0.0001)
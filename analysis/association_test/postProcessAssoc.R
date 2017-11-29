setwd("/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/association_test")
#fn = "test.fisher.tsv"  #"ExAC.r1.sites.vep.biallelic.combine.fisher.tsv"
# fn = "ExAC.r1.sites.vep.biallelic.combine.fisher.tsv"
fn = "assoc_results/ExAC.r1.sites.vep.biallelic.combine.fisher.NFE.152gene.tsv"
data = read.table( fn, sep="\t", fill=T, head=T)
title = "ExAC.r1.fisher.NFE.152g"

data_MAF0.01 = data[data$ExAC_AC/data$ExAC_AN < 0.01,]
# fn = "ExAC.r1.sites.vep.biallelic.combine.fisher.maf0.01.tsv"
# write.table(file=fn, data_MAF0.01, quote=F, sep="\t", row.names = F, col.names=F)

#Genomic correction
sink( sprintf( 'out/lambda.%s.MAF0.01.txt', title) )
data2=data_MAF0.01[!is.na(data_MAF0.01$P),]# & -log10(data$P)<20,]
ch = qchisq( data2$P, 1, lower.tail=F)
data2$ch=ch
theMedian = median(ch)
theLambda = median(ch)/0.456
cat( "median=", theMedian, "lambda=", theLambda, "\n")
sink()

source( 'qqman.R' )
jpeg( sprintf( 'out/qqplot.%s.MAF0.01.jpg', title), width=1200, height=1200)
qq( data2$P )
dev.off()

library(reshape2)

split_var = colsplit(string=data2$Var, pattern=":", names=c("CHR", "BP","ID","REF","ALT"))
data3 = cbind(data2,split_var)
data4 = data3[!(data3$CHR %in% c("X","Y")),]
data4$CHR = as.numeric(data4$CHR)
jpeg( sprintf( 'out/manhattan.%s.MAF0.01.jpg', title) , width=1200, height=1200)
manhattan( data4[,c("CHR","BP","P")], main=title)
dev.off()
##### plot_CharGer_summary.R #####
# Kuan-lin Huang @ WashU 2016 Jan
# plot results from CharGer

### dependencies ###
bdir = "/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/CharGer/PCGP/"
setwd(bdir)
source("../../global_aes_out.R")

### get input date and files ###

summary_AD_f = "2015_stJude_germline_nejm_S4_AD_charger.txt"
summary_AR_f = "2015_stJude_germline_nejm_S4_AR_charger.txt"

summary_AD = read.table(header=T, sep="\t",quote ="\"",stringsAsFactors = F,fill =T, file = summary_AD_f)
summary_AR = read.table(header=T, sep="\t",quote ="\"",stringsAsFactors = F,fill =T, file = summary_AR_f)

# missense_AD = summary_AD[summary_AD$Class == "missense",]
# missense_AR = summary_AR[summary_AR$Class == "missense",]
# missense = rbind(missense_AD, missense_AR)
# 
# truncation_AD = summary_AD[summary_AD$Class != "missense",]
# truncation_AR = summary_AR[summary_AR$Class != "missense",]
# truncation = rbind(truncation_AD, truncation_AR)

PCGP_all = rbind(summary_AD,summary_AR)

plot_barplot = function(matrix, x_string, fill_string=NULL, fileName="data.pdf"){
  
  fn = paste(pd, fileName,sep ="_")
  
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

##### PCGP #####
cat("##### PCGP #####")
### plotting bubble plots ###
#summary_AD$CharGer_Score = summary_AD$Positive_CharGer_Score-summary_AD$Negative_CharGer_Score 
PCGP_all = PCGP_all[PCGP_all$CharGer_Classification != "",]
PCGP_all$CharGer_Score = as.numeric(PCGP_all$CharGer_Score)

# jitter plot 
PCGP_all$Panel.Decision = as.character(PCGP_all$Panel.Decision)
PCGP_all$Panel.Decision[PCGP_all$Panel.Decision=="B"] = "Benign"
PCGP_all$Panel.Decision[PCGP_all$Panel.Decision=="PB"] = "Probably Benign"
PCGP_all$Panel.Decision[PCGP_all$Panel.Decision=="U"] = "Uncertain"
PCGP_all$Panel.Decision[PCGP_all$Panel.Decision=="PP"] = "Probably Pathogenic"
PCGP_all$Panel.Decision[PCGP_all$Panel.Decision=="P"] = "Pathogenic"
PCGP_all$Panel.Decision = factor(PCGP_all$Panel.Decision,levels=rev(c("Pathogenic","Probably Pathogenic","Uncertain","Probably Benign","Benign")))

p = ggplot(data=PCGP_all, aes(y=CharGer_Score, x=Panel.Decision))#, color = Panel.Decision))
p = p + geom_jitter(height = 0.2, alpha=0.2)
p = p + scale_colour_manual(values = rev(col_vector), guide=F) #+ scale_size_continuous(range = c(0, 20))
p = p + theme_bw()
p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_text(colour="black", size=14, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))
p = p + geom_vline(aes(xintercept=3.5), alpha=.8)
p = p + geom_hline(aes(yintercept=4.5), alpha=.8)
p
fn = "out/2015_stJude_germline_nejm_S4_PCGP_charger_concordance_jitter.pdf"
ggsave(file=fn, useDingbats=FALSE)

PCGP_all_noTP53 = PCGP_all[PCGP_all$Gene != "TP53",]
p = ggplot(data=PCGP_all_noTP53, aes(y=CharGer_Score, x=Panel.Decision))#, color = Panel.Decision))
p = p + geom_jitter(alpha=0.2)
p = p + scale_colour_manual(values = rev(col_vector), guide=F) #+ scale_size_continuous(range = c(0, 20))
p = p + theme_bw()
p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_text(colour="black", size=14, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))
p = p + geom_vline(aes(xintercept=3.5), alpha=.8)
p = p + geom_hline(aes(yintercept=5.5), alpha=.8)
p
fn = "out/2015_stJude_germline_nejm_S4_PCGP_charger_concordance_jitter_woTP53.pdf"
ggsave(file=fn, useDingbats=FALSE)

# stat table 
#table(PCGP_all$Panel.Decision, PCGP_all$CharGer_Classification)

cross_table = data.frame(table(PCGP_all$Panel.Decision, PCGP_all$CharGer_Classification))
colnames(cross_table) = c("Panel_Decision","CharGer_Classification","Freq")
outTable = dcast(cross_table, Panel_Decision ~ CharGer_Classification)
outTable = outTable[,c(1:3,6,4,5)]
#outTable[,c(1,2,5,3,4)]
# colnames(outTable) = gsub(" ","_",colnames(outTable))
# outTable$LikelyPath_or_Pathogenic = outTable$Pathogenic + outTable$Likely_Pathogenic
# 
# newrowP = outTable[outTable$Panel_Decision == "Pathogenic",2:6] + outTable[outTable$Panel_Decision == "Probably Pathogenic",2:6] #TODO
# newrowP$Panel_Decision = "ProbablyPathogenic_or_Pathogenic"
# newrowP = newrowP[,c(6,1:5)]
# newrowB = outTable[outTable$Panel_Decision == "Uncertain",2:6] + outTable[outTable$Panel_Decision == "Probably Benign",2:6] + outTable[outTable$Panel_Decision == "Benign",2:6]#TODO
# newrowB$Panel_Decision = "Uncertain_or_ProbablyBenign_or_Benign"
# newrowB = newrowB[,c(6,1:5)]

tn = "out/2015_stJude_germline_nejm_S4_PCGP_charger_concordance.tsv"
write.table(outTable, file=tn, quote=F, sep = '\t', col.names=NA)

outTable2 = data.frame()
outTable2[1,1] = sum(outTable[1:3,2:4])
outTable2[1,2] = sum(outTable[1:3,5:6])
outTable2[2,1] = sum(outTable[4:5,2:4])
outTable2[2,2] = sum(outTable[4:5,5:6])

tn = "out/2015_stJude_germline_nejm_S4_PCGP_charger_concordance_brief.tsv"
write.table(outTable2, file=tn, quote=F, sep = '\t', col.names=NA)

#outTableN = rbind(outTable,newrowP,newrowB)
# printTable = rbind(newrowP,newrowB)#outTableN[c(6,7),c(1,5,4)]
# data.frame(printTable)
cat("Total number of variants analyzed:", sum(outTable[,-1]),"\n")
cat("Sensitivity (TPR):", sum(outTable[4:5,5:6])/sum(outTable[4:5,-1]),"\n")
cat("False positive rate (FPR):", sum(outTable[1:3,5:6])/sum(outTable[1:3,-1]),"\n")

cross_table_noTP53 = data.frame(table(PCGP_all_noTP53$Panel.Decision, PCGP_all_noTP53$CharGer_Classification))
colnames(cross_table_noTP53) = c("Panel_Decision","CharGer_Classification","Freq")
outTable_noTP53 = dcast(cross_table_noTP53, Panel_Decision ~ CharGer_Classification)
outTable_noTP53 = outTable_noTP53[,c(1:3,6,4,5)]
cat("Total number of variants analyzed:", sum(outTable_noTP53[,-1]),"\n")
cat("Sensitivity (TPR):", sum(outTable_noTP53[4:5,5:6])/sum(outTable_noTP53[4:5,-1]),"\n")
cat("False positive rate (FPR):", sum(outTable_noTP53[1:3,5:6])/sum(outTable_noTP53[1:3,-1]),"\n")

# St_J_P =   sum(cross_table[cross_table$Panel.Decision=="P",]$Freq)
# both_P = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Pathogenic",]$Freq
# St_J_P_Charger_LP = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Likely Pathogenic",]$Freq
# TP_rate = both_P/St_J_P
# loose_TP_rate = (both_P+St_J_P_Charger_LP)/St_J_P
# cat("Number of St. Jude PCGP pathogenic variants: ",St_J_P, "\n")
# cat("Number of St. Jude PCGP pathogenic variants classified as Pathogenic by CharGer: ",both_P, "\n")
# cat("Number of St. Jude PCGP pathogenic variants classified as Likely Pathogenic by CharGer: ",St_J_P_Charger_LP, "\n")
# cat("True positive rate in PCGP genes: ",TP_rate, "\n")
# cat("Loose true positive rate in PCGP genes: ",loose_TP_rate, "\n")
# 
# # table for plotting
# class_table = table(PCGP_all$Panel.Decision, PCGP_all$CharGer_Score)
# class_table.m = melt(class_table)
# colnames(class_table.m) = c("PCGP_Panel_Decision","CharGer_Score","Count")
# class_table.m$ModCount = class_table.m$Count
# #class_table.m$ModCount[class_table.m$ModCount > 20] = 21
# class_table.m$PCGP_Panel_Decision = factor(class_table.m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# 
# # for grant
# if (FALSE){
#   class_table.m$PCGP_Panel_Decision = as.character(class_table.m$PCGP_Panel_Decision)
#   class_table.m$PCGP_Panel_Decision[class_table.m$PCGP_Panel_Decision=="B"] = "Benign"
#   class_table.m$PCGP_Panel_Decision[class_table.m$PCGP_Panel_Decision=="PB"] = "Probably Benign"
#   class_table.m$PCGP_Panel_Decision[class_table.m$PCGP_Panel_Decision=="U"] = "Uncertain"
#   class_table.m$PCGP_Panel_Decision[class_table.m$PCGP_Panel_Decision=="PP"] = "Probably Pathogenic"
#   class_table.m$PCGP_Panel_Decision[class_table.m$PCGP_Panel_Decision=="P"] = "Pathogenic"
#   class_table.m$PCGP_Panel_Decision = factor(class_table.m$PCGP_Panel_Decision, 
#                                              levels=rev(c("Pathogenic","Probably Pathogenic","Uncertain","Probably Benign","Benign")))
#   fn = paste(pd, "2015_stJude_germline_nejm_S4_PCGP_charger_concordance_grant.pdf",sep ="_")
#   p = ggplot(data=class_table.m, aes(y=CharGer_Score, x=PCGP_Panel_Decision, size=log10(Count), fill = PCGP_Panel_Decision))
#   p = p + geom_point(colour="black",pch=21)
#   p = p + scale_colour_manual(values = rev(brewer.pal(5,"Set1")), guide=F) + scale_size_continuous(range = c(0, 20))
#   p = p + theme_bw()
#   p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_text(colour="black", size=14, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))
#   p = p + geom_vline(aes(xintercept=3.5), alpha=.8)
#   p = p + geom_hline(aes(yintercept=5.5), alpha=.8)
#   p
#   ggsave(file=fn, useDingbats=FALSE)
# }
# 
# fn = "out/2015_stJude_germline_nejm_S4_PCGP_charger_concordance.pdf"
# p = ggplot(data=class_table.m, aes(x=CharGer_Score, y=PCGP_Panel_Decision, size=log10(Count), color = PCGP_Panel_Decision))
# #p = p + facet_grid(.~level,drop=T,scales = "free", space = "free")
# p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
# #p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
# p = p + scale_colour_manual(values = rev(brewer.pal(5,"Set1")), guide=F) + scale_size_continuous(range = c(0, 20))
# p = p + theme_bw() #+ labs( x = "Druggable Target", y = "") 
# p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_text(colour="black", size=14, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))
# p = p + geom_vline(aes(xintercept=5.5), alpha=.8)
# p = p + geom_hline(aes(yintercept=3.5), alpha=.8)
# p
# ggsave(file=fn, useDingbats=FALSE)
# 
# class_table_g = ftable(PCGP_all$Panel.Decision, PCGP_all$CharGer_Score, PCGP_all$HUGO_Symbol)
# class_table_g_m = melt(class_table_g)
# colnames(class_table_g_m) = c("PCGP_Panel_Decision","CharGer_Score","Gene","Count")
# #class_table_g_m$ModCount = class_table.m$Count
# class_table_g_m$Count[class_table_g_m$Count == 0] = NA
# class_table_g_m$PCGP_Panel_Decision = factor(class_table_g_m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# class_table_g_m$CharGer_Score = as.numeric(as.character(class_table_g_m$CharGer_Score))
# 
# fn = 'out/2015_stJude_germline_nejm_S4_PCGP_charger_gene_concordance.pdf'
# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
# p = ggplot(data=class_table_g_m)
# p = p + facet_grid(.~PCGP_Panel_Decision,drop=T,scales = "free", space = "free")
# p = p + geom_tile(aes(x=CharGer_Score, y=Gene, fill=Count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit =c(0,30))
# p = p + geom_text(aes(x=CharGer_Score, y=Gene, label = Count, stringsAsFactors=FALSE), color="black", size=3)
# p = p + geom_vline(aes(xintercept=5.5), alpha=.8,color="#fed976")
# p = p + geom_vline(aes(xintercept=8.5), alpha=.8,color="#e31a1c")
# p = p  + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=12, width=16, useDingbats=FALSE)
# 
# gene_table = table(PCGP_all[PCGP_all$Panel.Decision == "PP" | PCGP_all$Panel.Decision == "P",]$HUGO_Symbol)
# genes = names(gene_table[gene_table>=1])
# genes = genes[-1]
# class_table_g_m_s = class_table_g_m[!is.na(class_table_g_m$Gene) & class_table_g_m$Gene %in% genes & !is.na(class_table_g_m$Count),]
# 
# fn = 'out/2015_stJude_germline_nejm_S4_PCGP_charger_path_gene_concordance.pdf'
# p = ggplot(data=class_table_g_m_s)
# p = p + facet_grid(.~PCGP_Panel_Decision,drop=T)
# p = p + geom_tile(aes(x=CharGer_Score, y=Gene, fill=Count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit =c(0,30))
# p = p + geom_text(aes(x=CharGer_Score, y=Gene, label = Count, stringsAsFactors=FALSE), color="black", size=3)
# p = p + geom_vline(aes(xintercept=5.5), alpha=.8,color="#fed976")
# p = p + geom_vline(aes(xintercept=8.5), alpha=.8,color="#e31a1c")
# p = p  + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=8, width=16, useDingbats=FALSE)
# 
# # cancer classified tables
# class_table_g = ftable(PCGP_all$Panel.Decision, PCGP_all$CharGer_Score, PCGP_all$Subtype)
# class_table_g_m = melt(class_table_g)
# colnames(class_table_g_m) = c("PCGP_Panel_Decision","CharGer_Score","Subtype","Count")
# # class_table_g_m$Count[class_table_g_m$Count == 0] = NA
# class_table_g_m = class_table_g_m[class_table_g_m$Count != 0,]
# class_table_g_m$PCGP_Panel_Decision = factor(class_table_g_m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# class_table_g_m$CharGer_Score = as.numeric(as.character(class_table_g_m$CharGer_Score))
# #class_table_g_m$Count = as.numeric(as.character(class_table_g_m$Count))
# 
# fn = 'out/2015_stJude_germline_nejm_S4_PCGP_charger_cancer_concordance.pdf'
# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
# p = ggplot(data=class_table_g_m)
# p = p + facet_grid(.~PCGP_Panel_Decision)#,drop=T,scales = "free", space = "free")
# p = p + geom_tile(aes(x=CharGer_Score, y=Subtype, fill=Count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit =c(0,30))
# p = p + geom_text(aes(x=CharGer_Score, y=Subtype, label = Count, stringsAsFactors=FALSE), color="black", size=3)
# p = p + geom_vline(aes(xintercept=5.5), alpha=.8,color="#fed976")
# p = p + geom_vline(aes(xintercept=8.5), alpha=.8,color="#e31a1c")
# p = p  + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=6, width=18, useDingbats=FALSE)
# 
# ##### AD #####
# cat("##### AD #####")
# ### plotting bubble plots ###
# #summary_AD$CharGer_Score = summary_AD$Positive_CharGer_Score-summary_AD$Negative_CharGer_Score 
# summary_AD = summary_AD[summary_AD$CharGer_Classification != "",]
# summary_AD$CharGer_Score = as.numeric(summary_AD$CharGer_Score)
# 
# # stat table 
# table(summary_AD$Panel.Decision, summary_AD$CharGer_Classification)
# 
# cross_table = data.frame(table(summary_AD$Panel.Decision, summary_AD$CharGer_Classification))
# colnames(cross_table) = c("Panel.Decision","CharGer_Classification","Freq")
# St_J_P =   sum(cross_table[cross_table$Panel.Decision=="P",]$Freq)
# both_P = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Pathogenic",]$Freq
# St_J_P_Charger_LP = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Likely Pathogenic",]$Freq
# TP_rate = both_P/St_J_P
# loose_TP_rate = (both_P+St_J_P_Charger_LP)/St_J_P
# cat("Number of St. Jude AD pathogenic variants: ",St_J_P, "\n")
# cat("Number of St. Jude AD pathogenic variants classified as Pathogenic by CharGer: ",both_P, "\n")
# cat("Number of St. Jude AD pathogenic variants classified as Likely Pathogenic by CharGer: ",St_J_P_Charger_LP, "\n")
# cat("True positive rate in AD genes: ",TP_rate, "\n")
# cat("Loose true positive rate in AD genes: ",loose_TP_rate, "\n")
# 
# # table for plotting
# class_table = table(summary_AD$Panel.Decision, summary_AD$CharGer_Score)
# class_table.m = melt(class_table)
# colnames(class_table.m) = c("PCGP_Panel_Decision","CharGer_Score","Count")
# class_table.m$ModCount = class_table.m$Count
# class_table.m$ModCount[class_table.m$ModCount > 20] = 21
# class_table.m$PCGP_Panel_Decision = factor(class_table.m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# 
# fn = paste(pd, "2015_stJude_germline_nejm_S4_AD_charger_concordance.pdf",sep ="_")
# p = ggplot(data=class_table.m, aes(x=CharGer_Score, y=PCGP_Panel_Decision, size=ModCount, color = PCGP_Panel_Decision))
# #p = p + facet_grid(.~level,drop=T,scales = "free", space = "free")
# p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
# #p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
# p = p + scale_colour_manual(values = rev(brewer.pal(5,"Set1")), guide=F) + scale_size_continuous(range = c(0, 20))
# p = p + theme_bw() #+ labs( x = "Druggable Target", y = "") 
# p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_text(colour="black", size=14, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))
# p
# ggsave(file=fn, useDingbats=FALSE)
# 
# class_table_g = ftable(summary_AD$Panel.Decision, summary_AD$CharGer_Score, summary_AD$HUGO_Symbol)
# class_table_g_m = melt(class_table_g)
# colnames(class_table_g_m) = c("PCGP_Panel_Decision","CharGer_Score","Gene","Count")
# #class_table_g_m$ModCount = class_table.m$Count
# class_table_g_m$Count[class_table_g_m$Count == 0] = NA
# class_table_g_m$PCGP_Panel_Decision = factor(class_table_g_m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# class_table_g_m$CharGer_Score = as.numeric(as.character(class_table_g_m$CharGer_Score))
# 
# fn = paste(pd, '2015_stJude_germline_nejm_S4_AD_charger_gene_concordance.pdf',sep ="_")
# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
# p = ggplot(data=class_table_g_m)
# p = p + facet_grid(.~PCGP_Panel_Decision,drop=T,scales = "free", space = "free")
# p = p + geom_tile(aes(x=CharGer_Score, y=Gene, fill=Count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit =c(0,30))
# p = p + geom_text(aes(x=CharGer_Score, y=Gene, label = Count, stringsAsFactors=FALSE), color="black", size=3)
# p = p + geom_vline(aes(xintercept=5.5), alpha=.8,color="#fed976")
# p = p + geom_vline(aes(xintercept=8.5), alpha=.8,color="#e31a1c")
# p = p  + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=8, width=16, useDingbats=FALSE)
# 
# ##### AR #####
# cat("\n##### AR #####")
# summary_AR = summary_AR[summary_AR$CharGer_Classification != "",]
# summary_AR$CharGer_Score = as.numeric(summary_AR$CharGer_Score)
# 
# # stat table 
# table(summary_AR$Panel.Decision, summary_AR$CharGer_Classification)
# 
# cross_table = data.frame(table(summary_AR$Panel.Decision, summary_AR$CharGer_Classification))
# colnames(cross_table) = c("Panel.Decision","CharGer_Classification","Freq")
# cross_table
# 
# St_J_P =   sum(cross_table[cross_table$Panel.Decision=="P",]$Freq)
# both_P = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Pathogenic",]$Freq
# St_J_P_Charger_LP = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Likely Pathogenic",]$Freq
# TP_rate = both_P/St_J_P
# loose_TP_rate = (both_P+St_J_P_Charger_LP)/St_J_P
# cat("Number of St. Jude AR pathogenic variants: ",St_J_P, "\n")
# cat("Number of St. Jude AR pathogenic variants classified as Pathogenic by CharGer: ",both_P, "\n")
# cat("Number of St. Jude AR pathogenic variants classified as Likely Pathogenic by CharGer: ",St_J_P_Charger_LP, "\n")
# cat("True positive rate in AR genes: ",TP_rate, "\n")
# cat("Loose true positive rate in AR genes: ",loose_TP_rate, "\n")
# 
# ### plotting bubble plots ###
# class_table = table(summary_AR$Panel.Decision, summary_AR$CharGer_Score)
# class_table.m = melt(class_table)
# colnames(class_table.m) = c("PCGP_Panel_Decision","CharGer_Score","Count")
# class_table.m$ModCount = class_table.m$Count
# class_table.m$ModCount[class_table.m$ModCount > 20] = 21
# class_table.m$PCGP_Panel_Decision = factor(class_table.m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# 
# fn = paste(pd, "2015_stJude_germline_nejm_S4_AR_charger_concordance.pdf",sep ="_")
# p = ggplot(data=class_table.m, aes(x=CharGer_Score, y=PCGP_Panel_Decision, size=ModCount, color = PCGP_Panel_Decision))
# #p = p + facet_grid(.~level,drop=T,scales = "free", space = "free")
# p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
# #p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
# p = p + scale_colour_manual(values = rev(brewer.pal(5,"Set1")), guide=F) + scale_size_continuous(range = c(0, 20))
# p = p + theme_bw() #+ labs( x = "Druggable Target", y = "") 
# p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_text(colour="black", size=14, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))
# p
# ggsave(file=fn, useDingbats=FALSE)
# 
# PCGP_P = summary_AR[summary_AR$Panel.Decision=="P",]
# class_table_g = ftable(summary_AR$Panel.Decision, summary_AR$CharGer_Score, summary_AR$HUGO_Symbol)
# class_table_g_m = melt(class_table_g)
# colnames(class_table_g_m) = c("PCGP_Panel_Decision","CharGer_Score","Gene","Count")
# #class_table_g_m$ModCount = class_table.m$Count
# class_table_g_m$Count[class_table_g_m$Count == 0] = NA
# class_table_g_m$PCGP_Panel_Decision = factor(class_table_g_m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# class_table_g_m$CharGer_Score = as.numeric(as.character(class_table_g_m$CharGer_Score))
# 
# fn = paste(pd, '2015_stJude_germline_nejm_S4_AR_charger_gene_concordance.pdf',sep ="_")
# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
# p = ggplot(data=class_table_g_m)
# p = p + facet_grid(.~PCGP_Panel_Decision,drop=T,scales = "free", space = "free")
# p = p + geom_tile(aes(x=CharGer_Score, y=Gene, fill=Count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit =c(0,30))
# p = p + geom_text(aes(x=CharGer_Score, y=Gene, label = Count, stringsAsFactors=FALSE), color="black", size=3)
# p = p + geom_vline(aes(xintercept=5.5), alpha=.8,color="#fed976")
# p = p + geom_vline(aes(xintercept=8.5), alpha=.8,color="#e31a1c")
# p = p  + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=5, width=9, useDingbats=FALSE)
# 
# ##### missense #####
# cat("\n##### missense #####")
# missense = missense[missense$CharGer_Classification != "",]
# missense$CharGer_Score = as.numeric(missense$CharGer_Score)
# 
# # stat table 
# table(missense$Panel.Decision, missense$CharGer_Classification)
# 
# cross_table = data.frame(table(missense$Panel.Decision, missense$CharGer_Classification))
# colnames(cross_table) = c("Panel.Decision","CharGer_Classification","Freq")
# 
# St_J_P =   sum(cross_table[cross_table$Panel.Decision=="P",]$Freq)
# both_P = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Pathogenic",]$Freq
# St_J_P_Charger_LP = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Likely Pathogenic",]$Freq
# TP_rate = both_P/St_J_P
# loose_TP_rate = (both_P+St_J_P_Charger_LP)/St_J_P
# cat("Number of St. Jude missense pathogenic variants: ",St_J_P, "\n")
# cat("Number of St. Jude missense pathogenic variants classified as Pathogenic by CharGer: ",both_P, "\n")
# cat("Number of St. Jude missense pathogenic variants classified as Likely Pathogenic by CharGer: ",St_J_P_Charger_LP, "\n")
# cat("True positive rate in missense genes: ",TP_rate, "\n")
# cat("Loose true positive rate in missense genes: ",loose_TP_rate, "\n")
# 
# ### plotting bubble plots ###
# class_table = table(missense$Panel.Decision, missense$CharGer_Score)
# class_table.m = melt(class_table)
# colnames(class_table.m) = c("PCGP_Panel_Decision","CharGer_Score","Count")
# class_table.m$ModCount = class_table.m$Count
# class_table.m$ModCount[class_table.m$ModCount > 20] = 21
# class_table.m$PCGP_Panel_Decision = factor(class_table.m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# 
# fn = paste(pd, "2015_stJude_germline_nejm_S4_missense_charger_concordance.pdf",sep ="_")
# p = ggplot(data=class_table.m, aes(x=CharGer_Score, y=PCGP_Panel_Decision, size=ModCount, color = PCGP_Panel_Decision))
# #p = p + facet_grid(.~level,drop=T,scales = "free", space = "free")
# p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
# #p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
# p = p + scale_colour_manual(values = rev(brewer.pal(5,"Set1")), guide=F) + scale_size_continuous(range = c(0, 20))
# p = p + theme_bw() #+ labs( x = "Druggable Target", y = "") 
# p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_text(colour="black", size=14, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))
# p
# ggsave(file=fn, useDingbats=FALSE)
# 
# PCGP_P = missense[missense$Panel.Decision=="P",]
# class_table_g = ftable(missense$Panel.Decision, missense$CharGer_Score, missense$HUGO_Symbol)
# class_table_g_m = melt(class_table_g)
# colnames(class_table_g_m) = c("PCGP_Panel_Decision","CharGer_Score","Gene","Count")
# #class_table_g_m$ModCount = class_table.m$Count
# class_table_g_m$Count[class_table_g_m$Count == 0] = NA
# class_table_g_m$PCGP_Panel_Decision = factor(class_table_g_m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# class_table_g_m$CharGer_Score = as.numeric(as.character(class_table_g_m$CharGer_Score))
# 
# fn = paste(pd, '2015_stJude_germline_nejm_S4_missense_charger_gene_concordance.pdf',sep ="_")
# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
# p = ggplot(data=class_table_g_m)
# p = p + facet_grid(.~PCGP_Panel_Decision,drop=T,scales = "free", space = "free")
# p = p + geom_tile(aes(x=CharGer_Score, y=Gene, fill=Count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit =c(0,30))
# p = p + geom_text(aes(x=CharGer_Score, y=Gene, label = Count, stringsAsFactors=FALSE), color="black", size=3)
# p = p + geom_vline(aes(xintercept=5.5), alpha=.8,color="#fed976")
# p = p + geom_vline(aes(xintercept=8.5), alpha=.8,color="#e31a1c")
# p = p  + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=6, width=9, useDingbats=FALSE)
# 
# ##### truncation #####
# cat("\n##### truncation #####")
# truncation = truncation[truncation$CharGer_Classification != "",]
# truncation$CharGer_Score = as.numeric(truncation$CharGer_Score)
# 
# # stat table 
# table(truncation$Panel.Decision, truncation$CharGer_Classification)
# 
# cross_table = data.frame(table(truncation$Panel.Decision, truncation$CharGer_Classification))
# colnames(cross_table) = c("Panel.Decision","CharGer_Classification","Freq")
# 
# St_J_P =   sum(cross_table[cross_table$Panel.Decision=="P",]$Freq)
# both_P = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Pathogenic",]$Freq
# St_J_P_Charger_LP = cross_table[cross_table$Panel.Decision=="P" & cross_table$CharGer_Classification=="Likely Pathogenic",]$Freq
# TP_rate = both_P/St_J_P
# loose_TP_rate = (both_P+St_J_P_Charger_LP)/St_J_P
# cat("Number of St. Jude truncation pathogenic variants: ",St_J_P, "\n")
# cat("Number of St. Jude truncation pathogenic variants classified as Pathogenic by CharGer: ",both_P, "\n")
# cat("Number of St. Jude truncation pathogenic variants classified as Likely Pathogenic by CharGer: ",St_J_P_Charger_LP, "\n")
# cat("True positive rate in truncation genes: ",TP_rate, "\n")
# cat("Loose true positive rate in truncation genes: ",loose_TP_rate, "\n")
# 
# ### plotting bubble plots ###
# class_table = table(truncation$Panel.Decision, truncation$CharGer_Score)
# class_table.m = melt(class_table)
# colnames(class_table.m) = c("PCGP_Panel_Decision","CharGer_Score","Count")
# class_table.m$ModCount = class_table.m$Count
# class_table.m$ModCount[class_table.m$ModCount > 20] = 21
# class_table.m$PCGP_Panel_Decision = factor(class_table.m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# 
# fn = paste(pd, "2015_stJude_germline_nejm_S4_truncation_charger_concordance.pdf",sep ="_")
# p = ggplot(data=class_table.m, aes(x=CharGer_Score, y=PCGP_Panel_Decision, size=ModCount, color = PCGP_Panel_Decision))
# #p = p + facet_grid(.~level,drop=T,scales = "free", space = "free")
# p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
# #p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
# p = p + scale_colour_manual(values = rev(brewer.pal(5,"Set1")), guide=F) + scale_size_continuous(range = c(0, 20))
# p = p + theme_bw() #+ labs( x = "Druggable Target", y = "") 
# p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_text(colour="black", size=14, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))
# p
# ggsave(file=fn, useDingbats=FALSE)
# 
# PCGP_P = truncation[truncation$Panel.Decision=="P",]
# class_table_g = ftable(truncation$Panel.Decision, truncation$CharGer_Score, truncation$HUGO_Symbol)
# class_table_g_m = melt(class_table_g)
# colnames(class_table_g_m) = c("PCGP_Panel_Decision","CharGer_Score","Gene","Count")
# #class_table_g_m$ModCount = class_table.m$Count
# class_table_g_m$Count[class_table_g_m$Count == 0] = NA
# class_table_g_m$PCGP_Panel_Decision = factor(class_table_g_m$PCGP_Panel_Decision, levels=rev(c("P","PP","U","PB","B")))
# class_table_g_m$CharGer_Score = as.numeric(as.character(class_table_g_m$CharGer_Score))
# 
# fn = paste(pd, '2015_stJude_germline_nejm_S4_truncation_charger_gene_concordance.pdf',sep ="_")
# getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
# p = ggplot(data=class_table_g_m)
# p = p + facet_grid(.~PCGP_Panel_Decision,drop=T,scales = "free", space = "free")
# p = p + geom_tile(aes(x=CharGer_Score, y=Gene, fill=Count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit =c(0,30))
# p = p + geom_text(aes(x=CharGer_Score, y=Gene, label = Count, stringsAsFactors=FALSE), color="black", size=3)
# p = p + geom_vline(aes(xintercept=5.5), alpha=.8,color="#fed976")
# p = p + geom_vline(aes(xintercept=8.5), alpha=.8,color="#e31a1c")
# p = p  + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=5, width=9, useDingbats=FALSE)
# 
# 
# # 
# # gene_table = table(summary[summary$CharGer_Classification %in% c("Likely Pathogenic","Pathogenic"),]$HUGO_Symbol)
# # genes = names(gene_table[gene_table>=20])
# # selectedGenes = genes
# # plot_barplot(summary[summary$CharGer_Classification %in% c("Likely Pathogenic","Pathogenic") & summary$HUGO_Symbol %in% genes,], 
# #              x_string = "HUGO_Symbol", fill_string = "CharGer_Classification", 
# #              fileName="20160410_pan8000_truncation_CharGer_Classification_by_gene_bar.pdf")
# # 
# # gene_table = table(summary[summary$ClinVar_Pathogenicity %in% c("Likely Pathogenic","Pathogenic", "Pathogenic/Likely pathogenic"),]$HUGO_Symbol)
# # genes = names(gene_table[gene_table>=10])
# # plot_barplot(summary[summary$ClinVar_Pathogenicity %in% c("Likely Pathogenic","Pathogenic", "Pathogenic/Likely pathogenic") & summary$HUGO_Symbol %in% genes,], 
# #              x_string = "HUGO_Symbol", fill_string = "ClinVar_Pathogenicity", 
# #              fileName="20160410_pan8000_truncation_ClinVar_Pathogenicity_by_gene_bar.pdf")
# # 
# # fn = paste(pd, 'charger_pan8000_truncation_Positive_CharGer_Score_by_CharGer_Classification_bar.pdf',sep ="_")
# # summary$Positive_CharGer_Score = as.numeric(as.character(summary$Positive_CharGer_Score))
# # summary$score = as.numeric(as.character(summary$Positive_CharGer_Score)) - as.numeric(as.character(summary$Negative_CharGer_Score))
# # summary2 = summary[summary$CharGer_Classification %in% c("Likely Pathogenic","Pathogenic","Uncertain Significance"),]
# # p = ggplot(summary2,aes_string(x = "score", fill = "CharGer_Classification"))
# # p = p + geom_bar(position = "identity") + theme_bw() + theme_nogrid() + scale_y_log10(breaks=c(1,10,100,1000))
# # p = p + labs(x = "CharGer Score", y="counts")
# # p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
# # p
# # ggsave(file=fn, useDingbats=FALSE)
# # 
# # ##### heatmap #####
# # 
# # counts = read.table(header=F, quote = "", sep="\t", file = "Results/04_08_2016_vep.vcf_pan8000_all_truncation.tsv_count.tsv")
# # colnames(counts)= c("caller","cancer","gene","classification","count")
# # sum_c = counts[counts$gene %in% selectedGenes,]
# # 
# # ### charger
# # sum_c_char = sum_c[sum_c$caller == "CharGer",]
# # fn = paste(pd, 'charger_truncation_heatmap_summary.pdf',sep ="_")
# # plot_heatmap(sum_c_char)
# # ggsave(file=fn, height=6, width=10, useDingbats=FALSE)
# # 
# # ### clinvar
# # sum_c_clin = sum_c[sum_c$caller == "ClinVar",]
# # fn = paste(pd, 'clinvar_truncation_heatmap_summary.pdf',sep ="_")
# # plot_heatmap(sum_c_clin)
# # ggsave(file=fn, height=6, width=12, useDingbats=FALSE)
# # 
# # ##### MISSENSE #####
# # summary = read.table(header=T, quote = "", sep="\t", fill =T, file = "04_08_2016_vep.vcf_pan8000_all_missense.tsv")
# # 
# # ### plotting some bar plots ###
# # 
# # gene_table = table(summary[summary$CharGer_Classification %in% c("Likely Pathogenic","Pathogenic"),]$HUGO_Symbol)
# # genes = names(gene_table[gene_table>=3])
# # selectedGenes = genes
# # plot_barplot(summary[summary$CharGer_Classification %in% c("Likely Pathogenic","Pathogenic") & summary$HUGO_Symbol %in% genes,], 
# #              x_string = "HUGO_Symbol", fill_string = "CharGer_Classification", 
# #              fileName="20160410_pan8000_missense_CharGer_Classification_by_gene_bar.pdf")
# # 
# # gene_table = table(summary[summary$ClinVar_Pathogenicity %in% c("Likely Pathogenic","Pathogenic", "Pathogenic/Likely pathogenic"),]$HUGO_Symbol)
# # genes = names(gene_table[gene_table>=4])
# # plot_barplot(summary[summary$ClinVar_Pathogenicity %in% c("Likely Pathogenic","Pathogenic", "Pathogenic/Likely pathogenic") & summary$HUGO_Symbol %in% genes,], 
# #              x_string = "HUGO_Symbol", fill_string = "ClinVar_Pathogenicity", 
# #              fileName="20160410_pan8000_missense_ClinVar_Pathogenicity_by_gene_bar.pdf")
# # 
# # summary$Positive_CharGer_Score = as.numeric(as.character(summary$Positive_CharGer_Score))
# # plot_barplot(summary[summary$CharGer_Classification %in% c("Likely Pathogenic","Pathogenic"),],x_string = "Positive_CharGer_Score", fill_string = "CharGer_Classification",
# #              fileName="20160410_pan8000_missense_Positive_CharGer_Score_gt2_by_CharGer_Classification_bar.pdf")
# # 
# # fn = paste(pd, 'charger_pan8000_missense_Positive_CharGer_Score_by_CharGer_Classification_bar.pdf',sep ="_")
# # summary$Positive_CharGer_Score = as.numeric(as.character(summary$Positive_CharGer_Score))
# # summary$score = as.numeric(as.character(summary$Positive_CharGer_Score)) - as.numeric(as.character(summary$Negative_CharGer_Score))
# # summary2 = summary[summary$CharGer_Classification %in% c("Likely Pathogenic","Pathogenic","Uncertain Significance"),]
# # p = ggplot(summary2,aes_string(x = "score", fill = "CharGer_Classification"))
# # p = p + geom_bar(position = "identity") + theme_bw() + theme_nogrid() + scale_y_log10(breaks=c(1,10,100,1000))
# # p = p + labs(x = "CharGer Score", y="counts")
# # p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
# # p
# # ggsave(file=fn, useDingbats=FALSE)
# # 
# # ##### heatmap #####
# # counts = read.table(header=F, quote = "", sep="\t", file = "Results/04_08_2016_vep.vcf_pan8000_all_missense.tsv_count.tsv")
# # colnames(counts)= c("caller","cancer","gene","classification","count")
# # sum_c = counts[counts$gene %in% selectedGenes,]
# # 
# # ### charger
# # sum_c_char = sum_c[sum_c$caller == "CharGer" & sum_c$classification %in% c("Likely Pathogenic","Pathogenic","Uncertain Significance"),]
# # fn = paste(pd, 'charger_missense_heatmap_summary.pdf',sep ="_")
# # plot_heatmap(sum_c_char)
# # ggsave(file=fn, height=6, width=12, useDingbats=FALSE)
# # 
# # ### clinvar
# # sum_c_clin = sum_c[sum_c$caller == "ClinVar"& sum_c$classification %in% c("Likely Pathogenic","Pathogenic","Uncertain Significance"),]
# # fn = paste(pd, 'clinvar_missense_heatmap_summary.pdf',sep ="_")
# # plot_heatmap(sum_c_clin)
# # ggsave(file=fn, height=6, width=12, useDingbats=FALSE)

##### plotPathVarExpression.R #####
# Kuan-lin Huang @ WashU 201711
# plot assoc results for pathogenic variants

bdir = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/analysis/expression_effect"
setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")

p = ggplot(pathVarOT,aes(x=0,y=expressionQuantile, fill=binary_type))
p = p + facet_grid(.~Gene_Classification, scale = "free", space = "free", drop=T)
p = p + geom_dotplot(dotsize=0.4,binwidth=.02, binaxis= "y")#,stackdir ="centerwhole")
p = p + theme_bw() 
p = p + guides(fill=FALSE) + ylab("Expression Quantile")
p = p + scale_y_continuous(breaks = seq(0,1, by= 0.25))
p = p + getVarColorScale()
p = p + theme(axis.title = element_text(size=16), 
              axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),
              axis.text.y = element_text(colour="black", size=14),
              legend.position = "bottom")#element_text(colour="black", size=14))
p
fn = "out/pathVarExpression.pdf"
ggsave(file=fn, w=6, h=6,useDingbats=FALSE)

### Expression vs. CharGer ###
# first make a dataframe with frequencies

pathVarOT$exp_q.cut = cut(pathVarOT$expressionQuantile, breaks=seq(0,1,0.25))
df = as.data.frame(with(pathVarOT, table(exp_q.cut,Gene_Classification)))
# next: compute percentages per group
df = ddply(df, .(Gene_Classification), transform, p = Freq/sum(Freq))

tn = "out/pathVar_germline_onco.vs.tsg_exp.txt"
write.table(df, quote=F, sep="\t", file = tn, row.names = F)

### look into location effect of truncation on expression ###
cDNA_location_pre = gsub("-.*","",truncations$cDNA_position)
cDNA_location = as.numeric(cDNA_location_pre)

truncations$transcript_quantile = cDNA_location/truncations$transcript.length

truncations$bp_to_end = truncations$transcript.length - cDNA_location

truncations$bp_boolean = TRUE
truncations[!is.na(truncations$bp_to_end) & truncations$bp_to_end <= 50, ]$bp_boolean = FALSE

p = ggplot(truncations,aes(x=bp_to_end, y =expressionQuantile, color=CharGer_Classification))#, color=Variant_Classification))
p = p + facet_grid(Variant_Classification~.,drop=T)
p = p + geom_point(alpha=0.2) 
p = p + theme_bw() + theme_nogrid() #+ guides(color=FALSE)
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + geom_vline(xintercept=50)
p = p + scale_x_log10()
p
fn = "out/pathVar_trunc_exp_vs_transcript_bp.pdf"
ggsave(file=fn, useDingbats=FALSE)

truncations=truncations[!is.na(truncations$CharGer_Classification),]
p = ggplot(truncations,aes(x=transcript_quantile, y =expressionQuantile))#, color=Variant_Classification))
p = p + facet_grid(.~CharGer_Classification,drop=T)
#p = p + geom_point(alpha=0.1) 
p = p + theme_bw() + theme_nogrid() #+ guides(color=FALSE)
#p = p + geom_abline(intercept = 0, slope=1, alpha=0.2) #+ geom_density2d(alpha=0.5)
#p = p + stat_smooth(method=lm, fullrange=TRUE, alpha = 0.1) + ylim(0,1)
#p = p + geom_density2d(alpha=0.8)
#p = p + stat_density2d(aes(fill = ..level..),geom="tile")
p = p + stat_density2d(aes(fill = ..level..),size = 0.01, bins = 32,  geom="polygon")
p = p + scale_fill_gradientn(colours=getPalette(1000))
#p = p + scale_fill_manual(values =  brewer.pal(length(unique(truncations$CharGer_Classification)),"Set1"))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + coord_equal()
p
fn = "out/pathVar_trunc_exp_vs_transcript_quantile.pdf"
ggsave(file=fn, useDingbats=FALSE)

# look into both the truncation location and the AI relationship to expression 
truncations_exp = truncations[!is.na(truncations$expressionQuantile),]
fn = paste(pd, 'pan8000_MAF0.05_trunc_exp_AI_vs_transcript_location.pdf',sep ="_")
p = ggplot(truncations_exp,aes(x=transcript_quantile, y = tumor_by_normal_VAF, color=expressionQuantile))
p = p + facet_grid(.~CharGer_Classification,drop=T)
p = p + geom_point(alpha=0.3) + theme_bw() + theme_nogrid() 
p = p + scale_color_gradientn(colours=getPalette(1000))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + coord_equal()
p
ggsave(file=fn, width=16, useDingbats=FALSE)

p = ggplot(truncations[truncations$HUGO_Symbol == "BRCA2",],aes(x=bp_to_end, y =expressionQuantile, color=Variant_Classification))
p = p + facet_grid(.~CharGer_Classification,drop=T)
p = p + geom_point(alpha=0.8) + theme_bw() + theme_nogrid() #+ guides(color=FALSE)
p = p + scale_fill_manual(values =  brewer.pal(length(unique(truncations$CharGer_Classification)),"Set1"))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + geom_vline(xintercept=50)
p = p + scale_x_log10()
p
fn = 'out/trunc_exp_vs_transcript_location_BRCA2.pdf'
ggsave(file=fn, width=8, useDingbats=FALSE)

fn = paste(pd, 'pan8000_MAF0.05_trunc_exp_vs_transcript_location_BRCA1.pdf',sep ="_")
p = ggplot(truncations[truncations$HUGO_Symbol == "BRCA1",],aes(x=bp_to_end, y =expressionQuantile, color=Variant_Classification))
p = p + facet_grid(.~CharGer_Classification,drop=T)
p = p + geom_point(alpha=0.8) + theme_bw() + theme_nogrid() #+ guides(color=FALSE)
p = p + scale_fill_manual(values =  brewer.pal(length(unique(truncations$CharGer_Classification)),"Set1"))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14))#element_text(colour="black", size=14))
p = p + geom_vline(xintercept=50)
p = p + scale_x_log10()
p
fn = 'out/trunc_exp_vs_transcript_location_BRCA1.pdf'
ggsave(file=fn, width=8, useDingbats=FALSE)

# fit a linear model
# how does the association changes with different type of truncations
fit = lm(expressionQuantile ~ CharGer_Classification + transcript_quantile + Variant_Classification + LOH_Sig, data=truncations) 
summary(fit)
fit = lm(expressionQuantile ~ transcript_quantile, data=truncations) 
summary(fit)

fit = lm(expressionQuantile ~ CharGer_Classification + bp_to_end + Variant_Classification, data=truncations) 
summary(fit)

fit = lm(expressionQuantile ~ CharGer_Classification + bp_boolean + Variant_Classification, data=truncations) 
summary(fit)
for (type in unique(truncations$Variant_Classification)){
  truncations_type = truncations[truncations$Variant_Classification == type, ]
  fit = lm(expressionQuantile ~ CharGer_Classification + bp_boolean, data=truncations_type) 
  summary(fit)
}

fit = lm(expressionQuantile ~ CharGer_Classification + transcript_quantile + bp_boolean + Variant_Classification, data=truncations) 
summary(fit)
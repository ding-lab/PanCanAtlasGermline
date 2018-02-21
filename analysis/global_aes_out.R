### common libs and dependencies ### 
# especially for plotting #

# dependencies
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(UpSetR)
library(ggrepel)

# mis
system("mkdir out")
date=Sys.time()
date = sub(" .*","",date)
outPath=paste("out/", date, "/",sep="")

col_paletteB = colorRampPalette(brewer.pal(9,"Blues"))
col_paletteR = colorRampPalette(brewer.pal(9,"Reds"))
RdBu = brewer.pal(9, "RdBu") 
getPalette = colorRampPalette(rev(RdBu))
RdBu1024 = colorRampPalette(rev(RdBu))(1024)
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette2 = colorRampPalette(YlGnBu)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
set1 = brewer.pal(9,"Set1")
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals) ))

getPCACancerColor = function() {
  # according to google spreadsheet: https://docs.google.com/spreadsheets/d/1Nb9mMkonAhZR1_2OI9nv4ylCei0LZjAf-2vYTRQcXKw/edit#gid=1704872109
  colors = c(
    "#C1A72F",
    "#FAD2D9",
    "#ED2891",
    "#F6B667",
    "#104A7F",
    "#9EDDF9",
    "#3953A4",
    "#007EB5",
    "#B2509E",
    "#97D1A9",
    "#ED1C24",
    "#F8AFB3",
    "#EA7075",
    "#754C29",
    "#D49DC7",
    "#CACCDB",
    "#D3C3E0",
    "#A084BD",
    "#542C88",
    "#D97D25",
    "#6E7BA2",
    "#E8C51D",
    "#7E1918",
    "#DAF1FC",
    "#00A99D",
    "#BBD642",
    "#00AEEF",
    "#BE1E2D",
    "#F9ED32",
    "#CEAC8F",
    "#FBE3C7",
    "#F89420",
    "#009444")    
  color.names = c("ACC",
                  "BLCA",
                  "BRCA",
                  "CESC",
                  "CHOL",
                  "COAD",
                  "DLBC",
                  "ESCA",
                  "GBM",
                  "HNSC",
                  "KICH",
                  "KIRC",
                  "KIRP",
                  "LAML",
                  "LGG",
                  "LIHC",
                  "LUAD",
                  "LUSC",
                  "MESO",
                  "OV",
                  "PAAD",
                  "PCPG",
                  "PRAD",
                  "READ",
                  "SARC",
                  "SKCM",
                  "STAD",
                  "TGCT",
                  "THCA",
                  "THYM",
                  "UCEC",
                  "UCS",
                  "UVM")
  names(colors) = color.names
  color.scale = scale_color_manual(name="Cancer", values=colors)
  return(color.scale)
}

getPCACancerFill = function() {
  # according to google spreadsheet: https://docs.google.com/spreadsheets/d/1Nb9mMkonAhZR1_2OI9nv4ylCei0LZjAf-2vYTRQcXKw/edit#gid=1704872109
  colors = c(
    "#C1A72F",
    "#FAD2D9",
    "#ED2891",
    "#F6B667",
    "#104A7F",
    "#9EDDF9",
    "#3953A4",
    "#007EB5",
    "#B2509E",
    "#97D1A9",
    "#ED1C24",
    "#F8AFB3",
    "#EA7075",
    "#754C29",
    "#D49DC7",
    "#CACCDB",
    "#D3C3E0",
    "#A084BD",
    "#542C88",
    "#D97D25",
    "#6E7BA2",
    "#E8C51D",
    "#7E1918",
    "#DAF1FC",
    "#00A99D",
    "#BBD642",
    "#00AEEF",
    "#BE1E2D",
    "#F9ED32",
    "#CEAC8F",
    "#FBE3C7",
    "#F89420",
    "#009444")    
  color.names = c("ACC",
                  "BLCA",
                  "BRCA",
                  "CESC",
                  "CHOL",
                  "COAD",
                  "DLBC",
                  "ESCA",
                  "GBM",
                  "HNSC",
                  "KICH",
                  "KIRC",
                  "KIRP",
                  "LAML",
                  "LGG",
                  "LIHC",
                  "LUAD",
                  "LUSC",
                  "MESO",
                  "OV",
                  "PAAD",
                  "PCPG",
                  "PRAD",
                  "READ",
                  "SARC",
                  "SKCM",
                  "STAD",
                  "TGCT",
                  "THCA",
                  "THYM",
                  "UCEC",
                  "UCS",
                  "UVM")
  names(colors) = color.names
  color.scale = scale_fill_manual(name="Cancer", values=colors)
  return(color.scale)
}

getVarColorScale = function() {
  colors = c(NA, "#b2182b", "#2166ac") #positive is dark grey       
  color.names = c("Prioritized VUS","Pathogenic","Likely Pathogenic")
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="Variant Classification", values=colors)
  return(clinical.color.scale)
}

getVarColorScale2 = function() {
  colors = c("#252525", "#b2182b", "#2166ac") #positive is dark grey       
  color.names = c("Prioritized VUS","Pathogenic","Likely Pathogenic")
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="Variant Classification", values=colors)
  return(clinical.color.scale)
}

getVarFillScale = function() {
  colors = c(NA, "#b2182b", "#2166ac") #positive is dark grey       
  color.names = c("Prioritized VUS","Pathogenic","Likely Pathogenic")
  names(colors) = color.names
  clinical.color.scale = scale_fill_manual(name="Variant Classification", values=colors)
  return(clinical.color.scale)
}

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

theme0 = function(...) theme( legend.position = "none",
                               panel.background = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               #panel.margin = unit(0,"null"),
                               axis.ticks = element_blank(),
                               axis.text.x = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               axis.ticks.length = unit(0,"null"),
                               axis.ticks.margin = unit(0,"null"),
                               panel.border=element_rect(color=NA),...)

plot_top_counts = function(data, x_string, n=10, fill_string=NULL){
  top_feature = names(table(data[,x_string])[order(table(data[,x_string]),decreasing = T)][1:n])
  data_top_feature = data[data[,x_string] %in% top_feature,] 
  data_top_feature[,x_string] = factor(data_top_feature[,x_string],levels=top_feature)
  
  if (is.null(fill_string)){
    p = ggplot(data_top_feature,aes_string(x = x_string))
  } else{
    p = ggplot(data_top_feature,aes_string(x = x_string, fill = fill_string))
  }
  p = p + geom_bar() + theme_bw()
  p = p + labs(x = x_string, y="counts")
  p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
  return(p)
}

vep_truncations = c("transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","start_lost")
vep_inframe = c("inframe_insertion" , "inframe_deletion" ,"stop_lost" )

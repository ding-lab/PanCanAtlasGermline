### common libs and dependencies ### 
# especially for plotting #

# dependencies
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# mis
system("mkdir out")
date=Sys.time()
date = sub(" .*","",date)
outPath=paste("out/", date, "/",sep="")
# system(paste("mkdir ", outPath, sep=""))
# pd = paste(outPath, date, sep="")

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

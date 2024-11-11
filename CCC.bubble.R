CCC <- function(
  pfile, 
  mfile, 
  neg_log10_th= -log10(0.05), 
  means_exp_log2_th=1, 
  notused.cell=NULL, 
  used.cell=NULL, 
  neg_log10_th2=3, 
  means_exp_log2_th2=c(-4,6), 
  cell.pair=NULL, 
  gene.pair=NULL
){
 

  library(tidyverse)
  library(RColorBrewer)
  library(scales)
  library(reshape2)
  
  #-----------------------------------------------------------
  pvalues=read.table(pfile,header = T,sep = "\t",stringsAsFactors = F)
  pvalues=pvalues[,c(2,12:dim(pvalues)[2])]
  RMpairs=names(sort(table(pvalues$interacting_pair))[sort(table(pvalues$interacting_pair)) > 1])
  pvalues=pvalues[!(pvalues$interacting_pair %in% RMpairs),]
  pvalues.df1=melt(pvalues,id="interacting_pair")
  colnames(pvalues.df1)=c("geneA_geneB","cellA_cellB","pvalue")
  pvalues.df1$neg_log10=-log10(pvalues.df1$pvalue)
  pvalues.df1$geneA_geneB_cellA_cellB=paste(pvalues.df1$geneA_geneB,pvalues.df1$cellA_cellB,sep = ",")
  
  means=read.table(mfile,header = T,sep = "\t",stringsAsFactors = F)
  means=means[,c(2,12:dim(means)[2])]
  rmpairs=names(sort(table(means$interacting_pair))[sort(table(means$interacting_pair)) > 1])
  means=means[!(means$interacting_pair %in% rmpairs),]
  means.df1=melt(means,id="interacting_pair")
  colnames(means.df1)=c("geneA_geneB","cellA_cellB","means_exp")
  means.df1$geneA_geneB_cellA_cellB=paste(means.df1$geneA_geneB,means.df1$cellA_cellB,sep = ",")
  means.df1=means.df1[,c("geneA_geneB_cellA_cellB","means_exp")]
  
  raw.df=merge(pvalues.df1,means.df1,by="geneA_geneB_cellA_cellB")
  raw.df$means_exp_log2=log2(raw.df$means_exp)
  
  #-----------------------------------------------------------

  final.df=raw.df%>%filter(neg_log10 > neg_log10_th & means_exp_log2 > means_exp_log2_th)
  
  final.df$geneA=str_replace(final.df$geneA_geneB,"_.*$","") 
  final.df$geneB=str_replace(final.df$geneA_geneB,"^.*_","") 
  final.df$cellA=str_replace(final.df$cellA_cellB,"\\..*$","") 
  final.df$cellB=str_replace(final.df$cellA_cellB,"^.*\\.","")
  

  if (!is.null(notused.cell)) {
    final.df=final.df[!(final.df$cellA %in% notused.cell),]
    final.df=final.df[!(final.df$cellB %in% notused.cell),]
  }

  final.df=final.df[!(final.df$cellA==final.df$cellB),]
  
  #-----------------------------------------------------------

  final.df.gene=unique(final.df$geneA_geneB)
  final.df.cell=unique(final.df$cellA_cellB)

  if (!is.null(used.cell)){
    tmp_cell=c()
    for (i in used.cell) {
      tmp_cell=union(tmp_cell,final.df.cell[str_detect(final.df.cell,i)])
    }
    final.df.cell=tmp_cell
  }
  raw.df=raw.df[raw.df$geneA_geneB %in% final.df.gene ,]
  raw.df=raw.df[raw.df$cellA_cellB %in% final.df.cell ,]
  

  raw.df$neg_log10=ifelse(raw.df$neg_log10 > neg_log10_th2,neg_log10_th2,raw.df$neg_log10)
  raw.df$means_exp_log2=ifelse(raw.df$means_exp_log2 > means_exp_log2_th2[2],means_exp_log2_th2[2],
                               ifelse(raw.df$means_exp_log2 < means_exp_log2_th2[1],means_exp_log2_th2[1],raw.df$means_exp_log2))
  

  raw.df$cellA_cellB=as.character(raw.df$cellA_cellB)
  if (!is.null(cell.pair)) {
    tmp_pair=intersect(cell.pair,unique(raw.df$cellA_cellB))
    raw.df=raw.df[raw.df$cellA_cellB %in% tmp_pair,]
    raw.df$cellA_cellB=factor(raw.df$cellA_cellB,levels = tmp_pair)
  } else {
    tmp_pair=sort(unique(raw.df$cellA_cellB))
    raw.df$cellA_cellB=factor(raw.df$cellA_cellB,levels = tmp_pair)
  }
  #geneA-geneB的排列顺序
  raw.df$geneA_geneB=as.character(raw.df$geneA_geneB)
  if (!is.null(gene.pair)) {
    tmp_pair=intersect(gene.pair,unique(raw.df$geneA_geneB))
    raw.df=raw.df[raw.df$geneA_geneB %in% tmp_pair,]
    raw.df$geneA_geneB=factor(raw.df$geneA_geneB,levels = tmp_pair)
  } else {
    tmp_pair=sort(unique(raw.df$geneA_geneB))
    raw.df$geneA_geneB=factor(raw.df$geneA_geneB,levels = tmp_pair)
  }
  
  #-----------------------------------------------------------
  p=raw.df%>%ggplot(aes(geneA_geneB,cellA_cellB))+geom_point(aes(size=neg_log10,color=means_exp_log2))+
    scale_x_discrete("")+scale_y_discrete("")+
    scale_color_gradientn(colours = rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61","#FFFFB3","#ABD9E9", "#4575B4","#313695")))+
    theme_bw()+
    theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 0.5, size=8, angle = 90))+
    coord_flip()
  return(p)
}

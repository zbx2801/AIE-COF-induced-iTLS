sce2@meta.data$cellType1[sce2$RNA_snn_res.1.6%in%c("0","2","7")]<-"Cd8+Tcell"
sce2@meta.data$cellType1[sce2$RNA_snn_res.1.6%in%c("6","5")]<-"Treg"
sce2@meta.data$cellType1[sce2$RNA_snn_res.1.6%in%c("1","4","3","8")]<-"Cd4+Tcell"
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)  
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)



library(SingleR)


refdata <- ref_Human_all
load("sce2.Rdata")
testdata <- GetAssayData(scRNA_harmony, slot="data")
###把scRNA数据中的seurat_clusters提取出来，注意这里是因子类型的
clusters <- scRNA_harmony@meta.data$seurat_clusters
###开始用singler分析
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
scRNA_harmony@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
DimPlot(scRNA_harmony,group.by = "orig.ident")

DimPlot(scRNA_harmony,group.by = "celltype")


setwd("E:/super.lesson/lesson9/")
sc.combined=scRNA_harmony

table(sc.combined@meta.data$celltype)
table(sc.combined@meta.data$orig.ident)
type=c("B cells","Cd4+Tcell","Cd8+Tcell",
       "Epithelial cells","Myeloid cells", "NK cells","Plasma cells"
         ,"Stromal cells","Treg")


dir.create("deg/")
setwd("deg/")




library(future)
plan("multiprocess", workers =5)
options(future.globals.maxSize = 2000 * 1024^2)

r.deg=data.frame()
table(sce2@meta.data$orig.ident)

for (i in 1:9) {
  Idents(sce2)="cellType1"
  deg=FindMarkers(sce2,ident.1 = "con1",ident.2 = "ex",
                  group.by = "orig.ident",subset.ident =type[i]   )
  
  write.csv(deg,file = paste0( type[i],'deg.csv') )
  deg$gene=rownames(deg)
  deg$celltype=type[i]
  deg$unm=i-1
  r.deg=rbind(deg,r.deg)
  
}


table(r.deg$unm) 
#############################################################################
r.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.2)
r.deg$threshold <- as.factor(ifelse(r.deg$avg_log2FC > 0 , 'Up', 'Down'))
dim(r.deg)

r.deg$adj_p_signi <- as.factor(ifelse(r.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
r.deg$thr_signi <- paste0(r.deg$threshold, "_", r.deg$adj_p_signi)
r.deg$unm %<>% as.vector(.) %>% as.numeric(.)



top_up_label <- r.deg %>% 
  subset(., threshold%in%"Up") %>% 
  group_by(unm) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  as.data.frame()

top_down_label <- r.deg %>% 
  subset(., threshold %in% "Down") %>% 
  group_by(unm) %>% 
  top_n(n = -5, wt = avg_log2FC) %>% 
  as.data.frame()

top_label <- rbind(top_up_label,top_down_label)
top_label$thr_signi %<>% 
  factor(., levels = c("Up_Highly","Down_Highly","Up_Lowly","Down_Lowly"))


#write.csv(top_label, "easy_input_label.csv", quote = F)



colnames(r.deg)

background_position <- r.deg %>%
  dplyr::group_by(unm) %>%
  dplyr::summarise(Min = min(avg_log2FC) - 0.2, Max = max(avg_log2FC) + 0.2) %>%
  as.data.frame()
## `summarise()` ungrouping output (override with `.groups` argument)
background_position$unm %<>% as.vector(.) %>% as.numeric(.)
background_position$start <- background_position$unm - 0.4
background_position$end <- background_position$unm + 0.4


cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$unm - 0.5
cluster_bar_position$end <- cluster_bar_position$unm + 0.5
cluster_bar_position$unm %<>% 
  factor(., levels = c(0:max(as.vector(.))))


cols_thr_signi <- c("Up_Highly" = "#d7301f",
                    "Down_Highly" = "#225ea8",
                    "Up_Lowly" = "black",
                    "Down_Lowly" = "black")
cols_cluster <- c("0" = "#35978f",
                  "1" = "#8dd3c7",
                  "2" = "#ffffb3",
                  "3" = "#bebada",
                  "4" = "#fb8072",
                  "5" = "#80b1d3",
                  "6" = "#fdb462",
                  "7" = "#b3de69","8" = "#b3de89")

p= ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,
                                            ymax = Max),
            fill = "#525252", alpha = 0.1) +
  geom_jitter(data = r.deg, aes(x =unm, y = avg_log2FC, colour = thr_signi),
              size = 1,position = position_jitter(seed = 1)) +
  scale_color_manual(values = cols_thr_signi) +
  scale_x_continuous(limits = c(-0.5, max(r.deg$unm) + 0.5),
                     breaks = seq(0, max(r.deg$unm), 1),
                     label = seq(0, max(r.deg$unm),1)) 
  
  geom_text_repel(data = top_label, aes(x =unm, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                  box.padding = unit(0, "lines")) +
  
  geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.4,
                                             ymax = 0.4, fill = unm), color = "black", alpha = 1, show.legend = F) +
  scale_fill_manual(values = cols_cluster) +
  labs(x = "Cluster", y = "average log2FC") +
  theme_bw()

plot1 <- p + theme(panel.grid.minor = element_blank(), 
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(colour = 'black', size = 14),
                   axis.text.x = element_text(colour = 'black', size = 14, vjust = 73), 
                   panel.border = element_blank(), 
                   axis.ticks.x = element_blank(), 
                   axis.line.y = element_line(colour = "black")) 

plot1
ggsave(filename = "deg_pointplot.pdf", plot = plot1, width = 9, height = 6)

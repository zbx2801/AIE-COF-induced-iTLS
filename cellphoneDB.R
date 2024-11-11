setwd("C:/Users/jk/Desktop/COF/单细胞/细胞通讯/con1cellType1")
setwd("C:/Users/jk/Desktop/COF/单细胞/细胞通讯/excellType1")
setwd("C:/Users/jk/Desktop/COF/单细胞/细胞通讯/totalcellType1")
######################################################################################################
all_pval = read.table("out/pvalues.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
all_means = read.table("out/means.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)

# 提取要画的受体配体对和互作的细胞类型
intr_pairs = all_pval$interacting_pair
all_pval = all_pval[,-c(1:11)]
all_means = all_means[,-c(1:11)]

# 加载要画的受体配体对
selected_rows <- read.table("out/means.txt",header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)

selected_rows <- sample(c("PDCD1_CD274","CTLA4_CD86","CTLA4_CD80"
,"CSF1R_CSF1","CSF1_SIRPA"))



selected_rows <- sample(c("CSF1R_CSF1","SIRPA_CD47","CD80_CD274"))
"CSF1R_CSF1" "SIRPA_CD47" CD80_CD274

# 加载要画的互作的细胞类型
selected_columns <- sample(c("Treg|B cells","Treg|Cd4+Tcell","Treg|Cd8+Tcell","Treg|Epithelial cells",
"Treg|Myeloid cells","Treg|NK cells","Treg|Plasma cells","Treg|Stromal cells","Treg|Treg"))

selected_columns <- sample(c("B cells|Treg","Cd4+Tcell|Treg","Cd8+Tcell|Treg","Epithelial cells|Treg",
"Myeloid cells|Treg","NK cells|Treg","Plasma cells|Treg","Stromal cells|Treg","Treg|Treg"))


selected_columns <- sample(c("Epithelial cells|B cells","Epithelial cells|Cd4+Tcell","Epithelial cells|Cd8+Tcell","Epithelial cells|Epithelial cells",
"Epithelial cells|Myeloid cells","Epithelial cells|NK cells","Epithelial cells|Plasma cells","Epithelial cells|Stromal cells","Epithelial cells|Treg"))
selected_columns <- sample(c("Myeloid cells|B cells","Myeloid cells|Cd4+Tcell","Myeloid cells|Cd8+Tcell","Myeloid cells|Epithelial cells",
"Myeloid cells|Myeloid cells","Myeloid cells|NK cells","Myeloid cells|Plasma cells","Myeloid cells|Stromal cells","Myeloid cells|Treg"))

sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

df_names = expand.grid(selected_rows, selected_columns)
pval = unlist(sel_pval)
pval[pval==0] = 0.0009
plot.data = cbind(df_names,pval)
pr = unlist(as.data.frame(sel_means))
pr[pr==0] = 1
plot.data = cbind(plot.data,log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
head(plot.data)

my_palette <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(20))

pdf("receptorLigand.pdf", width = 15, height = 10)
ggplot(plot.data, aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', 
                        colors=my_palette) + # 用自定义颜色画点
  theme_bw() +
  theme(panel.grid.minor = element_blank(), #不画网格
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")) #边框
dev.off()

#####Figure 2#####
##analysis in SCP---
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(SCP)
library(BiocParallel)
register(MulticoreParam(workers = 8, progressbar = TRUE))

#load ST data
load('/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/data/HCC_region.rdata')

##Fig2A----
table(ST@meta.data$Region)
#exclude TR
subST<-subset(ST,Region %in% c('T','PT','TC','DPT'))
#add group
subST@meta.data$Group[subST@meta.data$orig.ident %in% c('P1','P2')] <- 'MVI-'
subST@meta.data$Group[subST@meta.data$orig.ident %in% c('P3','P4')] <- 'MVI+'
subST@meta.data$label<-paste(subST@meta.data$Region,subST@meta.data$Group,sep = '_')
subST@meta.data$label[subST@meta.data$orig.ident %in% c('N1','N3')] <- 'DPT'
table(subST$label)

Idents(subST)<-subST@meta.data$label

SpatialPlot(subST,stroke = 0,  pt.size.factor = 1.4,image.alpha = 0)


subST<-PrepSCTFindMarkers(subST)
#min.pct:0.6
subDEG <-FindAllMarkers(subST, only.pos = T, min.pct = 0.6, logfc.threshold = 0.25)

top100_genes = subDEG %>% group_by(cluster) %>% top_n(100, avg_log2FC)
top100_genes$cluster = factor(top100_genes$cluster, levels=c("T_MVI+", "T_MVI-","TC_MVI+","TC_MVI-", "PT_MVI+", "PT_MVI-", "DPT"))

ht <- FeatureHeatmap(
  srt = subST2, group.by = "label", split.by = "orig.ident",
  feature_split = top200_genes$cluster,
  features = top200_genes$gene,anno_terms = FALSE,
  group_palette = "Set1",
  cell_split_palcolor = c("#FB8072","#FB9A99","#55B1B1", "#8DD3C7","#A2CD5A","#CCEBC5"),
  feature_split_palette = "Set1",
)


##Fig2B----
library(clusterProfiler)
#SYMBOL to ENTREZID
gid <- bitr(unique(subDEG$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(subDEG, gid, by=c('gene' = 'SYMBOL'))

#enrichKEGG
KEGG = compareCluster(ENTREZID ~ cluster, data = markers, fun='enrichKEGG',pvalueCutoff=1)

#subset T_MVI+KEGG
T_MVI_pos<-subset(KEGG@compareClusterResult,cluster == 'T_MVI+')

df1 <- subset(T_MVI_pos, ID %in% c('hsa00190','hsa03030','hsa04110','hsa03410','hsa03430'))

df1$LogP<-log10(df1$pvalue)
df1$LogP <- -df1$LogP

df1$labelx=rep(0,nrow(df1))
df1$labely=seq(nrow(df1),1)

ggplot(data = df1, 
       aes(LogP, reorder(Description, LogP))) +
  geom_bar(stat = "identity",
           alpha = 0.5,
           fill = "#FE8D3C",
           width = 0.8) + 
  geom_text(aes(x = labelx,
                y = labely,
                label = Description),
            size = 5,  # 修改此处的 size 值以增大字体
            hjust = 0) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour = 'black', linewidth = 1),
    axis.text.x = element_text(colour = 'black', size = 12),
    axis.ticks.x = element_line(colour = 'black', linewidth = 1),
    axis.title.x = element_text(colour = 'black', size = 14),
    text = element_text(size = 16)  # 加大全局字体大小
  ) +
  xlab("-log10(pvalue)") +
  ggtitle("T_MVI+ KEGG Enrichment") +
  scale_x_continuous(expand = c(0, 0))


##Fig2C----
#subset T_MVI-KEGG
T_MVI_neg<-subset(KEGG@compareClusterResult,cluster == 'T_MVI-')
df2 <- subset(T_MVI_neg, ID %in% c('hsa04141','hsa01212','hsa04142','hsa04979','hsa04152'))

df2$LogP<-log10(df2$pvalue)
df2$LogP <- -df2$LogP

df2$labelx=rep(0,nrow(df2))
df2$labely=seq(nrow(df2),1)

ggplot(data = df2, 
       aes(LogP, reorder(Description, LogP))) +
  geom_bar(stat = "identity",
           alpha = 0.5,
           fill = "#FE8D3C",
           width = 0.8) + 
  geom_text(aes(x = labelx,
                y = labely,
                label = Description),
            size = 5,  # 修改此处的 size 值以增大字体
            hjust = 0) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour = 'black', linewidth = 1),
    axis.text.x = element_text(colour = 'black', size = 12),
    axis.ticks.x = element_line(colour = 'black', linewidth = 1),
    axis.title.x = element_text(colour = 'black', size = 14),
    text = element_text(size = 16)  # 加大全局字体大小
  ) +
  xlab("-log10(pvalue)") +
  ggtitle("T_MVI- KEGG Enrichment") +
  scale_x_continuous(expand = c(0, 0.2))




##Fig2D----
#TCGA-LIHC MVI+DEG
deg<-read.csv('~/DATA/luo/tongji/TCGA/1011/TCGA_LIHC_MVI_DEG.csv')
#up-regulated gene
UP_MVI<-subset(deg,Group %in% 'UP')

#AddModuleScore
ST<-AddModuleScore(ST,features =list(UP_MVI$X) ,name = 'MVI_score')

#Score in Public data
HCC1 <- readRDS("~/DATA/luo/pancancer/HCC/1L/HCC1.rds")
DefaultAssay(HCC1)<-'SCT'
HCC2 <- readRDS("~/DATA/luo/pancancer/HCC/2L/HCC2.rds")
DefaultAssay(HCC2)<-'SCT'

HCC1<-AddModuleScore(HCC1,features =list(UP_MVI$X) ,name = 'MVI_score')
HCC2<-AddModuleScore(HCC2,features =list(UP_MVI$X) ,name = 'MVI_score')


#MVI+ Score visualization
plot_gene = function (cluster){
  p1 = SpatialPlot(ST, features  = cluster,stroke = 0.15, pt.size.factor =1.4,images = 'C2',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p2 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 2.3,images = 'C4',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p3 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'C7',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p4 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'C8',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p5 = SpatialPlot(HCC1, features = cluster,stroke = 0.15,pt.size.factor = 2,image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p6 = SpatialPlot(HCC2, features = cluster,stroke = 0.15,pt.size.factor = 1.6,image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  ggarrange(p1,p2,p3,p4,p5,p6,ncol  = 2,nrow = 3)
}
plot_gene(cluster='MVI_score1')






##Fig2E----
##Scissor analysis
library(Scissor)

##HCC scRNA-seq data(HRA001748)
library(tidydr)
load('~/DATA/luo/tongji/HCC_nature/Hep_score.rdata')
DimPlot(Hep,group.by = 'seurat_clusters',pt.size = 0.05)

##TCGA-LIHC 
LIHC<-read.csv('/home/zhaojingwei/DATA/luo/tongji/TCGA/LIHC_TCGA.csv',header = T,row.names = 1)
table(LIHC$MVI)

#exclude Normal
LIHCl1<-subset(LIHC,MVI %in% c('MVI','non_MVI'))


#MVI Phenotype information
phenotype<-LIHCl1[,1, drop = FALSE]
table(phenotype$MVI)
phenotype$Group[phenotype$MVI %in% 'MVI']<-'1'
phenotype$Group[phenotype$MVI %in% 'non_MVI']<-'0'
table(phenotype$Group)
phenotype$MVI=NULL
table(phenotype)

##transposition
LIHCl1$MVI<-NULL
LIHCl1<-t(LIHCl1)

##random sampling 20%
unique(Hep$orig.ident)
Idents(Hep)<-Hep@meta.data$orig.ident
allCells=names(Idents(Hep))
allType = levels(Idents(Hep))


choose_Cells = unlist(lapply(allType, function(x) {
  cgCells = allCells[Idents(Hep) == x]
  
  numToSelect = ceiling(0.2 * length(cgCells))
  cg = sample(cgCells, numToSelect)
  cg
}))

#subset 
Hep_sub = Hep[, allCells %in% choose_Cells]
Hep_sub
DimPlot(Hep_sub,group.by = 'Label2')
Hep_sub<-FindNeighbors(Hep_sub)

##phenotype conver to numeric
class(phenotype)
phenotype <- as.numeric(unlist(phenotype))
tag <- c('MVI-', 'MVI+')

#Run Scissor
infos4 <- Scissor(LIHCl1, Hep, phenotype, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "~/DATA/luo/zhongnan/HCC_nature/Scissor_HCC_MVI_mutation.RData")

##visualization
Scissor_select <- rep(0, ncol(Hep_sub))
names(Scissor_select) <- colnames(Hep_sub)

#Scissor result add to metadata
Hep_sub <- AddMetaData(Hep_sub, metadata = Scissor_select, col.name = "Scissor")

#Scissor UMAP
DimPlot(Hep_sub, reduction = 'umap', group.by = 'Scissor', cols = c('grey','royalblue','indianred1'), pt.size = 1.2, order = c('Scissor_Pos','Scissor_Neg'))+theme_dr(xlength=0.22,ylength=0.22, arrow=grid::arrow(length=unit(0.15,"inches"),type="closed"))+theme(panel.grid=element_blank()) + theme(legend.text = element_text(size = 14))+ggtitle('Scissor')+ theme(plot.title = element_text(size = 16))


##Fig2F----
cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999")

DimPlot(Hep,group.by = 'seurat_clusters',cols = cols,pt.size = 0.05)+theme_dr(xlength=0.22,ylength=0.22, arrow=grid::arrow(length=unit(0.15,"inches"),type="closed"))+theme(panel.grid=element_blank())


##Fig2G----
#Scissor_Pos DEG
Scissor_deg2<-FindMarkers(Hep_sub2,ident.1 = c('6','7'),ident.2 = c('0','3','4','15'),logfc.threshold = 0.25,min.pct = 0.25)

Scissor_deg2$cluster<-'ns'
Scissor_deg2$cluster[Scissor_deg2$avg_log2FC >= 0.25 & Scissor_deg2$p_val_adj < 0.05]<-'Scissor_Pos'
Scissor_deg2$cluster[Scissor_deg2$avg_log2FC <= -0.25 & Scissor_deg2$p_val_adj < 0.05]<-'Scissor_Neg'
table(subDEG$cluster)

#subset
subDEG<-subset(Scissor_deg2,cluster %in% c('Scissor_Neg','Scissor_Pos'))
subDEG$gene=rownames(subDEG)

#enrichKEGG
library(clusterProfiler)
#SYMBOLto ENTREZID
gid <- bitr(unique(subDEG$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(subDEG, gid, by=c('gene' = 'SYMBOL'))

KEGG = compareCluster(ENTREZID ~ cluster, data = markers, fun='enrichKEGG')

#subset Scissor_Pos KEGG
PosKEGG<-subset(KEGG@compareClusterResult,cluster == 'Scissor_Pos')

df1 <- subset(PosKEGG, Description %in% c('Cell cycle','Oxidative phosphorylation','Proteasome','DNA replication','Spliceosome','Nucleotide excision repair','Mismatch repair','Base excision repair','Nucleotide metabolism','p53 signaling pathway'))

df1$LogP<-log10(df1$pvalue)
df1$LogP <- -df1$LogP

df1$labelx=rep(0,nrow(df1))
df1$labely=seq(nrow(df1),1)

ggplot(data = df1, 
       aes(LogP, reorder(Description, LogP))) +
  geom_bar(stat = "identity",
           alpha = 0.5,
           fill = "#FE8D3C",
           width = 0.8) + 
  geom_text(aes(x = labelx,
                y = labely,
                label = Description),
            size = 5,  # 修改此处的 size 值以增大字体
            hjust = 0) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour = 'black', linewidth = 1),
    axis.text.x = element_text(colour = 'black', size = 12),
    axis.ticks.x = element_line(colour = 'black', linewidth = 1),
    axis.title.x = element_text(colour = 'black', size = 14),
    text = element_text(size = 16)  # 加大全局字体大小
  ) +
  xlab("-log10(pvalue)") +
  ggtitle("Scissor Pos KEGG Enrichment") +
  scale_x_continuous(expand = c(0, 0))



##Fig2H----
#Cluster6,7 DEG
load( '~/DATA/luo/tongji/HCC_nature/Scissor/DEG/Scissor_deg.rdata')
table(Scissor_deg2$cluster)

#MVI scRNA-seq UP
scDEG<-read.csv( '~/DATA/luo/tongji/HCC2/GSE242889/scDEG.csv')

#ST MVI UP
stDEG<- readRDS("~/DATA/luo/tongji/ST/DEG/ST_MVI+_Tumor_DEG.rds")
stDEG$cluster<-'ns'
stDEG$cluster[stDEG$avg_log2FC >= 0.5 & stDEG$p_val_adj < 0.05 & stDEG$pct.1 >= 0.5]<-'UP'
stDEG$cluster[stDEG$avg_log2FC <= -0.5 & stDEG$p_val_adj < 0.05 & stDEG$pct.1 >= 0.5]<-'Down'
table(stDEG$cluster)


##Fig2I----
#load TCGA-LIHC MVI DEG
deg<-read.csv('~/DATA/luo/tongji/TCGA/1011/TCGA_LIHC_MVI_DEG.csv')
#subset up-regulated 
UP_MVI<-subset(deg,Group %in% 'UP')


##Fig2J----
load( '~/DATA/luo/tongji/HCC_nature/Scissor/DEG/Scissor_deg.rdata')
table(Scissor_deg2$cluster)
Pos<-subset(Scissor_deg2, cluster %in% 'Scissor_Pos')

#HCC AddModuleScore
ST<-AddModuleScore(ST,features = list(rownames(Pos)),name = 'Scissor_Pos')
#Score in Public data
HCC1<-AddModuleScore(HCC1,features = list(rownames(Pos)),name = 'Scissor_Pos')
HCC2<-AddModuleScore(HCC2,features = list(rownames(Pos)),name = 'Scissor_Pos')

#visualization
plot_gene = function (cluster){
  p1 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'C7',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p2 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'C8',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p3 = SpatialPlot(HCC1, features = cluster,stroke = 0.15,pt.size.factor = 2,image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p4 = SpatialPlot(HCC2, features = cluster,stroke = 0.15,pt.size.factor = 1.6,image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  ggarrange(p1,p2,p3,p4,ncol  = 2,nrow = 2)
}
plot_gene(cluster='Scissor_Pos1')

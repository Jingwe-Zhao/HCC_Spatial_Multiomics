######Figure6#####
library(RColorBrewer)
library(Seurat)
library(Matrix)
library(ggpubr)

#load ST data
load('/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/data/HCC_region.rdata')

SpatialPlot(ST,group.by = 'label',images = c('C2','C4','C7','C8'))

##Fig6A----
ST$label2<-'Other'
ST$label2[ST$clusters %in% '3']<-'C3'
ST$label2[ST$clusters %in% '5']<-'C5'
ST$label2[ST$label %in% 'Boundary']<-'Bdy'

col = c('C3'="#fdc086",'C5'="#f0027f", 'Bdy'="#1965B0", 'Other'="#FFFAFA")

plot_gene = function (cluster){
  p1 = SpatialPlot(ST, group.by = cluster,stroke = 0.15, pt.size.factor =1.4,images = 'C2',image.alpha = 0,cols = col)+ggtitle('P1')
  p2 = SpatialPlot(ST, group.by = cluster,stroke = 0.15,pt.size.factor = 2.3,images = 'C4',image.alpha = 0,cols = col)+ggtitle('P2')
  p3 = SpatialPlot(ST, group.by = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'C7',image.alpha = 0,cols = col)+ggtitle('P3')
  p4 = SpatialPlot(ST, group.by = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'C8',image.alpha = 0,cols = col)+ggtitle('P4')
  ggarrange(p1,p2,p3,p4,ncol  = 2,nrow =2,common.legend = T,legend = 'right')
}
plot_gene(cluster='label2')

#label meta.data
meta<-ST@meta.data[,'label2',drop=F]
#rownames
meta$CB<-rownames(meta)
meta$CB <- gsub("-", "_", meta$CB)
meta$CB<-paste0(meta$CB,'-1')
rownames(meta)<-meta$CB




##Fig6B----
#RCTD result
load("~/DATA/luo/zhongnan/ST/new/data/ST_region.rdata")

table(ST$orig.ident)
subST<-subset(ST, orig.ident %in% c('C2','C4','C7','C8'))

#AddMetaData
subST<-AddMetaData(subST,metadata =meta )

#subset
subST2<-subset(subST, label2 %in% c('Bdy','C3','C5'))

feature = c('Fibroblast','Endothelial','B','TNK','Myeloid','Hepatocyte','Tumor')

colors <- c("#A6761D",  "#1965B0",  "#FDB462", "#B2DF8A", "#55B1B1", "#E78AC3","#DC050C" )

#order
Idents(subST2) <-subST2$label2
My_levels <- c('C3','C5','Bdy')
Idents(subST2) <- factor(Idents(subST2), levels= My_levels)


VlnPlot(subST2,features=feature,stack=TRUE,flip=TRUE,cols=colors)+
  theme(legend.position="none")



##Fig6C----
load('/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/data/HCC_region.rdata')
SpatialPlot(ST,group.by = 'Region')

#subset TC
table(ST$Region)
TC<-subset(ST, Region %in%  c('TC') & orig.ident %in% c('C2','C4','C7','C8'))
SpatialPlot(TC,group.by = 'Region')


TC@meta.data$Group<-'MVI+'
TC@meta.data$Group[TC@meta.data$orig.ident %in% c('C2','C4')]<-'MVI-'
#switch Idents
Idents(TC)<-TC@meta.data$Group
SpatialPlot(TC)

#SCTransform
TC<-SCTransform(TC, assay = "Spatial", verbose = T)
#DEG
DEG<-FindAllMarkers(TC,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




library(clusterProfiler)
#SYMBOL to ENTREZID
gid <- bitr(unique(DEG$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(DEG, gid, by=c('gene' = 'SYMBOL'))

#enrichKEGG
KEGG = compareCluster(ENTREZID ~ cluster, data = markers, fun='enrichKEGG')

list=c('B cell receptor signaling pathway','NF−kappa B signaling pathway','Primary immunodeficiency','Endocytosis','Yersinia infection','Cell cycle','Chemokine signaling pathway','TNF signaling pathway','Complement and coagulation cascades','Alcoholic liver disease','PPAR signaling pathway','Fatty acid degradation','Primary bile acid biosynthesis','Sulfur metabolism')


#dotplot
dotplot(KEGG, label_format=40,showCategory=list) + theme(axis.text.x = element_text(angle=45, hjust=1))  +
  scale_color_distiller(palette = "Blues", direction = -1)

#save
write.csv(KEGG@compareClusterResult,file = '~/DATA/luo/zhongnan/ST/TC/TC_DEG_KEGG.csv')




##Fig6D----
#CAF signatures
CAF<-read.csv('/home/zhaojingwei/DATA/luo/zhongnan/ST/CAF/CAF_Marker.csv')
table(CAF$cluster)

iCAF<-subset(CAF ,cluster %in% 'iCAF')

TC<-AddModuleScore(TC,features = list(iCAF$gene),name = 'iCAF')

cols=c("MVI+" = "#DC050C", "MVI-" = "#386cb0")
my_comparisons<-list(c('MVI+','MVI-'))

##ST score--left
VlnPlot(TC, features = 'iCAF1', group.by = 'Group', pt.size = 0, cols = c('#386cb0','#DC050C')) +
  ylab('Score in ST') + xlab('') + ggtitle('iCAF signatures')+
  geom_boxplot(width = 0.2, col = "black", fill = NA, outlier.shape = NA) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 24), 
        text = element_text(size = 18)) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif', size = 14, label.y = 1.25) +
  ylim(-0.78, 1.5) + NoLegend()



##TCGA score--middle
LIHC<-read.csv('/home/zhaojingwei/DATA/luo/zhongnan/TCGA/LIHC_TCGA.csv',row.names = 1)
LIHCl1<-subset(LIHC,MVI %in% c('MVI','non_MVI'))

#MVI phenotype
meta<-LIHC[,'MVI', drop = FALSE]
meta$CB<-row.names(meta)
meta$Group<-meta$MVI
meta$Group[meta$Group %in% 'MVI']<-'MVI+'
meta$Group[meta$Group %in% 'non_MVI']<-'MVI-'

LIHCl1$MVI<-NULL
#transposed
LIHCl1<-t(LIHCl1)

##GSVA
ES <- gsva(LIHCl1, list(iCAF$gene))
ES<-as.data.frame(t(ES))
ES$CB<-row.names(ES)
colnames(ES)<-'iCAF'

#merge
score<-merge(ES,meta,by='CB')

my_comparisons <- list( c("MVI-", "MVI+"))
ggplot(score, aes(x = Group, y = iCAF, fill = Group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, color = "black", alpha = 0.6) +
  labs(title = "iCAF signatures", x = "", y = "Score in TCGA LIHC") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 24, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 24),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("MVI+" = "#DC050C", "MVI-" = "#386cb0")) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif', size = 14,label.y = 0.65)+
  ylim(-0.8, 0.9)+guides(fill = FALSE) 




##GSE10141--right
exp<-read.csv('~/DATA/luo/zhongnan/TCGA/CAF/GSE10141/exprSet.csv',row.names = 1)
exp2<-as.matrix(exp)
#read phenotype
meta<-read.csv('~/DATA/luo/zhongnan/TCGA/CAF/GSE10141/group.csv',row.names = 1)
rownames(meta)=meta$V1
meta$V1=NULL
meta$Group=meta$V2

meta$Group[meta$Group %in% 'positive']<-'MVI+'
meta$Group[meta$Group %in% 'negative']<-'MVI-'
table(meta$Group)
meta$V2=NULL

##GSVA 
ES <- gsva(exp2, list(iCAF$gene))
ES<-as.data.frame(t(ES))
colnames(ES)<-'iCAF'

score<-cbind(ES,meta)
table(score$Group)


my_comparisons <- list( c("MVI+", "MVI-"))

ggplot(score, aes(x = Group, y = iCAF, fill = Group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, color = "black", alpha = 0.6) +  # 显示离散点
  labs(title = "iCAF signatures", x = "", y = "Score in GSE10141") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 24, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 24),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("MVI+" = "#DC050C", "MVI-" = "#386cb0")) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 10,label.y = 0.48)+
  ylim(-0.57, 0.6)+
  guides(fill = FALSE) 






##Fig6E----
#TOP15 iCAF survival


##Fig6F----
table(CAF$cluster)

mCAF<-subset(CAF ,cluster %in% 'mCAF')
vCAF<-subset(CAF ,cluster %in% 'vCAF')
iCAF<-subset(CAF ,cluster %in% 'iCAF')
apCAF<-subset(CAF ,cluster %in% 'apCAF')

#CAF AddModuleScore
ST<-AddModuleScore(ST,features = list(mCAF$gene),name = 'mCAF')
ST<-AddModuleScore(ST,features = list(vCAF$gene),name = 'vCAF')
ST<-AddModuleScore(ST,features = list(iCAF$gene),name = 'iCAF')
ST<-AddModuleScore(ST,features = list(apCAF$gene),name = 'apCAF')

SpatialPlot(ST,images = 'C7',features = c('vCAF1','mCAF1','apCAF1','iCAF1'))




##Fig6G----
#load distance result
load('~/DATA/luo/zhongnan/ST/Region/distance_Bdy_TR.rdata')
distance<-distance[,8:9]

#subset TC
TC<-subset(ST, Region %in%  c('TC') & orig.ident %in% c('C2','C4','C7','C8'))
SpatialPlot(TC)

#AddMetaData
TC <- AddMetaData(object = TC,
                  metadata = distance) 

SpatialPlot(TC,images = 'C7',features = 'distance',image.alpha = 0.8)&
  theme_bw()&
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())



##Fig6H----
table(TC$orig.ident)
sample_list<-SplitObject(TC,split.by = 'orig.ident')


#plot_correlation
plot_correlation <- function(correlation_result, x, y, x_label, y_label) {
  ggplot(TC@meta.data, aes(x = x, y = y)) +
    geom_point(alpha = 1, col = '#7A67EE') +  
    stat_smooth(method = "lm", col = "red") +  
    labs(
      title = paste("Spearman's R = ", round(correlation_result$estimate, 2),  "\n", 
                    paste("P < ", formatC(correlation_result$p.value, format = "e", digits = 2))),
      x = x_label,
      y = y_label
    ) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title = element_text(size = 18, face = "bold", color = "black"),
      axis.text = element_text(size = 15, face = "bold", color = "black"),
      plot.title = element_text(size = 20, face = "bold", color = "black")
    )
}

#Distance-CAF Corr
Cor_vCAF <- cor.test(sample_list[["C7"]]$distance, sample_list[["C7"]]$vCAF1, method = 'spearman')
Cor_mCAF <- cor.test(sample_list[["C7"]]$distance, sample_list[["C7"]]$mCAF1, method = 'spearman')
Cor_apCAF <- cor.test(sample_list[["C7"]]$distance, sample_list[["C7"]]$apCAF1, method = 'spearman')
Cor_iCAF <- cor.test(sample_list[["C7"]]$distance, sample_list[["C7"]]$iCAF1, method = 'spearman')

#plot
plot_correlation(Cor_vCAF, TC@meta.data$distance, TC@meta.data$vCAF1, "Distance to Boundary", "vCAF")
plot_correlation(Cor_mCAF, TC@meta.data$distance, TC@meta.data$mCAF1, "Distance to Boundary", "mCAF")
plot_correlation(Cor_apCAF, TC@meta.data$distance, TC@meta.data$apCAF1, "Distance to Boundary", "apCAF")
plot_correlation(Cor_iCAF, TC@meta.data$distance, TC@meta.data$iCAF1, "Distance to Boundary", "iCAF")


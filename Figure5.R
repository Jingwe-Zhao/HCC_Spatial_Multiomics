######Figure5######
library(RColorBrewer)
library(Seurat)
library(Matrix)
library(ggpubr)


#载入样本
load('/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/data/HCC_region.rdata')
load("~/DATA/Pathway/Glycer_list.rdata")
Dia<-Glycer_list[[7]][["diacylglycerol biosynthetic process"]]

##Fig5A----
#HCC public scRNA-srq
load( '~/DATA/luo/tongji/HCC_nature/Scissor/Hep_sub_Scisso.rdata')
DimPlot(Hep_sub, reduction = 'umap', group.by = 'Scissor', cols = c('grey','royalblue','indianred1'), pt.size = 1.2, order = c('Scissor_Pos','Scissor_Neg'))+ggtitle('Scissor')+ theme(plot.title = element_text(size = 16))

DimPlot(Hep_sub,group.by = 'RNA_snn_res.0.4')

#Grouping
Hep_sub@meta.data$Label3<-'Background'
Hep_sub@meta.data$Label3[Hep_sub@meta.data$RNA_snn_res.0.4 %in% c('6','7')]<-'STMN1+HMGN2+GPC3+ cell'
Hep_sub@meta.data$Label3[Hep_sub@meta.data$RNA_snn_res.0.4 %in% c('0','3','4','15')]<-'Other tumor cell'
DimPlot(Hep_sub,group.by = 'Label3')

Hep_sub<-AddModuleScore(Hep_sub,features = list(Dia),name = 'Diacylglycerol' )

#exclude Background
Hep_sub2<-subset(Hep_sub, Label3 %in% c('STMN1+HMGN2+GPC3+ cell','Other tumor cell'))


my_comparisons <- list(c("STMN1+HMGN2+GPC3+ cell", "Other tumor cell"))

VlnPlot(Hep_sub2,features = c('Diacylglycerol'),group.by = 'Label3',pt.size = 0,cols = c("#386cb0","#DC050C"))+ggtitle('HRA001748')+NoLegend()+geom_boxplot(width=.2,col="black",fill="white",outlier.shape = NA)+xlab('')+ylab('DG biosynthetic process')+
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif', size=12,label.y = 0.7) +  #
  ylim(-0.4, 0.8)+ theme(text = element_text(size = 20)) + theme(axis.text = element_text(size = 18))+ theme(axis.text.x=element_text(angle=0,hjust=0.5)) #5.70  7.66




#Fig5B---- 
#DG AddModuleScore
ST<-AddModuleScore(ST,features =  list(Dia),name = 'Diacylglycerol')
stTumor<-subset(ST, region %in% 'T' & orig.ident %in%  c('C2','C4','C7','C8'))
#Grouping
stTumor$Group<-'STMN1+HMGN2+GPC3+ cell'
stTumor$Group[stTumor$orig.ident %in%  c('C2','C4')]<-'Other tumor cell'

#left
my_comparisons<-list(c('T_MVI+','T_MVI-'))

VlnPlot(stTumor, features = 'Diacylglycerol1', group.by = 'Group', pt.size = 0, cols = c('#386cb0','#DC050C')) +
  ylab('Score in ST DATA') + xlab('') + ggtitle('Diacylglycerol biosynthetic process')+
  geom_boxplot(width = 0.2, col = "black", fill = NA, outlier.shape = NA) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 24), 
        text = element_text(size = 18)) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif', size = 14, label.y = 0.75) +
  ylim(-0.78, 1) + NoLegend()


#right
LIHC<-read.csv('/home/zhaojingwei/DATA/luo/tongji/TCGA/LIHC_TCGA.csv',row.names = 1)
LIHCl1<-subset(LIHC,MVI %in% c('MVI','non_MVI'))

#MVI phenotype
meta<-LIHC[,'MVI', drop = FALSE]
meta$CB<-row.names(meta)
meta$Group<-meta$MVI
meta$Group[meta$Group %in% 'MVI']<-'MVI+'
meta$Group[meta$Group %in% 'non_MVI']<-'MVI-'
LIHCl1$MVI<-NULL
#Transpose
LIHCl1<-t(LIHCl1)

##GSVA
ES <- gsva(LIHCl1, list(Tri))
ES<-as.data.frame(t(ES))
ES$CB<-row.names(ES)

#merge GSVA and phenotype
score<-merge(ES,meta,by='CB')

my_comparisons <- list( c("MVI-", "MVI+"))
ggplot(score, aes(x = Group, y = V1, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, color = "black", alpha = 0.6) +
  labs(title = "Diacylglycerol biosynthetic process", x = "", y = "Score in TCGA LIHC") +
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


#Fig5C---- 
limit=c(-0.7, 0.7)
breaks = c(-0.7,0, 0.7)
colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100)

plot_gene = function (cluster){
  p1 = SpatialPlot(ST, features = cluster,stroke = 0.15, pt.size.factor =1.4,images = 'C2',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks) 
  p2 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 2.3,images = 'C4',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks) 
  p3 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'C7',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks) 
  p4 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'C8',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks) 
  ggarrange(p1,p2,p3,p4,ncol  = 2,nrow = 2)
}
plot_gene(cluster='Diacylglycerol')



##Fig5D----
p1<-SpatialPlot(ST,features = c('AGPAT2','LPIN3','DGAT1'),images = 'C2',pt.size.factor = 1.6,ncol = 1)
p2<-SpatialPlot(ST,features = c('AGPAT2','LPIN3','DGAT1'),images = 'C8',pt.size.factor = 1.7,ncol = 1)
ggarrange(p1,p2)



#Fig5E----

##Fig5F----
#DG survival analysis
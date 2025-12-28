###Figure3: Integrated Multi-omics analysis--------------
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(devil)


#load ST data
load('/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/data/HCC_region.rdata')
able(ST@meta.data$Region)

ST = subset(ST, sampleid %in% c("P3", "P4", "P1", "P2"))
ST = subset(ST, Region %in% c("T"))

ST@meta.data$group <- ifelse(ST@meta.data$sampleid %in% c("P1", "P2"), "MVI_Neg",
                             ifelse(ST@meta.data$sampleid %in% c("P3", "P4"), "MVI_Pos",
                                    ST@meta.data$sampleid))
table(ST$group)

counts = as.matrix(GetAssayData(ST[["Spatial"]], slot = "counts"))
# counts <- as.matrix(ST@assays$Spatial)
counts[1:5,1:5]

c = "MVI_Pos"
idx_cluster <- which(ST$group == c)
idx_others <- which(!(ST$group == c))

design_matrix <- model.matrix(~group, dplyr::tibble(group = ST$group == c))

cell_idx <- c(idx_cluster, idx_others)
# First filter
dm <- design_matrix[cell_idx,]
table(dm[,2])

fit <- devil::fit_devil(input_matrix=counts, design_matrix=dm, verbose=T, size_factors=F, overdispersion = F)
summary(fit)
colnames(fit$beta)

library(dplyr)
table(ST$sampleid)
ST$sampleid = factor(ST$sampleid, levels = c("P3", "P4", "P1", "P2"))
patient_ids <- as.numeric(as.factor(ST$sampleid))
patient_ids <- patient_ids[cell_idx]

res <- devil::test_de(fit, contrast = c(0,1), clusters = patient_ids, max_lfc = 50) %>% dplyr::mutate(cluster = c)
res

filtered_res <- res %>%
  filter(adj_pval < 0.05, abs(lfc) > 1)
dim(filtered_res)

write.csv(res, "devil_T_MVI_degs.csv")

sig_genes=read.csv("~/Project/MVI/devil_T_MVI_degs.csv")

##Fig3A1----
sig_genes=sig_genes[,-1]
library(ggplot2)
library(ggforce)
library(ggrastr)
library(ggrepel)

## define Significance
sig_genes$Significance <- ifelse(sig_genes$adj_pval<0.05 & abs(sig_genes$lfc) > 1, 
                                 ifelse(sig_genes$lfc > 1, "Up", "Down"), "NS")
table(sig_genes$Significance)
top_genes_up <- sig_genes %>%
  arrange(desc(lfc)) %>%
  head(10)

# top 10 genes
top_genes_down <- sig_genes %>%
  arrange(lfc) %>%
  head(10)

range(-log10(sig_genes$adj_pval))

p <- ggplot(sig_genes, aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(aes(color = Significance), size = 1.5, alpha = 0.8) + 
  scale_color_manual(values = c("Up" = "#FA7E4D", "Down" = "#377EB8", "NS" = "grey")) +
  theme_classic(base_size = 14) +
  ggtitle("T_MVI+ vs T_MVI-") +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 13.5, color = "black"),
    axis.title = element_text(size = 15.5, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5,size = 17.5)
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") + 
  ylab("-log10(FDR)")+
  xlab("log2(Fold Change)")+
  xlim(-40,20)+
  ylim(0,300) 
p


load('/home/zhaojingwei/DATA/luo/tongji/ST/Mycode/data/HCC_region.rdata')
table(ST@meta.data$Region)
table(ST@meta.data$sampleid)
table(ST@meta.data$Region, ST@meta.data$sampleid)

SpatialDimPlot(ST, group.by = "Region")

ST = subset(ST, sampleid %in% c("C7", "C8", "C2", "C4"))
ST = subset(ST, Region %in% c("PT", "T", "TC"))

ST@meta.data$sample_Region = paste0(ST@meta.data$sampleid, "_", ST@meta.data$Region)
table(ST@meta.data$sample_Region)

###read DEGs result
degs_mvi = read.csv("~/Project/MVI/devil_T_MVI_degs.csv")
head(degs_mvi)

# degs_mvi = degs_mvi[degs_mvi$Significance == "Up", ]
degs_mvi = degs_mvi[degs_mvi$lfc > 0, ]

top10_degs <- degs_mvi[order(degs_mvi$pval, -abs(degs_mvi$lfc)), ][1:10, ]
# top10_degs <- degs_mvi[order(degs_mvi$pval), ][1:10, ]
top10_degs
top10_degs$name


avg_expr1 <- AverageExpression(ST, 
                               assays = "Spatial", 
                               features = top10_degs$name,  
                               group.by = "sample_Region")  
avg_expr = avg_expr1$SCT

class(avg_expr)

avg_expr = avg_expr[,c("P3_T","P4_T","P1_T","P2_T","P3_TC","P4_TC","P1_TC","P2_TC","P3_PT","P4_PT","P1_PT","P2_PT")]


library(pheatmap)
# anno data
annotation_df <- data.frame(
  MVI_Status = c("MVI+", "MVI+", "MVI-", "MVI-",
                 "MVI+", "MVI+", "MVI-", "MVI-", 
                 "MVI+", "MVI+", "MVI-", "MVI-")
)
rownames(annotation_df) <- colnames(avg_expr)  

# anno col
annotation_colors <- list(
  MVI_Status = c("MVI+" = "#FA7E4D", "MVI-" = "#377EB8")
)

##Fig3A2----
p = pheatmap(avg_expr,
             color = colorRampPalette(c("navy", "white", "#FA7E4D"))(100),
             scale = "row",
             cluster_rows = TRUE, 
             treeheight_row = 0,  
             cluster_cols = FALSE,
             annotation_col = annotation_df,
             annotation_colors = annotation_colors,
             angle_col = 45,  
             border_color = NA,
             show_rownames = TRUE,
             show_colnames = TRUE,
             main = "Average Expression Heatmap",
             fontsize_row = 11,
             fontsize_col = 11)
p
pdf("~/Project/MVI/heatmap_UP_top10.pdf", width=7, height=4.5)
p
dev.off()


###Fig3B####
library(clusterProfiler)
library(org.Hs.eg.db)

sig_genes=subset(sig_genes,  adj_pval < 0.05)
diff=subset(sig_genes, abs(lfc)> 1)
write.csv(diff,file="/home/yangxinmeng/Project/MVI/T_MVI_padj0.05_lfc1_devil.csv",row.names = F)

diff_up=subset(sig_genes, lfc> 1)

sig_gene_symbols <- diff_up$name

# convert Entrez ID to SYMBOL
sig_entrez <- bitr(sig_gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID",                            
                   OrgDb = org.Hs.eg.db)


library(KEGG.db)
kk<-enrichKEGG(sig_entrez$ENTREZID,organism="hsa",pvalueCutoff=1,qvalueCutoff=1,use_internal_data = T)
# convert Entrez ID to SYMBOL
KEGG_id = setReadable(kk, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID") 


write.csv(sig_genes,file = "/home/yangxinmeng/Project/MVI/sig_genes.csv")
write.csv(KEGG_id@result,file = "/home/yangxinmeng/Project/MVI/T_MVI+_KEGG_padj0.05_lfc1_devil.csv",row.names = F)

df=read.csv("/home/yangxinmeng/Project/MVI/T_MVI+_KEGG_padj0.05_lfc1_devil.csv")
#df <- subset(KEGG_id@result, pvalue<0.05)
df1 <- subset(df, ID %in% c('hsa00190','hsa04714','hsa04110','hsa03030','hsa03430'))
df1$LogP<-log10(df1$pvalue)
df1$LogP <- -df1$LogP

df1$labelx=rep(0,nrow(df1))
df1$labely=seq(nrow(df1),1)

ggplot(data = df1, 
       aes(LogP, reorder(Description, LogP))) +
  geom_bar(stat = "identity",
           alpha = 0.5,
           fill = "#fbb1a2",
           width = 0.8) + 
  geom_text(aes(x = labelx,
                y = labely,
                label = Description),
            size = 5, 
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
    text = element_text(size = 16) 
  ) +
  xlab("-log10(pvalue)") +
  ggtitle("T_MVI+ KEGG Enrichment") + 
  scale_x_continuous(expand = c(0, 0))

###Fig3C####
sig_genes=read.csv("/home/yangxinmeng/Project/MVI/T_MVI_padj0.05_lfc1_devil.csv")
diff_down=subset(sig_genes, lfc< -1)
sig_gene_symbols <- diff_down$name

sig_entrez <- bitr(sig_gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID",                            
                   OrgDb = org.Hs.eg.db)


library(KEGG.db)
kk<-enrichKEGG(sig_entrez$ENTREZID,organism="hsa",pvalueCutoff=1,qvalueCutoff=1,use_internal_data = T)

KEGG_id = setReadable(kk, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID") 


all_path = as.data.frame(KEGG.db::KEGGPATHID2NAME)
KEGG_id@result$Description = all_path$path_name[match(KEGG_id@result$ID, all_path$path_id)]
head(KEGG_id@result$Description)


write.csv(KEGG_id@result,file = "/home/yangxinmeng/Project/MVI/T_MVI-_KEGG_padj0.05_lfc1_devil.csv",row.names = F)

df=read.csv("/home/yangxinmeng/Project/MVI/T_MVI-_KEGG_padj0.05_lfc1_devil.csv")
df2 <- subset(df, ID %in% c('hsa00982','hsa04820','hsa04080','hsa00830','hsa04974'))
df2$LogP<-log10(df2$pvalue)
df2$LogP <- -df2$LogP

df2$labelx=rep(0,nrow(df2))
df2$labely=seq(nrow(df2),1)

ggplot(data = df2, 
       aes(LogP, reorder(Description, LogP))) +
  geom_bar(stat = "identity",
           alpha = 0.5,
           fill = "#fbb1a2",
           width = 0.8) + 
  geom_text(aes(x = labelx,
                y = labely,
                label = Description),
            size = 5,  
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
    text = element_text(size = 16)  
  ) +
  xlab("-log10(pvalue)") +
  ggtitle("T_MVI- KEGG Enrichment") +
  scale_x_continuous(expand = c(0, 0.02))




##Fig3D----
#TCGA-LIHC MVI+DEG
deg<-read.csv('~/DATA/luo/zhongnan/TCGA/1011/TCGA_LIHC_MVI_DEG.csv')
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
  p1 = SpatialPlot(ST, features  = cluster,stroke = 0.15, pt.size.factor =1.4,images = 'P1',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p2 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 2.3,images = 'P2',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p3 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'P3',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p4 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'P4',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p5 = SpatialPlot(HCC1, features = cluster,stroke = 0.15,pt.size.factor = 2,image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p6 = SpatialPlot(HCC2, features = cluster,stroke = 0.15,pt.size.factor = 1.6,image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  ggarrange(p1,p2,p3,p4,p5,p6,ncol  = 2,nrow = 3)
}
plot_gene(cluster='MVI_score1')






##Fig3E----
##Scissor analysis
library(Scissor)

##HCC scRNA-seq data(HRA001748)
library(tidydr)
load('~/DATA/luo/zhongnan/HCC_nature/Hep_score.rdata')
DimPlot(Hep,group.by = 'seurat_clusters',pt.size = 0.05)

##TCGA-LIHC 
LIHC<-read.csv('/home/zhaojingwei/DATA/luo/zhongnan/TCGA/LIHC_TCGA.csv',header = T,row.names = 1)
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


##Fig3F----
cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999")

DimPlot(Hep,group.by = 'seurat_clusters',cols = cols,pt.size = 0.05)+theme_dr(xlength=0.22,ylength=0.22, arrow=grid::arrow(length=unit(0.15,"inches"),type="closed"))+theme(panel.grid=element_blank())


##Fig3G----
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



##Fig3H----
#Cluster6,7 DEG
load( '~/DATA/luo/zhongnan/HCC_nature/Scissor/DEG/Scissor_deg.rdata')
table(Scissor_deg2$cluster)

#MVI scRNA-seq UP
scDEG<-read.csv( '~/DATA/luo/zhongnan/HCC2/GSE242889/scDEG.csv')

#ST MVI UP
stDEG<- readRDS("~/DATA/luo/zhongnan/ST/DEG/ST_MVI+_Tumor_DEG.rds")
stDEG$cluster<-'ns'
stDEG$cluster[stDEG$avg_log2FC >= 0.5 & stDEG$p_val_adj < 0.05 & stDEG$pct.1 >= 0.5]<-'UP'
stDEG$cluster[stDEG$avg_log2FC <= -0.5 & stDEG$p_val_adj < 0.05 & stDEG$pct.1 >= 0.5]<-'Down'
table(stDEG$cluster)


##Fig3I----
#load TCGA-LIHC MVI DEG
deg<-read.csv('~/DATA/luo/zhongnan/TCGA/1011/TCGA_LIHC_MVI_DEG.csv')
#subset up-regulated 
UP_MVI<-subset(deg,Group %in% 'UP')


##Fig3J----
load( '~/DATA/luo/zhongnan/HCC_nature/Scissor/DEG/Scissor_deg.rdata')
table(Scissor_deg2$cluster)
Pos<-subset(Scissor_deg2, cluster %in% 'Scissor_Pos')

#HCC AddModuleScore
ST<-AddModuleScore(ST,features = list(rownames(Pos)),name = 'Scissor_Pos')
#Score in Public data
HCC1<-AddModuleScore(HCC1,features = list(rownames(Pos)),name = 'Scissor_Pos')
HCC2<-AddModuleScore(HCC2,features = list(rownames(Pos)),name = 'Scissor_Pos')

#visualization
plot_gene = function (cluster){
  p1 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'P1',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p2 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'P2',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p3 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'P3',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  p4 = SpatialPlot(ST, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'P4',image.alpha = 0)+
    theme(legend.key.size = unit(0.7, "lines"))
  ggarrange(p1,p2,p3,p4,ncol  = 2,nrow = 2)
}
plot_gene(cluster='Scissor_Pos1')

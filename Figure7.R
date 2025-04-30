######Figure7######
library(MetaboDiff)
library(Seurat)
library(Matrix)
library(ggpubr)
library(MetaboDiff)

#load MSI-ST-Seurat
load('/home/zhaojingwei/DATA/luo/zhongnan/ST/ST_MSI_combine.rdata')

#subset sample
HCC<-subset(MSI,orig.ident %in% c('P1','P2','P3','P4'))

HCC@meta.data$REGION <- ifelse(HCC@meta.data$clusters %in% c('3','5'), "C3C5", "Other")
SpatialPlot(HCC,group.by = 'REGION',cols = c('C3C5'='royalblue','Other'='grey90'))

msi_count<-as.data.frame(HCC@assays$MSI@data)
coldata<-HCC@meta.data

#read MSI metadata
Metabolic_metadata<-read.table('/home/zhaojingwei/DATA/luo/zhongnan/oebiotech/DZLM2023120285/DZLM2023120285/2.uniondata/C-2-1/pos/metadata.txt',header = T,row.names = 1)
Metabolic_metadata <- Metabolic_metadata[, 4:12]
rownames(Metabolic_metadata)<-Metabolic_metadata$Formula

#DEM by MetaboDiff
(met <- create_mae(msi_count,Metabolic_metadata,coldata))

#_SMPDBanno
met <- get_SMPDBanno(met,
                     column_kegg_id=8,
                     column_hmdb_id=5,
                     column_chebi_id=9)

(met = knn_impute(met,cutoff=0.4))

##normalize_met (vsn) 
(met <- normalize_met(met))

met = diff_test(met,
                group_factors = "REGION")


mydem=met@metadata[["ttest_REGION_C3C5_vs_Other"]]
save(mydem,file = '~/DATA/luo/zhongnan/SM/DEM/mydeg.rdata')


##Fig7A-----
load('~/DATA/luo/zhongnan/SM/DEM/mydeg.rdata')
table(mydeg$describe)

#VolcanoPlot script
source('/home/zhaojingwei/DATA/luo/zhongnan/HCC3/MetaboDiff/VolcanoPlot2.R')

dif=data.frame(
  symbol=mydeg$metabolite,
  log2FoldChange=mydeg$dm,
  padj=mydeg$adj_pval)
head(dif, n=3)

#VolcanoPlot
VolcanoPlot2(dif, padj=0.05, log2FC=0.15,title="TC VS Other region",label.symbols = c('C2H7NS','C2H7NO3S'))
#6.2,5.4



##Fig7B-----
#TC DEM Corr
library(corrplot)
library(RColorBrewer)

##读入MSI的metadata
C2Pos_metadata<-read.table('/home/zhaojingwei/DATA/luo/zhongnan/oebiotech/DZLM2023120285/DZLM2023120285/2.uniondata/C-2-1/pos/metadata.txt',header = T,row.names = 1)

sel<-C2Pos_metadata[,c(1,4,6,8,11,18)]

mydeg$Formula=mydeg$metabolite

df<-merge(mydeg,sel,by='Formula')
#write.csv(df,file = '~/DATA/luo/zhongnan/SM/DEM/HCC_TC_DEM_result.csv')
table(df$describe)

##Grouping
mydeg$describe<-'ns'
mydeg$describe[mydeg$pval<0.05 & mydeg$dm>0.15]<-'up'
mydeg$describe[mydeg$pval<0.05 & mydeg$dm<0.15]<-'down'
table(mydeg$describe)

##subset up_DEM
up_DEM<-subset(mydeg,describe %in% 'up')

#Image registration results
msi<-readRDS('/home/zhaojingwei/DATA/luo/zhongnan/HCC3/MetaboDiff/MSI_Exp.rds')

##subset up_DEM Formula
upMSI_list<-up_DEM$Formula

##subset DEM intensity
upMSI<-subset(msi,rownames(msi) %in% upMSI_list)
upMSI_transposed <- t(upMSI)

#Corr
cor_matrix <- cor(upMSI_transposed)
color2 = colorRampPalette(c("#1E90FF","white","firebrick3"))(100)
#corrplot
p<-corrplot(cor_matrix, method = 'color',tl.col = 'white',tl.cex = 0.5,col = color2,order = 'AOE')



##Fig7C----
DEM_corr<-p[["corr"]]

#Set2 DEM
mz<-DEM_corr[131:166,131:166]

#et2 DEM KEGG
#Metaboanalyst website

#read result
df<-read.csv('~/DATA/luo/zhongnan/SM/Corr/DEM/TC_msea_ora_result.csv',row.names = 1,header = 1)
df$LogP<-log10(df$Raw.p)
df$LogP <- -df$LogP
df$Description<-rownames(df)


df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)

ggplot(data = df, 
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
    axis.text.x = element_text(colour = 'black', size = 14),
    axis.ticks.x = element_line(colour = 'black', linewidth = 1),
    axis.title.x = element_text(colour = 'black', size = 16),
    text = element_text(size = 16)  # 加大全局字体大小
  ) +
  xlab("-log10(pvalue)") +
  ggtitle("TC metabolits enrichment") +
  scale_x_continuous(expand = c(0, 0))




##Fig7D----
#RCTD-MSI Corr
load('~/DATA/luo/zhongnan/SM/Corr/C2_cell_msi_cor.rdata')
c2Corr<-cor_results
load('~/DATA/luo/zhongnan/SM/Corr/C4_cell_msi_cor.rdata')
c4Corr<-cor_results
load('~/DATA/luo/zhongnan/SM/Corr/C7_cell_msi_cor.rdata')
c7Corr<-cor_results
load('~/DATA/luo/zhongnan/SM/Corr/C8_cell_msi_cor.rdata')
c8Corr<-cor_results


#correlation_list
correlation_list <- list(c2Corr, c4Corr, c7Corr, c8Corr)
#Calculate the mean
mean_correlation_matrix <- Reduce(`+`, correlation_list) / length(correlation_list)


###load TC DEM_corr
load('~/DATA/luo/zhongnan/SM/Corr/DEM_corr.rdata')

##Set2 DEM
mz<-DEM_corr[131:166,131:166]
mz<-rownames(mz)

##subset Set2 DEM
sub_cor<-subset(mean_correlation_matrix,rownames(mean_correlation_matrix) %in% mz)
sub_cor_transposed <- t(sub_cor)


library(reshape2)
sub_cor_melted <- melt(sub_cor_transposed, id.vars = "Cell_Type")

color2 = colorRampPalette(c("#1E90FF","white","firebrick3"))(100)

#heatmap
ggplot(sub_cor_melted, aes(x = Var2, y = Var1, fill = value, label = round(value, 2))) +
  geom_tile(color = "white") +  
  geom_text(color = "black", size = 2.8) +  
  scale_fill_gradientn(colors = color2, na.value = "grey50",
                       name = "Correlation") +  
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),  
        axis.text.y = element_text(color = "black"),  
        axis.title = element_text(color = "black"),   
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank()) + 
  labs(title = NULL, x = NULL, y = NULL) +
  coord_fixed() 
#12 x 6




##Fig7E----
#read CAF signatures
CAF<-read.csv('/home/zhaojingwei/DATA/luo/zhongnan/ST/CAF/CAF_Marker.csv')
table(CAF$cluster)

#four CAF signatures
caf_list <- list(
  mCAF = CAF$gene[CAF$cluster == "mCAF"],
  iCAF = CAF$gene[CAF$cluster == "iCAF"], 
  apCAF = CAF$gene[CAF$cluster == "apCAF"],
  vCAF = CAF$gene[CAF$cluster == "vCAF"]
)

#AddModuleScore
CAF_type <- c("mCAF", "iCAF", "apCAF", "vCAF")

for (caf in CAF_type) {
  for (i in 1:4) {
    msi_obj[[i]] <- AddModuleScore(
      msi_obj[[i]],
      features = list(caf_list[[caf]]),
      name = caf
    )
  }
}

#rename
for (i in 1:4) {
  colnames(msi_obj[[i]]@meta.data)[colnames(msi_obj[[i]]@meta.data) == "mCAF1"] <- "mCAF"
}

for (i in 1:4) {
  colnames(msi_obj[[i]]@meta.data)[colnames(msi_obj[[i]]@meta.data) == "iCAF1"] <- "iCAF"
}

for (i in 1:4) {
  colnames(msi_obj[[i]]@meta.data)[colnames(msi_obj[[i]]@meta.data) == "apCAF1"] <- "apCAF"
}

for (i in 1:4) {
  colnames(msi_obj[[i]]@meta.data)[colnames(msi_obj[[i]]@meta.data) == "vCAF1"] <- "vCAF"
}

#subset CAF score
C2_CAF <- msi_obj[[1]]@meta.data[,  c("mCAF", "iCAF", "apCAF", "vCAF")]
C4_CAF <- msi_obj[[2]]@meta.data[,  c("mCAF", "iCAF", "apCAF", "vCAF")]
C7_CAF <- msi_obj[[3]]@meta.data[,  c("mCAF", "iCAF", "apCAF", "vCAF")]
C8_CAF <- msi_obj[[4]]@meta.data[,  c("mCAF", "iCAF", "apCAF", "vCAF")]


#get_tau_data
get_tau_data <- function(obj) {
  msi_data <- as.data.frame(t(obj@assays[["MSI"]]@data))
  return(as.data.frame(msi_data$C2H7NO3S))
}


C2_Tau <- get_tau_data(msi_obj[[1]])
colnames(C2_Tau)<-'C2H7NO3S'
C4_Tau <- get_tau_data(msi_obj[[2]])
colnames(C4_Tau)<-'C2H7NO3S'
C7_Tau <- get_tau_data(msi_obj[[3]])
colnames(C7_Tau)<-'C2H7NO3S'
C8_Tau <- get_tau_data(msi_obj[[4]])
colnames(C8_Tau)<-'C2H7NO3S'


#merge CAF-Tau
C2_df<-cbind(C2_CAF,C2_Tau)
C4_df<-cbind(C4_CAF,C4_Tau)
C7_df<-cbind(C7_CAF,C7_Tau)
C8_df<-cbind(C8_CAF,C8_Tau)

# calculate_mean_correlation
calculate_mean_correlation <- function(df_list) {
  correlations <- lapply(df_list, cor)  # 计算每个数据框的相关性
  return(Reduce(`+`, correlations) / length(correlations))  # 计算均值
}

#mean_correlation_matrix
mean_correlation_matrix <- calculate_mean_correlation(list(C2_df, C4_df, C7_df, C8_df))

library(RColorBrewer)

colors =colorRampPalette(brewer.pal(n = 8, name ="YlOrBr"))(100)
colors = colorRampPalette(c("white", "firebrick3"))(100)

pheatmap::pheatmap(mean_correlation_matrix,color = colors,show_rownames  = T, display_numbers = T,cluster_cols = F,cluster_rows = F,
                   fontsize = 14,angle_col = 45,
                   main = "Taurine CAF Correlation",
                   legend_breaks = c(0.25,0.50,0.75,1.0),
                   number_color = "black",
                   border_color ='grey20')
#save(mean_correlation_matrix,file = '~/DATA/luo/zhongnan/ST/Fig/250108/Fig7E-Taurine-CAF-Cor.rdata')


##Fig7F----
##Fig7G----
#wet lab

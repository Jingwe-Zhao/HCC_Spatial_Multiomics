###Figure4: MVI associatecd analysis--------------
library(MetaboDiff)
library(Seurat)
library(Matrix)
library(ggpubr)


# load ST data
load('/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/data/HCC_region.rdata')

Mal<-subset(ST, Region %in%  c('T') & orig.ident %in% c('P1','P2','P3','P4'))
Mal@meta.data$Group<-'T_MVI+'
Mal@meta.data$Group[Mal@meta.data$orig.ident %in% c('P1','P2')]<-'T_MVI-'
Mal<-SCTransform(Mal, assay = "Spatial", verbose = T)

# load pathway signatures
signatures<-read.gmt('~/DATA/luo/KEGGPathway/20230205_kegg_hsa.gmt')

CellCycle<-subset(signatures,term %in% 'hsa04110_Cell_cycle')
DNArep<-subset(signatures,term %in% 'hsa03030_DNA_replication')
MisRep<-subset(signatures,term %in% 'hsa03430_Mismatch_repair')

#signatures scoring
Mal<-AddModuleScore(Mal,features = list(CellCycle),name = 'Cell cycle')
Mal<-AddModuleScore(Mal,features = list(DNArep),name = 'DNA replication')
Mal<-AddModuleScore(Mal,features = list(MisRep),name = 'Mismatch repair')

cols=c("T_MVI+" = "#DC050C", "T_MVI-" = "#386cb0")
my_comparisons<-list(c('T_MVI+','T_MVI-'))

##Fig4A----
#Cell cycle
VlnPlot(TC, features = 'Cell cycle', group.by = 'Group', pt.size = 0, cols = c('#386cb0','#DC050C')) +
  ylab('Score in ST') + xlab('') + ggtitle('Cell cycle')+
  geom_boxplot(width = 0.2, col = "black", fill = NA, outlier.shape = NA) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 24), 
        text = element_text(size = 18)) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 14, label.y = 1.25) +
  ylim(-0.3, 0.6) + NoLegend()

#DNA replication
VlnPlot(TC, features = 'DNA replication', group.by = 'Group', pt.size = 0, cols = c('#386cb0','#DC050C')) +
  ylab('Score in ST') + xlab('') + ggtitle('DNA replication')+
  geom_boxplot(width = 0.2, col = "black", fill = NA, outlier.shape = NA) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 24), 
        text = element_text(size = 18)) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 14, label.y = 1.25) +
  ylim(-0.3, 0.6) + NoLegend()

#Mismatch repair
VlnPlot(TC, features = 'Mismatch repair', group.by = 'Group', pt.size = 0, cols = c('#386cb0','#DC050C')) +
  ylab('Score in ST') + xlab('') + ggtitle('Mismatch repair')+
  geom_boxplot(width = 0.2, col = "black", fill = NA, outlier.shape = NA) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 24), 
        text = element_text(size = 18)) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 14, label.y = 1.25) +
  ylim(-0.3, 0.6) + NoLegend()


##Fig4B----
##Scissor+ geneset----
load('~/DATA/luo/tongji/HCC_nature/Scissor/DEG/Scissor_deg.rdata')
table(Scissor_deg2$cluster)
Scissor_pos<-subset(Scissor_deg2,cluster %in% 'Scissor_Pos' )

library(tidyverse)
Scissor_pos  %>%
  top_n(n = 100, wt = avg_log2FC) -> top

#load TCGA-LIHC
LIHC<-read.csv('/home/zhaojingwei/DATA/luo/tongji/TCGA/LIHC_TCGA.csv',header = T,row.names = 1)

##Scissor+ GSVA
library(GSVA)
EXP<-as.matrix(t(LIHC))

ES = gsva(EXP, list(top$gene))
ES<-as.data.frame(t(ES))
colnames(ES)<-'Scissor_Pos'
ES$CB=rownames(ES)


#Scissor AFP Corr
cli3<- read_tsv('/home/zhaojingwei/DATA/luo/tongji/TCGA/lihc_tcga_clinical_data.tsv')
rownames(cli3)<-cli3$`Sample ID`
cli3<-cli3[,3:7]
cli3$CB=paste0(cli3$`Sample ID`,'A')

cli3<-merge(cli3,ES,by='CB')
cli3$AFP=log10(cli3$`AFP At Procurement`)

cor_test <- cor.test(cli3$AFP, cli3$Scissor_Pos, method = 'spearman')
print(cor_test)


ggplot(cli3, aes(x = AFP, y = Scissor_Pos)) +
  geom_point() +
  geom_smooth(method = 'lm', color = 'blue', se = FALSE) +
  labs(title = "R=0.31 P<0.0001",
       x = "Log10 AFP At Procurement",
       y = "MVI-associated HCC cell subtype")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+ theme(text = element_text(size = 16))  # 添加边框



##Fig4C----
TCGA_LIHC_clin <- read_tsv('/home/zhaojingwei/DATA/luo/tongji/TCGA/TCGA-LIHC.GDC_phenotype.tsv.gz')
rownames(clin)<-clin$submitter_id.samples
clin$CB=clin$submitter_id.samples
clin2<-subset(clin, sample_type.samples %in% c('Primary Tumor', 'Recurrent Tumor'))
table(clin2$sample_type.samples)


ES2 <- merge(ES, clin2, by = 'CB')

rownames(ES2)<-ES2$CB
table(ES2$neoplasm_histologic_grade)

library(ggpubr)
my_comparisons <- list( c("G1", "G2"), c("G1","G3"), c( "G1","G4"), c("G2","G3"), c("G2","G4"),c( "G3","G4"))

ggplot(ES2, aes(x = neoplasm_histologic_grade, y = Scissor_Pos, fill = neoplasm_histologic_grade)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, color = "black", alpha = 0.6) +  # 显示离散点
  labs(title = "MVI-associated HCC cell subtype", x = "", y = "Scissor_Pos Score") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+NoLegend()+
  scale_fill_manual(values = c('G4'="#DC050C",  'G1'="#386cb0",'G2'="#FDB462",'G3'='#B17BA6'))+
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 9)+
  ylim(-0.6, 1.5)+RotatedAxis()



##Fig4D----
load('gse109211.RData')
exp2<-as.matrix(gse109211)
meta<-gse109211_cli
rownames(meta)=meta$ID
meta$ID=NULL

##GSVA--
library(GSVA)
ES2 <- gsva(exp2, list(rownames(top)))
ES2<-as.data.frame(t(ES2))
colnames(ES2)<-'Scissor_Pos_score'

score2<-cbind(ES2,meta)
table(score2$status)
score2$Group[score2$status %in% 'non-responder']<-'NR'
score2$Group[score2$status %in% 'responder']<-'R'


table(score2$group)
score_Sor<-subset(score2, group %in% 'treatment: Sor')


my_comparisons <- list( c("NR", "R"))

ggplot(score_Sor, aes(x = Group, y = Scissor_Pos_score, fill = Group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, color = "black", alpha = 0.6) + 
  labs(title = "MVI-associated HCC cell subtype", x = "", y = "Score in GSE109211(Sor)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 24),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("R" = "#DC050C", "NR" = "#386cb0")) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 10,label.y = 0.63)+
  ylim(-0.5, 0.8)+
  guides(fill = FALSE) 



##Fig4E----
load('GSE104580.RData')
exp2<-as.matrix(gse104580)
meta<-gse104580_cli
rownames(meta)=meta$id
meta$id=NULL

##GSVA--
library(GSVA)
ES <- gsva(exp2, list(rownames(top)))
ES<-as.data.frame(t(ES))
colnames(ES)<-'Scissor_Pos_score'

score<-cbind(ES,meta)
unique(score$group)

score$Group <- sub("(^[^_]+_[^_]+).*", "\\1", score$group)

table(score$Group)
score$Group2[score$Group %in% 'TACE_Non-Response']<-'NR'
score$Group2[score$Group %in% 'TACE_Response']<-'R'
table(score$Group2)

##box plot
library(ggpubr)
my_comparisons <- list( c("NR", "R"))

ggplot(score, aes(x = Group2, y = Scissor_Pos_score, fill = Group2)) +
  geom_boxplot(outlier.shape = NA) +  # 不显示boxplot的离群点
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, color = "black", alpha = 0.6) +  # 显示离散点
  labs(title = "MVI-associated HCC cell subtype", x = "", y = "Score in GSE104580(TACE)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 24),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("R" = "#DC050C", "NR" = "#386cb0")) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 10,label.y = 0.63)+
  ylim(-0.5, 0.8)+
  guides(fill = FALSE)
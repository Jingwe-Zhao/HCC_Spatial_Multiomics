####Figure3####
##Fig3A----
load('~/DATA/luo/tongji/HCC_nature/Scissor/Hep_sub_Scisso.rdata')

Hep_sub@meta.data$Label3<-'Background'
Hep_sub@meta.data$Label3[Hep_sub@meta.data$RNA_snn_res.0.4 %in% c('6','7')]<-'Scissor_Pos'
Hep_sub@meta.data$Label3[Hep_sub@meta.data$RNA_snn_res.0.4 %in% c('0','3','4','15')]<-'Scissor_Neg'
DimPlot(Hep_sub,group.by = 'Label3')
#remove Background
Hep_sub2<-subset(Hep_sub, Label3 %in% c('Scissor_Pos','Scissor_Neg'))

#up-regulated
TOP<-c('STMN1','H2AFZ','HMGB2' ,'PTTG1' , 'TUBA1B','TOP2A', 'JPT1',  'BIRC5','UBE2C','CKS2' ,'TUBB','UBE2S','FABP5','CDKN3','HMGN2','NUSAP1','KPNA2','HIST1H4C','GSTP1','CCNB1')

DotPlot(Hep_sub2,features = TOP,group.by = 'Label3')+RotatedAxis()+  
  scale_color_gradientn(colours = colorRampPalette(colors =rev( brewer.pal(11, name = "RdBu")))(11)) + theme(legend.text = element_text(size = 10))+ theme(legend.title = element_text(size = 12))+ylab('seurat_clusters')



##Fig3B----
p1<-FeaturePlot(Hep_sub,features = c('STMN1'),ncol = 3)+  scale_color_gradientn(colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100))+theme_dr(xlength=0.22,ylength=0.22, arrow=grid::arrow(length=unit(0.15,"inches"),type="closed"))+theme(panel.grid=element_blank())
p2<-FeaturePlot(Hep_sub,features = c('HMGN2'),ncol = 3)+  scale_color_gradientn(colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100))+theme_dr(xlength=0.22,ylength=0.22, arrow=grid::arrow(length=unit(0.15,"inches"),type="closed"))+theme(panel.grid=element_blank())
p3<-FeaturePlot(Hep_sub,features = c('TUBB'),ncol = 3)+  scale_color_gradientn(colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100))+theme_dr(xlength=0.22,ylength=0.22, arrow=grid::arrow(length=unit(0.15,"inches"),type="closed"))+theme(panel.grid=element_blank())
ggarrange(p1,p2,p3,ncol = 3)




##Fig3C----
#Prediction
expr1 <- as.data.frame(Hep_sub@assays[["RNA"]]@data)

TOP2<-c('HMGN2','STMN1')
expr <- expr1[TOP2,]

score_df <- expr %>% t() %>% as.data.frame() %>%
  mutate(CB=rownames(.))

##cluster
meta<-Hep_sub@meta.data
meta$CB<-rownames(meta)

score_df <- merge(score_df,meta[,c("CB","RNA_snn_res.0.4")],by="CB")


##cluter6,7:1,Other:0，sampling
score_df <- score_df %>%
  mutate(group=as.factor(ifelse(score_df$RNA_snn_res.0.4 %in% c(6,7),1,0)))%>%
  group_by(group) %>% 
  sample_n(1000)

table(score_df$group)


train_sub <- sample(nrow(score_df),0.8*nrow(score_df))
train_data <- score_df[train_sub,]
test_data <- score_df[-train_sub,]


#2. randomForest
gene_rf <- randomForest(group ~ HMGN2+STMN1, 
                        data=train_data, 
                        importance=TRUE,
                        proximity=TRUE)
print(gene_rf)
varImpPlot(gene_rf, main="variable importance") 
prob <- as.data.frame(gene_rf[["votes"]]) 


#ROC on the training set
ROC <- roc(train_data$group, prob$`1`,ci=T,auc=T)
plot(ROC, legacy.axes = TRUE,
     print.auc=T,
     max.auc.polygon=T,
     print.thres=T,
     auc.polygon=T,
     auc.polygon.col = "#fca311")

#Predict
pred <- predict(gene_rf,newdata = test_data)
table(test_data$group, pred)

pred <- predict(gene_rf,newdata = test_data,type="prob")
ROC <- roc(test_data$group, as.data.frame(pred)[,2],ci=T,auc=T)
plot(ROC, legacy.axes = TRUE,
     print.auc=T,
     max.auc.polygon=T,
     print.thres=T,
     auc.polygon=T,
     auc.polygon.col = "#fca311"
)

##TOP2  HMGN2+STMN1,
gene_rf1 <- randomForest(group ~ HMGN2, 
                         data=train_data,proximity=TRUE)
gene_rf2 <- randomForest(group ~ STMN1, 
                         data=train_data,proximity=TRUE)
gene_rf3 <- randomForest(group ~ HMGN2+STMN1, 
                         data=train_data,proximity=TRUE)
# ROC training set
ROC1 <- roc(train_data$group, 
            as.data.frame(gene_rf1[["votes"]])[,2],
            ci=T,auc=T)
ROC2 <- roc(train_data$group, 
            as.data.frame(gene_rf2[["votes"]])[,2],
            ci=T,auc=T)
ROC3 <- roc(train_data$group, 
            as.data.frame(gene_rf3[["votes"]])[,2],
            ci=T,auc=T)

roc_obj3 <- plot(ROC3,col ="#E71D36",print.auc=T,print.auc.x=0.8,print.auc.y=0.8)
roc_obj1 <- lines.roc(ROC2, col ="#2EC4B6") 
roc_obj2 <- lines.roc(ROC1,col ="#ff9f1c")


legend("bottomright",legend = c("HMGN2+STMN1","STMN1","HMGN2"),
       col = c("#E71D36","#2EC4B6","#ff9f1c"),
       lwd=2)
title(main = 'Classifier of Scissor Pos', line = 2.2, cex.main = 1) 



##Fig3D-E----
#wet lab

##Fig3F----
df<-read.csv('~/DATA/luo/tongji/bulkRNA/verify/gene_verify.csv',row.names = 1)
table(df$MVI)


library(ggpubr)
my_comparisons <- list( c("MVI+", "MVI-"))


p1<-ggplot(df, aes(x = MVI, y = STMN1_HMGN2, fill = MVI)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, color = "black", alpha = 0.6) +
  labs(title = "STMN1_HMGN2", x = "Group", y = "Cell ratio") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),  
    axis.text.y = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5, color = "black"),  
    plot.subtitle = element_text(color = "black"),  
    axis.title.x = element_text(color = "black"),  
    axis.title.y = element_text(color = "black"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.text = element_text(color = "black"), 
    legend.title = element_text(color = "black"), 
    plot.caption = element_text(color = "black"), 
    plot.tag = element_text(color = "black"),  
    plot.tag.position = "top" 
  ) +
  scale_fill_manual(values = c("MVI+" = "#DC050C", "MVI-" = "#386cb0")) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 8, paired = F,label.y = 20,method='t.test') +
  theme(legend.position = "none")+ theme(text = element_text(size = 24))  


p2<-ggplot(df, aes(x = MVI, y = STMN1, fill = MVI)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, color = "black", alpha = 0.6) +  # 显示离散点
  labs(title = "STMN1", x = "Group", y = "Cell ratio") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"), 
    axis.text.y = element_text(color = "black"), 
    plot.title = element_text(hjust = 0.5, color = "black"),
    plot.subtitle = element_text(color = "black"), 
    axis.title.x = element_text(color = "black"),  
    axis.title.y = element_text(color = "black"), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.text = element_text(color = "black"),  
    legend.title = element_text(color = "black"), 
    plot.caption = element_text(color = "black"),  
    plot.tag = element_text(color = "black"),  
    plot.tag.position = "top"  
  ) +
  scale_fill_manual(values = c("MVI+" = "#DC050C", "MVI-" = "#386cb0")) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 8, paired = F,label.y = 40,method='t.test') +
  theme(legend.position = "none")+ theme(text = element_text(size = 24))  


p3<-ggplot(df, aes(x = MVI, y = HMGN2, fill = MVI)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, color = "black", alpha = 0.6) +
  labs(title = "HMGN2", x = "Group", y = "Cell ratio") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),  
    axis.text.y = element_text(color = "black"),  
    plot.title = element_text(hjust = 0.5, color = "black"),  
    plot.subtitle = element_text(color = "black"), 
    axis.title.x = element_text(color = "black"), 
    axis.title.y = element_text(color = "black"), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.text = element_text(color = "black"), 
    legend.title = element_text(color = "black"), 
    plot.caption = element_text(color = "black"),  
    plot.tag = element_text(color = "black"),  
    plot.tag.position = "top"  
  ) +
  scale_fill_manual(values = c("MVI+" = "#DC050C", "MVI-" = "#386cb0")) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.label', size = 8, paired = F,label.y = 66,method='t.test') +
  theme(legend.position = "none")+ theme(text = element_text(size = 24)) 
ggarrange(p1,p2,p3,ncol = 3)



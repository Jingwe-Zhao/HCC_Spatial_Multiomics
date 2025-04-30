######Figure4######
library(RColorBrewer)

##Fig4A-D----
#oebiotech


##Fig4E-F----
#MSI_Cluster
msi<-read.csv('/home/chenweiming/Project/renwu/MSI_Cluster_character_value.csv')
msi<-msi[,1:16]
dim(msi)

library(dplyr)

column_names <- colnames(msi)[2:16]
max_values <- apply(msi[, 2:16], 1, max)
max_columns <- apply(msi[, 2:16], 1, function(x) {
  names(x)[which.max(x)]
})

msi$label <- max_columns
table(msi$label)

MSI_11.12_mz = msi[msi$label %in% c("Cluster.11","Cluster.9"), "mz"]

sel_clus = c('Cluster.2','Cluster.10', "Cluster.13",'Cluster.15')

sub_msi<-subset(msi, label %in% sel_clus)  # "Cluster.11","Cluster.12",'Cluster.2','Cluster.10', "Cluster.13",'Cluster.15'
sub_msi<-sub_msi[,c('mz',sel_clus, "label")]  # "Cluster.11","Cluster.12",'Cluster.2','Cluster.10', "Cluster.13",'Cluster.15'


#load TC DEM
load('/home/zhaojingwei/DATA/luo/zhongnan/SM/Corr/DEM_corr.rdata')
##Set1
mz1<-colnames(DEM_corr[1:130,1:130])
##Set2
mz2<-colnames(DEM_corr[131:166,131:166])
mz1
mz2

if (length(sel_clus) == 4){
  lg_label = "Tumor"
  clun = "Cluster(2.10.13.15)"
}else{
  lg_label = "TC"
  clun = "Cluster(9.11)"
}

##quantitative --
quantitative = read.csv("/home/chenweiming/Project/renwu/quantitative.csv")
quantitative[1:5, 1:5]
quantitative$Sub.Class[1:3]
colnames(quantitative)
dim(quantitative)

Super_Class_data0 = quantitative[, c("mz", "Metabolites", "Super.Class", "Class", "Sub.Class")]
head(Super_Class_data0)
dim(Super_Class_data0)

Super_Class_data0 <- Super_Class_data0 %>%
  mutate(Super.Class = ifelse(grepl(";", Super.Class), str_split_fixed(Super.Class, ";", 2)[, 1], Super.Class))
Super_Class_data0 <- Super_Class_data0 %>%
  mutate(Class = ifelse(grepl(";", Class), str_split_fixed(Class, ";", 2)[, 1], Class))
Super_Class_data0 <- Super_Class_data0 %>%
  mutate(Sub.Class = ifelse(grepl(";", Sub.Class), str_split_fixed(Sub.Class, ";", 2)[, 1], Sub.Class))
head(Super_Class_data0)
table(Super_Class_data0$Super.Class)


####### visualization
sel_cla = "Super.Class"  # Super.Class  Sub.Class  Class

if (sel_cla == "Super.Class"){
  Super_Class_data = Super_Class_data0
}else{
  Super_Class_data = subset(Super_Class_data0, subset = Super.Class == "Lipids and lipid-like molecules")
}
table(Super_Class_data$Class)
head(Super_Class_data)

Super_Class_data$class_data = Super_Class_data[[sel_cla]]

table(Super_Class_data$class_data)
library(dplyr)
if (sel_cla == "Super.Class"){
  class_counts <- Super_Class_data %>%
    dplyr::count(class_data, sort = TRUE)
  least_common_classes <- class_counts %>%  
    top_n(-6, n) %>%
    pull(class_data)
  Super_Class_data <- Super_Class_data %>%  
    mutate(class_data = ifelse(class_data %in% least_common_classes, "others", class_data))
}

table(Super_Class_data$class_data)

total_num = as.data.frame(table(Super_Class_data$class_data))
colnames(total_num) <- c("Class", "Count")
total_num$Region = rep("Total", dim(total_num)[1])
dim(total_num)
total_num

# cluster2,10,13,15
all(sub_msi$mz %in%  Super_Class_data$mz)
clusters_num = data.frame(Class = total_num$Class)
clusters_num$Region = rep("Sel_Clusters", dim(clusters_num)[1])
clusters_num$Count = rep(0, dim(clusters_num)[1])
rownames(clusters_num) = clusters_num$Class
dim(clusters_num)
head(clusters_num)

#Count the frequency
merge_data = left_join(sub_msi,Super_Class_data, by = "mz")
merge_data = na.omit(merge_data)
#Count the occurrence frequency
count_data <- merge_data %>%
  group_by(class_data) %>%
  summarise(Count = n())
clusters_num[count_data$class_data, "Count"] =  count_data$Count
dim(clusters_num)

#rbind
SuperClass_count = rbind(total_num, clusters_num)
rownames(SuperClass_count) = NULL
head(SuperClass_count)
SuperClass_count <- SuperClass_count %>%
  arrange(desc(Count))

#delete 0
`%!in%` <- Negate(`%in%`)
if (sel_cla == "Sub.Class"){
  names = SuperClass_count[which(SuperClass_count$Count == 0), "Class"]
  SuperClass_count = SuperClass_count[which(SuperClass_count$Count %!in% names),]
}

SuperClass_count$Region <- factor(SuperClass_count$Region, levels = c("Total", "TC")) # order

if (sel_cla == "Super.Class"){
  p <- ggplot(data=SuperClass_count, aes(x=reorder(Class, -Count), y=Count, fill=Region)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    theme_minimal() +
    scale_fill_manual(values=c('#E69F00','#999999'), labels=c('Total',lg_label)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 17.8, color = "black"),  
      axis.text.y = element_text(size = 17, color = "black"),  
      axis.title.x = element_text(size = 18), 
      axis.title.y = element_text(size = 18), 
      legend.text = element_text(size = 17.8), 
      legend.title = element_text(size = 17.8),    
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),  
      panel.border = element_blank(), 
      axis.line = element_line(color = "black", size = 0.5),  
      plot.margin = unit(c(1, 1, 1, 3), "cm"), 
      legend.position = "top" 
    ) +
    labs(x = "Super Class") +  
    ylim(0, max(SuperClass_count$Count) * 1.2) 
  p

  ggsave(paste0("/home/chenweiming/Project/renwu/Metabolites_Class_Analyse/", clun,"Super.Class_barplot2.pdf"), p, width = 11.6, height = 7.6, units = "in", dpi = 300, limitsize = FALSE)  # Cluster(11.12)
}

if (sel_cla == "Class"){
  p <- ggplot(data=SuperClass_count, aes(x=reorder(Class, -Count), y=Count, fill=Region)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=0.67) +
    theme_minimal() +
    scale_fill_manual(values=c('#E69F00','#999999'), labels=c('Total',lg_label)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15.3, color = "black"),  # 调整x轴字体大小和颜色
      axis.text.y = element_text(size = 14.6, color = "black"), 
      axis.title.x = element_text(size = 15.5),  
      axis.title.y = element_text(size = 15.5), 
      legend.text = element_text(size = 15), 
      legend.title = element_text(size = 15),  
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line = element_line(color = "black", size = 0.5),  
      plot.margin = unit(c(1, 1, 1, 2), "cm") ,  
      legend.position = "top" 
    ) +
    labs(x = "Class")+   # 设置y轴限制，例如将上限设置为数据最大值的115%
  
  ggsave(paste0("/home/chenweiming/Project/renwu/Metabolites_Class_Analyse/", clun,"Class_barplot.pdf"), p, width = 9, height = 7, units = "in", dpi = 300, limitsize = FALSE)   # Cluster(11.12)
}

if (sel_cla == "Sub.Class"){
  p <- ggplot(data=SuperClass_count, aes(x=reorder(Class, -Count), y=Count, fill=Region)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    theme_minimal() +
    scale_fill_manual(values=c('#E69F00','#999999'), labels=c('Total',lg_label)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), 
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 14), 
      axis.title.y = element_text(size = 14), 
      legend.text = element_text(size = 12), 
      legend.title = element_text(size = 14),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_blank(), 
      axis.line = element_line(color = "black", size = 0.5), 
      plot.margin = unit(c(1, 1, 1, 2), "cm"), 
      legend.position = "top" 
    ) +
    labs(x = "Sub Class")+ 
    ylim(0, max(SuperClass_count$Count) * 1.15)  
  

  ggsave(paste0("/home/chenweiming/Project/renwu/Metabolites_Class_Analyse/", clun,"Sub.Class_barplot.pdf"), p, width = 11, height = 6, units = "in", dpi = 300, limitsize = FALSE)  # Cluster(11.12)
}


###Fisher test----
setwd("/home/chenweiming/Project/renwu/Metabolites_Class_Analyse/")

#all_mz
all_mz = Super_Class_data0$mz
all_mz_num = length(all_mz)

#DEM
sel_clus = c('Cluster.2','Cluster.10', "Cluster.13",'Cluster.15')
# sel_clus = c('Cluster.9','Cluster.11')
sel_clus_mz = msi[msi$label %in% sel_clus, "mz"]
sel_clus_mz_num = length(sel_clus_mz)
head(sel_clus_mz)
all(sel_clus_mz %in% all_mz)

if (length(sel_clus) == 4){
  lg_label = "Tumor"
  clun = "Cluster(2.10.13.15)"
}else{
  lg_label = "TC"
  clun = "Cluster(9.11)"
}


class_type = "Class"  # Super.Class  Sub.Class  Class

if (class_type == "Super.Class"){
  Super_Class_data = Super_Class_data0
  all_mz = Super_Class_data$mz
  all_mz_num = length(all_mz)
}else{
  Super_Class_data = subset(Super_Class_data0, subset = Super.Class == "Lipids and lipid-like molecules")
  
  if (class_type == "Sub.Class"){
    tb = table(Super_Class_data$Sub.Class)
    Super_Class_data = subset(Super_Class_data, Sub.Class %in% names(tb[tb>=10]))  
  }
}
head(Super_Class_data)
table(Super_Class_data[[class_type]])

sel_class = unique(Super_Class_data[[class_type]])
sel_class

fhtest_results <- data.frame(
  Description = character(),
  p.value = numeric(),
  estimate = numeric(),
  conf.int_lower = numeric(),
  conf.int_upper = numeric(),
  MCRatio = character(), 
  BgRatio = character(),
  Count = numeric(), 
  MZ = character()
)

for (i in seq_along(sel_class)){
  print(paste0("------------ ", sel_class[i], " ------------------------"))
  class_mz = Super_Class_data[Super_Class_data[[class_type]] == sel_class[i], "mz"]
  Class_mz_num = length(class_mz)
  
  inter_mz = intersect(sel_clus_mz,class_mz)
  inter_num = length(inter_mz)
  num_12 = Class_mz_num-inter_num
  num_21 = sel_clus_mz_num-inter_num
  num_22 = all_mz_num-sel_clus_mz_num-num_12
  
  test_df = data.frame(In_cluster = c(inter_num, num_12),
                       Not_in_cluster = c(num_21, num_22))
  rownames(test_df) = c("In_Class", "Not_in_Class")
  
  test_result = fisher.test(test_df)
  
  result <- data.frame(
    Description = sel_class[i],
    p.value = test_result$p.value,
    estimate = test_result$estimate,
    conf.int_lower = test_result$conf.int[1],
    conf.int_upper = test_result$conf.int[2],
    MCRatio = paste0(inter_num, "/", sel_clus_mz_num),  
    BgRatio = paste0(Class_mz_num, "/", all_mz_num),
    Count = inter_num,
    MZ = ifelse(all(is.na(inter_mz)), "", paste(na.omit(inter_mz), collapse = "/"))
  )
  rownames(result) = NULL
  fhtest_results = rbind(fhtest_results, result) 
}

fhtest_results$ClassType = class_type
fhtest_results = fhtest_results[, c("ClassType", colnames(fhtest_results)[-dim(fhtest_results)[2]])]
# odds ratio forest
data_sig = fhtest_results
class(data_sig)

library(stringr) 
data_sig$Description
head(data_sig)
#order
data_sig$Description <- reorder(data_sig$Description, data_sig$estimate, FUN = function(x) x)

p = ggplot(data = data_sig, aes(x = estimate, y = Description,
                                xmin = conf.int_lower, xmax = conf.int_upper)) +
  geom_pointrange(color = "#999999", fatten = 1.5, size = 1, linewidth = 1) + 
  geom_point(color = "#97a0c8", size = 3) +  
  theme_bw() +
  labs(x = NULL, y = "Odds ratio", title = class_type) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45, color = "black", size = 16.5),
        axis.text.y = element_text(hjust = 1, color = "black", size = 20),
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 22))
p

ggsave(paste0(clun,"_", class_type,"_odds_ratio_forestplot.pdf"), plot = p, dpi = 400, width = 7.5, height = 7.5)



data_sig$NumBeforeSlash <- sapply(strsplit(data_sig$BgRatio, "/"), function(x) as.integer(x[1]))
rows_with_num_gt_10 <- which(data_sig$NumBeforeSlash > 10)
data_sig = data_sig[rows_with_num_gt_10, ]
data_sig

library(stringr) 
data_sig$Description <- str_wrap(data_sig$Description, width = 20) 
data_sig$Description
head(data_sig)
data_sig$Description <- reorder(data_sig$Description, -data_sig$estimate, FUN = function(x) x)


data_sig$point_color <- ifelse(data_sig$p.value < 0.05, "#e97785", "#97a0c8")
p = ggplot(data = data_sig, aes(x = Description, y = estimate,
                                ymin = conf.int_lower, ymax = conf.int_upper))+
  geom_pointrange(color = "#999999", fatten = 1.5, size = 1, linewidth = 1) +
  geom_point(color = data_sig$point_color, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#999999", size = 1) + 
  theme_bw() +
  labs(x = NULL, y = "Odds ratio", title = "Distribution of Metabolites in Clusters 9 and 11") +
  theme(axis.text.x = element_text(hjust = 1, angle = 45, color = "black", size = 16.5),
        axis.text.y = element_text(hjust = 1, color = "black", size = 20),
        text = element_text(size = 20),
        panel.grid = element_blank(), #
        plot.title = element_text(hjust = 0.5, face = "bold", size = 22)) +
  ylim(0, max(data_sig$conf.int_upper))+  
  scale_y_continuous(breaks = function(x) c(1, pretty(x)))  
p

ggsave(paste0(clun,"_", class_type,"_odds_ratio_forestplot2.pdf"), plot = p, dpi = 400, width = 10.5, height = 7.5)

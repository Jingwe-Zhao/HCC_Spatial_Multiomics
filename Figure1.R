###HCC Spatial Multiomics####
library(Matrix, lib.loc = "/public/software/miniconda3/envs/r4.1.3/lib/R/library")
library(Seurat)
library(tidyverse)
library(ggpubr)
library(tidydr)


#SpaceRanger-out
input_data_dir <-"/home/bioinformatics/Data_download/zhaojw/Spatial_data/WFB_HCC/1.SpaceRanger"
sample_list <- list.files(input_data_dir)
sample_list

sample_list <- sample_list[c(2:7)]
sample_list

## file path
samples_dir<- sample_list %>% file.path(input_data_dir, .)
samples_dir
##slice id
sample_names <- c("C2", "C4", "C7", "C8", "P2", "P7")

## read samples
sample_objects <- purrr::map(1:length(sample_list), function(x) {
  ## read data
  one_dir <- samples_dir[x]
  sample_id <- sample_list[x]
  slice_id <- sample_names[x]
  sample_object <- Load10X_Spatial(
    data.dir = one_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = slice_id,
    filter.matrix = TRUE
  )
  sample_object@project.name <- sample_id
  sample_object@meta.data$orig.ident <- slice_id
  sample_object <- RenameCells(object = sample_object, add.cell.id = slice_id)
  
  return(sample_object)
})

#SCT transform
sample_objects <- lapply(sample_objects, 
                         SCTransform, 
                         assay = "Spatial", 
                         method = "poisson",variable.features.n = 3000)
#merge
ST<-merge(sample_objects[[1]], y = sample_objects[2:6])
DefaultAssay(ST) <- "SCT"

#VariableFeatures duplicate removal
uniFeature <- unique(unlist(lapply(sample_objects, VariableFeatures)))
VariableFeatures(ST) <- uniFeature

#PCA
ST<-RunPCA(ST, assay = "SCT",features = VariableFeatures(object = ST), verbose = FALSE)

library(harmony)
ST <- RunHarmony(ST,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")

#FindClusters
ST <- FindNeighbors(ST, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4)
ST <- RunUMAP(ST, dims = 1:20, reduction = "harmony")

#rename Cluster+1
ST@meta.data$clusters <- as.numeric(as.character(ST@meta.data$SCT_snn_res.0.4)) + 1
ST@meta.data$clusters <- as.factor(ST@meta.data$clusters)



##Fig1A----
col = c('1'="#7fc97f", '2'="#beaed4", '3'="#fdc086", '4'="#386cb0", '5'="#f0027f", '6'="#a34e3b", '7'="#666666", '8'="#1b9e77", '9'="#d95f02", '10'="#6c67a5", '11'="#d01b2a", '12'="#43acde", '13'="#efbd25")

umap<-DimPlot(ST,cols = col,group.by = 'clusters')+theme_dr(xlength=0.22,ylength=0.22, arrow=grid::arrow(length=unit(0.15,"inches"),type="closed"))+theme(panel.grid=element_blank())
umap
ggsave(plot = umap,filename = "/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/Fig/Fig1/Fig1A_ST_UMAP.pdf",width = 4.6,height = 3.8)




##Fig1B----
list<-c('HAMP','MT1X','MT1G','KIAA1522','IL11RA','CD3D','CD3E','IGHG1','IGHA1','COL1A1','COL1A2','ACTA2','MYH11','SMO','IFI6','TMEM45B','LAMB3','PKLR','SYT7','SLC22A7','SERPINA7','TUBB','GPC3','LGALS4')
#AverageHeatmap script
source('/home/zhaojingwei/DATA/code/AverageHeatmap.R')
Heatmap<-AverageHeatmap(object = ST,
               group.by = 'clusters',
               assays = "SCT",
               markerGene = list,
               htCol = c("#0099CC", "white", "#CC0033"),
               cluster.order = c(12,10,1,6,8,7,5,3,11,13,9,4,2),
               column_names_rot = -90,
               row_title ='',
               annoCol=T,
               clusterAnnoName = F,
               myanCol=c("#43acde","#6c67a5","#7fc97f","#a34e3b",
                         "#1b9e77","#666666","#f0027f","#fdc086",
                         "#d01b2a","#efbd25","#d95f02","#386cb0","#beaed4"))
pdf("/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/Fig/Fig1/Fig1B_Heatmap.pdf", width = 3.9, height = 5.8)

Heatmap
dev.off()



##Fig1C-E----
#H&E + Clusters
h1<-SpatialPlot(ST,  pt.size.factor =0,images = 'C2')+NoLegend()
c1<-SpatialPlot(ST, group.by = 'clusters',stroke = 0.15, pt.size.factor =1.4,images = 'C2',image.alpha = 0,cols = col)+NoLegend()

h2<-SpatialPlot(ST,  pt.size.factor =0,images = 'C4')+NoLegend()
c2<-SpatialPlot(ST, group.by = 'clusters',stroke = 0.15, pt.size.factor =2.3,images = 'C4',image.alpha = 0,cols = col)+NoLegend()

h3<-SpatialPlot(ST,  pt.size.factor =0,images = 'C7')+NoLegend()
c3<-SpatialPlot(ST, group.by = 'clusters',stroke = 0.15, pt.size.factor =1.6,images = 'C7',image.alpha = 0,cols = col)+NoLegend()

h4<-SpatialPlot(ST,  pt.size.factor =0,images = 'C8')+NoLegend()
c4<-SpatialPlot(ST, group.by = 'clusters',stroke = 0.15, pt.size.factor =1.5,images = 'C8',image.alpha = 0,cols = col)+NoLegend()

h5<-SpatialPlot(ST,  pt.size.factor =0,images = 'P2')+NoLegend()
c5<-SpatialPlot(ST, group.by = 'clusters',stroke = 0.15, pt.size.factor =1.6,images = 'P2',image.alpha = 0,cols = col)+NoLegend()

h6<-SpatialPlot(ST, pt.size.factor =0,images = 'P7')+NoLegend()
c6<-SpatialPlot(ST, group.by = 'clusters',stroke = 0.15, pt.size.factor =1.3,images = 'P7',image.alpha = 0,cols = col)+NoLegend()

p<-ggarrange(h1,c1,h3,c3,h5,c5,h2,c2,h4,c4,h6,c6,ncol = 6,nrow = 2)
p
ggsave(plot = p,filename = "/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/Fig/Fig1/Fig1_C-E.pdf",width = 11.8,height = 5.3)


##Fig1 F----
#cellRatio script
source('/home/zhaojingwei/DATA/code/cellRatioPlot.R')
cellRatio<-cellRatioPlot(object = ST,
                         sample.name = "orig.ident",
                         celltype.name = "clusters",
                         flow.curve = 0.5,fill.col = col)+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5))+xlab('Sample')
cellRatio 
ggsave("/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/Fig/Fig1/Fig1F_Cluster_Ratio.pdf", plot = cellRatio, width = 7.2, height = 5.1)


###Fig1G----
#load RCTD result
HCC<-readRDS('/home/bioinformatics/Project/zhluo/HCC_wangfubing/HCC_RCTD_SeuratObject.rds')

cell_prop = as.data.frame(HCC@assays[["RCTD"]]@counts) %>% t() %>% 
  as.data.frame() %>%
  rownames_to_column("spot_id") %>%
  pivot_longer(-spot_id)
cell_prop$id<-cell_prop$spot_id

#head(cell_prop)
niche_info = HCC@meta.data %>% as.data.frame() %>% 
  rownames_to_column("spot_id") %>%
  select_at(c("spot_id", "orig.ident", "clusters")) %>%
  dplyr::rename("Niche" = "clusters") %>%
  mutate(Niche = paste0("Cluster_", Niche))


cellprops_info =  cell_prop  %>%
  left_join(niche_info, c("spot_id" = "spot_id")) %>%
  na.omit()

#head(cellprops_info)
cell_props_summary_CT_pat <- cellprops_info %>%
  group_by(name, Niche) %>%
  summarize(median_CT = median(value))
head(cell_props_summary_CT_pat)

# Check niches that are unique to some patients
cell_prop <-  cell_prop %>%
  left_join(niche_info, c("spot_id" = "spot_id")) %>%
  na.omit() %>%
  dplyr::select(-spot_id) %>%
  group_by(name) %>%
  nest() %>%
  mutate(wres = map(data, function(dat) {
    
    niches <- dat$Niche %>%
      unique() %>%
      set_names()
    
    map(niches, function(g) {
      
      test_data <- dat %>%
        mutate(test_group = ifelse(.data[["Niche"]] == g,
                                   "target", "rest")) %>%
        mutate(test_group = factor(test_group,
                                   levels = c("target", "rest")))
      
      wilcox.test(value ~ test_group, 
                  data = test_data,
                  alternative = "greater") %>%
        broom::tidy()
    }) %>% enframe("Niche") %>%
      unnest()
    
  }))

wilcox_types <- cell_prop %>%
  dplyr::select(wres) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(adj_pval = p.adjust(p.value)) %>%
  dplyr::mutate(log_adj_pval = -log10(adj_pval)) %>%
  dplyr::mutate(sign = ifelse(adj_pval < 0.005, "*", ""))

#head(wilcox_types)
ct_median_desc <- cell_prop %>%
  dplyr::select(data) %>%
  unnest() %>%
  group_by(name, Niche) %>%
  summarise(median_prop = median(value)) %>%
  mutate(scaled_median_prop = (median_prop - mean(median_prop))/sd(median_prop))
dim(ct_median_desc)

niche_car_df <- left_join(wilcox_types, ct_median_desc) %>% na.omit()
dim(niche_car_df)

ct_median_desc_mat <- niche_car_df %>% dplyr::select(name, Niche, scaled_median_prop) %>%
  pivot_wider(names_from = Niche, values_from = scaled_median_prop) %>%
  column_to_rownames("name") %>%
  as.matrix()

ct_sign_desc_mat <- niche_car_df %>% dplyr::select(name, Niche, sign) %>%
  pivot_wider(names_from = Niche, values_from = sign) %>%
  column_to_rownames("name") %>%
  as.matrix()

ct_median_desc_mat
ct_sign_desc_mat

library("ComplexHeatmap")
Heatmap(ct_median_desc_mat, name = "scaled comp", 
        rect_gp = gpar(col = "black", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf(ct_sign_desc_mat[i, j]), x, y, gp = gpar(fontsize = 14))
        })




###Fig1H----
library(mistyR)
source("./misty_utilities.R")

run_colocalization <- function(slide, 
                               assay, 
                               useful_features, 
                               out_label, 
                               misty_out_alias = "./results/tissue_structure/misty/cell_map/cm_") {
  
  # Define assay of each view ---------------
  view_assays <- list("main" = assay,
                      "juxta" = assay,
                      "para" = assay)
  # Define features of each view ------------
  view_features <- list("main" = useful_features, 
                        "juxta" = useful_features,
                        "para" = useful_features)
  # Define spatial context of each view -----
  view_types <- list("main" = "intra", 
                     "juxta" = "juxta",
                     "para" = "para")
  # Define additional parameters (l in case of paraview,
  # n of neighbors in case of juxta) --------
  view_params <- list("main" = NULL, 
                      "juxta" = 5,
                      "para" = 15)
  
  misty_out <- paste0(misty_out_alias, 
                      out_label, "_", assay)
  
  run_misty_seurat(visium.slide = slide,
                   view.assays = view_assays,
                   view.features = view_features,
                   view.types = view_types,
                   view.params = view_params,
                   spot.ids = NULL,
                   out.alias = misty_out)
  
  return(misty_out)
}

pathway_obj = readRDS(file="/home/bioinformatics/Project/zhluo/HCC_wangfubing/HCC_RCTD_SeuratObject.rds")

sample_obs_list = list()
samplenames = c("P1", "P2", "P3", "P4", "N1", "N3")
for (i in c(1:length(samplenames))){
  samplename = samplenames[i]
  print(paste0("-------------------------------", samplename, "----------------------------------"))
  ST_obs_sub = subset(pathway_obj, subset = sampleid == samplename)
  ST_obs_sub@images <- ST_obs_sub@images[c(samplename)]
  sample_obs_list[[i]] <- ST_obs_sub}
pathway_obj = sample_obs_list


for (i in 1:6){
  assay_label <- "RCTD"
  slide = pathway_obj[[i]]
  out_dir_1 = file.path(getwd(), "cell_type_data","misty_")
  sample_name = names(slide@images)
  
  assay <- assay_label
  DefaultAssay(slide) <- assay
  useful_features <- rownames(slide)
  useful_features <- useful_features[! useful_features %in% "prolif"]
  
  mout <- run_colocalization(slide = slide,
                             useful_features = useful_features,
                             out_label = sample_name,
                             assay = assay,
                             misty_out_alias = out_dir_1)
}

misty_outs <- list.dirs(file.path(getwd(), "cell_type_data"),full.names = T,recursive = F)
misty_outs
misty_outs <- misty_outs[grepl("RCTD", misty_outs)]
misty_res <- collect_results(misty_outs)

pdf("/home/chenweiming/Jupyter/renwu/6.03rw/cell_type_data/MistyR_ct_interaction_heatmap.pdf", height = 5, width = 5)
mistyR::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
dev.off()




###Fig1I----
#region defination
ST@meta.data$Region<-NA
ST@meta.data$Region[ST@meta.data$cluster %in% c('2','4','9','11')] <- "PT"
ST@meta.data$Region[ST@meta.data$cluster %in% c('1','6','7','8','10','12','13')] <- "T"
##C2样本的Cluster4为癌灶
ST@meta.data$Region[ST@meta.data$cluster %in% c('4') & ST@meta.data$orig.ident %in% c('C2')] <- "T"
ST@meta.data$Region[ST@meta.data$cluster %in% c('3','5')] <- "TC"
##C2样本的Cluster7为过渡区TR
ST@meta.data$Region[ST@meta.data$cluster %in% c('7')& ST@meta.data$orig.ident %in% c('C7')] <- "TR"

col2<-c('PT'='#FDB462','TC'='#1965B0','T'='#DC050C','TR'='#B17BA6')

Idents(ST)<-ST@meta.data$Region
plot_gene = function (cluster){
  p1 = SpatialPlot(ST, group.by = cluster,stroke = 0.15, pt.size.factor =1.4,images = 'C2',image.alpha = 0,cols = col2)
  p2 = SpatialPlot(ST, group.by = cluster,stroke = 0.15,pt.size.factor = 2.3,images = 'C4',image.alpha = 0,cols = col2)
  p3 = SpatialPlot(ST, group.by = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'C7',image.alpha = 0,cols = col2)
  p4 = SpatialPlot(ST, group.by = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'C8',image.alpha = 0,cols = col2)
  ggarrange(p1,p2,p3,p4,ncol  = 4,nrow = 1)
}
Region<-plot_gene(cluster='Region')
Region

ggsave("/home/zhaojingwei/DATA/luo/zhongnan/ST/Mycode/Fig/Fig1/Fig1I_Region.pdf", plot = Region, width = 10.8, height = 3)


#save
saveRDS(ST,file = '/home/zhaojingwei/DATA/luo/zhongnan/HCC2/HCC6-defined.rds')

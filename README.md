## Project Overview
### In this project, we integrated *spatial transcriptomics*, *spatial metabolomics*, *single-cell transcriptomics* and *bulk RNA-seq* data to identify the metabolic-transcriptomic features on tumor foci and tumor capsule in Microvascular Invasion HCC
<img width="2614" height="1481" alt="image" src="https://github.com/user-attachments/assets/c8c3e1e5-3aa5-4722-92ef-f2c9d1ff67ac" />


## Analysis workflow
#### :pencil:Figure1: Spatial transcriptomics analysis
- [ ] Unsupervised clustering, Differentially Expressed Genes(DEGs) , Pathological annotation, Cell percent ratio, RCTD cell type deconvolution, Cell type spatial co-localization, Spatial region division

#### :pencil:Figure2: Spatial metabolomics analysis
- [ ] Unsupervised clustering, Pathological annotation, Differentially Expressed Metabolites(DEMs), Metabolic cluster correlation, Distribution of Metabolites

#### :pencil:Figure3: Integrated Multi-omics analysis
- [ ] DEGs,  Marker gene expression, KEGG pathway enrichment, Spatial gene expression scoring , Phenotypic association analysis, Unsupervised clustering,  Comparison of DEGs
#### :pencil:Figure4: Clinical relevance analysis
- [ ] Marker gene expression, Random Forest model, Clinical indicator distribution, Cell subtype proportion, Disease-free survival (DFS) and overall survival (OS) analysis
#### :pencil:Figure5: Metabolic reprogramming analysis
- [ ] Metabolic Score, Spatial distribution of metabolic scoring, Spatial expression of enzyme genes, Spatial intensity mapping of metabolites, Overall survival (OS) analysis
#### :pencil:Figure6: Spatial microenvironment analysis
- [ ] Identification of tumor capsule and tumor boundary, RCTD cell type deconvolution,  KEGG pathway enrichment, Cancer-associated fibroblast (CAF) subtype scoring, Overall survival (OS) analysis,  Spatial distribution of CAF subtypes, Correlation between CAF subtypes and distance
#### :pencil:Figure7: CAF subtypes specific metabolite analysis
- [ ] Differentially Expressed Metabolites(DEMs), Metabolites correlation, KEGG pathway enrichment, CAF-Metabolites correlation, iCAF-Taurine correlation

## SessionInfo
### R version 4.1.3 (2022-03-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)
### R software version
Seurat_4.3.0 , SeuratObject_4.1.3 , Matrix_1.5-4.1, sctransform_0.4.1, harmony_1.2.0, spacexr_2.2.1, mistyR_1.10.1  , Scissor_2.0.0, clusterProfiler_4.8.3, GSVA_1.48.3 , ggplot2_3.4.4 , ggpubr_0.6.0  , pheatmap_1.0.12 , infercnv_1.19.1 , GSEABase_1.62.0 ,  patchwork_1.1.3 ,  randomForest_4.7-1.1, Cardinal_3.0.1, survival_3.8-3, MetaboDiff_0.9.5, ComplexHeatmap_2.15.2, CellChat_2.1.1, devil_0.1.0

## Data access 
### Generated in this study
#### Spatial transcriptomics data: 
:link: https://ngdc.cncb.ac.cn/gsa-human/browse/HRA011344
#### Spatial metabolomics data: 
:link: https://ngdc.cncb.ac.cn/omix/release/OMIX006935

### Public data
#### scRNA-Seq data: 
:white_check_mark: Gene Expression Omnibus (GSE149614 and GSE242889)

:white_check_mark: National Genomics Data Center (HRA001748)
#### Spatial transcriptomics data: 
:white_check_mark: National Genomics Data Center (HRA000437)
#### bulk RNA-seq data: 
:white_check_mark: UCSC Xena (the TCGA-LIHC, IGCG-LIHC)

:white_check_mark: Gene Expression Omnibus (GSE14520, GSE40873 & GSE54236)




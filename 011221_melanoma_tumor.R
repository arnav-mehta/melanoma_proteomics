# analysis of melanoma tumor data obtained from Keren and Moshe
# samples were processed using the SmartSeq2 platform

rm(list=ls())

library(dplyr)
library(Seurat)
library(cowplot)
library(xlsx)
library(ggplot2)

data_folder <- "/Users/arnavmehta/Dropbox (Personal)/Manuscripts/Cytokine, ctDNA, exosome ML immuntherapy responders/Analysis/scRNAseq data/melanoma-cell-keren/counts/"
setwd(data_folder)
save_folder <- "/Users/arnavmehta/Dropbox (Personal)/Manuscripts/Cytokine, ctDNA, exosome ML immuntherapy responders/Analysis/scRNAseq data/melanoma-cell-keren/"

# # load Seurat objects for each sample
# P2.counts <- read.csv("./sc_tumor_Data_post_P2.csv", row.names = 1, check.names = FALSE)
# P2 <- CreateSeuratObject(counts=P2.counts, project="melanoma_scrnaseq", 
#                               min.cells=0, min.features=0)
# P2$sample.ident = "P2"
# 
# P5.counts <- read.csv("./sc_tumor_Data_post_P5.csv", row.names = 1, check.names = FALSE)
# P5 <- CreateSeuratObject(counts=P5.counts, project="melanoma_scrnaseq", 
#                          min.cells=0, min.features=0)
# P5$sample.ident = "P5"
# 
# P28.counts <- read.csv("./sc_tumor_Data_post_P28.csv", row.names = 1, check.names = FALSE)
# P28 <- CreateSeuratObject(counts=P28.counts, project="melanoma_scrnaseq", 
#                          min.cells=0, min.features=0)
# P28$sample.ident = "P28"
# 
# P28.2.counts <- read.csv("./sc_tumor_Data_post_P28_2.csv", row.names = 1, check.names = FALSE)
# P28.2 <- CreateSeuratObject(counts=P28.2.counts, project="melanoma_scrnaseq", 
#                           min.cells=0, min.features=0)
# P28.2$sample.ident = "P28.2"
# 
# P30.counts <- read.csv("./sc_tumor_Data_post_P30.csv", row.names = 1, check.names = FALSE)
# P30 <- CreateSeuratObject(counts=P30.counts, project="melanoma_scrnaseq", 
#                          min.cells=0, min.features=0)
# P30$sample.ident = "P30"
# 
# P34.counts <- read.csv("./sc_tumor_Data_pre_P34.csv", row.names = 1, check.names = FALSE)
# P34 <- CreateSeuratObject(counts=P34.counts, project="melanoma_scrnaseq", 
#                           min.cells=0, min.features=0)
# P34$sample.ident = "P34"
# 
# P38.counts <- read.csv("./sc_tumor_Data_pre_P38.csv", row.names = 1, check.names = FALSE)
# P38 <- CreateSeuratObject(counts=P38.counts, project="melanoma_scrnaseq", 
#                             min.cells=0, min.features=0)
# P38$sample.ident = "P38"
# 
# # include the percent of mitochondrial genes
# P2[["percent.mt"]] <- PercentageFeatureSet(P2, pattern = "^MT-")
# P5[["percent.mt"]] <- PercentageFeatureSet(P5, pattern = "^MT-")
# P28[["percent.mt"]] <- PercentageFeatureSet(P28, pattern = "^MT-")
# P28.2[["percent.mt"]] <- PercentageFeatureSet(P28.2, pattern = "^MT-")
# P30[["percent.mt"]] <- PercentageFeatureSet(P30, pattern = "^MT-")
# P34[["percent.mt"]] <- PercentageFeatureSet(P34, pattern = "^MT-")
# P38[["percent.mt"]] <- PercentageFeatureSet(P38, pattern = "^MT-")
# 
# # Vvsualize QC metrics as a violin and scatterplots
# sample <- P38
# pdf(paste0(save_folder,"results/P38_QC.pdf"))
# VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3, group.by = "sample.ident")
# plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample.ident") + NoLegend() +
#   FontSize(x.text=10, y.text=10)
# plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample.ident") + NoLegend()+
#   FontSize(x.text=10, y.text=10)
# plot1 + plot2
# dev.off()
# 
# # filtering steps
# P2 <- subset(P2, subset = nFeature_RNA > 500 & percent.mt < 20)
# P5 <- subset(P5, subset = nFeature_RNA > 500 & percent.mt < 20)
# P28 <- subset(P28, subset = nFeature_RNA > 500 & percent.mt < 20)
# P28.2 <- subset(P28.2, subset = nFeature_RNA > 500 & percent.mt < 20)
# P30 <- subset(P30, subset = nFeature_RNA > 500 & percent.mt < 20)
# P34 <- subset(P34, subset = nFeature_RNA > 500 & percent.mt < 20)
# P38 <- subset(P38, subset = nFeature_RNA > 500 & percent.mt < 20)
# 
# # normalize data
# P2 <- NormalizeData(P2, normalization.method = "LogNormalize", scale.factor = 10000)
# P5 <- NormalizeData(P5, normalization.method = "LogNormalize", scale.factor = 10000)
# P28 <- NormalizeData(P28, normalization.method = "LogNormalize", scale.factor = 10000)
# P28.2 <- NormalizeData(P28.2, normalization.method = "LogNormalize", scale.factor = 10000)
# P30 <- NormalizeData(P30, normalization.method = "LogNormalize", scale.factor = 10000)
# P34 <- NormalizeData(P34, normalization.method = "LogNormalize", scale.factor = 10000)
# P38 <- NormalizeData(P38, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# # merge objects
# tumor_cells <- merge(P2, y=c(P5, P28, P28.2, P30, P34, P38),
#                      add.cell.ids = c("P2", "P5", "P28", "P28.2", "P30", "P34", "P38"),
#                      project = "melanoma_scrnaseq")
# 
# # find highly variable features
# tumor_cells <- FindVariableFeatures(tumor_cells, selection.method = "vst", nfeatures = 2000)
# 
# # identify the 10 most highly variable genes
# # plot variable features with and without labels
# top10 <- head(VariableFeatures(tumor_cells), 10)
# plot1 <- VariableFeaturePlot(tumor_cells)
# plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# pdf(paste0(save_folder,"results/merged_variable_features.pdf"))
# plot1
# dev.off()
# 
# # scale the data
# all.genes <- rownames(tumor_cells)
# tumor_cells <- ScaleData(tumor_cells, features = all.genes)
# 
# # PCA
# tumor_cells <- RunPCA(tumor_cells, features = VariableFeatures(object = tumor_cells))
# print(tumor_cells[["pca"]], dims = 1:5, nfeatures = 5)
# 
# pdf(paste0(save_folder,"results/merged_PCA_loadings.pdf"))
# VizDimLoadings(tumor_cells, dims = 1:2, reduction = "pca")
# dev.off()
# 
# pdf(paste0(save_folder,"results/merged_PCA.pdf"))
# DimPlot(tumor_cells, reduction = "pca", group.by = "sample.ident")
# dev.off()
# 
# pdf(paste0(save_folder,"results/merged_PCA_loadings_heatmap.pdf"))
# DimHeatmap(tumor_cells, dims = 1:9, cells = 500, balanced = TRUE)
# dev.off()
# 
# # assessing dimensionality
# tumor_cells <- JackStraw(tumor_cells, num.replicate = 100)
# tumor_cells <- ScoreJackStraw(tumor_cells, dims = 1:20)
# 
# pdf(paste0(save_folder,"results/merged_PCA_jackstraw.pdf"))
# JackStrawPlot(tumor_cells, dims = 1:15)
# dev.off()
# 
# ElbowPlot(tumor_cells)
# 
# # perform clustering
# tumor_cells <- FindNeighbors(tumor_cells, dims = 1:10)
# tumor_cells <- FindClusters(tumor_cells, resolution = 0.5)
# 
# # visualize clusters on UMAP
# tumor_cells <- RunUMAP(tumor_cells, dims = 1:10)
# 
# pdf(paste0(save_folder,"results/merged_UMAP.pdf"))
# DimPlot(tumor_cells, reduction = "umap", group.by = "sample.ident")
# dev.off()
# 
# # save object 
# saveRDS(tumor_cells, file = "./tumor_cells_merge.RDS")

# load seurat object created in script above
tumor_cells <- readRDS("./tumor_cells_merge.RDS")

# now upload Olink DE data
olink_data <- "/Users/arnavmehta/Dropbox (Personal)/Manuscripts/Cytokine, ctDNA, exosome ML immuntherapy responders/Analysis/All data from Olink/Cohort 2 (Validation)"
lme.all <- read.csv(paste0(olink_data, "/Differential expression/LME_all_results.csv"))
lme.sig <- filter(lme.all, Threshold == "Significant")
lme.timepoint <- filter(lme.sig, term=="Timepoint")
lme.response <- filter(lme.sig, term=="grouped_response")
lme.interaction <- filter(lme.sig, term=="Timepoint:grouped_response")

geneNames_olink_time <- lme.timepoint$Assay
geneNames_olink_response <- lme.response$Assay
geneNames_olink_interaction <- lme.interaction$Assay

geneNames_olink_time<-recode(geneNames_olink_time, "LAIR-2"="LAIR2", "ADAM 8"="ADAM8", 
                             "MMP-2"="MMP2", "CLM-1"="CD300LF", "CNTN-1"="F3", "CRTAC1"="CEP68", "FcRL2"="FCRL2", 
                             "Gal-1"="LGALS1", "GCP5"="TUBGCP5", "GDNF"="ATF1", "IFN-gamma"="IFNG", 
                             "IL-15RA"="IL15RA", "IL-17A"="IL17A", "IL-18BP"="IL18BP", "IL-1RT2"="IL1R2", 
                             "IL-6RA"="IL6R", "IL2-RA"="IL2RA", "MAD homolog 5"="SMAD5", "MCP-3"="CCL7",
                             "N-CDase"="ASAH2", "NELL1"="NRP1", "PD-L1"="CD274", "SMOC2"="SMAP2", "ST2"="IL1RL1", 
                             "TLT-2"="TREML2", "TNF-R1"="TNFRSF1A", "TNF-R2"="TNFRSF1B", "TRAIL"="TNFSF10",
                             "TWEAK"="TNFSF12", "U-PAR"="PLAUR", "uPA"="PLAU")

geneNames_olink_response<-recode(geneNames_olink_response, "ADGRG1"="GPR56", 
                                 "CLM-1"="CD300LF", "CNTN-1"="F3", "EN-RAGE"="S100A12", 
                                 "GDF-15"="GDF15", "GCP5"="TUBGCP5", "IGFBP-2"="IGFBP2",
                                 "IL-18BP"="IL18BP", "IL-18R1"="IL18R1", "MAD homolog 5"="SMAD5",
                                 "MB"="PVALB", "MCP-4"="CCL13", "MK"="MDK", "N2DL-2"="ULBP2",
                                 "N-CDase"="ASAH2", "OPN"="SPP1", "PROC"="APC", "ST2"="IL1RL1",
                                 "TGF-alpha"="TGFA", "TNF-R1"="TNFRSF1A", "TNF-R2"="TNFRSF1B",
                                 "TNFB"="LTA", "TRAIL"="TNFSF10", "TWEAK"="TNFSF12", "U-PAR"="PLAUR")

geneNames_olink_interaction<-recode(geneNames_olink_interaction, "CLM-1"="CD300LF", 
                                    "BLM hydrolase"="BLMH", "ADGRG1"="GPR56", "CNTN-1"="F3", "CRTAC1"="CEP68", 
                                    "GM-CSF-R-alpha"="CSF2RA", "EN-RAGE"="S100A12", "FcRL2"="FCRL2", "Gal-1"="LGALS1",
                                    "GCP5"="TUBGCP5", "GDF-15"="GDF15", "IFN-gamma"="IFNG", "IGFBP-2"="IGFBP2",
                                    "IL-6RA"="IL6R", "IL2-RA"="IL2RA", "MAD homolog 5"="SMAD5", "MK"="MDK", "MMP-2"="MMP2",
                                    "N-CDase"="ASAH2", "SYND1"="SDC1", "SMOC2"="SMAP2", "ST2"="IL1RL1", "TLT-2"="TREML2", 
                                    "TRAIL"="TNFSF10", "U-PAR"="PLAUR")

# lets now visualize these time and response related genes in dotplots from scrnaseq data
assays <- intersect(unique(geneNames_olink_response[1:50]), rownames(tumor_cells))

pdf(paste0(save_folder,"./results/merged_response-significant_heatmap.pdf")) 
DoHeatmap(tumor_cells, features = assays, group.by = "sample.ident", size = 4)
dev.off()

pdf(paste0(save_folder,"./results/merged_timepoint-significant_dotplot.pdf"), 10, 10) 
DotPlot(tumor_cells, features = assays, group.by = "sample.ident") + RotatedAxis()
dev.off()

RidgePlot(tumor_cells, features = c("IL8"), group.by = "sample.ident")

FeaturePlot(tumor_cells, features = c("IL8"))

# END OF FILE
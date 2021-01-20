# analysis of melanoma immune cells obtained from Keren and Moshe
# samples were processed using the SmartSeq2 platform

rm(list=ls())

library(dplyr)
library(Seurat)
library(cowplot)
library(xlsx)
library(ggplot2)

data_folder <- "/Users/arnavmehta/Dropbox (Personal)/Manuscripts/Cytokine, ctDNA, exosome ML immuntherapy responders/Analysis/scRNAseq data/Nir Cell paper/Count matrix/"
setwd(data_folder)
save_folder <- "../results/"

# load counts and import into Seurat object
LOCATION <- "/Volumes/seq_hacohenlab1/david"
rnaseqFolder <- file.path(LOCATION, "rna-seq")
load("./2017-10-16_Melanoma_Keren_CELL.RData")
source(file.path(rnaseqFolder, "R_scripts/scRNAseq.R"))
moshe_useAllCells(TRUE, includeDNDP = FALSE)

immune_cells <- CreateSeuratObject(counts = plotDefault_exprDF,
                         project = project, meta.data = plotDefault_samplesDF)
Idents(immune_cells) <- if (USE_ALL_CELLS) "cluster" else "cluster_cd8"

# scale data
all.genes <- rownames(immune_cells)
immune_cells <- ScaleData(immune_cells, features = all.genes)

# find variable features
immune_cells <- FindVariableFeatures(immune_cells, selection.method = "vst", nfeatures = 2000)

# PCA
immune_cells <- RunPCA(immune_cells, features = VariableFeatures(object = immune_cells))
print(immune_cells[["pca"]], dims = 1:5, nfeatures = 5)

pdf(paste0(save_folder,"merged_PCA_loadings.pdf"))
VizDimLoadings(immune_cells, dims = 1:2, reduction = "pca")
dev.off()

pdf(paste0(save_folder,"merged_PCA.pdf"))
DimPlot(immune_cells, reduction = "pca", group.by = "cluster")
dev.off()

pdf(paste0(save_folder,"merged_PCA_loadings_heatmap.pdf"))
DimHeatmap(immune_cells, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

# assessing dimensionality
immune_cells <- JackStraw(immune_cells, num.replicate = 100)
immune_cells <- ScoreJackStraw(immune_cells, dims = 1:20)

pdf(paste0(save_folder,"results/merged_PCA_jackstraw.pdf"))
JackStrawPlot(immune_cells, dims = 1:15)
dev.off()

ElbowPlot(immune_cells)

# perform clustering
immune_cells <- FindNeighbors(immune_cells, dims = 1:10)
immune_cells <- FindClusters(immune_cells, resolution = 0.5)

# visualize clusters on UMAP
immune_cells <- RunUMAP(immune_cells, dims = 1:10)

pdf(paste0(save_folder,"merged_UMAP_cluster.pdf"))
DimPlot(immune_cells, reduction = "umap", group.by = "cluster")
dev.off()

# save object 
# saveRDS(immune_cells, file = "./immune_cells_merge.RDS")

# load seurat object created in script above
immune_cells <- readRDS("./immune_cells_merge.RDS")

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
assays <- intersect(unique(geneNames_olink_response[1:100]), rownames(immune_cells))

pdf(paste0(save_folder,"merged_response-significant_heatmap.pdf")) 
DoHeatmap(immune_cells, features = assays, group.by = "cluster", size = 4)
dev.off()

pdf(paste0(save_folder,"merged_timepoint-significant_dotplot.pdf")) 
DotPlot(immune_cells, features = assays, group.by = "cluster") + RotatedAxis()
dev.off()

RidgePlot(immune_cells, features = c("IL8"), group.by = "cluster")

pdf(paste0(save_folder,"UMAP_CLEC4D.pdf"))
FeaturePlot(immune_cells, features = c("CLEC4D"))
dev.off()

# END OF FILE
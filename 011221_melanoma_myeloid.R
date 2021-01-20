# analysis of melanoma myeloid cells obtained from Keren and Moshe
# samples were processed using the SmartSeq2 platform

rm(list=ls())

library(dplyr)
library(Seurat)
library(cowplot)
library(xlsx)
library(readxl)
library(ggplot2)
library(patchwork)
library(future)
library(SeuratWrappers)
library(magrittr)
library(monocle3)

data_folder <- "/Users/arnavmehta/Dropbox (Personal)/Manuscripts/Cytokine, ctDNA, exosome ML immuntherapy responders/Analysis/scRNAseq data/myeloid-subsets-sherry/"
setwd(data_folder)
save_folder <- "./results/"

# load Seurat objects for sample
myeloid.counts <- read.delim("./Myeloid_cells_TPM_more_data.txt", header = FALSE)
rownames(myeloid.counts) <- myeloid.counts[,1]
myeloid.counts <- myeloid.counts[,-1]
colnames(myeloid.counts) <- myeloid.counts[2,]
myeloid.counts <- myeloid.counts[-2,]
original.cluster <- data.frame(t(rbind(colnames(myeloid.counts), as.integer(myeloid.counts[1,]))))
colnames(original.cluster) <- c("cell", "cluster")
rownames(original.cluster) <- original.cluster$cell
myeloid.counts <- myeloid.counts[-c(1,2),]

myeloid <- CreateSeuratObject(counts=myeloid.counts, project="melanoma_myeloid",
                              min.cells=0, min.features=0)
myeloid <- AddMetaData(myeloid, original.cluster)

# find highly variable features
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)

# identify the 10 most highly variable genes
# plot variable features with and without labels
top10 <- head(VariableFeatures(myeloid), 10)
plot1 <- VariableFeaturePlot(myeloid)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(paste0(save_folder,"results/variable_features.pdf"))
plot1
dev.off()

# scale the data
all.genes <- rownames(myeloid)
myeloid <- ScaleData(myeloid, features = all.genes)

# PCA
myeloid <- RunPCA(myeloid, features = VariableFeatures(object = myeloid))
print(myeloid[["pca"]], dims = 1:5, nfeatures = 5)

pdf(paste0(save_folder,"merged_PCA_loadings.pdf"))
VizDimLoadings(myeloid, dims = 1:2, reduction = "pca")
dev.off()

pdf(paste0(save_folder,"merged_PCA.pdf"))
DimPlot(myeloid, reduction = "pca")
dev.off()

pdf(paste0(save_folder,"merged_PCA_loadings_heatmap.pdf"))
DimHeatmap(myeloid, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

# assessing dimensionality
myeloid <- JackStraw(myeloid, num.replicate = 100)
myeloid <- ScoreJackStraw(myeloid, dims = 1:20)

pdf(paste0(save_folder,"merged_PCA_jackstraw.pdf"))
JackStrawPlot(myeloid, dims = 1:15)
dev.off()

ElbowPlot(myeloid)

# perform clustering
myeloid <- FindNeighbors(myeloid, dims = 1:15)
myeloid <- FindClusters(myeloid, resolution = 0.3)

# visualize clusters on UMAP
myeloid <- RunUMAP(myeloid, dims = 1:10)

# this is from the immune cell object that includes all cells
# import response data
response <- data.frame(immune_cells$response)
sample_name <- immune_cells$sampleName
rownames(response) <- sample_name

myeloid <- AddMetaData(myeloid, response)

pdf(paste0(save_folder,"merged_UMAP.pdf"))
DimPlot(myeloid, reduction = "umap", group.by ="cluster")
DimPlot(myeloid, reduction = "umap", group.by ="immune_cells.response")
DimPlot(myeloid, reduction = "umap")
dev.off()

# save object 
# saveRDS(myeloid, file = "./myeloid_cells.RDS")
myeloid <- readRDS("./myeloid_cells.RDS")

# lets find marker genes for each of the identified clusters
myeloid.markers <- FindAllMarkers(myeloid, only.pos = TRUE, min.pct = 0.25, 
                                             logfc.threshold = 0.25)
myeloid.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(myeloid.markers, "./results/myeloid-markers.csv")

top5 <- myeloid.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf("./marker-gene-heatmap.pdf")
DoHeatmap(myeloid, features = top5$gene) + NoLegend()
dev.off()

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
assays <- intersect(unique(geneNames_olink_response[1:50]), rownames(myeloid))

pdf(paste0(save_folder,"merged_response-significant_heatmap.pdf")) 
DoHeatmap(myeloid, features = assays, size = 4)
dev.off()

pdf(paste0(save_folder,"merged_timepoint-significant_dotplot.pdf")) 
DotPlot(myeloid, features = assays) + RotatedAxis()
dev.off()

pdf(paste0(save_folder,"IL8-ridge-plot.pdf")) 
RidgePlot(myeloid, features = c("IL8"))
dev.off()

pdf(paste0(save_folder,"OSM-umap-plot.pdf")) 
FeaturePlot(myeloid, features = c("OSM"))
dev.off()

# now consider proteins higher in non-responders
response.posthoc <- read_excel("/Users/arnavmehta/Dropbox (Personal)/Manuscripts/Cytokine, ctDNA, exosome ML immuntherapy responders/Analysis/All data from Olink/Cohort 2 (Validation)/Differential expression/021720_final-LME-groups/MGH_Mehta-LME.xlsx", sheet=3)
interact.posthoc <- read_excel("/Users/arnavmehta/Dropbox (Personal)/Manuscripts/Cytokine, ctDNA, exosome ML immuntherapy responders/Analysis/All data from Olink/Cohort 2 (Validation)/Differential expression/021720_final-LME-groups/MGH_Mehta-LME.xlsx", sheet=4)
response.posthoc.upR <- filter(response.posthoc, estimate < 0)
response.posthoc.upNR <- filter(response.posthoc, estimate > 0)

assays <- intersect(unique(response.posthoc.upNR$Assay[1:100]), rownames(myeloid))

pdf(paste0(save_folder,"merged_upNR_heatmap.pdf")) 
DoHeatmap(myeloid, features = assays, size = 4)
dev.off()

pdf(paste0(save_folder,"merged_upNR_dotplot.pdf")) 
DotPlot(myeloid, features = assays) + RotatedAxis()
dev.off()

# lets now look at these cells in pseudotime
cds <- as.cell_data_set(myeloid)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

myeloid.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(myeloid.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

cluster2.cells <- WhichCells(object = myeloid, idents = 2)
cds <- order_cells(cds, root_cells = cluster2.cells)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

# END OF FILE
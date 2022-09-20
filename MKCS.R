# Packages for Spatial Transcriptomics
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
BiocManager::install('limma')
library(limma)
library(dplyr)
library(hdf5r)
library(rhdf5)
install.packages("BiocManager")
library(Matrix)
library(ggpubr)
library(tidyverse)
library(sctransform)

# Specifying work directory
root_dir<- "C:/Users/Owner/Downloads/MKCS"
setwd(root_dir)

# File names
spatial_folder_name<-"spatial.tar.gz"
untar(spatial_folder_name)
untar(spatial_folder_name,list=TRUE)
list.files(spatial_folder_name) # files in spatial subdirectory
h5_mat_name<-"filtered_feature_bc_matrix.h5"
h5ls("filtered_feature_bc_matrix.h5") # open data in h5 file
# Load 10x genomics into seurat object
MKCS_data<-Seurat::Load10X_Spatial(
  data.dir=root_dir,
  filename=h5_mat_name,
  assay="Spatial",
  slice="slice1",
  filter.matrix=TRUE,
  to.upper=FALSE)
MKCS_data
# Exploring Seurat Object
# a) Dimension of active assay
dim(x=MKCS_data)
nrow(x=MKCS_data)
ncol(x=MKCS_data)
# b) Features and sample name
head(x=rownames(MKCS_data),n=5)
tail(x=colnames(MKCS_data),n=5)
class(MKCS_data[[]])
colnames(MKCS_data[[]])
head(MKCS_data)
head(MKCS_data@meta.data)
MKCS_data$nCount_Spatial[1:3]
sum(MKCS_data$nFeature_Spatial ==  colSums(MKCS_data@assays$Spatial@counts > 0))
sum(MKCS_data$nCount_Spatial ==  colSums(MKCS_data@assays$Spatial@counts))

names(x=MKCS_data)
MKCS_data[["Spatial"]]
MKCS_data[["slice1"]]
MKCS_data@assays$Spatial@counts[]
MKCS_data[["Spatial"]]@meta.features
head(MKCS_data[["Spatial"]][[]])

MKCS_data@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(MKCS_data@images[["slice1"]]@coordinates[["tissue"]])
MKCS_data@images[["slice1"]]@coordinates[["row"]] <- as.integer(MKCS_data@images[["slice1"]]@coordinates[["row"]])
MKCS_data@images[["slice1"]]@coordinates[["col"]] <- as.integer(MKCS_data@images[["slice1"]]@coordinates[["col"]])
MKCS_data@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(MKCS_data@images[["slice1"]]@coordinates[["imagerow"]])
MKCS_data@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(MKCS_data@images[["slice1"]]@coordinates[["imagecol"]])
# Quality Control 1
plotA <- VlnPlot(MKCS_data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plotA

plotB <- SpatialFeaturePlot(MKCS_data,features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plotA, plotB)
# Quality Control 2
MKCS_data[["percent.mt"]] <- PercentageFeatureSet(MKCS_data,pattern = "^mt-")
VlnPlot(
  MKCS_data, features = c("nFeature_Spatial", "nCount_Spatial","percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
# Jointly (rather than separately) consider the QC metrics when filtering
plot1 <- FeatureScatter(
  MKCS_data, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(
  MKCS_data, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
  NoLegend()
plot1 + plot2
# Filtering if applicable
MKCS_dataA1 <- subset(MKCS_data, subset = nFeature_Spatial < 9500 & nFeature_Spatial > 2000 & nCount_Spatial < 100000)
print(paste("Filter out", ncol(MKCS_data) - ncol(MKCS_dataA1), 
            "samples because of the outlier QC metrics, with", ncol(MKCS_dataA1),
            "samples left."))

# Normalization
SpatialFeaturePlot(
  MKCS_dataA1,features = c("nFeature_Spatial", "nCount_Spatial")) &
  theme(legend.position = "bottom")
MKCS_dataA1norm <- SCTransform(MKCS_dataA1, assay = "Spatial", verbose = FALSE)
names(MKCS_dataA1norm)
#dim(MKCS_data@assays$SCT@counts)
#dim(MKCS_data@assays$SCT@data)
#dim(MKCS_data@assays$SCT@scale.data)
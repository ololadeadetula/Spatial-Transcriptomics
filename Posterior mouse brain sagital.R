---
  title: "Analysis, visualization, and integration of spatial datasets with Seurat"
output:
  html_document:
  theme: united
df_print: kable
pdf_document: default
date: 'Compiled: `r Sys.Date()`'
---
  
  ```{r setup, include=FALSE}
all_times <- list()  # store the time for each chunk
knitr::knit_hooks$set(time_it = local({
  now <- NULL
  function(before, options) {
    if (before) {
      now <<- Sys.time()
    } else {
      res <- difftime(Sys.time(), now, units = "secs")
      all_times[[options$label]] <<- res
    }
  }
}))
knitr::opts_chunk$set(
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 95),
  message = FALSE,
  warning = FALSE,
  fig.width = 10,
  time_it = TRUE
)
```

# Overview

This tutorial demonstrates how to use Seurat (>=3.2) to analyze spatially-resolved RNA-seq data. While the analytical pipelines are similar to the Seurat workflow for [single-cell RNA-seq analysis](pbmc3k_tutorial.html), we introduce updated interaction and visualization tools, with a particular emphasis on the integration of spatial and molecular information. This tutorial will cover the following tasks, which we believe will be common for many spatial analyses:
  
  * Normalization 
* Dimensional reduction and clustering
* Detecting spatially-variable features
* Interactive visualization
* Integration with single-cell RNA-seq data
* Working with multiple slices

For our first vignette, we analyze a dataset generated with the [Visium technology](https://www.10xgenomics.com/spatial-transcriptomics/) from 10x Genomics. We will be extending Seurat to work with additional data types in the near-future, including [SLIDE-Seq](https://science.sciencemag.org/content/363/6434/1463), [STARmap](https://science.sciencemag.org/content/361/6400/eaat5691), and [MERFISH](https://science.sciencemag.org/content/362/6416/eaau5324).

First, we load Seurat and the other packages necessary for this vignette.

```{r install}
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
```

```{r libraries.for.rmd, echo = FALSE}
library("htmltools")
library("vembedr")
```


# 10x Visium

## Dataset

Here, we will be using a recently released dataset of sagital mouse brain slices generated using the Visium v1 chemistry. There are two serial anterior sections, and two (matched) serial posterior sections. 

You can download the data [here](https://support.10xgenomics.com/spatial-gene-expression/datasets), and load it into Seurat using the `Load10X_Spatial()` function. This reads in the output of the [spaceranger](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger) pipeline, and returns a Seurat object that contains both the spot-level expression data along with the associated image of the tissue slice. You can also use our [SeuratData package](https://github.com/satijalab/seurat-data) for easy data access, as demonstrated below. After installing the dataset, you can type `?stxBrain` to learn more.

```{r data.install, eval = FALSE}
InstallData("stxBrain")
```

```{r data}
pbrain <- LoadData('stxBrain', type = 'posterior1')
```
pbrain

```{r qc, fig.height=5}
plot1 <- VlnPlot(pbrain, features = 'nCount_Spatial', pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(pbrain, features = 'nCount_Spatial') + theme(legend.position = "right")
wrap_plots(plot1, plot2)
```

```{r preprocess}
pbrain <- SCTransform(pbrain, assay = "Spatial", verbose = FALSE)
```
## Gene expression visualization 
```{r featureplot}
SpatialFeaturePlot(pbrain, features = c("Hpca", "Ttr"))

```{r fpe1}
p1 <- SpatialFeaturePlot(pbrain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(pbrain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2
```
# Dimensionality reduction, clustering, and visulization
```{r dim.cluster}
pbrain <- RunPCA(pbrain, assay = "SCT", verbose = FALSE)
pbrain <- FindNeighbors(pbrain, reduction = "pca", dims = 1:30)
pbrain <- FindClusters(pbrain, verbose = FALSE)
pbrain <- RunUMAP(pbrain, reduction = "pca", dims = 1:30)
```
```{r dim.plots,fig.height=5}
p1 <- DimPlot(pbrain, reduction = "umap", label = TRUE) 
p2 <- SpatialDimPlot(pbrain, label = TRUE, label.size = 3) 
p1 + p2
```
```{r facetdim}
SpatialDimPlot(pbrain, cells.highlight = CellsByIdentities(object = pbrain,idents = c(2, 1, 4, 3, 5, 8)), facet.highlight = TRUE, ncol = 3)

#Interactive plotting
```{r ispatialdimplot, eval = FALSE}
SpatialDimPlot(pbrain, interactive = TRUE)
``
```{r ispatialfeatureplot, eval = FALSE}
SpatialFeaturePlot(pbrain, features = "Ttr", interactive = TRUE)
```
```{r linkedplot, eval=FALSE}
LinkedDimPlot(pbrain)
```
#Identification of Spatially Variable Features
```{r de, fig.height = 4}
de_markers <- FindMarkers(pbrain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = pbrain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
```
```{r spatial.vf}
pbrain <- FindSpatiallyVariableFeatures(pbrain, assay = 'SCT', features = VariableFeatures(pbrain)[1:1000], selection.method = 'markvariogram')
```
```{r spatial.vf.plot, fig.height=8}
top.features <- head(SpatiallyVariableFeatures(pbrain, selection.method = 'markvariogram'),6)
SpatialFeaturePlot(pbrain, features = top.features, ncol = 3, alpha = c(0.1, 1))
```
#Subset out anatomical regions
```{r subset1}
cortex <- subset(pbrain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400 | image_imagecol < 150))
cortex <- subset(cortex, posterior1_imagerow > 400 | posterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, posterior1_imagerow > 275 & posterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, posterior1_imagerow > 250 & posterior1_imagecol > 440, invert = TRUE)
```

```{r subset1.plot, fig.height = 4}
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
```
#Integration with single-cell data

```{r sc.data}
allen_reference <- readRDS("../data/allen_cortex.rds")
```

```{r sc.data2}
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k cells
# this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
```

```{r sc.data3, fig.width=8, fig.align="center"}
# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = 'Spatial', verbose = FALSE) %>% RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = 'subclass', label = TRUE)
```
```{r sc.data5}
anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE, weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
```
  
  ```{r sc.data7}
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
```
```{r sc.data8, fig.height = 10}
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "markvariogram", features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
```

```{r sc.data9,fig.height=20,fig.width=10}
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT", "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))
```
#Working the multiple slices in Seurat
```{r brain2data}
pbrain2 <- LoadData('stxBrain', type = 'anterior1')
pbrain2 <- SCTransform(pbrain2, assay = "Spatial", verbose = FALSE)
```
```{r merge}
pbrain.merge <- merge(pbrain, pbrain2)
```
```{r joint.analysis}
DefaultAssay(pbrain.merge) <- "SCT"
VariableFeatures(pbrain.merge) <- c(VariableFeatures(pbrain), VariableFeatures(pbrain2))
pbrain.merge <- RunPCA(pbrain.merge, verbose = FALSE)
pbrain.merge <- FindNeighbors(pbrain.merge, dims = 1:30)
pbrain.merge <- FindClusters(pbrain.merge, verbose = FALSE)
pbrain.merge <- RunUMAP(pbrain.merge, dims = 1:30)
```
```{r joint.viz, fig.height = 4}
DimPlot(pbrain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
```
```{r joint.viz2}
SpatialDimPlot(pbrain.merge)
```

```{r joint.viz3, fig.height = 10}
SpatialFeaturePlot(pbrain.merge, features = c('Hpca', 'Plp1'))
```
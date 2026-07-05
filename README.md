
# scBaptism

<!-- badges: start -->
<!-- badges: end -->


## What is scBaptism and what can it do?

scBaptism is a tool that enables the user to use multiple existing single cell annotation methods at once (SCINA, CIA, SingleR, Seurat, scPred, scClassify, clustifyr, scmap, CelliD marker based, CelliD reference based), as well as methods to integrate the obtained annotations and options for visualizations.

Single Cell transcriptomic data can be used to identify different cells inside a sample, the so called annotation. Many different tools are available to automatically perform this process. For researchers this poses the challenge of being familiar with many different tools and how to use them at once. scBaptism offers the possibility to run 10 different annotation methods from 9 existing R packages at once and with the same input data. This makes the process of annotating faster and more accessible. Additionally the resulting annotations can be collapsed into one annotation by majority voting and through a meta cell approach.

## Install the package

# You can install scBaptism directly from github:

``` r
#install scBaptism
remotes::install_github("alinajenn/scBaptism")
```

## Load the package

``` r
library("scBaptism")
```

## Example

This is a basic example were we use scBaptism to perform annotations with SingleR (reference-based) and CIA (marker-based), using a dataset from the iUSEiSEE package

``` r
# load neccessary packages
library(iUSEiSEE)
library(Seurat)
library(dplyr)

# load example SCE from iUSEiSEE package (also used as reference for this example)

sce_annotated <- readRDS(
 file = system.file("datasets", 
                    "sce_pbmc3k.RDS", 
                    package = "iUSEiSEE"))

# get the markers list using Seurat

myseu <- Seurat::as.Seurat(sce_annotated)
myseu <- Seurat::ScaleData(myseu)
Seurat::Idents(myseu) <- "labels_main"
seu_all_markers <- Seurat::FindAllMarkers(myseu, 
                                          test.use = "wilcox", 
                                          only.pos = TRUE,
                                          min.pct = 0.25, 
                                          logfc.threshold = 0.25)
top_k_markers <- 50
markers_lists <- seu_all_markers %>%
 dplyr::group_by(cluster) %>%
 dplyr::top_n(n = top_k_markers, wt = avg_log2FC)
markers_lists <- split(markers_lists$gene, 
                      markers_lists$cluster)

#run scBaptism to annotate with SingleR and CIA

anno_result <- run_scBaptism(sce_query = sce_annotated,
                             anno_methods = c("SingleR", "CIA"),
                             markers_list = markers_lists,
                             reference = sce_annotated,
                             ref_labs = "labels_main")

```


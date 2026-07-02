library("iUSEiSEE")
library("Seurat")
library("dplyr")

# load example SCE from iUSEiSEE package (also used as reference)

sce_annotated <- readRDS(
 file = system.file("datasets", "sce_pbmc3k.RDS", package = "iUSEiSEE"))

# get the markers list using Seurat
suppressWarnings({
  myseu <- Seurat::as.Seurat(sce_annotated)
})
myseu <- Seurat::ScaleData(myseu)
Seurat::Idents(myseu) <- "labels_main"

seu_all_markers <- Seurat::FindAllMarkers(myseu, test.use = "wilcox", only.pos = TRUE,
                                          min.pct = 0.25, logfc.threshold = 0.25)

top_k_markers <- 50

markers_lists <- seu_all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = top_k_markers, wt = avg_log2FC)
markers_lists <- split(markers_lists$gene, markers_lists$cluster)


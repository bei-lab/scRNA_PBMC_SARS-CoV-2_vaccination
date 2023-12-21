packages <- list("Seurat", "Matrix", "stringr", "stringi", "ggplot2", "plyr", "ggthemes", "cowplot", "data.table", "SeuratData", "SeuratDisk", "future", "ggsci", "harmony", "ggpubr")
lapply(packages, library, character.only = TRUE)

options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multiprocess", workers = 20)
plan()

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]

###------------------------------------------------------batch effect correction-----------------------------------------------------------###
subset_cells <- get(load("Seurat.object.RData"))
subset_cells <- NormalizeData(subset_cells, verbose = TRUE)
subset_cells <- FindVariableFeatures(subset_cells, nfeatures = 2000, selection.method = 'vst', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf))
subset_cells <- ScaleData(subset_cells, verbose = TRUE,  vars.to.regress = c("nCount_RNA", "percent.mito" ) )
subset_cells <- RunPCA(subset_cells, features = setdiff(subset_cells@assays$RNA@var.features, row.names(subset_cells) %>% grep(pattern = "^MT-|^RPL|^RPS|^IGK|^IGV|^IGL|^IGH", v = T)) , npcs = 50, verbose = TRUE)

### correct batch effect wiht harmony
subset_cells <- RunHarmony(subset_cells, group.by.vars = c("Samples_new", "Batch_new"))

### run UMAP and find clusters
subset_cells <- RunUMAP(subset_cells, reduction = "harmony", dims = 1:30)
subset_cells <- FindNeighbors(subset_cells, reduction = "harmony", dims = 1:30)
subset_cells <- FindClusters(subset_cells, resolution = 1)

png("All_cell_major_cluster_by_cluster_with_legend.png",
    width = 15, height = 15, units = "in", res = 300)
p5 <- DimPlot(object = subset_cells[, WhichCells(subset_cells, downsample = 30000)], reduction = 'umap', label = F, cols = color_used, pt.size = 0.0001, label.size = 5, raster = F) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5))) +
  NoLegend()
print(p5)
dev.off()

###---------------------DEG calculation-------------------###
plan("multiprocess", workers = 1)

result <- mclapply(as.numeric(levels(subset_cells@active.ident)),
                   FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL, max.cells.per.ident = 500)},
                   mc.cores = 36)
RESULT <- result

roundN <- 1
while(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
  if(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
    recalculate_clusters <- which(mapply(length, result, SIMPLIFY = TRUE)!=5)-1
    print(recalculate_clusters)
    result1 <- mclapply(recalculate_clusters,
                        FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL, max.cells.per.ident = 500)},
                        mc.cores = 35)
  }
  print(roundN + 1)
  for(i in 1:length(recalculate_clusters)){
    result[c(recalculate_clusters+1)[i]] <- result1[i]
  }
}

all_markers <- do.call(rbind, result)
all_markers$gene <- unlist(mapply(rownames, result))
all_markers$cluster <- rep(levels(subset_cells@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
subset_cells.markers <- all_markers

subset_cells.markers %>% TOP_N(50, pct.1 = 0.2) -> top50 ## TOP_N: please ref to AHCA project
subset_cells.markers <- subset_cells.markers %>% TOP_N(5000)

pdf("Major_cluster_Vlnplot_of_markers.pdf", height = 15, width = 15)
VlnPlot(subset_cells,
        features = c("CD3E", "CD3D" , "CD3G", "CCR7", "SELL", "TCF7", "LEF1", "GNLY", "PRF1", "GZMA", "GZMB", "FOXP3", "TRDV2", "TRGV9",
                     "MKI67", "SLC4A10", "KLRF1", "NKG7", "MS4A1", "CD79A", "TCL1A", "CD27",
                     "XBP1", "CD14", "FCGR3A", "CLEC9A", "CLEC10A", "IL3RA", "GATA2", "CYTL1", "PPBP"),
        pt.size = 0,
        log = F,
        ncol = 1,
        stack = TRUE,
        flip = TRUE,
        cols = color_used) #+
  #theme(legend.position = "none") + ggtitle("identities on y-axis")
dev.off()

write.table(top50,
            file = "All_cell_major_cluster_top50_DEGs.csv",
            sep = ",",
            row.names = T,
            quote = F)

write.table(subset_cells.markers,
            file = "All_cell_major_cluster_all_DEGs.csv",
            sep = ",",
            row.names = T,
            quote = F)

png(paste0("featureplot_cells.png"), height = 30, width = 30, units = "in", res = 600)
FeaturePlot(subset_cells[, WhichCells(subset_cells, downsample = 30000)],
            features =  c("CD3E", "CD3D" , "CD3G", "CCR7", "SELL", "TCF7", "LEF1", "GNLY", "PRF1", "GZMA", "GZMB", "FOXP3", "TRDV2", "TRGV9",
                          "MKI67", "SLC4A10", "KLRF1", "NKG7", "MS4A1", "CD79A", "TCL1A", "CD27",
                          "XBP1", "CD14", "FCGR3A", "CLEC9A", "CLEC10A", "IL3RA", "GATA2", "CYTL1", "PPBP"), 
            pt.size = 0.0001,
            ncol = 6,
            # cols = c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12))[],
            reduction = "umap",
            raster = F)  
dev.off()

### Below is the demo script
library(Seurat)
library(Matrix)
library(stringr)
library(gridExtra)
library(ggplot2)
library(stringi)
library(plyr)
library(ggthemes)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(reshape2)
library(scales)
library(rlang)
library(future)
library(parallel)
library(ggrepel)
library(ggsci)
library(ape)
library(dplyr)
library(sigminer)
library(ggpubr)
library(rstatix)

options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multiprocess", workers = 20)
plan()

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]

###------------------------------load data------------------------------###
load("/data4/heshuai/RAW_data/1-SingleCell/21.COV19_vaccine/1.House_in_dat/3.Analysis/Reanalysis_COVID19_Virus/Annotation.all_cells.meta.data_788534_COVID19_virus.RData")
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features = VariableFeatures(subset_cells))

####---------------------------------------------for COVID19---------------------------------------------###
subset_cells <- subset_cells[1:50, ]
Idents(subset_cells) <- subset_cells$Final_Major_cluster
meta.data <- read.table("/data4/heshuai/RAW_data/1-SingleCell/21.COV19_vaccine/1.House_in_dat/3.Analysis/All/meta.data_hallmark.genesets.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
selected_data <- meta.data[subset_cells %>% colnames(), ]

subset_cells@assays$RNA@data <- selected_data[, 30:79] %>% t

dat_names <- data.frame(orig_name = c(subset_cells@assays$RNA@counts %>% row.names)[1:50], fixed_names = subset_cells@assays$RNA@data %>% row.names %>% gsub(pattern = "1$|HALLMARK_", replacement = ""))
row.names(subset_cells@assays$RNA@data) <- dat_names$orig_name

result <- lapply(X = c("mild/moderate-convalescence-", "mild/moderate-progression-", "severe/critical-convalescence-", "severe/critical-progression-"), FUN = function(y){
  Vaccine_control <- lapply(subset_cells@active.ident %>% levels %>% grep(pattern = "HSC|Pro |cDC1", invert = T, value = T), FUN = function(x){
    lis <-  FindMarkers(subset_cells, ident.2 = "control-control-", ident.1 = y,  group.by = 'Condition', subset.ident = x, logfc.threshold = 0)
    lis$cluster <- x
    lis$gene <- row.names(lis)
    lis$group <- y
    print(c(x, y))
    return(lis)
  })
  results <- do.call(rbind, Vaccine_control)
}
)

COVID_hallmark <- do.call(rbind, result) %>% as.data.frame()
COVID_hallmark$gene <- mapvalues(COVID_hallmark$gene, from = dat_names$orig_name %>% as.character, to = dat_names$fixed_names %>% as.character())

COVID_hallmark_up <- subset(COVID_hallmark, p_val_adj < 0.05 & avg_logFC > 0)
COVID_hallmark_up <- COVID_hallmark_up %>% subset(grepl(group, pattern = "progression"))
COVID_hallmark_up$gene %>% table %>% sort(decreasing = T)
selected_TFs_up <- names(COVID_hallmark_up$gene %>% table %>% sort(decreasing = T))[c(COVID_hallmark_up$gene %>% table %>% sort(decreasing = T)) > 10]
selected_TFs_up_table <- c(COVID_hallmark_up$gene %>% table %>% sort(decreasing = T))[c(COVID_hallmark_up$gene %>% table %>% sort(decreasing = T)) > 10]


COVID_hallmark_down <- subset(COVID_hallmark, p_val_adj < 0.05 & avg_logFC < 0)
COVID_hallmark_down <- COVID_hallmark_down %>% subset(grepl(group, pattern = "progression"))
COVID_hallmark_down$gene %>% table %>% sort(decreasing = T)
selected_TFs_down <- names(COVID_hallmark_down$gene %>% table %>% sort(decreasing = T))[c(COVID_hallmark_down$gene %>% table %>% sort(decreasing = T)) > 10]
selected_TFs_down_table <- c(COVID_hallmark_down$gene %>% table %>% sort(decreasing = T))[c(COVID_hallmark_down$gene %>% table %>% sort(decreasing = T)) > 10]

COVID_hallmark %>% write.table("COVID_hallmark_DEGs.txt", col.names = T, row.names = T, sep = "\t")

COVID_hallmark$Groups <- paste0(COVID_hallmark$cluster, "_", COVID_hallmark$group)
selected_dat <- COVID_hallmark %>% select(p_val_adj, gene, Groups)
selected_dat$p_val_adj <- ifelse(selected_dat$p_val_adj < 0.001, "***",
                                 ifelse(selected_dat$p_val_adj > 0.001 & selected_dat$p_val_adj < 0.01, "**",
                                        ifelse(selected_dat$p_val_adj > 0.01 & selected_dat$p_val_adj <= 0.05, "*", "")))

control_group <- selected_dat %>% subset(grepl(Groups, pattern = "severe/critical-progression-"))
control_group$Groups <- control_group$Groups %>% gsub(pattern = "severe/critical-progression-", replacement = "control-control-")
selected_dat <- rbind(control_group, selected_dat) %>% as.data.frame()
selected_dat$p_val_adj[grepl(selected_dat$Groups, pattern = "control-control-")] <- ""
dcast_mat <- dcast(selected_dat, gene ~ Groups, value.var = "p_val_adj")

# save(subset_cells, file = "data_for_TF_analysis.RData")## save seurat data for TF anlysis, especifically for TF DEGs calculation.
###--------------------------------------normalized module score--------------------------------###
cellInfo <- subset_cells@meta.data[, c("Final_Major_cluster", "Sampleinfo", "Groups", "Condition")]
cellInfo$seuratCluster <- paste0(cellInfo$Final_Major_cluster, "_", cellInfo$Condition)

cellInfo <- cellInfo %>% subset(!seuratCluster %>% grepl(pattern = "Pro |HSC|Virus|cDC1"))

subset_cells <- subset_cells[, cellInfo %>% row.names()]

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(subset_cells@assays$RNA@data[, cells])) %>% as.data.frame()

regulonActivity_byCellType_Scaled <- dplyr::select(regulonActivity_byCellType, !contains(c("mRNA", "Virus", "HP", "HSC", "cDC1")))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[, COVID_cluster_order]
row.names(regulonActivity_byCellType_Scaled) <- mapvalues(regulonActivity_byCellType_Scaled %>% row.names(), from = dat_names$orig_name %>% as.character, to = dat_names$fixed_names %>% as.character())

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType_Scaled[, ]), center = T, scale=T))
col <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# anto <- data.frame(cluster = colnames(regulonActivity_byCellType_Scaled))
# row.names(anto) <- colnames(regulonActivity_byCellType_Scaled)
# 
# names(cluster) <- anto$cluster
# ancols <- list(cluster = cluster)

set.seed(123)
mat <-  regulonActivity_byCellType_Scaled[c(selected_TFs_up, selected_TFs_down), COVID_cluster_order]
row.names(mat) <- mat %>% row.names %>% gsub(pattern = "HALLMARK_", replacement = "")

Group <- regulonActivity_byCellType_Scaled %>% colnames() %>% as.character

cols <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8][1:43]
cols <- data.frame(X1 = cols[seq(1, 42, 3)],
                   x2 = alpha(cols[seq(1, 42, 3) + 1], 0.1),
                   x3 = alpha(cols[seq(1, 42, 3) + 1], 0.5),
                   x4 = alpha(cols[seq(1, 42, 3) + 2], 0.1),
                   X5 = alpha(cols[seq(1, 42, 3) + 2], 0.5))
cols <- cols%>%as.matrix%>%t%>%c

names(cols) <- Group
column_ha <-  HeatmapAnnotation(Group = regulonActivity_byCellType_Scaled %>% colnames() %>% as.character,
                                col = list(Group = cols))

row_ha <- rowAnnotation(Up_down = c(rep("UP", time = length(selected_TFs_up_table)), rep("Down", time = length(selected_TFs_down))),
                        Freq = anno_barplot(c(selected_TFs_up_table, selected_TFs_down_table))
)

dcast_mat[, "gene"] <- dcast_mat[, "gene"] %>% gsub(pattern = "HALLMARK_", replacement = "")
row.names(dcast_mat) <- dcast_mat[, "gene"]
dcast_mat[, "gene"] <- NULL
selected_padj_mat <- dcast_mat[mat %>% row.names, mat %>% colnames]

pdf("COVID19_hallmark.pheatmap_large.pdf", height = 14, width = 25)
Heatmap(mat,
        name = "mat",
        col = colorRampPalette(c("#181495", "blue", "white", "red", "#930107"))(100),
        top_annotation = column_ha,
        right_annotation = row_ha,
        cluster_rows = F,
        cluster_columns = F,
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", selected_padj_mat[i, j]), x, y)
        }
)
dev.off()

##-----------------------------------------order of cluster-----------------------------------###
COVID_cluster_order = c(
  'Naive T_control-control-',
  'Naive T_mild/moderate-progression-',
  'Naive T_mild/moderate-convalescence-',
  'Naive T_severe/critical-progression-',
  'Naive T_severe/critical-convalescence-',
  'Active T_control-control-',
  'Active T_mild/moderate-progression-',
  'Active T_mild/moderate-convalescence-',
  'Active T_severe/critical-progression-',
  'Active T_severe/critical-convalescence-',
  'Treg_control-control-',
  'Treg_mild/moderate-progression-',
  'Treg_mild/moderate-convalescence-',
  'Treg_severe/critical-progression-',
  'Treg_severe/critical-convalescence-',
  'MAIT_control-control-',
  'MAIT_mild/moderate-progression-',
  'MAIT_mild/moderate-convalescence-',
  'MAIT_severe/critical-progression-',
  'MAIT_severe/critical-convalescence-',
  'GD T_control-control-',
  'GD T_mild/moderate-progression-',
  'GD T_mild/moderate-convalescence-',
  'GD T_severe/critical-progression-',
  'GD T_severe/critical-convalescence-',
  'NK_control-control-',
  'NK_mild/moderate-progression-',
  'NK_mild/moderate-convalescence-',
  'NK_severe/critical-progression-',
  'NK_severe/critical-convalescence-',
  'Naive B_control-control-',
  'Naive B_mild/moderate-progression-',
  'Naive B_mild/moderate-convalescence-',
  'Naive B_severe/critical-progression-',
  'Naive B_severe/critical-convalescence-',
  'Memony B_control-control-',
  'Memony B_mild/moderate-progression-',
  'Memony B_mild/moderate-convalescence-',
  'Memony B_severe/critical-progression-',
  'Memony B_severe/critical-convalescence-',
  'Plasma_control-control-',
  'Plasma_mild/moderate-progression-',
  'Plasma_mild/moderate-convalescence-',
  'Plasma_severe/critical-progression-',
  'Plasma_severe/critical-convalescence-',
  'CD14 monocyte_control-control-',
  'CD14 monocyte_mild/moderate-progression-',
  'CD14 monocyte_mild/moderate-convalescence-',
  'CD14 monocyte_severe/critical-progression-',
  'CD14 monocyte_severe/critical-convalescence-',
  'CD16 monocyte_control-control-',
  'CD16 monocyte_mild/moderate-progression-',
  'CD16 monocyte_mild/moderate-convalescence-',
  'CD16 monocyte_severe/critical-progression-',
  'CD16 monocyte_severe/critical-convalescence-',
  'cDC2_control-control-',
  'cDC2_mild/moderate-progression-',
  'cDC2_mild/moderate-convalescence-',
  'cDC2_severe/critical-progression-',
  'cDC2_severe/critical-convalescence-',
  'pDC_control-control-',
  'pDC_mild/moderate-progression-',
  'pDC_mild/moderate-convalescence-',
  'pDC_severe/critical-progression-',
  'pDC_severe/critical-convalescence-',
  'Platelet_control-control-',
  'Platelet_mild/moderate-progression-',
  'Platelet_mild/moderate-convalescence-',
  'Platelet_severe/critical-progression-',
  'Platelet_severe/critical-convalescence-'
)


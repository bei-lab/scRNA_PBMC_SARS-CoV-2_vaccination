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
library(EnhancedVolcano)

options(future.globals.maxSize = 400*1000 * 1024^2)
plan("multiprocess", workers = 40)
plan()

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10), pal_jama()(7))[-8]
condition_colors <- c("#E64B35FF", "#00A087FF", "#ADE2D0FF", "#3C5488FF", "#8491B4FF", "#D6D6CEFF", "#B09C85FF", "#7E6148FF")[c(6:8, 1:5)]

###--------------------------------calculation of cell proportion----------------------------###
###--------------merge all the singlet cell identifed by hashtag----------------#####
orig.data <- read.table("meta.data.TNK.cell.txt", header = T, row.names = 1, stringsAsFactors = F, sep = "\t")
meta.data <- orig.data

meta.data$Condition <- meta.data$Condition %>%
  as.character %>%
  gsub(pattern = "control-control-", replacement = "HD") %>% 
  gsub(pattern = "mild/moderate-progression-", replacement = "MP") %>% 
  gsub(pattern = "mild/moderate-convalescence-", replacement = "MC") %>% 
  gsub(pattern = "severe/critical-progression-", replacement = "SP") %>% 
  gsub(pattern = "severe/critical-convalescence-", replacement = "SC") %>% 
  gsub(pattern = "Virus-A0", replacement = "IVV-D0") %>% 
  gsub(pattern = "Virus-A28", replacement = "IVV-D28") %>% 
  gsub(pattern = "Virus-A42", replacement = "IVV-D42") %>%
  gsub(pattern = "A42-A28", replacement = "D42-D28")

df1 <- meta.data %>% select(Patient_name, Sex, Vaccine_type, IgG, Batch_new, Samples_new, Major_celltype, Sampleinfo, Groups, Condition)

df1$Condition <- factor(df1$Condition %>% as.character(), levels = c("IVV-D0",
                                                                     "IVV-D28",
                                                                     "IVV-D42",
                                                                     "HD",
                                                                     "MP",
                                                                     "MC",
                                                                     "SP",
                                                                     "SC"))

result <- xtabs(data = df1, ~ Sampleinfo + Major_celltype) %>% prop.table(., margin = 1) %>% round(5) %>% as.data.frame()
setnames(result, old = "Major_celltype", new = "Celltype")

Pmapvalues <- function(dat1, dat2, x, from, to){
  if(length(from) != length(to)){
    stop("provide from and to with sampe length")
  } else {
    for (i in 1:length(from)) {
      print(c(to[i], from[i]))
      dat1[, to[i]] <- plyr::mapvalues(x = dat1[, x] %>% as.character(),
                                       from = dat2[, from[i]] %>% as.character(),
                                       to = dat2[, to[i]] %>% as.character(),
                                       warn_missing = F)
    }
    return(dat1)
  }
}

results <- Pmapvalues(result, meta.data, "Sampleinfo", from = rep("Sampleinfo", 6), to = c("Condition", "IgG", "Sex", "Patient_name", "Groups", "Sampleinfo"))
results <- setnames(results, old = "Patient_name", "Name")

results$Celltype <- factor(results$Celltype,
                                      levels = c('CD4 TN',
                                                 'CD4 TCM',
                                                 "CD4 TEM",
                                                 'CD4 TEFF',
                                                 'Treg',
                                                 'CD8 TN',
                                                 'CD8 TCM',
                                                 "CD8 TEM",
                                                 'CD8 TEFF',
                                                 'MAIT', 
                                                 'GD T1',
                                                 'GD T2',
                                                 'CD16 NK',
                                                 'CD56 NK',
                                                 "Pro lymphocyte"))

celltypes <- c('CD4 TN',
               'CD4 TCM',
               "CD4 TEM",
               'CD4 TEFF',
               'Treg',
               'CD8 TN',
               'CD8 TCM',
               "CD8 TEM",
               'CD8 TEFF',
               'MAIT', 
               'GD T1',
               'GD T2',
               'CD16 NK',
               'CD56 NK',
               "Pro lymphocyte")

results$Condition <- factor(results$Condition,
                                       levels = c(
                                         "IVV-D0",
                                         "IVV-D28",
                                         "IVV-D42",
                                         "HD",
                                         "MP",
                                         "MC",
                                         "SP",
                                         "SC"))

df1 <- results %>% subset(Groups == "Virus" & !is.na(Celltype)) 
df1 <- df1[df1$Celltype %in% c(df1$Celltype %>% unique() %>% grep(pattern = "*", value = T, invert = F)), ]
df1 <- df1 %>% arrange(Condition, Name) 

stat.test <- df1 %>%
  group_by(IgG, Celltype) %>%
  wilcox_test(Freq ~ Condition, p.adjust.method = "fdr", paired = T) %>%
  add_significance() %>% 
  add_xy_position(x = 'Condition')

stat.test$xmin <- stat.test$xmin

ggplot(df1, aes(x = Condition, y = Freq, color = Condition)) +
  scale_color_manual(values = c(alpha(condition_colors, 0.8)[])) +
  geom_boxplot() + 
  geom_point(aes(color = Condition), size = 1) +
  facet_grid(IgG ~ Celltype, scales = "free_y") +
  geom_line(data = df1,
            aes(x = Condition, y = Freq, group = Name), linetype = "dashed", color =  condition_colors[7], size = 0.5) +
  theme_classic2() +
  labs(title  = "Virus",
       y = "cell proportation") +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    # legend.position = "none",
    axis.text.x = element_blank()
    # axis.text.y = element_blank(),
    # axis.title = element_blank()
  ) + 
  theme(strip.background = element_rect(fill = alpha("grey", 0.4), colour = NA)) +
  stat_pvalue_manual(stat.test, label = "{p.adj}", hide.ns = T)

###---------------------------------------Vlnplot---------------------------###
migration <- read.table("migration.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
meta.data$migration_score <- mapvalues(meta.data %>% row.names(), from = migration %>% row.names(), to = migration$modulescore_migration1)
meta.data$Condition <- meta.data$Condition %>%
  as.character %>%
  gsub(pattern = "control-control-", replacement = "HD") %>% 
  gsub(pattern = "mild/moderate-progression-", replacement = "MP") %>% 
  gsub(pattern = "mild/moderate-convalescence-", replacement = "MC") %>% 
  gsub(pattern = "severe/critical-progression-", replacement = "SP") %>% 
  gsub(pattern = "severe/critical-convalescence-", replacement = "SC") %>% 
  gsub(pattern = "Virus-A0", replacement = "IVV-D0") %>% 
  gsub(pattern = "Virus-A28", replacement = "IVV-D28") %>% 
  gsub(pattern = "Virus-A42", replacement = "IVV-D42")

meta.data$Condition <- factor(meta.data$Condition %>% as.character(), levels = c("IVV-D0",
                                                                                 "IVV-D28",
                                                                                 "IVV-D42",
                                                                                 "HD",
                                                                                 "MP",
                                                                                 "MC",
                                                                                 "SP",
                                                                                 "SC"))

df1 <- meta.data %>% subset(Groups == "Virus")
df1$migration_score <- df1$migration_score %>% as.numeric()

stat.test <- df1 %>%
  group_by(IgG) %>%
  wilcox_test(modulescore_apoptosis1 ~ Condition, p.adjust.method = "fdr", paired = F) %>%
  add_significance() %>% 
  add_xy_position(x = 'Condition')

ggviolin(df1, "Condition", "modulescore_apoptosis1", fill = "Condition",
         # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         palette = c(alpha(condition_colors, 0.8)),
         facet.by = "IgG",
         add = "boxplot",
         add.params = list(fill = "white"),
         color = "Condition",
         width = 1
) + 
  stat_pvalue_manual(stat.test, label = "{p.adj.signif}", hide.ns = T) +
  theme(
    axis.title.y.right = element_blank(),  #
    # axis.ticks.y = element_blank(),      #
    axis.text.y = element_text(margin = margin(r = 0),colour = 'black',size = 10), #
    strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(colour = 'black',size = 10, hjust = 1, vjust = 0.5, angle = 90),
    panel.spacing = unit(-0.3, "mm"),
    strip.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")
    # strip.background = element_rect(size = 1, colour = "black",fill = NA, linetype = "solid")
  )

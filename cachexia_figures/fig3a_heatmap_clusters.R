# Author: Sonia Boscenco
# Create heatmap of body composition changes stratified by cluster 

rm(list = ls())
gc()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
bodycomp_deltas                        <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0302.csv")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clusters                               <- as.data.frame(as.factor(bodycomp_deltas$cluster)) 

rownames(clusters)                     <- bodycomp_deltas$MRN

delta_values                           <- bodycomp_deltas %>% 
                                          dplyr::select(starts_with(("delta_"))) %>%
                                          as.matrix()
rownames(delta_values)                 <- bodycomp_deltas$MRN

f1                                     <- colorRamp2(seq(-50, 50, length = 3), c("#4E79A7FF","white","#E15759FF"))

clean_labels                           <- gsub("delta_", "", colnames(delta_values))
clean_labels                           <- gsub("(?<=\\w)(Density)", " \\1", clean_labels, perl = TRUE)
clean_labels                           <- gsub("(?<=\\w)(Volume)", " \\1", clean_labels, perl = TRUE)
clean_labels                           <- gsub("Area", " Volume", clean_labels, perl = TRUE)
clean_labels                           <- gsub("Muscle Volume", "SKM", clean_labels, perl = TRUE)
clean_labels                           <- gsub("BMDLStandard", " Bone Mineral Density", clean_labels, perl = TRUE)

mat                 <- t(delta_values)
row_labels          <- clean_labels  
column_labels       <- rownames(delta_values) 

ht <- Heatmap(
  (mat),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,     
  show_column_names = FALSE, 
  row_names_side = "left",
  cluster_column_slices = FALSE, 
  row_labels = row_labels,
  column_labels = column_labels,
  col = f1,
  row_dend_gp = gpar(lwd = 0.1),
  row_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  row_title_gp = gpar(fontsize = 6, fontfamily = "ArialMT", fontface = "bold"),
  column_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontfamily = "ArialMT"),
    labels_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
    legend_height = unit(5, "mm"), 
    legend_width = unit(20, "mm"),     
    grid_height = unit(2, "mm"),      
    grid_width  = unit(2, "mm"),
    direction = 'horizontal',
    title = expression(Delta~Value),
    title_position = "topcenter"
  ),
  show_column_dend = FALSE,
  show_parent_dend_line = FALSE,
  show_row_dend = FALSE,
  # row split must now be indexed by patients (columns)
  column_split = factor(clusters$`as.factor(bodycomp_deltas$cluster)`, labels = (c("Type A", "Type B", "Type C"))),
  column_title_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  column_gap = unit(2, "mm"),
  height = unit(nrow(mat) * 3, "mm"),
  width = unit(60, "mm"),
  border = "black",
  border_gp = gpar(lwd = 0.1),
  use_raster = FALSE
)

pdf(file = "~/Desktop/reznik/bodycomp_main/results/cachexia//heatmap_clusters_0415.pdf", width = 4, height = 5)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Sonia Boscenco
# Validation cohort 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

uw_bodycomp                         <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/UW_CachexiaCohort_NoPHI_NewMarkers.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# choose only first event, and non pediatric patients 
uw_bodycomp$delta_bmi               <- uw_bodycomp$end_bmi - uw_bodycomp$start_bmi
uw_bodycomp                         <- uw_bodycomp %>% 
                                       filter(delta_bmi < 0) %>% 
                                       filter(PatientAge >= 18) %>% 
                                       group_by(coded_id) %>% 
                                       slice_min(start_day) %>%
                                       filter(StartEvent_IVContrast_IVContrastPresent == "TRUE" & EndEvent_IVContrast_IVContrastPresent == "TRUE")
 
# make sure numeric and QC using physiological ranges                     
to_numeric                           <- function(df, cols) {
  for (col in cols) {
    df[[col]] <- as.numeric(df[[col]])
    }
  return(df)
  }
num_cols                             <- colnames(uw_bodycomp)[grepl("Start|End", colnames(uw_bodycomp))]
uw_bodycomp                          <- to_numeric(uw_bodycomp, num_cols)
qc_rules                             <- list(
  "StartEvent_MuscleValues_L3TotalMuscleMedianHU"= c(-50, 200),
  "EndEvent_MuscleValues_L3TotalMuscleMedianHU"= c(-50, 200),
  "StartEvent_MuscleValues_L3TotalMuscleArea"    = c(25, 500),
  "EndEvent_MuscleValues_L3TotalMuscleArea"    = c(25, 500),
  "StartEvent_L3FatValues_L3VATArea"             = c(0.1, 1200),
  "EndEvent_L3FatValues_L3VATArea"             = c(0.1, 1200),
  "StartEvent_L3FatValues_L3SATArea"             = c(0.1, 1000),
  "EndEvent_L3FatValues_L3SATArea"             = c(0.1, 1000),
  "StartEvent_L3FatValues_L3SATMedian"           = c(-120, -30),
  "EndEvent_L3FatValues_L3SATMedian"           = c(-120, -30),
  "StartEvent_L3FatValues_L3VATMedian"           = c(-120, -60),
  "EndEventL3FatValues_L3VATMedian"           = c(-120, -60),
  
  "StartEvent_LiverValues_LiverMedianHU"         = c(-50, 180),
  "StartEvent_LiverValues_LiverVolume"           = c(100, 5000),
  "EndEvent_LiverValues_LiverMedianHU"         = c(-50, 180),
  "EndEvent_LiverValues_LiverVolume"           = c(100, 5000),
  
  "StartEvent_SpleenValues_SpleenMedianHU"       = c(10, 500),
  "StartEvent_SpleenValues_SpleenVolume"         = c(50, 6000),
  "EndEvent_SpleenValues_SpleenMedianHU"       = c(10, 500),
  "EndEvent_SpleenValues_SpleenVolume"         = c(50, 6000),
  
  "StartEvent_KidneyValues_KidneyMedianHU"       = c(5, 300),
  "StartEvent_KidneyValues_KidneyVolume"         = c(50, 750),
  "EndEvent_KidneyValues_KidneyMedianHU"       = c(5, 300),
  "EndEvent_KidneyValues_KidneyVolume"         = c(50, 750),
  
  "StartEvent_PancreasValues_PancreasMedianHU"   = c(-10, 260),
  "StartEvent_PancreasValues_PancreasVolume"     = c(20, 140),
  "EndEvent_PancreasValues_PancreasMedianHU"   = c(-10, 260),
  "EndEvent_PancreasValues_PancreasVolume"     = c(20, 140),
  
  "StartEvent_BMDL3Values_BMDL3HighSensitivityHU"= c(-50, 1200),
  "StartEvent_BMDL3Values_BMDL3StandardHU"       = c(-50, 1200),
  "EndEvent_BMDL3Values_BMDL3HighSensitivityHU"= c(-50, 1200),
  "EndEvent_BMDL3Values_BMDL3StandardHU"       = c(-50, 1200)
)

for (col in names(qc_rules)) {
  if (col %in% names(uw_bodycomp)) {
    print(col)
    range <- qc_rules[[col]]
    uw_bodycomp[[col]] <- ifelse(uw_bodycomp[[col]] < range[1] | uw_bodycomp[[col]] > range[2], NA, uw_bodycomp[[col]])
  }
}

uw_bodycomp                        <- uw_bodycomp  %>% 
                                      filter(StartEvent_L3FatValues_L3VATArea > 35 & EndEvent_L3FatValues_L3TATArea > 35)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate deltas
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bodycomp_deltas                   <- uw_bodycomp %>% 
                                     group_by(coded_id) %>%
                                     reframe(delta_MuscleDensity = EndEvent_MuscleValues_L3TotalMuscleMedianHU - StartEvent_MuscleValues_L3TotalMuscleMedianHU,
                                             delta_SATDensity    = EndEvent_L3FatValues_L3SATMedian - StartEvent_L3FatValues_L3SATMedian,
                                             delta_VATDensity    = EndEvent_L3FatValues_L3VATMedian - StartEvent_L3FatValues_L3VATMedian,
                                             delta_KidneyDensity = EndEvent_KidneyValues_KidneyMedianHU - StartEvent_KidneyValues_KidneyMedianHU,
                                             delta_PancreasDensity = EndEvent_PancreasValues_PancreasMedianHU - StartEvent_PancreasValues_PancreasMedianHU,
                                             delta_SpleenDensity   = EndEvent_SpleenValues_SpleenMedianHU - StartEvent_SpleenValues_SpleenMedianHU,
                                             delta_LiverDensity    = EndEvent_LiverValues_LiverMedianHU - StartEvent_LiverValues_LiverMedianHU,
    
                                             delta_LiverArea    = (EndEvent_LiverValues_LiverVolume - StartEvent_LiverValues_LiverVolume) / StartEvent_LiverValues_LiverVolume * 100,
                                             delta_KidneyVolume = (EndEvent_KidneyValues_KidneyVolume - StartEvent_KidneyValues_KidneyVolume) / StartEvent_KidneyValues_KidneyVolume * 100,
                                             delta_PancreasVolume = (EndEvent_PancreasValues_PancreasVolume - StartEvent_PancreasValues_PancreasVolume) / StartEvent_PancreasValues_PancreasVolume * 100,
                                             delta_SpleenVolume   = (EndEvent_SpleenValues_SpleenVolume - StartEvent_SpleenValues_SpleenVolume) / StartEvent_SpleenValues_SpleenVolume * 100,
    
                                             delta_BMDLStandard  = EndEvent_BMDL3Values_BMDL3StandardHU - StartEvent_BMDL3Values_BMDL3StandardHU,
                                             
                                             delta_MuscleArea = (EndEvent_MuscleValues_L3TotalMuscleArea - StartEvent_MuscleValues_L3TotalMuscleArea) / StartEvent_MuscleValues_L3TotalMuscleArea * 100,
                                             delta_SAT        = (EndEvent_L3FatValues_L3SATArea - StartEvent_L3FatValues_L3SATArea) / StartEvent_L3FatValues_L3SATArea * 100,
                                             delta_VAT        = (EndEvent_L3FatValues_L3VATArea- StartEvent_L3FatValues_L3VATArea) / StartEvent_L3FatValues_L3VATArea * 100,
                                             delta_IMAT       = (EndEvent_MuscleValues_L3IMATArea - StartEvent_MuscleValues_L3IMATArea) / (StartEvent_MuscleValues_L3IMATArea+ 1) * 100)

# keep just one instance per patient 
bodycomp_deltas                         <- bodycomp_deltas %>% 
                                           distinct(coded_id, .keep_all = TRUE) %>% 
                                           as.data.frame()
# some patients had same scan for pre and post cachexia so need to remove those  
bodycomp_deltas                         <- bodycomp_deltas %>% 
                                           filter(rowSums(across(all_of((starts_with("delta_"))))) != 0) %>%
                                           filter(delta_VAT < 500) %>%
                                           filter(delta_SAT < 500)
rownames(bodycomp_deltas)               <- bodycomp_deltas$coded_id

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# repeat subtyping as before 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

delta_cols                              <- colnames(bodycomp_deltas)[grepl("delta_", colnames(bodycomp_deltas))]
bodycomp_deltas_knn                     <- VIM::kNN(bodycomp_deltas, k = 10, imp_var = FALSE, variable = delta_cols)
rownames(bodycomp_deltas_knn)           <- bodycomp_deltas_knn$coded_id
bodycomp_deltas_knn                     <- bodycomp_deltas_knn[grepl("delta_", colnames(bodycomp_deltas_knn))]
bodycomp_deltas_knn_scaled              <- scale(bodycomp_deltas_knn, center = TRUE, scale = TRUE)
# find optimal number of clusters
wss                                      <- sapply(1:10, function(k) {kmeans(bodycomp_deltas_knn_scaled, centers = k, nstart = 25)$tot.withinss})

plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares")

km                                        <- kmeans(bodycomp_deltas_knn_scaled, centers = 3, nstart = 50)

cluster_df                                <- as.data.frame(km$cluster)
colnames(cluster_df)                      <- c("cluster")
cluster_df$MRN                            <- as.numeric((rownames(cluster_df)))


bodycomp_deltas$cluster                   <- cluster_df$cluster[match(bodycomp_deltas$coded_id, cluster_df$MRN)]
bodycomp_deltas$cluster_name              <- ifelse(bodycomp_deltas$cluster == 1, "Type A", ifelse(bodycomp_deltas$cluster == 2, "Type B", "Type C"))
clusters                                  <- as.data.frame(as.factor(bodycomp_deltas$cluster)) 
rownames(clusters)                        <- bodycomp_deltas$coded_id
rownames(bodycomp_deltas)                 <- bodycomp_deltas$coded_id
delta_values                              <- bodycomp_deltas %>% 
                                             dplyr::select(-c("coded_id", "cluster", "cluster_name")) %>%
                                             as.matrix()
rownames(delta_values)                    <- bodycomp_deltas$coded_id
f1                                        <- colorRamp2(seq(-50, 50, length = 3), c("#4E79A7FF", "white", "#E15759FF"))

clean_labels                              <- gsub("delta_", "", colnames(delta_values))
clean_labels                              <- gsub("(?<=\\w)(Density)", " \\1", clean_labels, perl = TRUE)
clean_labels                              <- gsub("(?<=\\w)(Volume)", " \\1", clean_labels, perl = TRUE)
clean_labels                              <- gsub("Area", " Volume", clean_labels, perl = TRUE)
clean_labels                              <- gsub("BMDLStandard", " Bone Mineral Density", clean_labels, perl = TRUE)

mat                                       <- t(delta_values) 
row_labels                                <- clean_labels  
column_labels                             <- rownames(delta_values) 

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
  border_gp = gpar(lwd = 0.1)
)

pdf(file = "~/Desktop/reznik/bodycomp_main/results/cachexia/clusters/heatmap_clusters_0227_uw.pdf", width = 4, height = 5)
draw(ht, heatmap_legend_side = "bottom")
dev.off()
bodycomp_deltas$cluster_name      <- factor(bodycomp_deltas$cluster_name, levels = c("Type A", "Type B", "Type C"))

pairwise.t.test(bodycomp_deltas$delta_SAT, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_SAT ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_SAT, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"SAT (%)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_SAT_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_VAT, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_VAT ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_VAT, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"VAT (%)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_VAT_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_MuscleArea, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_MuscleArea ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_MuscleArea, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"SKM (%)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_skm_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_LiverArea, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_LiverArea ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_LiverArea, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"Liver Volume (%)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_liver_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_LiverDensity, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_LiverDensity ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_LiverDensity, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"Liver Density (HU)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_liverdensity_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_SpleenVolume, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_SpleenVolume ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_SpleenVolume, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"Spleen Volume (%)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_spleen_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_SpleenDensity, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_SpleenDensity ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_SpleenDensity, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"Spleen Density (HU)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_spleendensity_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_PancreasVolume, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_PancreasVolume ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_PancreasVolume, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"Pancreas Volume (%)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_pancreas_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_PancreasDensity, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_PancreasDensity ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_PancreasDensity, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"Pancreas Density (HU)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 
ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_pancreasdensity_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_KidneyVolume, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_KidneyVolume ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_KidneyVolume, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"Kidney Volume (%)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_kidney_wastednonwasted.pdf", width =2, height =2, units = "in")

pairwise.t.test(bodycomp_deltas$delta_KidneyDensity, bodycomp_deltas$cluster_name, p.adjust.method = "BH")
bodycomp_deltas %>% cohens_d(delta_KidneyDensity ~ cluster_name)

p <- ggplot(bodycomp_deltas, aes(x = cluster_name, y = pmin(delta_KidneyDensity, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"Kidney Density (HU)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank()) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  scale_x_discrete(labels = c("Type A", "Type B", "Type C"), position = "bottom") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_boxplot_kidneydensity_wastednonwasted.pdf", width =2, height =2, units = "in")

bodycomp_deltas$cancertype   <- uw_bodycomp$Cancer.Primary.Site.Detail[match(bodycomp_deltas$coded_id, uw_bodycomp$coded_id)]
bodycomp_deltas$primary_site <- uw_bodycomp$primary_site_description[match(bodycomp_deltas$coded_id, uw_bodycomp$coded_id)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# histrogram distibutions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


bodycomp_changes_long <- pivot_longer(bodycomp_deltas, cols = starts_with("delta"), 
                                      names_to = "loss_organ", 
                                      values_to = "loss")

bodycomp_changes_long <- bodycomp_changes_long %>% filter(!loss_organ %in% c("delta_TotalAortic"))

bodycomp_changes_long$loss_organ <- factor(bodycomp_changes_long$loss_organ, levels = c("delta_SAT", "delta_VAT", "delta_IMAT",
                                                                                        "delta_MuscleArea", 
                                                                                        "delta_SATDensity", "delta_VATDensity",
                                                                                        "delta_MuscleDensity","delta_BMDLStandard",
                                                                                        "delta_KidneyVolume", "delta_SpleenVolume",
                                                                                        "delta_LiverArea", "delta_PancreasVolume",
                                                                                        "delta_KidneyDensity", "delta_SpleenDensity",
                                                                                        "delta_LiverDensity", "delta_PancreasDensity"))

props                           <- bodycomp_changes_long %>%
                                   group_by(loss_organ) %>%
                                   summarise(prop_neg = mean(loss < 0, na.rm = TRUE))

unique_bodycomp                <- unique(bodycomp_changes_long$loss_organ)
densities                      <- unique_bodycomp[grepl("Density|BMD", unique_bodycomp)]
non_densities                  <- unique_bodycomp[!unique_bodycomp %in% densities]

p <- ggplot(filter(bodycomp_changes_long, loss_organ %in% non_densities), aes(x = loss, fill = loss < 0)) + 
  geom_histogram(aes(fill = after_stat(x) < 0),
                 color = "white", linewidth = 0.1, bins = 20, stat = "bin") +
  scale_fill_manual(values = c("TRUE" = "#4E79A7FF", "FALSE" = "#E15759FF"), guide = "none") +
  theme_std() + 
  labs(x = "", y = "% change across cachexia") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(n.breaks = 4) +
  facet_wrap(~loss_organ, scales = "free", nrow = 2, 
             labeller = labeller(loss_organ = c(
               "delta_SAT" = "SAT",
               "delta_SATDensity" = "SAT Density",
               "delta_VAT" = "VAT",
               "delta_VATDensity" = "VAT Density",
               "delta_IMAT" = "IMAT",
               "delta_MuscleArea" = "SKM",
               "delta_MuscleDensity" = "SKM Density",
               "delta_SpleenVolume" = "Spleen",
               "delta_SpleenDensity" = "Spleen Density",
               "delta_LiverArea" = "Liver",
               "delta_LiverDensity" = "Liver Density",
               "delta_KidneyVolume" = "Kidney",
               "delta_KidneyDensity" = "Kidney Density",
               "delta_PancreasVolume" = "Pancreas",
               "delta_PancreasDensity" = "Pancreas Density",
               "delta_BMDLStandard" = "BMD"
             ))) + 
  theme(strip.background = element_blank())


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_histograms_bodycompchanges_allorgans.pdf", width = 5.5, height = 2)

p <- ggplot(filter(bodycomp_changes_long, loss_organ %in% densities), aes(x = loss)) + 
  geom_histogram(fill = "black",
                 color = "white", linewidth = 0.1, bins = 20, stat = "bin") +
  #scale_fill_manual(values = c("TRUE" = "#E15759", "FALSE" = "black"), guide = "none") +
  theme_std() + 
  labs(x = "", y = "raw change across cachexia (HU)") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(n.breaks = 4) +
  facet_wrap(~loss_organ, scales = "free", nrow = 2, 
             labeller = labeller(loss_organ = c(
               "delta_SAT" = "SAT",
               "delta_SATDensity" = "SAT Density",
               "delta_VAT" = "VAT",
               "delta_VATDensity" = "VAT Density",
               "delta_IMAT" = "IMAT",
               "delta_MuscleArea" = "SKM",
               "delta_MuscleDensity" = "SKM Density",
               "delta_SpleenVolume" = "Spleen",
               "delta_SpleenDensity" = "Spleen Density",
               "delta_LiverArea" = "Liver",
               "delta_LiverDensity" = "Liver Density",
               "delta_KidneyVolume" = "Kidney",
               "delta_KidneyDensity" = "Kidney Density",
               "delta_PancreasVolume" = "Pancreas",
               "delta_PancreasDensity" = "Pancreas Density",
               "delta_BMDLStandard" = "BMD"
             ))) + 
  theme(strip.background = element_blank())

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/uw_histograms_bodycompchanges_allorgansdensities.pdf",width = 6.25, height = 2.25)

n_patients <- nrow(bodycomp_deltas)
n_lost_sat <- sum(bodycomp_deltas$delta_SAT < 0)
paste0("Prop patients lost SAT:", n_lost_sat / n_patients)

n_lost_vat <- sum(bodycomp_deltas$delta_VAT < 0)
paste0("Prop patients lost VAT:", n_lost_vat / n_patients)

n_lost_muscle <- sum(bodycomp_deltas$delta_MuscleArea < 0)
paste0("Prop patients lost SKM:", n_lost_muscle / n_patients)

n_lost_imat <- sum(bodycomp_deltas$delta_IMAT < 0)
paste0("Prop patients lost IMAT:", n_lost_imat / n_patients)

n_lost_kidney <- sum(bodycomp_deltas$delta_KidneyVolume < 0, na.rm = TRUE)
paste0("Prop patients lost kidney:", n_lost_kidney / n_patients)

n_lost_spleen <- sum(bodycomp_deltas$delta_SpleenVolume < 0, na.rm = TRUE)
paste0("Prop patients lost spleen:", n_lost_spleen / n_patients)

n_lost_liver <- sum(bodycomp_deltas$delta_LiverArea < 0, na.rm = TRUE)
paste0("Prop patients lost liver:", n_lost_liver / n_patients)

n_lost_pancreas <- sum(bodycomp_deltas$delta_PancreasVolume < 0, na.rm = TRUE)
paste0("Prop patients lost pancreas:", n_lost_pancreas / n_patients)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# linear model for pancreatic adenocarcinoma 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bodycomp_deltas$contrast <- ifelse(bodycomp_deltas$cancertype %in% c("C25.9-Acinar cell carcinoma of pancreas (HCC)", "C25.1-Cancer of pancreas, body (HCC)", "C25.0-Adenocarcinoma of head of pancreas (HCC)", "Pancreas"), 1, 0)
bodycomp_deltas$contrast <- as.factor(bodycomp_deltas$contrast)

delta_cols <- colnames(bodycomp_deltas)[grepl("delta_", colnames(bodycomp_deltas))]

results_lm <- lapply(delta_cols, function(col) {
  
  model <- glm(
    as.formula(paste0(col, " ~ contrast")),
    data = bodycomp_deltas)
  tidy(model, conf.int = TRUE) |>
    filter(!grepl('Intercept', term)) |>
    mutate(bodycomp = col)}) |>
  bind_rows()

results_lm$p.adj <- p.adjust(results_lm$p.value, method = "BH")
results_lm$sig   <- ifelse(results_lm$p.adj < 0.05, T, F)

means <- bodycomp_deltas %>%
  summarise(across(all_of(delta_cols), ~mean(., na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "bodycomp", values_to = "mean_delta") %>%
  mutate(cleaned_bodycomp = str_replace(bodycomp, "delta_", "")) %>%
  mutate(cleaned_bodycomp = str_replace(cleaned_bodycomp, "Volume", " Volume")) %>%
  mutate(cleaned_bodycomp = str_replace(cleaned_bodycomp, "LiverArea", "Liver Volume")) %>%
  mutate(cleaned_bodycomp = str_replace(cleaned_bodycomp, "Density", " Density")) %>%
  mutate(cleaned_bodycomp = str_replace(cleaned_bodycomp, "Area", "")) %>%
  mutate(cleaned_bodycomp = str_replace(cleaned_bodycomp, "Muscle", "SKM")) %>%
  mutate(cleaned_bodycomp = str_replace(cleaned_bodycomp, "BMDLStandard", "BMD"))

results_lm$cleaned_bodycomp <- str_replace(results_lm$bodycomp, "delta_", "")
results_lm$cleaned_bodycomp <- str_replace(results_lm$cleaned_bodycomp, "Volume", " Volume")
results_lm$cleaned_bodycomp <- str_replace(results_lm$cleaned_bodycomp, "LiverArea", "Liver Volume")
results_lm$cleaned_bodycomp <- str_replace(results_lm$cleaned_bodycomp, "Density", " Density")
results_lm$cleaned_bodycomp <- str_replace(results_lm$cleaned_bodycomp, "Area", "")
results_lm$cleaned_bodycomp <- str_replace(results_lm$cleaned_bodycomp, "Muscle", "SKM")
results_lm$cleaned_bodycomp <- str_replace(results_lm$cleaned_bodycomp, "BMDLStandard", "BMD")

row_order <- c("Spleen Volume",
               "SAT Density",
               "VAT Density",
               "Pancreas Density",
               "Spleen Density",
               "Liver Density",
               "Kidney Density",
               "SKM Density",
               "BMD",
               "Liver Volume",
               "Kidney Volume",
               "SKM",
               "Pancreas Volume",
               "IMAT",
               "VAT",
               "SAT")


results_lm_with_means <- results_lm %>%
  left_join(means, by = "cleaned_bodycomp")
results_lm_with_means$cleaned_bodycomp <- factor(
  results_lm_with_means$cleaned_bodycomp, 
  levels = row_order
)

mat <- results_lm_with_means %>%
  dplyr::select(cleaned_bodycomp, mean_delta) %>%
  distinct() %>%
  arrange(cleaned_bodycomp) %>%
  column_to_rownames("cleaned_bodycomp") %>%
  as.matrix()

sig_mat <- results_lm_with_means %>%
  dplyr::select(cleaned_bodycomp, sig) %>%
  distinct() %>%
  column_to_rownames("cleaned_bodycomp") %>%
  as.matrix()


colnames(mat) <- "Pancreatic Cancer"
colnames(sig_mat) <- "Pancreatic Cancer"

col_fun <- colorRamp2(
  breaks = c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)), 
  colors = c("#59A14FFF", "white", "#D37295FF")
)


sig_layer <- function(j, i, x, y, w, h, fill) {
  if (!is.na(sig_mat[i, j]) && sig_mat[i, j]) {
    grid.points(x, y, pch = 8, size = unit(1, "mm"), gp = gpar(lwd = 0.3))
  }
}

p <- Heatmap(
  mat,
  col = col_fun,
  cluster_columns = FALSE,  
  cluster_rows = FALSE,      
  cell_fun = sig_layer,
  border = "white",
  row_names_side = "left",
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  column_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontfamily = "ArialMT"),
    labels_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
    legend_height = unit(5, "mm"), 
    legend_width = unit(20, "mm"),     
    grid_height = unit(2, "mm"),      
    grid_width  = unit(2, "mm"),
    title_position = "topcenter",
    title = expression("Mean\nChange"),
    color_space = "RGB"
  ),
  show_column_dend = FALSE,
  column_names_side = "bottom",
  width = unit(0.4, "cm"),
)

pdf(file = "~/Desktop/reznik/bodycomp_main/results/cachexia/UW_pancreatic_heatmap_cancertypeassociations_mean_value.pdf", width = 1.25, height = 3.1)
draw(p, heatmap_legend_side = "right")
dev.off()

mat <- results_means %>%
  dplyr::select(cancer_type, cleaned_bodycomp, mean_delta) %>%
  pivot_wider(names_from = cancer_type, values_from = mean_delta) %>%
  column_to_rownames("cleaned_bodycomp") %>%
  as.matrix()

sig_mat <- results_means %>%
  dplyr::select(cancer_type, cleaned_bodycomp, sig) %>%
  pivot_wider(names_from = cancer_type, values_from = sig) %>%
  column_to_rownames("cleaned_bodycomp") %>%
  as.matrix()

# Rest of your heatmap code
col_fun <- colorRamp2(
  breaks = c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)), 
  colors = c("#59A14FFF", "white", "#D37295FF")
)

sig_layer <- function(j, i, x, y, w, h, fill) {
  if (!is.na(sig_mat[i, j]) && sig_mat[i, j]) {
    grid.points(x, y, pch = 8, size = unit(1, "mm"), gp = gpar(lwd = 0.3))
  }
}


p <- Heatmap(
  mat,
  col = col_fun,
  cluster_columns = TRUE, 
  cluster_rows = TRUE,    
  cell_fun = sig_layer,
  border = "white",
  row_names_side = "left",
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  column_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontfamily = "ArialMT"),
    labels_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
    legend_height = unit(5, "mm"), 
    legend_width = unit(20, "mm"),     
    grid_height = unit(2, "mm"),      
    grid_width  = unit(2, "mm"),
    # direction = 'horizontal',
    title_position = "topcenter",
    title = expression("Mean\nChange"),
    color_space = "RGB"
  ),
  show_column_dend = FALSE,
  column_names_side = "bottom"
)

p <- ggplot(results_lm, aes(x = reorder(cleaned_bodycomp, OR) , y = log2(OR))) + 
  geom_stripped_cols() +
  ggtitle("UW Pancreatic Cancer") +
  geom_point(stroke = NA, aes(colour = sig)) + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  geom_errorbar(aes(ymin = log2(CI_low_OR), ymax = log2(CI_high_OR), colour = sig), linewidth = 0.1, width = 0.1) +
  coord_flip() +
  theme_std() +
  theme(legend.position = "none",
        plot.title = element_text(family = "ArialMT", size = 7, hjust = 0.5)) +
  xlab("") +
  ylab(expression(log[2](OR))) + 
  scale_colour_manual(values = c("black", "#E15759FF"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/UW_pancreatic_vs_rest_forestplot.pdf", width = 1.75, height = 2.25)


p <- ggplot(bodycomp_deltas, aes(x = factor(contrast, levels = c(1,0)), y = pmin(delta_SAT, 100))) + 
  geom_quasirandom(alpha = 0.5, aes(colour = contrast), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  labs(x = "") +
  ylab(expression(Delta~"(%)")) + 
  scale_x_discrete(labels = c("0" = "Other", "1" = "Pancreatic\nCancer")) +
  ggtitle("UW SAT") +
  theme_std() +
  theme(plot.title = element_text(family = "ArialMT", size = 7, face = "bold", hjust = 0.5),
        legend.position = "none") +
  scale_colour_manual(values = c("#F28E2BFF", "#D4A6C8FF"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/UW_pancreatic_vs_rest_SAT_boxplot.pdf", width = 1.5, height = 1.75)

p <- ggplot(bodycomp_deltas, aes(x = factor(contrast, levels = c(1,0)), y = delta_SATDensity)) + 
  geom_quasirandom(alpha = 0.5, aes(colour = contrast), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  labs(x = "") +
  ylab(expression(Delta~"(HU)")) + 
  scale_x_discrete(labels = c("0" = "Other", "1" = "Pancreatic\nCancer")) +
  ggtitle("UW SAT Density") +
  theme_std() +
  theme(plot.title = element_text(family = "ArialMT", size = 7, hjust = 0.5),
        legend.position = "none") +
  scale_colour_manual(values = c("#F28E2BFF", "#D4A6C8FF"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/UW_pancreatic_vs_rest_SATDensity_boxplot.pdf", width = 1.5, height = 1.75)

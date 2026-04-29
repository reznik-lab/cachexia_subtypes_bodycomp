# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Sonia Boscenco 
# Association of body composition changes within cancer types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
gc()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
bodycomp_deltas                        <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0302.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate means instead and do 1 vs. the rest for linear model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
delta_cols     <- colnames(bodycomp_deltas)[grepl("delta_", colnames(bodycomp_deltas))]
results_means  <- lapply(delta_cols, function(col) {
  cancer_types <- unique(bodycomp_deltas$CANCER_TYPE_DETAILED)
   pval_list   <- lapply(cancer_types, function(ct) {bodycomp_deltas$contrast <- ifelse(bodycomp_deltas$CANCER_TYPE_DETAILED == ct, 1, 0)

                                                     bodycomp_deltas$contrast <- as.factor(bodycomp_deltas$contrast)
                                                     model <- glm(as.formula(paste0(col, " ~ contrast + GENDER + AGE_CCX + START_CCX_BMI")), data = bodycomp_deltas)
                                                     tidy(model) %>%
                                                     filter(term == "contrast1") %>%
                                                     mutate(cancer_type = ct, bodycomp = col) %>%
                                                     dplyr::select(cancer_type, p.value, bodycomp, estimate)}) %>%
                                                     bind_rows()

  means        <- bodycomp_deltas %>%
                  group_by(CANCER_TYPE_DETAILED) %>%
                  summarise(mean_delta = mean(!!sym(col), na.rm = TRUE), .groups = 'drop', ) %>%
                  mutate(cancer_type = CANCER_TYPE_DETAILED) %>%
                  dplyr::select(-c(CANCER_TYPE_DETAILED)) %>%
                  mutate(bodycomp = col)
  
  means %>% left_join(pval_list, by = c('cancer_type', 'bodycomp')) %>%
            select(cancer_type, mean_delta, p.value, bodycomp, estimate)}) %>%
            bind_rows()


results_means$p.adj            <- p.adjust(results_means$p.value, method = "BH")
results_means$sig              <- ifelse(results_means$p.adj < 0.05, TRUE, FALSE)
results_means$cleaned_bodycomp <- str_replace(results_means$bodycomp, "delta_", "")
results_means$cleaned_bodycomp <- str_replace(results_means$cleaned_bodycomp, "Volume", " Volume")
results_means$cleaned_bodycomp <- str_replace(results_means$cleaned_bodycomp, "LiverArea", "Liver Volume")
results_means$cleaned_bodycomp <- str_replace(results_means$cleaned_bodycomp, "Density", " Density")
results_means$cleaned_bodycomp <- str_replace(results_means$cleaned_bodycomp, "Area", "")
results_means$cleaned_bodycomp <- str_replace(results_means$cleaned_bodycomp, "Muscle", "SKM")
results_means$cleaned_bodycomp <- str_replace(results_means$cleaned_bodycomp, "BMDLStandard", "BMD")

write.csv(results_means, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/tables/supplementary_table_2_multivariate_cancer_type_associations_means.csv", row.names = FALSE)

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

draw(p, heatmap_legend_side = "right")
pdf(file = "~/Desktop/reznik/bodycomp_main/results/cachexia/heatmap_cancertypeassociations_mean_values.pdf", width = 3.25, height = 3.75)
draw(p, heatmap_legend_side = "right")
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pancreatic cancer vignette
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bodycomp_deltas$class                 <- ifelse(bodycomp_deltas$CANCER_TYPE_DETAILED == "Pancreatic Adenocarcinoma", "Pancreatic Adenocarcinoma", "Other")

p <- ggplot(bodycomp_deltas, aes(x = factor(class, levels = c("Pancreatic Adenocarcinoma", "Other")), y = pmin(delta_SAT, 100))) +
  geom_quasirandom(alpha = 0.5, aes(colour = class), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  labs(x = "") +
  ylab(expression(Delta~"(%)")) + 
  scale_x_discrete(labels = c("Pancreatic Adenocarcinoma" = "Pancreatic\nAdenocarcinoma")) +
  ggtitle("MSKCC SAT") +
  theme_std() +
  theme(plot.title = element_text(family = "ArialMT", size = 7, hjust = 0.5),
        legend.position = "none") +
  scale_colour_manual(values = c("#59A14FFF", "#D37295FF"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/pancreatic_vs_rest_SAT_boxplot.pdf", width = 1.5, height = 1.75)

p <- ggplot(bodycomp_deltas, aes(x = factor(class, levels = c("Pancreatic Adenocarcinoma", "Other")), y = delta_SATDensity)) +
  geom_quasirandom(alpha = 0.5, aes(colour = class), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  labs(x = "") +
  ylab(expression(Delta~"(HU)")) + 
  scale_x_discrete(labels = c("Pancreatic Adenocarcinoma" = "Pancreatic\nAdenocarcinoma")) +
  ggtitle("MSKCC SAT Density") +
  theme_std() +
  theme(plot.title = element_text(family = "ArialMT", size = 7, hjust = 0.5),
        legend.position = "none") +
  scale_colour_manual(values = c("#59A14FFF", "#D37295FF"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/pancreatic_vs_rest_SATDensity_boxplot.pdf", width = 1.5, height = 1.75)



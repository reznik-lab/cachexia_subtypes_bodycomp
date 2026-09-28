# Author: Sonia Boscenco
# Cox proportional hazard model to identify which body composition changes are most important for patient outcomes

rm(list = ls())
gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data and scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
master_file                   <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0828.csv")
longitudinal_data             <- read.csv("~/Desktop/reznik/bodycomp_main/data/master_processed/processed_bodycomp_w_metadata_0828.csv")
cachexia_labels               <- read.csv("~/Desktop/reznik/impact/data/cachexia/cac_episodes_all_frozen_0930.csv")
mets                          <- read.csv("~/Desktop/reznik/bodycomp/data/metpred_out_venise_080425.csv")

# flag if wanting to save files
save_files                    <- F

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

master_file$END_CCX_SCAN_DATE   <- as.Date(master_file$END_CCX_SCAN_DATE)
master_file$START_CCX_SCAN_DATE <- as.Date(master_file$START_CCX_SCAN_DATE)

# if they passed, thats OS, otherwise its the last follow-up
master_file$OS_DATE             <- ifelse(is.na(master_file$PT_DEATH_DTE) | master_file$PT_DEATH_DTE == "", master_file$PLA_LAST_CONTACT_DTE, master_file$PT_DEATH_DTE)

# calculate time from end of cachexia as OS 
master_file$OS_CCX              <- (as.numeric(as.Date(master_file$OS_DATE) - as.Date(master_file$END_CCX_SCAN_DATE))) / 30.44

# annotate age as time for birth to start of ccx
master_file$AGE_CCX             <- (as.numeric(as.Date(master_file$START_CCX_SCAN_DATE) - as.Date(master_file$PT_BIRTH_DTE))) / 365.25

# get bmi at start of cachexia (only looking at the first episode always)
cachexia_labels                 <- cachexia_labels %>% 
                                   group_by(MRN) %>% 
                                   slice_min(order_by = start_day)

mets                            <- mets %>% 
  filter(DMP_ID %in% master_file$DMP_ID)
mets$END_DATE                   <- master_file$END_CCX_SCAN_DATE[match(mets$DMP_ID, master_file$DMP_ID)]

# count number of instance liver mets recovered before the end of cachexia
mets                            <- mets %>% 
  filter(as.Date(RADIOLOGY_PERFORMED_DATE) <= as.Date(END_DATE)) %>%
  group_by(DMP_ID) %>% 
  summarise(n_liver = sum(Liver), 
            n_bone = sum(Bone))

# keep only patients without any liver mets 
no_liver_mets                   <- mets %>% 
  filter(n_liver == 0)

bodycomp_no_liver_mets          <- master_file %>% 
  filter(DMP_ID %in% no_liver_mets$DMP_ID)

master_file$HAS_LIVER_MET <- ifelse(master_file$DMP_ID %in% bodycomp_no_liver_mets$DMP_ID, FALSE, TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make piecewise because u-shaped vars and divide by 10 to make more interpretable
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

master_file$PancreasGain        <- (ifelse(master_file$delta_PancreasVolume > 0, master_file$delta_PancreasVolume/10, 0))
master_file$PancreasLoss        <- (ifelse(master_file$delta_PancreasVolume < 0, abs(master_file$delta_PancreasVolume/10), 0))

master_file$KidneyGain          <- (ifelse(master_file$delta_KidneyVolume > 0, master_file$delta_KidneyVolume/10, 0))
master_file$KidneyLoss          <- (ifelse(master_file$delta_KidneyVolume < 0, abs(master_file$delta_KidneyVolume/10), 0))

master_file$LiverGain           <- (ifelse(master_file$delta_LiverArea > 0, master_file$delta_LiverArea/10, 0))
master_file$LiverLoss           <- (ifelse(master_file$delta_LiverArea < 0, abs(master_file$delta_LiverArea/10), 0))

master_file$SpleenGain          <- (ifelse(master_file$delta_SpleenVolume > 0, master_file$delta_SpleenVolume/10, 0))
master_file$SpleenLoss          <- (ifelse(master_file$delta_SpleenVolume < 0, abs(master_file$delta_SpleenVolume/10), 0))

master_file$MuscleGain          <- (ifelse(master_file$delta_MuscleArea > 0, master_file$delta_MuscleArea/10, 0))
master_file$MuscleLoss          <- (ifelse(master_file$delta_MuscleArea < 0, abs(master_file$delta_MuscleArea/10), 0))

master_file$SATGain            <- (ifelse(master_file$delta_SAT > 0, master_file$delta_SAT/10, 0))
master_file$SATLoss            <- (ifelse(master_file$delta_SAT < 0, abs(master_file$delta_SAT/10), 0))

master_file$VATGain            <- (ifelse(master_file$delta_VAT > 0, master_file$delta_VAT/10, 0))
master_file$VATLoss            <- (ifelse(master_file$delta_VAT < 0, abs(master_file$delta_VAT/10), 0))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set factors correctly
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
master_file$GENDER            <- as.factor(master_file$GENDER)
master_file$CANCER            <- as.factor(master_file$CANCER_TYPE_DETAILED)
master_file$Status            <- ifelse(master_file$OS_STATUS == "1:DECEASED", 1, 0)
master_file$STAGE             <- as.factor(master_file$STAGE_CCX)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# do it across all cancers 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cancers                <- unique(master_file$CANCER_TYPE_DETAILED)
results                <- data.frame(
                          cancer_type = character(),
                          bodycomp = character(),
                          p.val = numeric(),
                          HR = numeric())

for(ct in cancers){
  master_subset        <- master_file %>% filter(CANCER_TYPE_DETAILED == ct)
  model                <- coxph(Surv(OS_CCX, Status) ~ VATLoss + VATGain + SATGain + SATLoss + MuscleGain + MuscleLoss + SpleenGain + SpleenLoss + LiverGain + LiverLoss + KidneyGain + KidneyLoss + PancreasGain + PancreasLoss + GENDER  + AGE_CCX + START_CCX_BMI + STAGE, data = master_subset)
  
  results              <- rbind(results, data.frame(cancer_type = ct,
                                                    bodycomp = rownames(summary(model)$coefficients),
                                                    HR = summary(model)$coefficients[, "exp(coef)"],
                                                    p.value = summary(model)$coefficients[, "Pr(>|z|)"])) 
}

results$p.adj           <- p.adjust(results$p.value, method = "BH")
results$sign            <- ifelse(results$p.adj < 0.05, TRUE, FALSE)

# only plot bodycomp metrics 
results_sub             <- results %>% filter(grepl("Loss|Gain", bodycomp))
results_sub             <- results_sub %>% filter(p.adj < 0.05)
demographics            <- c("GENDERMALE", "AGE_CCX", "START_CCX_BMI", "STAGEStage 1-3", "STAGEStage 4")
# find number of cancer types that are significant per cancer types 
results_sub_tbl         <- results %>%
                           filter(!bodycomp %in% demographics) %>%
                           group_by(bodycomp) %>% 
                           summarise(n_sign = sum(sign == TRUE)) %>%
                           arrange((n_sign))

bodycomp_features       <- c("VATLoss", "VATGain", "SATGain", "SATLoss", "MuscleGain", "MuscleLoss", "SpleenGain", "SpleenLoss", "LiverGain", "LiverLoss", "KidneyGain", "KidneyLoss", "PancreasGain", "PancreasLoss")

# set up complex heatmap
hr_mat                  <- results %>% 
                           filter(!bodycomp %in% demographics & !is.na(bodycomp)) %>%
                           mutate(bodycomp = factor(bodycomp, levels = results_sub_tbl$bodycomp)) %>%
                           dplyr::select(bodycomp, cancer_type, HR) %>%
                           pivot_wider(names_from = cancer_type, values_from = HR) %>%
                           column_to_rownames("bodycomp") %>%
                           as.matrix()

hr_mat                  <- log2(hr_mat)

col_fun                 <- colorRamp2(c(-0.5, 0, 1.8),c("#4E79A7FF", "#F7F7F7", "#E15759FF"))
bodycomp_labels         <- c("MuscleLoss"   = "SKM loss",
                             "SATLoss"      = "SAT loss",
                             "LiverGain"    = "Liver gain",
                             "PancreasLoss" = "Pancreas loss",
                             "MuscleGain"   = "SKM gain",
                             "SpleenLoss"   = "Spleen loss",
                             "KidneyLoss"   = "Kidney loss",
                             "KidneyGain"   = "Kidney gain",
                             "PancreasGain" = "Pancreas gain",
                             "VATGain"      = "VAT gain",
                             "SpleenGain"   = "Spleen gain",
                             "VATLoss"      = "VAT loss",
                             "LiverLoss"    = "Liver loss",
                             "SATGain"      = "SAT gain")

sig_mat                  <- results %>%
                            mutate(sig = p.adj < 0.05) %>%
                            filter(!bodycomp %in% demographics) %>%
                            mutate(bodycomp = factor(bodycomp, levels = results_sub_tbl$bodycomp)) %>%
                            dplyr::select(bodycomp, cancer_type, sig) %>%
                            pivot_wider(names_from = cancer_type, values_from = sig) %>%
                            column_to_rownames("bodycomp") %>%
                            as.matrix()
row_order                <- rev(results_sub_tbl$bodycomp)

hr_mat                   <- hr_mat[row_order, , drop = FALSE]
sig_mat                  <- sig_mat[row_order, , drop = FALSE]

# to annotate the stars 
sig_layer                <- function(j, i, x, y, w, h, fill) {
                            if (!is.na(sig_mat[i, j]) && sig_mat[i, j]) {
                            grid.points(x, y, pch = 8, size = unit(1, "mm"), gp = gpar(lwd = 0.3))}}
ht <- Heatmap(
  hr_mat,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_side = "top",
  row_names_side = "left",
  row_labels = bodycomp_labels[rownames(hr_mat)],
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
    title = expression(log[2](HR)),
    title_position = "topcenter"
  ),
  column_names_rot = 90,
  rect_gp = gpar(col = "white", lwd = 0.1),
  cell_fun = sig_layer
)

draw(ht, heatmap_legend_side = "bottom")

p2 <- ggplot(results_sub_tbl, aes(x = factor(bodycomp, levels = bodycomp), y = n_sign)) + 
      geom_bar(stat = "identity", fill = "black", width = 0.8) + 
      coord_flip() + 
     geom_text(aes(label = n_sign), hjust = 2, size = 2, colour = "white") +
      theme_std() + 
      scale_y_continuous(expand = c(0,0), position = "right", n.breaks = 6) + 
      scale_x_discrete(expand = c(0,0)) + 
      labs(x = "",  y = "") + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            text = element_text(family = "ArialMT", size = 6))

# save files
pdf(file = "~/Desktop/reznik/bodycomp_main/revision/figures//heatmap_cox_percancertype.pdf", width = 3, height = 4)
draw(ht, heatmap_legend_side = "bottom")
dev.off()
ggsave(p2, file = "~/Desktop/reznik/bodycomp_main/revision/figures//barplot_numcancertypes_sign_heatmap.pdf", height = 2.25, width = 1)
write.csv(results, file = "~/Desktop/reznik/bodycomp_main/revision//tables/supplementary_table_cox_proportional_across_cancers.csv", row.names = FALSE)

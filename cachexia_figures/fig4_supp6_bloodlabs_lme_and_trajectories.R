# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Sonia Boscenco
# Identify changes in blood labs between Type A and Type C, and Type B and Type C
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab                                 <- read.csv("~/Desktop/reznik/bodycomp_main/data/clinical/lab_cleaned_0909.csv")
bodycomp_metadata                   <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0302.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# clean data  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab                                <- lab %>%
                                      dplyr::filter(MRN %in% bodycomp_metadata$MRN) 
lab$PRE_CCX_DATE                   <- bodycomp_metadata$START_CCX_SCAN_DATE[match(lab$MRN, bodycomp_metadata$MRN)]
lab$POST_CCX_DATE                  <- bodycomp_metadata$END_CCX_SCAN_DATE[match(lab$MRN, bodycomp_metadata$MRN)]

# get only labs that were drawn during the cachectic window
lab_during_ccx                     <- lab %>% 
                                      filter(as.Date(Date) >= as.Date(PRE_CCX_DATE) & as.Date(Date) <= as.Date(POST_CCX_DATE)) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# derived blood lab markers  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab_during_ccx$NLR                              <- lab_during_ccx$Neut/lab_during_ccx$Lymph
lab_during_ccx$`Prognostic Nutritional Index`   <- lab_during_ccx$Albumin + 5*lab_during_ccx$Lymph

lab_during_ccx                       <- lab_during_ccx %>%
                                        mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA, .x)))

lab_names                            <- colnames(lab_during_ccx) 
lab_names                            <- lab_names[!lab_names %in% c("MRN", "Date", "PRE_CCX_DATE", "POST_CCX_DATE")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add metadata
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab_during_ccx$SEX                   <- as.factor(bodycomp_metadata$GENDER[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)])
lab_during_ccx$CANCERTYPE            <- as.factor(bodycomp_metadata$CANCER_TYPE_DETAILED[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)])
lab_during_ccx$AGE                   <- bodycomp_metadata$AGE_CCX[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)]
lab_during_ccx$BMI                   <- bodycomp_metadata$CCX_START_BMI[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)]
lab_during_ccx$STAGE                 <- as.factor(bodycomp_metadata$STAGE_CDM_DERIVED_GRANULAR[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)])
lab_during_ccx$CLUSTER               <- factor(bodycomp_metadata$cluster_name[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)], levels = c("Type C", "Type B", "Type A"))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up mixed linear effects model (across cancer types)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

labs_long                            <- lab_during_ccx %>% 
                                        pivot_longer(cols = all_of(lab_names), 
                                                     names_to = "lab_name", 
                                                     values_to = "lab_value")

options(emmeans = list(lmerTest.limit = 10000000, pbkrtest.limit = 10000000))

results_allcancers                   <- labs_long %>%
                                        group_by(lab_name) %>%
                                        nest() %>%
                                        mutate(model = map(data, ~ lmer(log(lab_value+0.00001) ~ CLUSTER  + CANCERTYPE + SEX + AGE + STAGE + BMI + (1 | MRN),, data = .x, REML = TRUE)),
                                               contrasts = map(model, ~ pairs(emmeans(.x, ~ CLUSTER), adjust = "none", reverse = TRUE)),
                                               tidied = map(contrasts, broom::tidy)) %>%
                                        unnest(tidied) %>%
                                        ungroup() %>% 
                                        filter(contrast %in% c("Type A - Type C", "Type B - Type C")) %>%
                                        mutate(adj.p.value = p.adjust(p.value, method = "BH")) %>%
                                        dplyr::select(lab_name, contrast, estimate, statistic,df, adj.p.value, std.error) %>%
                                        arrange(adj.p.value)

results_allcancers                   <- results_allcancers %>%
                                        left_join(labs_long %>%
                                        group_by(lab_name) %>%
                                        summarise(sd_lab = sd(lab_value, na.rm = TRUE), .groups = "drop"), by = "lab_name") %>%
                                        mutate(std_estimate = estimate / sd_lab) %>%
                                        dplyr::select(lab_name, contrast, estimate, std_estimate, std.error, adj.p.value)


results_allcancers$fold_change      <- exp(results_allcancers$estimate)
results_allcancers$colour           <- ifelse(results_allcancers$adj.p.value < 0.05, "sign", "no")

# make some lab names nicer when plotting 
results_allcancers$lab_name         <- gsub("..Total", "", results_allcancers$lab_name)
results_allcancers$lab_name         <- gsub("\\.", " ", results_allcancers$lab_name)
results_allcancers$lab_name         <- gsub("Eos", "Eosinophils", results_allcancers$lab_name)
results_allcancers$lab_name         <- gsub("Baso", "Basophil", results_allcancers$lab_name)
results_allcancers$lab_name         <- gsub("Neut", "Neutrophil", results_allcancers$lab_name)
results_allcancers$lab_name         <- gsub("Lymph", "Lymphocyte", results_allcancers$lab_name)
results_allcancers$lab_name         <- gsub("Mono", "Monocyte", results_allcancers$lab_name)
results_allcancers$lab_name         <- gsub("ALK", "Alk. Phos.", results_allcancers$lab_name)

# filtering for dual forest plot 
sub_results                         <- results_allcancers %>%
                                       filter(contrast %in% c("Type A - Type C", "Type B - Type C"))

list_bloodlabs                      <- sub_results %>%
                                       filter(contrast == "Type A - Type C") %>%
                                       arrange(estimate) %>%
                                       pull(lab_name)

p <- ggplot(sub_results, aes(x = factor(lab_name, levels = list_bloodlabs), y = (estimate), colour = contrast, shape = colour)) + 
  geom_stripped_cols(colour = NA) +
  geom_point(alpha = 1) + 
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.1, linetype = "dashed") +
  theme_std() + 
  labs(x = "", shape = "", colour = "") +
  ylab(expression(log[2](Fold-Change))) +
  scale_shape_manual(values = c("no" = 2, "sign" = 8),
                     labels = c("no" = "not-significant", "sign" = "significant")) +
  scale_colour_manual(labels = c("Type A - Type C" = "Type A vs. Type C",
                                 "Type B - Type C" = "Type B vs. Type C"),
                                 values = rev(c("#499894", "#B07AA1")))+
 
  theme(axis.ticks.y = element_blank()) 


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/lme_blood_labs.pdf", width = 4, height = 3.25)
write.csv(results_allcancers, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/tables/supp_all_labs_lme.csv", row.names = FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up mixed linear effects model (within each cancer type)
# sanitry check to make sure that trend holds 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unique_cancers              <- unique(lab_during_ccx$CANCERTYPE)

for(cancer in unique_cancers){
  formula <- if(cancer %in% c("Breast Invasive Ductal Carcinoma","Prostate","Uterine Endometrioid Carcinoma",
                              "Prostate Adenocarcinoma","High-Grade Serous Ovarian Cancer",
                              "Uterine Serous Carcinoma")){
    log(lab_value+0.000001) ~ CLUSTER + AGE + BMI + STAGE + (1 | MRN)
  } else {
    log(lab_value+0.000001) ~ CLUSTER + SEX + AGE + BMI + STAGE+ (1 | MRN)
  }
  
  
  df_cancer        <- labs_long %>% filter(CANCERTYPE == cancer)
  results <- df_cancer %>%
    group_by(lab_name) %>%
    nest() %>%
    mutate(model = map(data, ~ lmer(formula, data = .x, REML = TRUE)),
           contrasts = map(model, ~ pairs(emmeans(.x, ~ CLUSTER), adjust = "none", reverse = TRUE)),
           tidied = map(contrasts, broom::tidy)) %>%
    unnest(tidied) %>%
    ungroup() %>% 
    filter(contrast %in% c("Type A - Type C", "Type B - Type C")) %>%
    mutate(adj.p.value = p.adjust(p.value, method = "BH")) %>%
    dplyr::select(lab_name, contrast, estimate, statistic,df, adj.p.value, std.error) %>%
    arrange(adj.p.value)
  
  results_allcancers                   <- results %>%
    left_join(labs_long %>%
                group_by(lab_name) %>%
                summarise(sd_lab = sd(lab_value, na.rm = TRUE), .groups = "drop"), by = "lab_name") %>%
    mutate(std_estimate = estimate / sd_lab) %>%
    mutate(fold_change = exp(estimate)) %>%
    dplyr::select(lab_name, contrast, estimate, std_estimate, std.error, adj.p.value, fold_change)
  
  results_allcancers$colour           <- ifelse(results_allcancers$adj.p.value < 0.05, "sign", "no")
  
  # make some lab names nicer when plotting 
  results_allcancers$lab_name         <- gsub("..Total", "", results_allcancers$lab_name)
  results_allcancers$lab_name         <- gsub("\\.", " ", results_allcancers$lab_name)
  results_allcancers$lab_name         <- gsub("Eos", "Eosinophils", results_allcancers$lab_name)
  results_allcancers$lab_name         <- gsub("Baso", "Basophil", results_allcancers$lab_name)
  results_allcancers$lab_name         <- gsub("Neut", "Neutrophil", results_allcancers$lab_name)
  results_allcancers$lab_name         <- gsub("Lymph", "Lymphocyte", results_allcancers$lab_name)
  results_allcancers$lab_name         <- gsub("Mono", "Monocyte", results_allcancers$lab_name)
  results_allcancers$lab_name         <- gsub("ALK", "Alk. Phos.", results_allcancers$lab_name)
  
  list_bloodlabs                      <- results_allcancers %>%
    filter(contrast == "Type A - Type C") %>%
    arrange(fold_change) %>%
    pull(lab_name)
  
  p <- ggplot(results_allcancers, aes(x = factor(lab_name, levels = list_bloodlabs), y = (fold_change), colour = contrast, shape = colour)) + 
    geom_stripped_cols(colour = NA) +
    geom_point(alpha = 1) + 
    coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.1, linetype = "dashed") +
    theme_std() + 
    labs(x = "", shape = "", colour = "") +
    ylab(expression(~beta/SD)) +
    ggtitle(cancer) +
    scale_shape_manual(values = c("no" = 2, "sign" = 8),
                       labels = c("no" = "not-significant", "sign" = "significant")) +
    scale_colour_manual(labels = c("Type A - Type C" = "Type A vs. Type C",
                                   "Type B - Type C" = "Type B vs. Type C"),
                        values = c("#499894", "#B07AA1"))+
    
    theme(axis.ticks.y = element_blank(),
          plot.title = element_text(family = "ArialMT", size = 7, hjust = 0.5, face = "bold")) 
  
  # Save plot
  ggsave(filename = paste0("~/Desktop/reznik/bodycomp/results/cancertypemixedlmblood/volcano_", cancer, ".pdf"),
         plot = p, width = 6, height = 4, dpi = 300)
  message("Saved plot for ", cancer)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# trends of select labs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
copy_of_labs                     <- copy(lab_during_ccx)
copy_of_labs$Date                <- as.Date(copy_of_labs$Date)

# get normalized time
copy_of_labs                     <- copy_of_labs %>% 
                                    group_by(MRN) %>% 
                                    arrange(as.Date(Date)) %>% 
                                    mutate(days_since_first = as.numeric(Date - min(Date))) %>% 
                                    mutate(norm_time = days_since_first /( as.numeric(max(Date) - min(Date)))) %>% 
                                    ungroup()

# find mean within bins
mean_points                      <- copy_of_labs %>%
                                    mutate(time_bin = round(norm_time / 0.05) * 0.05) %>%   
                                    group_by(CLUSTER, time_bin) %>%
                                    summarise(mean_albumin = mean(Albumin, na.rm = TRUE),
                                              ci_albumin = qt(0.975, sum(!is.na(Albumin)) - 1)*(sd(Albumin, na.rm = TRUE) / sqrt(sum(!is.na(Albumin)))),
                                              mean_ALK = mean(ALK, na.rm = TRUE),
                                              mean_NLR = mean(NLR, na.rm = TRUE), 
                                              mean_HGB = mean(HGB, na.rm = TRUE),
                                              mean_ast = mean(AST, na.rm = TRUE),
                                              mean_billi = mean(Bilirubin..Total, na.rm = TRUE),
                                              ci_ast = qt(0.975, sum(!is.na(AST)))*(sd(AST, na.rm = TRUE) / sqrt(sum(!is.na(AST)))),
                                              ci_hgb = qt(0.975, sum(!is.na(HGB)))*(sd(HGB, na.rm = TRUE) / sqrt(sum(!is.na(HGB)))),
                                              ci_bill = qt(0.975, sum(!is.na(Bilirubin..Total)))*(sd(Bilirubin..Total, na.rm = TRUE) / sqrt(sum(!is.na(Bilirubin..Total)))),
                                              ci_NLR = qt(0.975, sum(!is.na(NLR)))*(sd(NLR, na.rm = TRUE) / sqrt(sum(!is.na(NLR)))),
                                              ci_ALK = qt(0.975, sum(!is.na(ALK)) - 1)*(sd(ALK, na.rm = TRUE) / sqrt(sum(!is.na(ALK)))),.groups = "drop")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# albumin
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p <- ggplot(copy_of_labs, aes(x = norm_time, y = Albumin)) + 
  geom_line(alpha = 0.05, linewidth = 0.05, mapping = aes(group = MRN), colour = "#BAB0ACFF") + 
  geom_smooth(data = mean_points, aes(y = mean_albumin, colour = CLUSTER, x = time_bin), method = "loess", se = FALSE, linewidth = 0.25, inherit.aes = FALSE) +
  geom_hline(yintercept = 3.5, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_hline(yintercept = 5.4, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_point(data = mean_points, aes(x = time_bin, y = mean_albumin, colour = CLUSTER), size = 0.5, stroke = NA) +
  geom_errorbar(data = mean_points, aes(x = time_bin, ymin = mean_albumin - ci_albumin, ymax = mean_albumin + ci_albumin, colour = CLUSTER), inherit.aes = FALSE, width = 0.01, linewidth = 0.1) +
  scale_colour_manual(values = rev(c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))) +
  #facet_wrap(~CLUSTER) + 
  scale_y_continuous(expand = c(0,0), limits = c(3, 4)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_std() +
  theme(
        panel.spacing = unit(0.5, "cm"),
        strip.background = element_rect(size = 0.1)) +
  labs(x = "Normalized time", y = "Albumin", colour = "")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/longitudinal_albumin_all.pdf", width = 2.5, height = 1.5, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# alk phos
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p <- ggplot(copy_of_labs, aes(x = norm_time, y = ALK)) + 
     geom_line(alpha = 0.05, linewidth = 0.05, mapping = aes(group = MRN), colour = "#BAB0ACFF") + 
     geom_smooth(data = mean_points, aes(y = mean_ALK, colour = CLUSTER, x = time_bin), method = "loess", se = FALSE, linewidth = 0.25, inherit.aes = FALSE) +
     geom_hline(yintercept = 44, colour = "black", linetype = "dashed", linewidth = 0.1) +
     geom_hline(yintercept = 147, colour = "black", linetype = "dashed", linewidth = 0.1) +
     geom_point(data = mean_points, aes(x = time_bin, y = mean_ALK, colour = CLUSTER), size = 0.5, stroke = NA) +
     geom_errorbar(data = mean_points, aes(x = time_bin, ymin = mean_ALK - ci_ALK, ymax = mean_ALK + ci_ALK, colour = CLUSTER), inherit.aes = FALSE, width = 0.01, linewidth = 0.1) +
     scale_colour_manual(values = rev(c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))) +
     #facet_wrap(~CLUSTER) + 
     scale_y_continuous(expand = c(0,0), limits = c(50, 300)) +
     scale_x_continuous(expand = c(0,0)) +
     theme_std() +
     theme(
           panel.spacing = unit(0.5, "cm"),
           strip.background = element_rect(size = 0.1)) +
     labs(x = "Normalized time", y = "Alk. Phos.", colour = "")


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/longitudinal_alkphos_all.pdf", width = 2.5, height = 1.5, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# billirubin
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p <- ggplot(copy_of_labs, aes(x = norm_time, y = Bilirubin..Total)) + 
  geom_line(alpha = 0.05, linewidth = 0.05, mapping = aes(group = MRN), colour = "#BAB0ACFF") + 
  geom_smooth(data = mean_points, aes(y = mean_billi, colour = CLUSTER, x = time_bin), method = "loess", se = FALSE, linewidth = 0.25, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.2, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_hline(yintercept = 1.2, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_point(data = mean_points, aes(x = time_bin, y = mean_billi, colour = CLUSTER), size = 0.5, stroke = NA) +
  geom_errorbar(data = mean_points, aes(x = time_bin, ymin = mean_billi - ci_bill, ymax = mean_billi + ci_bill, colour = CLUSTER), inherit.aes = FALSE, width = 0.01, linewidth = 0.1) +
  scale_colour_manual(values = rev(c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))) +
  #facet_wrap(~CLUSTER) + 
  scale_y_continuous(expand = c(0,0), limits = c(0.4, 1.25)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_std() +
  theme(
    panel.spacing = unit(0.5, "cm"),
    strip.background = element_rect(size = 0.1)) +
  labs(x = "Normalized time", y = "Bilirubin", colour = "")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/longitudinal_bilirubin_all.pdf", width = 2.5, height = 1.5, dpi = 300)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NLR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p <- ggplot(copy_of_labs, aes(x = norm_time, y = NLR)) + 
  geom_line(alpha = 0.05, linewidth = 0.05, mapping = aes(group = MRN), colour = "#BAB0ACFF") + 
  geom_smooth(data = mean_points, aes(y = mean_NLR, colour = CLUSTER, x = time_bin), method = "loess", se = FALSE, linewidth = 0.25, inherit.aes = FALSE) +
  geom_hline(yintercept = 1, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_hline(yintercept = 3, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_point(data = mean_points, aes(x = time_bin, y = mean_NLR, colour = CLUSTER), size = 0.5, stroke = NA) +
  geom_errorbar(data = mean_points, aes(x = time_bin, ymin = mean_NLR - ci_NLR, ymax = mean_NLR + ci_NLR, colour = CLUSTER), inherit.aes = FALSE, width = 0.01, linewidth = 0.1) +
  scale_colour_manual(values = rev(c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))) +
  #facet_wrap(~CLUSTER) + 
  scale_y_continuous(expand = c(0,0), limits = c(4, 12)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_std() +
  theme(
        panel.spacing = unit(0.5, "cm"),
        strip.background = element_rect(size = 0.1)) +
  labs(x = "Normalized time", y = "NLR", colour = "")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/longitudinal_nlr_all.pdf", width = 2.5, height = 1.5, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# hgb
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p <- ggplot(copy_of_labs, aes(x = norm_time, y = HGB)) + 
  geom_line(alpha = 0.05, linewidth = 0.05, mapping = aes(group = MRN), colour = "#BAB0ACFF") + 
  geom_smooth(data = mean_points, aes(y = mean_HGB, colour = CLUSTER, x = time_bin), method = "loess", se = FALSE, linewidth = 0.25, inherit.aes = FALSE) +
  geom_hline(yintercept = 12, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_hline(yintercept = 16, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_point(data = mean_points, aes(x = time_bin, y = mean_HGB, colour = CLUSTER), size = 0.5, stroke = NA) +
  geom_errorbar(data = mean_points, aes(x = time_bin, ymin = mean_HGB - ci_hgb, ymax = mean_HGB + ci_hgb, colour = CLUSTER), inherit.aes = FALSE, width = 0.01, linewidth = 0.1) +
  scale_colour_manual(values = rev(c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))) +
  #facet_wrap(~CLUSTER) + 
  scale_y_continuous(expand = c(0,0), limits = c(9, 12)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_std() +
  theme(
    panel.spacing = unit(0.5, "cm"),
    strip.background = element_rect(size = 0.1)) +
  labs(x = "Normalized time", y = "HGB", colour = "")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/longitudinal_hgp_all.pdf", width = 2.5, height = 1.5, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AST
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p <- ggplot(copy_of_labs, aes(x = norm_time, y = AST)) + 
  geom_line(alpha = 0.05, linewidth = 0.05, mapping = aes(group = MRN), colour = "#BAB0ACFF") + 
  geom_smooth(data = mean_points, aes(y = mean_ast, colour = CLUSTER, x = time_bin), method = "loess", se = FALSE, linewidth = 0.25, inherit.aes = FALSE) +
  geom_hline(yintercept = 10, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_hline(yintercept = 40, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_point(data = mean_points, aes(x = time_bin, y = mean_ast, colour = CLUSTER), size = 0.5, stroke = NA) +
  geom_errorbar(data = mean_points, aes(x = time_bin, ymin = mean_ast - ci_ast, ymax = mean_ast + ci_ast, colour = CLUSTER), inherit.aes = FALSE, width = 0.01, linewidth = 0.1) +
  scale_colour_manual(values = rev(c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))) +
  #facet_wrap(~CLUSTER) + 
  scale_y_continuous(expand = c(0,0), limits = c(20, 69)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_std() +
  theme(
    panel.spacing = unit(0.5, "cm"),
    strip.background = element_rect(size = 0.1)) +
  labs(x = "Normalized time", y = "AST", colour = "")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/longitudinal_ast_all.pdf", width = 2.5, height = 1.5, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# big matrix comparisons within each cancer type
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
results_allcancers            <- labs_long %>%
                                 group_by(lab_name, CANCERTYPE) %>%
                                 nest() %>%
                                 mutate(has_sex_variation = map_lgl(data,~ dplyr::n_distinct(.x$SEX[!is.na(.x$SEX)]) > 1),
                                        formula = map(has_sex_variation, ~ if (.x) {lab_value ~ CLUSTER + SEX + AGE + STAGE + BMI + (1 | MRN)} 
                                                      else {log(lab_value+0.000001) ~ CLUSTER + AGE + STAGE + BMI + (1 | MRN)}),
                                        model = map2(formula, data, ~ lmer(.x, data = .y, REML = TRUE)),
                                        contrasts = map(model,~ pairs(emmeans(.x, ~ CLUSTER), adjust = "none", reverse = TRUE)),
                                        tidied = map(contrasts, broom::tidy)) %>%
                                  unnest(tidied) %>%
                                  ungroup() %>%
                                  filter(contrast %in% c("Type A - Type C", "Type B - Type C")) %>%
                                  mutate(adj.p.value = p.adjust(p.value, method = "BH")) %>%
                                  mutate(fold_change = exp(estimate)) %>%
                                 dplyr::select(lab_name, CANCERTYPE, contrast, estimate, std.error, statistic, df, p.value, adj.p.value, fold_change) %>%
                                 arrange(adj.p.value)

results_allcancers$sign      <- ifelse(results_allcancers$p.value < 0.05, TRUE, FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compare slopes across lab values 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get first or labst labs 
first_last_labs               <- lab_during_ccx %>%
                                 group_by(MRN) %>%
                                 filter(Date == min(Date, na.rm = TRUE) | Date == max(Date, na.rm = TRUE))

# calculate delta 
first_last_labs               <- first_last_labs %>%
                                 group_by(MRN) %>%
                                 arrange(Date) %>%
                                 summarise(across(all_of(lab_names), list(delta = ~ dplyr::last(.) - dplyr::first(.)), .names = "{.col}_delta"))
                                  

first_last_labs$CLUSTER       <- bodycomp_metadata$cluster_name[match(first_last_labs$MRN, bodycomp_metadata$MRN)]
pairwise.t.test(first_last_labs$NLR_delta, first_last_labs$CLUSTER)
delta_cols <- names(first_last_labs)[grep("_delta$", names(first_last_labs))]

pairwise_results <- map_df(delta_cols, function(col) {
  test_result    <- pairwise.t.test(first_last_labs[[col]], first_last_labs$CLUSTER)
  cluster_means  <- first_last_labs %>% 
                    group_by(CLUSTER) %>%
                    summarise(mean_val = mean(!!sym(col), na.rm = TRUE), .groups = 'drop')
  
  p_values      <- test_result$p.value
  
  if (!is.null(p_values)) {as.data.frame(as.table(p_values)) %>%
                           filter(!is.na(Freq)) %>%
                           mutate(lab_value = col, cluster_1 = as.character(Var1), cluster_2 = as.character(Var2), pvalue = Freq) %>%
                           select(-Var1, -Var2, -Freq) %>%
                           left_join(cluster_means %>% dplyr::rename(cluster_1 = CLUSTER, mean_1 = mean_val), by = "cluster_1") %>%
                           left_join(cluster_means %>% dplyr::rename(cluster_2 = CLUSTER, mean_2 = mean_val), by = "cluster_2") %>%
                           mutate(FC = mean_1 / mean_2, log2FC = log2(FC)) %>%
                           select(lab_value, cluster_1, cluster_2, mean_1, mean_2, FC, log2FC, pvalue)}}) 
pairwise_results         <- pairwise_results %>% filter(cluster_1 == "Type C")
pairwise_results$q.value <- p.adjust(pairwise_results$pvalue, method = "BH")

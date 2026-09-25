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
bodycomp_metadata                   <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0828.csv")

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
lab_during_ccx$STAGE                 <- as.factor(bodycomp_metadata$STAGE_CCX[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)])
lab_during_ccx$CLUSTER               <- factor(bodycomp_metadata$cluster_name[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)], levels = c("Type C", "Type B", "Type A"))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up mixed linear effects model (across cancer types)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

labs_long                            <- lab_during_ccx %>% 
                                        pivot_longer(cols = all_of(lab_names), 
                                                     names_to = "lab_name", 
                                                     values_to = "lab_value") %>%
                                        group_by(lab_name) %>%
                                        mutate(lab_value_z = as.numeric(scale(lab_value))) %>%
                                        ungroup()

options(emmeans = list(lmerTest.limit = 10000000, pbkrtest.limit = 10000000))
plan(multisession, workers = future::availableCores() - 1)

fit_lab <- function(df) {
  
  m     <- glmmTMB(lab_value ~ CLUSTER + CANCERTYPE + SEX + AGE + STAGE + BMI + (1 | MRN),data = df,family = tweedie(link = "log"))
  emm   <- emmeans(m, ~ CLUSTER, type = "link")
  broom::tidy(pairs(emm, adjust = "none", reverse = TRUE))
}
results_allcancers <- labs_long %>%
                      filter(!is.na(lab_value)) %>%
                      group_by(lab_name) %>%
                      nest() %>%
                      mutate(tidied = future_map(data, fit_lab, .options = furrr_options(seed = TRUE))) %>%
                      select(-data) %>%
                      unnest(tidied) %>%
                      ungroup() %>%
                      filter(contrast %in% c("Type A - Type C", "Type B - Type C")) %>%
                      mutate(log2FC      = estimate / log(2),
                             log2FC_SE   = std.error / log(2),
                             adj.p.value = p.adjust(p.value, method = "BH")) %>%
                      dplyr::select(lab_name, contrast, log2FC, log2FC_SE, statistic, adj.p.value) %>%
                      arrange(adj.p.value)

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
                                       arrange(log2FC) %>%
                                       pull(lab_name)

p <- ggplot(sub_results, aes(x = factor(lab_name, levels = list_bloodlabs), y = (log2FC), colour = contrast, shape = colour)) + 
  geom_stripped_cols(colour = NA) +
  geom_point(alpha = 1) + 
  geom_errorbar(aes(ymin = log2FC - log2FC_SE, ymax = log2FC + log2FC_SE),
                width = 0.1, alpha = 0.6, linewidth = 0.1) +
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.1, linetype = "dashed") +
  theme_std() + 
  labs(x = "", shape = "", colour = "") +
  scale_y_continuous(limits = c(-0.1, 0.8)) +
  ylab(expression(log[2](Fold-Change))) +
  scale_shape_manual(values = c("no" = 2, "sign" = 8),
                     labels = c("no" = "not-significant", "sign" = "significant")) +
  scale_colour_manual(labels = c("Type A - Type C" = "Type A vs. Type C",
                                 "Type B - Type C" = "Type B vs. Type C"),
                                 values = rev(c("#499894", "#B07AA1")))+
 
  theme(axis.ticks.y = element_blank()) 


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision/results///lme_blood_labs.pdf", width = 4, height = 3.25)
write.csv(results_allcancers, file = "~/Desktop/reznik/bodycomp_main/revision///tables/supp_all_labs_lme_09252026.csv", row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# trends of select labs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
copy_of_labs                     <- copy(lab_during_ccx)
copy_of_labs$Date                <- as.Date(copy_of_labs$Date)

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
  #geom_hline(yintercept = 3.5, colour = "black", linetype = "dashed", linewidth = 0.1) +
  #geom_hline(yintercept = 5.4, colour = "black", linetype = "dashed", linewidth = 0.1) +
  geom_point(data = mean_points, aes(x = time_bin, y = mean_albumin, colour = CLUSTER), size = 0.5, stroke = NA) +
  geom_errorbar(data = mean_points, aes(x = time_bin, ymin = mean_albumin - ci_albumin, ymax = mean_albumin + ci_albumin, colour = CLUSTER), inherit.aes = FALSE, width = 0.01, linewidth = 0.1) +
  scale_colour_manual(values = rev(c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))) +
  #facet_wrap(~CLUSTER) + 
  scale_y_continuous(expand = c(0,0), limits = c(-2.5, 2.5)) +
  #scale_x_continuous(expand = c(0,0)) +
  theme_std() +
  theme(
        panel.spacing = unit(0.5, "cm"),
        strip.background = element_rect(size = 0.1)) +
  labs(x = "Normalized time", y = "Albumin (g/dL)", colour = "")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision//results//longitudinal_albumin_all.pdf", width = 2, height = 1.5, dpi = 300)

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
     labs(x = "Normalized time", y = "Alk. Phos. (IU/L)", colour = "")


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision//results//longitudinal_alkphos_all.pdf", width = 2, height = 1.5, dpi = 300)

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
  labs(x = "Normalized time", y = "Bilirubin (mg/dL)", colour = "")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision//results//longitudinal_bilirubin_all.pdf", width = 2, height = 1.5, dpi = 300)


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

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision//results//longitudinal_nlr_all.pdf", width = 2, height = 1.5, dpi = 300)

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
  labs(x = "Normalized time", y = "HGB (g/dL)", colour = "")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision//results//longitudinal_hgp_all.pdf", width = 2, height = 1.5, dpi = 300)

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

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision//results//longitudinal_ast_all.pdf", width = 2.5, height = 1.5, dpi = 300)

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

lab_long <- lab_during_ccx %>%
            select(MRN, Date, all_of(lab_names)) %>%
            pivot_longer(cols = all_of(lab_names), names_to = "lab_name", values_to = "value") %>%
            filter(!is.na(value), !is.na(Date))

# at least two time points
lfc_long <- lab_long %>%
            group_by(MRN, lab_name) %>%
            filter(n_distinct(Date) >= 2) %>%  
            filter(Date == min(Date) | Date == max(Date)) %>%
            arrange(Date, .by_group = TRUE) %>%
            summarise(n_obs      = n(),
                      first_val  = dplyr::first(value),
                      last_val   = dplyr::last(value),
                      log2fc     = if_else(first_val > 0 & last_val > 0,
                      log2(last_val / first_val), NA_real_), .groups = "drop")

lfc_wide         <- lfc_long %>%
                    select(MRN, lab_name, log2fc) %>%
                    pivot_wider(names_from = lab_name, values_from = log2fc, names_glue = "{lab_name}_log2fc")

lfc_wide$CLUSTER <- bodycomp_metadata$cluster_name[match(lfc_wide$MRN, bodycomp_metadata$MRN)]
lfc_cols         <- names(lfc_wide)[grep("_log2fc$", names(lfc_wide))]

pairwise_results <- map_df(lfc_cols, function(col) {
  test_result    <- pairwise.t.test(lfc_wide[[col]], lfc_wide$CLUSTER, pool.sd = FALSE)
  cluster_means  <- lfc_wide %>%
                    group_by(CLUSTER) %>%
                    summarise(mean_val = mean(!!sym(col), na.rm = TRUE), .groups = 'drop')
  
  p_values       <- test_result$p.value
  
  if (!is.null(p_values)) {
    as.data.frame(as.table(p_values)) %>%
      filter(!is.na(Freq)) %>%
      mutate(lab_value = col,
             cluster_1 = as.character(Var1),
             cluster_2 = as.character(Var2),
             pvalue    = Freq) %>%
      select(-Var1, -Var2, -Freq) %>%
      left_join(cluster_means %>% dplyr::rename(cluster_1 = CLUSTER, mean_1 = mean_val), by = "cluster_1") %>%
      left_join(cluster_means %>% dplyr::rename(cluster_2 = CLUSTER, mean_2 = mean_val), by = "cluster_2") %>%
      select(lab_value, cluster_1, cluster_2, mean_1, mean_2, pvalue)}})

pairwise_results <- pairwise_results %>%
                    filter((cluster_1 == "Type A" & cluster_2 == "Type C") |
                    (cluster_1 == "Type C" & cluster_2 == "Type A") |
                    (cluster_1 == "Type B" & cluster_2 == "Type C") |
                    (cluster_1 == "Type C" & cluster_2 == "Type B")) %>%
                    mutate(swap = cluster_2 == "Type C", 
                           comparison_cluster = if_else(swap, cluster_1, cluster_2),  
                           mean_log2fc_ref    = if_else(swap, mean_2, mean_1),
                           mean_log2fc_comp   = if_else(swap, mean_1, mean_2),
                           log2FC_diff        = mean_log2fc_comp - mean_log2fc_ref,
                           FC_diff            = 2^log2FC_diff) %>%
                     select(lab_value, comparison_cluster, mean_log2fc_ref, mean_log2fc_comp, log2FC_diff, FC_diff, pvalue)

pairwise_results$q.value <- p.adjust(pairwise_results$pvalue, method = "BH")

# Author: Sonia Boscenco
# Compare lipid panels wasted vs. Type C patients 

rm(list = ls())
gc()

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab                                 <- read.csv("~/Desktop/reznik/bodycomp_main/data/clinical/lipidpanels.csv")
bodycomp_metadata                   <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0302.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# clean data  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab                                <- lab %>%
                                      dplyr::filter(DMP_ID %in% bodycomp_metadata$DMP_ID)
lab$MRN                            <- bodycomp_metadata$MRN[match(lab$DMP_ID, bodycomp_metadata$DMP_ID)]
lab$PRE_CCX_DATE                   <- bodycomp_metadata$START_CCX_SCAN_DATE[match(lab$MRN, bodycomp_metadata$MRN)]
lab$POST_CCX_DATE                  <- bodycomp_metadata$END_CCX_SCAN_DATE[match(lab$MRN, bodycomp_metadata$MRN)]

# get only labs that were drawn during the cachectic window
lab_during_ccx                     <- lab %>% 
                                      filter(as.Date(LR_PERFORMED_DTE) >= as.Date(PRE_CCX_DATE) & as.Date(LR_PERFORMED_DTE) <= as.Date(POST_CCX_DATE)) 
lab_during_ccx$LR_RESULT_VALUE     <- as.numeric(lab_during_ccx$LR_RESULT_VALUE)
lab_during_ccx                     <- lab_during_ccx %>%
                                      mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA, .x))) 
lab_names                            <- unique(lab_during_ccx$CLEANED_TEST_NAME)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add metadata
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab_during_ccx$SEX                   <- as.factor(bodycomp_metadata$GENDER[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)])
lab_during_ccx$CANCERTYPE            <- as.factor(bodycomp_metadata$CANCER_TYPE_DETAILED[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)])
lab_during_ccx$AGE                   <- bodycomp_metadata$AGE_CCX[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)]
lab_during_ccx$CLUSTER               <- factor(bodycomp_metadata$cluster_name[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)], levels = c("Type C", "Type B", "Type A"))
lab_during_ccx$BMI                   <- bodycomp_metadata$CCX_START_BMI[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)]
lab_during_ccx$STAGE                 <- as.factor(bodycomp_metadata$STAGE_CDM_DERIVED_GRANULAR[match(lab_during_ccx$MRN, bodycomp_metadata$MRN)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up mixed linear effects model (across cancer types)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

results_allcancers                   <- lab_during_ccx %>%
                                        group_by(CLEANED_TEST_NAME) %>%
                                        nest() %>%
                                        mutate(model = map(data, ~ lmer(log(LR_RESULT_VALUE+0.00001) ~ CLUSTER  + CANCERTYPE + SEX + AGE + BMI + STAGE + (1 | MRN), data = .x, REML = TRUE)),
                                               contrasts = map(model, ~ pairs(emmeans(.x, ~ CLUSTER), adjust = "none", reverse = TRUE)),
                                               tidied = map(contrasts, broom::tidy)) %>%
                                        unnest(tidied) %>%
                                        ungroup() %>% 
                                        filter(contrast %in% c("Type A - Type C", "Type B - Type C")) %>%
                                        mutate(adj.p.value = p.adjust(p.value, method = "none")) %>%
                                        dplyr::select(CLEANED_TEST_NAME, contrast, estimate, statistic,df, adj.p.value, std.error) %>%
                                        arrange(adj.p.value)

results_allcancers$fold_change       <- exp(results_allcancers$estimate)
results_allcancers <- results_allcancers %>%
                      left_join(lab_during_ccx %>%
                                group_by(CLEANED_TEST_NAME) %>%
                                summarise(sd_lab = sd(LR_RESULT_VALUE, na.rm = TRUE), .groups = "drop"), by = "CLEANED_TEST_NAME") %>%
                      mutate(std_estimate = estimate / sd_lab) %>%
                      mutate(adj.p.value = p.adjust())
                      dplyr::select(CLEANED_TEST_NAME, contrast, estimate, std_estimate, std.error, adj.p.value, fold_change)

list_bloodlabs                      <- results_allcancers %>%
                                       filter(contrast == "Type A - Type C") %>%
                                       arrange(estimate) %>%
                                       pull(CLEANED_TEST_NAME)
results_allcancers$colour           <- ifelse(results_allcancers$adj.p.value < 0.05, "sign", "no")

p <- ggplot(results_allcancers, aes(x = factor(CLEANED_TEST_NAME, levels = list_bloodlabs), y = (estimate), colour = contrast, shape = colour)) + 
  geom_stripped_cols(colour = NA) +
  geom_point(alpha = 1) + 
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.1, linetype = "dashed") +
  theme_std() + 
  labs(x = "", shape = "", colour = "") +
  ylab(expression(log[2](Fold-Change))) +
  scale_shape_manual(values = c("no" = 2, "sign" = 8),
                     labels = c("no" = "not-significant", "sign" = "significant")) +
  scale_colour_manual(labels = c("(Type A) - (Type C)" = "Type A vs. Type C",
                                 "(Type B) - (Type C)" = "Type B vs. Type C"),
                      values = rev(c("#499894", "#B07AA1")))+
  
  theme(axis.ticks.y = element_blank()) 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/bloodlabs/lipids_forestplot_labs_linear_mixed_effects_wastedvsnonwasted.pdf", width = 3.25, height = 1.25)
write.csv(results_allcancers, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/tables/supp_lipids_lme.csv", row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# boxplot of just CRP to show 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

crp                         <- lab_during_ccx %>% 
                               filter(CLEANED_TEST_NAME == "C-Reactive Protein") %>%
                               group_by(MRN) %>%
                               summarise(mean(LR_RESULT_VALUE, na.rm = TRUE))

crp$CLUSTER                 <- bodycomp_metadata$cluster_name[match(crp$MRN, bodycomp_metadata$MRN)]
pairwise.t.test(crp$`mean(LR_RESULT_VALUE, na.rm = TRUE)`, crp$CLUSTER, "BH", pool.sd = FALSE)
p <- ggplot(crp, aes(x = CLUSTER, y = pmin(`mean(LR_RESULT_VALUE, na.rm = TRUE)`, 10), fill = CLUSTER)) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5) +
  geom_jitter(alpha = 0.5, width = 0.1, stroke = NA) +
  labs(x = "", y = "mean CRP (mg/dl)", fill = "") +
  theme_std() + 
  scale_fill_manual(values = rev(c("#A0CBE8FF", "#499894", "#B07AA1"))) +
  theme(legend.position = "none") +
  scale_x_discrete(expand = c(0,0.5))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/crp_boxplot_clusters.pdf", width = 1.75, height = 2)


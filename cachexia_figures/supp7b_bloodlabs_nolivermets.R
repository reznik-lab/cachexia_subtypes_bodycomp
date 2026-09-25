
rm(list = ls())
gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data and scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
all_bodycomp                    <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0828.csv")
mets                            <- read.csv("~/Desktop/reznik/bodycomp/data/metpred_out_venise_080425.csv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mets                            <- mets %>% 
                                   filter(DMP_ID %in% all_bodycomp$DMP_ID)
mets$END_DATE                   <- all_bodycomp$END_CCX_SCAN_DATE[match(mets$DMP_ID, all_bodycomp$DMP_ID)]

# count number of instance liver mets recovered before the end of cachexia
mets                            <- mets %>% 
                                   filter(as.Date(RADIOLOGY_PERFORMED_DATE) <= as.Date(END_DATE)) %>%
                                   group_by(DMP_ID) %>% 
                                   summarise(n_liver = sum(Liver),
                                   n_bone = sum(Bone))

# keep only patients without any liver or bone mets 
no_liver_mets                       <- mets %>% 
                                       filter(n_liver == 0 & n_bone == 0)

bodycomp_no_liver_mets              <- all_bodycomp %>% 
                                       filter(DMP_ID %in% no_liver_mets$DMP_ID)

bodycomp_no_liver_mets$cluster_name <- factor(bodycomp_no_liver_mets$cluster_name, levels = c("Type C", "Type B", "Type A"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# all other blood labs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab                                <- read.csv("~/Desktop/reznik/bodycomp_main/data/clinical/lab_cleaned_0909.csv") %>% 
                                      filter(MRN %in% bodycomp_no_liver_mets$MRN)

lab$PRE_CCX_DATE                   <- bodycomp_no_liver_mets$START_CCX_SCAN_DATE[match(lab$MRN, bodycomp_no_liver_mets$MRN)]
lab$POST_CCX_DATE                  <- bodycomp_no_liver_mets$END_CCX_SCAN_DATE[match(lab$MRN, bodycomp_no_liver_mets$MRN)]

# get only labs that were drawn during the cachectic window
lab_during_ccx                     <- lab %>% 
                                      filter(as.Date(Date) >= as.Date(PRE_CCX_DATE) & as.Date(Date) <= as.Date(POST_CCX_DATE)) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# derived blood lab markers  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab_during_ccx$NLR                            <- lab_during_ccx$Neut/lab_during_ccx$Lymph
lab_during_ccx$`Prognostic Nutritional Index` <- lab_during_ccx$Albumin + 5*lab_during_ccx$Lymph

lab_during_ccx                                <- lab_during_ccx %>%
                                                 mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA, .x)))

lab_names                            <- colnames(lab_during_ccx) 
lab_names                            <- lab_names[!lab_names %in% c("MRN", "Date", "PRE_CCX_DATE", "POST_CCX_DATE")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add metadata
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lab_during_ccx$SEX                   <- as.factor(bodycomp_no_liver_mets$GENDER[match(lab_during_ccx$MRN, bodycomp_no_liver_mets$MRN)])
lab_during_ccx$CANCERTYPE            <- as.factor(bodycomp_no_liver_mets$CANCER_TYPE_DETAILED[match(lab_during_ccx$MRN, bodycomp_no_liver_mets$MRN)])
lab_during_ccx$AGE                   <- bodycomp_no_liver_mets$AGE_CCX[match(lab_during_ccx$MRN, bodycomp_no_liver_mets$MRN)]
lab_during_ccx$BMI                   <- bodycomp_no_liver_mets$CCX_START_BMI[match(lab_during_ccx$MRN, bodycomp_no_liver_mets$MRN)]
lab_during_ccx$STAGE                 <- as.factor(bodycomp_no_liver_mets$STAGE_CCX[match(lab_during_ccx$MRN, bodycomp_no_liver_mets$MRN)])
lab_during_ccx$CLUSTER               <- factor(bodycomp_no_liver_mets$cluster_name[match(lab_during_ccx$MRN, bodycomp_no_liver_mets$MRN)], levels = c("Type C", "Type B", "Type A"))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up mixed linear effects model (across cancer types)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

labs_long                            <- lab_during_ccx %>% 
                                        pivot_longer(cols = all_of(lab_names), names_to = "lab_name", values_to = "lab_value") %>%
                                        ungroup()

options(emmeans = list(lmerTest.limit = 10000000, pbkrtest.limit = 10000000))
plan(multisession, workers = future::availableCores() - 1)

fit_lab <- function(df) {
  
  m       <- glmmTMB(lab_value ~ CLUSTER + CANCERTYPE + SEX + AGE + STAGE + BMI + (1 | MRN), data = df, family = tweedie(link = "log"))
  emm     <- emmeans(m, ~ CLUSTER, type = "link")
  broom::tidy(pairs(emm, adjust = "none", reverse = TRUE))}

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
                      dplyr::select(lab_name, contrast, log2FC, log2FC_SE, statistic, adj.p.value, p.value) %>%
                      arrange(p.value)

results_allcancers$colour           <- ifelse(results_allcancers$p.value < 0.05, "sign", "no")

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


p <- ggplot(sub_results, aes(x = factor(lab_name, levels = list_bloodlabs), y = log2FC, colour = contrast, shape = colour)) + 
  geom_stripped_cols(colour = NA) +
  geom_point(alpha = 1) + 
  geom_errorbar(aes(ymin = log2FC - log2FC_SE, ymax = log2FC + log2FC_SE),
                width = 0.1, alpha = 0.6, linewidth = 0.1) +
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.1, linetype = "dashed") +
  theme_std() + 
  labs(x = "", shape = "", colour = "") +
  ylab(expression(log[2](Fold-Change))) +
  #ggtitle("No Liver Metastasis (n = 1,177)") +
  scale_shape_manual(values = c("no" = 2, "sign" = 8),
                     labels = c("no" = "not significant", "sign" = "nominal P < 0.05")) +
  scale_colour_manual(labels = c("(Inflammatory-Wasted) - (Non-Wasted)" = "Inflammatory-Wasted vs. Non-Wasted",
                                 "(Atrophy-Wasted) - (Non-Wasted)" = "Atrophy-Wasted vs. Non-Wasted"),
                      values = rev(c("#499894", "#B07AA1")))+
  
  theme(axis.ticks.y = element_blank(),
        plot.title = element_text(family = "ArialMT", size = 7, hjust = 0.5, face = "bold")) 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision/figures/lme_blood_labs_no_liver_mets.pdf", width = 4, height = 3.25)


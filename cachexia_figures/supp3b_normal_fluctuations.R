# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Sonia Boscenco
# Normal fluctuations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
bodycomp                               <- read.csv("~/Desktop/reznik/bodycomp_main/data/main_processed/processed_bodycomp_w_metadata_0220.csv")
bodycomp_deltas                        <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0302.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# make sure cohort includes only patients without cachexia 
no_cachexia                           <- bodycomp %>% 
                                         filter(HAS_CCX == 0)

# at least 4 time points required per patient
no_cachexia                           <- no_cachexia %>% 
                                         group_by(MRN) %>% 
                                         mutate(n = n()) %>% 
                                         filter(n > 4) %>% 
                                         dplyr::select(-c(n)) %>% 
                                         ungroup()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to calculate changes in body composition from one time point to the next
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

calculate_body_composition_deltas <- function(df) {
                                     df %>% 
                                     arrange(MRN, SCAN_ORDER) %>%
                                     group_by(MRN) %>%
    mutate(
      delta_MuscleDensity   = lead(MuscleValues.L3TotalMuscleMedianHU) - MuscleValues.L3TotalMuscleMedianHU,
      delta_SATDensity      = lead(L3FatValues.L3SATMedian) - L3FatValues.L3SATMedian,
      delta_VATDensity      = lead(L3FatValues.L3VATMedian) - L3FatValues.L3VATMedian,
      delta_KidneyDensity   = lead(KidneyValues.KidneyMedianHU) - KidneyValues.KidneyMedianHU,
      delta_PancreasDensity = lead(PancreasValues.PancreasMedianHU) - PancreasValues.PancreasMedianHU,
      delta_SpleenDensity   = lead(SpleenValues.SpleenMedianHU) - SpleenValues.SpleenMedianHU,
      delta_LiverDensity    = lead(LiverValues.LiverMedianHU) - LiverValues.LiverMedianHU,
      
      delta_LiverArea =
        (lead(LiverValues.LiverVolume) - LiverValues.LiverVolume) /
        LiverValues.LiverVolume * 100,
      
      delta_KidneyVolume =
        (lead(KidneyValues.KidneyVolume) - KidneyValues.KidneyVolume) /
        KidneyValues.KidneyVolume * 100,
      
      delta_PancreasVolume =
        (lead(PancreasValues.PancreasVolume) - PancreasValues.PancreasVolume) /
        PancreasValues.PancreasVolume * 100,
      
      delta_SpleenVolume =
        (lead(SpleenValues.SpleenVolume) - SpleenValues.SpleenVolume) /
        SpleenValues.SpleenVolume * 100,
      
      delta_BMDLStandard =
        lead(BMDL3Values.BMDL3StandardHU) - BMDL3Values.BMDL3StandardHU,
      
      delta_MuscleArea =
        (lead(MuscleValues.L3TotalMuscleArea) - MuscleValues.L3TotalMuscleArea) /
        MuscleValues.L3TotalMuscleArea * 100,
      
      delta_SAT =
        (lead(L3FatValues.L3SATArea) - L3FatValues.L3SATArea) /
        L3FatValues.L3SATArea * 100,
      
      delta_VAT =
        (lead(L3FatValues.L3VATArea) - L3FatValues.L3VATArea) /
        L3FatValues.L3VATArea * 100,
      
      delta_IMAT =
        (lead(MuscleValues.L3IMATArea) - MuscleValues.L3IMATArea) /
        (MuscleValues.L3IMATArea) * 100,
    
      scans_days =
        as.numeric(as.Date(lead(SCAN_DATE)) - as.Date(SCAN_DATE)),
      
      scans_months = scans_days * 30.44
    ) %>%
    ungroup() %>%
    filter(!is.na(scans_days))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate normal fluctuations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

longitudinal_deltas                    <- calculate_body_composition_deltas(no_cachexia)

# normalize by time 
longitudinal_deltas                    <- longitudinal_deltas %>% 
                                          mutate(across(starts_with("delta_"), ~.x / (scans_days/30.44)))

# normalize cachexia by time as well so that they are comparable
bodycomp_deltas_time                   <- bodycomp_deltas %>% 
                                          mutate(across(starts_with("delta_"), ~.x / (scans_days/30.44)))

# take the median fluctuation of each non-cachectic control
longitudinal_deltas_per_patients       <- longitudinal_deltas %>% group_by(MRN) %>% 
                                          summarise(across(starts_with("delta_"), ~median(.x, na.rm = TRUE)),.groups = "drop")

# combine dfs for plotting
ccx                                    <- bodycomp_deltas_time %>% 
                                          dplyr::select(starts_with("delta_"))
no_ccx                                 <- longitudinal_deltas_per_patients %>% 
                                          dplyr::select(starts_with("delta_"))


ccx$type                               <- 'ccx'
no_ccx$type                            <- 'no_ccx'
combined_df                            <- rbind(ccx, no_ccx)
combined_df$type                       <- factor(combined_df$type, levels = c("ccx", "no_ccx"))
combined_df$class                      <- ifelse(combined_df$type == "ccx", 1, 0)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# statistics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

t.test(delta_SAT ~ type, data = combined_df, exact = FALSE, )
cohens_d(delta_SAT ~ type, data = combined_df)

t.test(delta_VAT ~ type, data = combined_df, exact = FALSE)
cohens_d(delta_VAT ~ type, data = combined_df)
t.test(delta_MuscleArea ~ type, data = combined_df, exact = FALSE)
cohens_d(delta_MuscleArea ~ type, data = combined_df)

t.test(delta_LiverArea ~ type, data = combined_df, exact = FALSE)
cohens_d(delta_LiverArea ~ type, data = combined_df)

t.test(delta_PancreasVolume ~ type, data = combined_df, exact = FALSE)
cohens_d(delta_PancreasVolume ~ type, data = combined_df)

t.test(delta_KidneyVolume ~ type, data = combined_df, exact = FALSE)
cohens_d(delta_KidneyVolume ~ type, data = combined_df)

t.test(delta_SpleenVolume ~ type, data = combined_df, exact = FALSE)
cohens_d(delta_SpleenVolume ~ type, data = combined_df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p <- ggplot(combined_df, aes(x = pmax(pmin(delta_SAT, 30), -30), fill = type)) + 
  geom_density(alpha = 0.8, linewidth = 0.1) + 
  ylab("Density") + 
  xlab(expression("%"~Delta~SAT/month)) + 
  labs(fill = "") + 
  scale_fill_manual(values = c("#D4A6C8FF", "#8CD17DFF"),
                    labels = c("no_ccx" = "No Cachexia",
                               "ccx" = "Cachexia")) + 
  theme_std() + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.key.size = unit(0.25, "cm"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/normal_fluct_density_SAT.pdf", width = 2.5, height = 1.5)

p <- ggplot(combined_df, aes(x = pmax(pmin(delta_VAT, 50), -50), fill = type)) + 
  geom_density(alpha = 0.8, linewidth = 0.1) + 
  ylab("Density") + 
  xlab(expression("%"~Delta~VAT/month)) + 
  labs(fill = "") + 
  scale_fill_manual(values = c("#D4A6C8FF", "#8CD17DFF"),
                    labels = c("no_ccx" = "No Cachexia",
                               "ccx" = "Cachexia")) + 
  theme_std() + 
  scale_y_continuous(expand = c(0,0)) + 
  #scale_x_continuous(limits = c(-50, 50)) + 
  theme(legend.key.size = unit(0.25, "cm"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/normal_fluct_density_VAT.pdf", width = 2.5, height = 1.5)


p <- ggplot(combined_df, aes(x = pmax(pmin(delta_MuscleArea, 15), -15), fill = type)) + 
  geom_density(alpha = 0.8, linewidth = 0.1) + 
  ylab("Density") + 
  xlab(expression("%"~Delta~SKM/month)) + 
  labs(fill = "") + 
  scale_fill_manual(values = c("#D4A6C8FF", "#8CD17DFF"),
                    labels = c("no_ccx" = "No Cachexia",
                               "ccx" = "Cachexia")) + 
  theme_std() + 
  scale_y_continuous(expand = c(0,0)) + 
  #scale_x_continuous(limits = c(-0.015, 0.015)) + 
  theme(legend.key.size = unit(0.25, "cm"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/normal_fluct_density_SKM.pdf", width = 2.5, height = 1.5)

p <- ggplot(combined_df, aes(x = pmax(pmin(delta_SpleenVolume, 40), -40), fill = type)) + 
  geom_density(alpha = 0.8, linewidth = 0.1) + 
  ylab("Density") + 
  xlab(expression("%"~Delta~Spleen/month)) + 
  labs(fill = "") + 
  scale_fill_manual(values = c("#D4A6C8FF", "#8CD17DFF"),
                    labels = c("no_ccx" = "No Cachexia",
                               "ccx" = "Cachexia")) + 
  theme_std() + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.key.size = unit(0.25, "cm"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/normal_fluct_density_SpleenVolume.pdf", width = 2.5, height = 1.5)

p <- ggplot(combined_df, aes(x = pmax(pmin(delta_KidneyVolume, 20), -20), fill = type)) + 
  geom_density(alpha = 0.8, linewidth = 0.1) + 
  ylab("Density") + 
  xlab(expression("%"~Delta~Kidney/month)) + 
  labs(fill = "") + 
  scale_fill_manual(values = c("#D4A6C8FF", "#8CD17DFF"),
                    labels = c("no_ccx" = "No Cachexia",
                               "ccx" = "Cachexia")) + 
  theme_std() + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.key.size = unit(0.25, "cm"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/normal_fluct_density_KidneyVolume.pdf", width = 2.5, height = 1.5)

p <- ggplot(combined_df, aes(x = pmax(pmin(delta_LiverArea, 30), -30), fill = type)) + 
  geom_density(alpha = 0.8, linewidth = 0.1) + 
  ylab("Density") + 
  xlab(expression("%"~Delta~Liver/month)) + 
  labs(fill = "") + 
  scale_fill_manual(values = c("#D4A6C8FF", "#8CD17DFF"),
                    labels = c("no_ccx" = "No Cachexia",
                               "ccx" = "Cachexia")) + 
  theme_std() + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.key.size = unit(0.25, "cm"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/normal_fluct_density_LiverVolume.pdf", width = 2.5, height = 1.5)

p <- ggplot(combined_df, aes(x = pmax(pmin(delta_PancreasVolume, 30), -30), fill = type)) + 
  geom_density(alpha = 0.8, linewidth = 0.1) + 
  ylab("Density") + 
  xlab(expression("%"~Delta~Pancreas/month)) + 
  labs(fill = "") + 
  scale_fill_manual(values = c("#D4A6C8FF", "#8CD17DFF"),
                    labels = c("no_ccx" = "No Cachexia",
                               "ccx" = "Cachexia")) + 
  theme_std() + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.key.size = unit(0.25, "cm"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/normal_fluct_density_PancreasVolume.pdf", width = 2.5, height = 1.5)


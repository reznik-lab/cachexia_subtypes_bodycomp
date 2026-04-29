# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Sonia Boscenco
# Create histogram of all body composition changes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# filtering outliers for visualization purposes only 
upper                              <- quantile(bodycomp_deltas$delta_KidneyVolume, 0.99, na.rm = TRUE)
bodycomp_deltas$delta_KidneyVolume <- pmin(bodycomp_deltas$delta_KidneyVolume,upper)

upper                              <- quantile(bodycomp_deltas$delta_BMDLStandard, 0.99, na.rm = TRUE)
bodycomp_deltas$delta_BMDLStandard <- pmin(bodycomp_deltas$delta_BMDLStandard,upper)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data for potting
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
                 color = "black", linewidth = 0.1, bins = 20, stat = "bin") +
  scale_fill_manual(values = c("FALSE" = "#E15759", "TRUE" = "#4E79A7FF"), guide = "none") +
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

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/histograms_bodycompchanges_allorgans.pdf", width = 7, height = 2.25)

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

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/histograms_bodycompchanges_allorgansdensities.pdf", width = 5.5, height = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Report percentages for text
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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




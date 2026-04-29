# Author: Sonia Boscenco
# Create delta file from processed body composition data
# Must have a column name [SCAN_TYPE] that denotes Start-Cachexia and Post-Cachexia

rm(list = ls())

source("~/Desktop/reznik/bodycomp_main/analysis/processing/calculate_deltas_func.R")

bodycomp                      <- read.csv("~/Desktop/reznik/bodycomp_main/data/main_processed/processed_bodycomp_w_metadata_0220.csv")
bodycomp$BMI_CHANGE           <- bodycomp$CCX_END_BMI - bodycomp$CCX_START_BMI
bodycomp$SCAN_DATE            <- as.Date(bodycomp$SCAN_DATE)

# filters to make sure no pediatric cases
# TAT filter for anasarca
# ensuring that all patients lose weight
bodycomp_sub                  <- bodycomp %>% 
                                 dplyr::filter(SCAN_TYPE %in% c("Start-Cachexia", "Post-Cachexia")) %>%
                                 dplyr::filter(BMI_CHANGE < 0) %>%
                                 dplyr::filter(CURRENT_AGE_DEID >= 18) %>%
                                 dplyr::filter(L3FatValues.L3TATArea > 35)

# to run need wide format with SCAN_TYPE column denoting "Start-Cachexia" or "Post-Cachexia"
bodycomp_deltas               <- calculate_body_composition_deltas(bodycomp_sub)

# keep only one row per patient
bodycomp_deltas               <- bodycomp_deltas %>% distinct(MRN, .keep_all = TRUE)

# removing abnormally large values 
bodycomp_deltas               <- bodycomp_deltas %>% subset(delta_SAT < 500 & delta_VAT < 500)

curr_date                     <- format(Sys.Date(), "%m%d")

write.csv(bodycomp_deltas, file = paste0("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_", curr_date, ".csv"), row.names = FALSE)

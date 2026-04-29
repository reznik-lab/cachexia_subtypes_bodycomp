
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AUTHOR: Sonia Boscenco
# process bodycomposition raw file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
gc()

source("~/Desktop/reznik/bodycomp_main/analysis/processing/clean_bodycomp_data_func.R")
source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

longitidunal_key            <- read.csv("~/Desktop/reznik/bodycomp_main/data/master_processed/scan_level_master_key_w_allmetadata_0220.csv")
bodycomp                    <- read.csv("~/Desktop/reznik/bodycomp_main/data/raw_bodycomp_files//PANCAN_20-320P_bodycomp_20260220.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make sure everything is numeric and within physiological range
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bodycomp                   <- clean_bodycomp_data(bodycomp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pick the axial scan by ensuring the level exists
# and its not IV constrast 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bodycomp_processed         <- bodycomp %>%
                              dplyr::filter(Levels.L3Level > 0) %>%
                              dplyr::filter(grepl("PRIMARY_AXIAL", DICOMHeader.ImageType)) %>% 
                              dplyr::filter(IVContrast.IVContrastPresent == "True" | grepl("WITH CONTRAST|WITH IV CONTRAST| W/ CONTRAST | W/CONTRAST", DICOMHeader.ProtocolName, ignore.case = TRUE))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# if multiple axial series choose one with least missingness in key tissues
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bodycomp_processed          <- bodycomp_processed %>%
                               mutate(non_missing_count = rowSums(!is.na(across(c("MuscleValues.L3TotalMuscleMeanHU", 
                                                                                  "MuscleValues.L3TotalMuscleArea",
                                                                                  "MuscleValues.L3IMATArea",
                                                                                  "MuscleValues.L3IMATMeanHU",
                                                                                  "L3FatValues.L3SATArea",
                                                                                  "L3FatValues.L3VATArea"))))) %>%
                               group_by(StudyInfo.Accession) %>%
                               slice_max(non_missing_count, n = 1, with_ties = FALSE) %>%
                               ungroup() %>%
                               dplyr::select(-non_missing_count)

bodycomp_processed_wmetadata  <- merge(bodycomp_processed, longitidunal_key, by.y = "CASE_ID", by.x = "StudyInfo.Accession", all.x = TRUE)

write.csv(bodycomp_processed, file = "~/Desktop/reznik/bodycomp_main/data/master_processed//processed_bodycomp_0220.csv", row.names = FALSE)
write.csv(bodycomp_processed_wmetadata, file = "~/Desktop/reznik/bodycomp_main/data/master_processed/processed_bodycomp_w_metadata_0220.csv", row.names = FALSE)

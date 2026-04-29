# Identify progression windows during cachexia across clusters
rm(list = ls())
gc()

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bodycomp_metadata                   <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0302.csv")
progression                         <- read.csv("~/Downloads/table_timeline_radiology_cancer_progression_predictions.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# keep only progression is true events 
progression                   <- progression %>% 
                                      dplyr::filter(MRN %in% bodycomp_metadata$MRN) 

progression_true <- progression %>%
                                      filter(PROGRESSION == "Yes")

progression_true$Pre_CCX_Scan_Date <- bodycomp_metadata$START_CCX_SCAN_DATE[match(progression_true$MRN, bodycomp_metadata$MRN)]
progression_true$Post_CCX_Scan_Date<- bodycomp_metadata$END_CCX_SCAN_DATE[match(progression_true$MRN, bodycomp_metadata$MRN)]

# choose only events that happen during the cachexia window 
progression_true_during_ccx        <- progression_true %>% 
                                      filter(as.Date(START_DATE) >= as.Date(Pre_CCX_Scan_Date) & as.Date(START_DATE) <= as.Date(Post_CCX_Scan_Date)) 

# get cachexia length to normalize windows 
progression_true_during_ccx$length_ccx <- as.numeric(as.Date(progression_true_during_ccx$Post_CCX_Scan_Date) - as.Date(progression_true_during_ccx$Pre_CCX_Scan_Date)) / 30.44


# calculate number of events per month 
progression_true_during_ccx            <- progression_true_during_ccx %>% 
                                          group_by(MRN) %>% 
                                          mutate(n_events = n()) %>% 
                                          mutate(n_events_permonth = n_events/as.numeric(length_ccx)) %>% 
                                          distinct(MRN, .keep_all = TRUE)

progression_true_during_ccx$cluster    <- bodycomp_metadata$cluster_name[match(progression_true_during_ccx$MRN, bodycomp_metadata$MRN)]
progression_true_during_ccx$CancerType <- bodycomp_metadata$CANCER_TYPE_DETAILED[match(progression_true_during_ccx$MRN, bodycomp_metadata$MRN)]

progression_true_during_ccx$cluster    <- factor(progression_true_during_ccx$cluster, levels = c("Type C", "Type B", "Type A"))

progression_true_during_ccx$time       <- as.Date(progression_true_during_ccx$START_DATE) - as.Date(progression_true_during_ccx$Pre_CCX_Scan_Date)

# get only first event (15 day buffer because otherwise could coincide w/ start of cachexia )
first_event                            <- progression_true_during_ccx %>% 
                                          filter(time >= 15) %>% group_by(MRN) %>% 
                                          slice_min(time) %>% 
                                          distinct(MRN, .keep_all = TRUE)
bodycomp_w_event                       <- bodycomp_metadata %>% filter(MRN %in% progression$MRN)
bodycomp_w_event$Status                <- ifelse(bodycomp_w_event$MRN %in% first_event$MRN, 1, 0)
bodycomp_w_event$Time                  <- first_event$time[match(bodycomp_w_event$MRN, first_event$MRN)]
bodycomp_w_event$Time                  <- ifelse(is.na(bodycomp_w_event$Time), bodycomp_w_event$scans_days, bodycomp_w_event$Time)
bodycomp_w_event$Time                  <- as.numeric(bodycomp_w_event$Time) / 30.44
bodycomp_w_event$cluster_name          <- factor(bodycomp_w_event$cluster_name, levels = c("Type A", "Type B", "Type C"))

# within each cancer type 
for(ct in unique(bodycomp_w_event$CANCER_TYPE_DETAILED)){
  ct_df <- bodycomp_w_event %>% filter(CANCER_TYPE_DETAILED == ct)
  
  fit_km <- survfit(Surv(Time, Status) ~ cluster_name, data = ct_df)
  names(fit_km$strata)                    <- gsub("cluster_name=", "", names(fit_km$strata))
  p <- ggsurvplot(
    fit_km,
    fun = "event",
    pval = TRUE,
    pval.size = 2,
    size = 0.5,
    censor.size = 2,
    conf.int = TRUE,
    axes.offset = FALSE,
    ggtheme = theme_std(base_size = 7),
    xlab = "Time from cachexia onset (months)",
    ylab = "Cumulative incidence",
    legend.title = "",
    legend = "right",
    palette = (c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))
  )  
  p$plot <- p$plot + 
    ggtitle(ct) +
    theme(legend.key.size = unit(0.25, "cm"),
          plot.title = element_text(family = "ArialMT", face = "bold", size = 6, hjust = 0.5))
  
  outfile <- paste0("~/Desktop/reznik/bodycomp_main/results/cachexia/pancan_progression_km_0317/", ct, ".pdf")
  ggsave(outfile, p$plot, width = 3, height = 1.5, dpi = 300)
  print(ct)
}

  
  
  
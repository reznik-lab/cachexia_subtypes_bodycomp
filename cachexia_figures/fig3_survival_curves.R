# Author: Sonia Boscenco
# KM curves for cluster assignments 

rm(list = ls())
gc()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data and scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("~/Desktop/reznik/bodycomp_master/analysis/prerequisites.R")
master_file                          <- read.csv("~/Desktop/reznik/bodycomp_master/data/cachexia/cachexia_deltas_w_metdata_0302.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
master_file$Status                   <- ifelse(master_file$OS_STATUS == "1:DECEASED", 1, 0)
master_file$cluster_name             <- factor(master_file$cluster_name, levels = c("Type A", "Type B", "Type C"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# do it for every cancertype individually
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(cancertype in unique(master_file$CANCER_TYPE_DETAILED)){
  df <- master_file %>% filter(CANCER_TYPE_DETAILED == cancertype)
  fit <- survfit(Surv(OS_CCX, Status) ~ cluster_name, data = df)
  names(fit$strata)                    <- gsub("cluster_name=", "", names(fit$strata))
  p <- ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    pval.size = 2,
    size = 0.5,
    censor.size = 2,
    conf.int = TRUE,
    ggtheme = theme_std(base_size = 7),
    xlab = "Time (months)",
    legend = "right",
    legend.title = "",
    palette = (c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))
  )
  p$plot <- p$plot + 
    ggtitle(cancertype) +
    theme(legend.key.size = unit(0.25, "cm"),
                           plot.title = element_text(family = "ArialMT", face = "bold", size = 6, hjust = 0.5))
  outfile <- paste0("~/Desktop/reznik/bodycomp_master/results/cachexia/pancan_km_0311/", cancertype, ".pdf")
  ggsave(outfile, p$plot, width = 3, height = 1.5, dpi = 300)
  print(cancertype)
}
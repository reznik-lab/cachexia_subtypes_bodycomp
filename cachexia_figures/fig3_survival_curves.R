# Author: Sonia Boscenco
# KM curves for cluster assignments 
library(ggrastr)
rm(list = ls())
gc()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data and scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
master_file                          <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0828.csv")

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
    legend = "bottom",
    legend.title = "",
    palette = (c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))
  )
  p$plot <- p$plot + 
    ggtitle(cancertype) +
    theme(legend.key.size = unit(0.25, "cm"),
                           plot.title = element_text(family = "ArialMT", face = "bold", size = 6, hjust = 0.5))
  outfile <- paste0("~/Desktop/reznik/bodycomp_main/revision/results/pancan_km_0828/", cancertype, ".pdf")
  ggsave(outfile, p$plot, width = 1.75, height = 1.5, dpi = 300)
  print(cancertype)
}

pval_results <- data.frame(
  CANCER_TYPE_DETAILED = character(),
  n_patients = integer(),
  n_clusters = integer(),
  pval = numeric(),
  stringsAsFactors = FALSE
)

for(cancertype in unique(master_file$CANCER_TYPE_DETAILED)){
  df <- master_file %>% filter(CANCER_TYPE_DETAILED == cancertype)
  
    fit <- survfit(Surv(OS_CCX, Status) ~ cluster_name, data = df)
  
  sdiff <- survdiff(Surv(OS_CCX, Status) ~ cluster_name, data = df)
  pval <- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  
  pval_results <- rbind(pval_results, data.frame(
    CANCER_TYPE_DETAILED = cancertype,
    n_patients = nrow(df),
    n_clusters = nlevels(df$cluster_name),
    pval = pval,
    stringsAsFactors = FALSE
  ))
}

pval_results$pval_BH <- p.adjust(pval_results$pval, method = "BH")
pval_results <- pval_results %>% arrange(pval_BH)

format_qval <- function(q){
  if (is.na(q)) return("q = NA")
  if (q < 1e-4) return("q < 0.0001")
  sprintf("q = %.4f", q)
}


plot_list <- list()

for(cancertype in unique(master_file$CANCER_TYPE_DETAILED)){
  
  df <- master_file %>% filter(CANCER_TYPE_DETAILED == cancertype)

  qval <- pval_results$pval_BH[pval_results$CANCER_TYPE_DETAILED == cancertype]
  qlabel <- format_qval(qval)
  
  fit <- survfit(Surv(OS_CCX, Status) ~ cluster_name, data = df)
  names(fit$strata) <- gsub("cluster_name=", "", names(fit$strata))
  
  p <- ggsurvplot(
    fit,
    data = df,
    pval = qlabel,        
    pval.size = 2,
    size = 0.5,
    censor.size = 2,
    conf.int = TRUE,
    ggtheme = theme_std(),
    xlab = "Time (months)",
    legend = "bottom",
    legend.title = "",
    palette = c("#B07AA1FF", "#499894FF", "#A0CBE8FF")
  )
  
  p$plot <- p$plot +
    ggtitle(cancertype) +
    theme(legend.key.size = unit(0.25, "cm"),
          plot.title = element_text(family = "ArialMT", face = "bold", size = 6, hjust = 0.5))
  
  # Rasterize the line/CI/censor layers (keeps text, axes, legend as vector)
  p$plot$layers <- lapply(p$plot$layers, function(l) rasterise(l, dpi = 1200))
  
  plot_list[[cancertype]] <- p
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Save each plot at 1200 dpi
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(cancertype in names(plot_list)){
  outfile <- paste0("~/Desktop/reznik/bodycomp_main/revision/results/pancan_km_0920/", cancertype, ".pdf")
  ggsave(outfile, plot = plot_list[[cancertype]]$plot,
         width = 1.75, height = 1.5, dpi = 1200)
}

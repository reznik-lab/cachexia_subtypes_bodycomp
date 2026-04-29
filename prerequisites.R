## set seed for reproducibility of figures
set.seed(123)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load required packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

required.packages <- c('data.table','ggplot2', 'dplyr', 'tidyr', 'stringr', 'RColorBrewer', 'ggsignif', 'ggrepel', 'binom', 'reshape2', 'ggpubr',
                       'survival', 'survminer', 'speedglm', 'ComplexUpset', 'pls', "ggbeeswarm", 'broom', 'purrr', 'pheatmap', 'paletteer',
                       'tibble', 'pls', 'glmnet', 'ggstats', "DESeq2", "biomaRt", "ComplexHeatmap", "fgsea", "qusage", "circlize", "lmerTest", "tximport",
                       "lme4", "broom.mixed", "cmprsk", "patchwork", "chisq.posthoc.test", "emmeans", "multcomp", "ggbreak", "limma", "GSEABase", "edgeR",
                       "ggfortify", "ggbiplot", "VIM", "rstatix", "rsq", "org.Hs.eg.db")
hide <- suppressMessages(lapply(required.packages, require, character.only = TRUE))
missing.packages <- required.packages[!required.packages %in% (.packages())]
if(length(missing.packages)>0) stop(paste('Could not load required packages:',paste(missing.packages,collapse=', ')))

cols <- paletteer::paletteer_d("ggthemes::Tableau_20")
curr_date                     <- format(Sys.Date(), "%m%d")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ggplot theme
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_std <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = 'ArialMT')  %+replace%
    theme(
      line = element_line(colour = "black", linewidth = 0.1),
      axis.ticks = element_line(linewidth = 0.1),
      axis.ticks.length = unit(1, "pt"))
}


calc_pairwise_fc <- function(data, group_col, value_col, log2 = FALSE) {
  
    summary_df <- data %>%
    group_by(.data[[group_col]]) %>%
    summarise(mean_value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  
  pairs <- combn(summary_df[[group_col]], 2, simplify = FALSE)
  
  out <- lapply(pairs, function(p) {
    a <- p[1]
    b <- p[2]
    va <- summary_df$mean_value[summary_df[[group_col]] == a]
    vb <- summary_df$mean_value[summary_df[[group_col]] == b]
    
    fc_ab <- va / vb
    fc_ba <- vb / va
    
    tibble(
      group1 = a,
      group2 = b,
      mean1 = va,
      mean2 = vb,
      FC_1_vs_2 = if (log2) log2(fc_ab) else fc_ab,
      FC_2_vs_1 = if (log2) log2(fc_ba) else fc_ba
    )
  })
  
  bind_rows(out)
}

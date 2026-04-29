## what metrics best explain BMI 
rm(list = ls())
gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("~/Desktop/reznik/bodycomp_master/analysis/prerequisites.R")
bodycomp_deltas                        <- read.csv("~/Desktop/reznik/bodycomp_master/data/cachexia/cachexia_deltas_w_metdata_0302.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bodycomp_deltas$bmi_change             <- bodycomp_deltas$CCX_END_BMI - bodycomp_deltas$CCX_START_BMI
bodycomp_deltas$percent_bmi_change     <- (bodycomp_deltas$CCX_END_BMI - bodycomp_deltas$CCX_START_BMI) / bodycomp_deltas$CCX_START_BMI * 100
bodycomp_deltas                        <- bodycomp_deltas %>%
                                          mutate(fill_color = ifelse(percent_bmi_change < -10, "#E15759FF", "black"))
# get mean and sd values
mean(bodycomp_deltas$percent_bmi_change)
sd(bodycomp_deltas$percent_bmi_change)

bodycomp_deltas$cluster_name           <- factor(bodycomp_deltas$cluster_name, levels = c("Inflammatory-Wasted", "Atrophy-Wasted", "Non-Wasted"))
prop_df                                <- bodycomp_deltas %>%
                                          group_by(cluster_name) %>%
                                          summarise(prop_below_neg10 = mean(percent_bmi_change < -10), .groups = "drop")
# to colour the bars below -10 red
breaks                                 <- c(seq(0, -10, length.out = 5),  seq(-12.5, -60, length.out = 20))

p <- ggplot(bodycomp_deltas, aes(x = percent_bmi_change, fill = fill_color)) + 
  geom_histogram(breaks = breaks, colour = "white", linewidth = 0.1) + 
  facet_wrap(~cluster_name) + 
  theme_std() + 
  labs(x = "BMI change (%)", y = "Number of patients") + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_identity() + theme(
    panel.border = element_rect(color = "black", size = 0.1),
    strip.background = element_rect(linewidth  = 0.1),
    strip.text = element_text(size = 5)
  )

ggsave(p, file = "~/Desktop/reznik/bodycomp_master/results/BMI_change_props_hist_clusters.pdf", width = 3, height =1.5 )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# association of body comp changes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p <- ggplot(bodycomp_deltas_std, aes(x = delta_SAT, y = delta_MuscleArea)) + 
  geom_point(alpha = 0.5, stroke = NA, aes(colour = cluster_name)) + 
  facet_wrap(~cluster_name, scale = "free") + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) + 
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.1) + 
  stat_cor(size = 2, font="ArialMT", label.x.npc = "center") +
  xlab(expression(Delta*"SAT(%)")) +
  ylab(expression(Delta*"SKM(%)")) +
  theme(
    panel.border = element_rect(color = "black", size = 0.1),
    strip.background = element_rect(linewidth  = 0.1),
    strip.text = element_text(size = 5),
    legend.position = "none"
  ) + 
  scale_colour_manual(values = c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))

ggsave(p, file = "~/Desktop/reznik/bodycomp_master/results/cachexia/satvsmuscle_clusters.pdf", width = 4, height = 2)
p <- ggplot(bodycomp_deltas_std, aes(x = delta_SAT, y = delta_VAT)) + 
  geom_point(alpha = 0.5, stroke = NA, aes(colour = cluster_name)) + 
  facet_wrap(~cluster_name, scale = "free") + 
  theme_std() + 
  xlab(expression(Delta*"SAT(%)")) +
  ylab(expression(Delta*"VAT(%)")) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) + 
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.1) + 
  stat_cor(size = 2, font="ArialMT", label.x.npc = "center") +
  theme(
    panel.border = element_rect(color = "black", size = 0.1),
    strip.background = element_rect(linewidth  = 0.1),
    strip.text = element_text(size = 5),
    legend.position = "none"
  ) + 
  scale_colour_manual(values = c("#B07AA1FF", "#499894FF", "#A0CBE8FF"))
ggsave(p, file = "~/Desktop/reznik/bodycomp_master/results/cachexia/satvsvat_clusters.pdf", width = 4, height = 2)





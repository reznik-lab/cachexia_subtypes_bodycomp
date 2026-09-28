# Author: Sonia Boscenco
# Boxplots of wasted vs. non wasted patients

rm(list = ls())
gc()

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bodycomp_metadata                   <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0828.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bodycomp_metadata$cluster_name      <- factor(bodycomp_metadata$cluster_name, levels = c("Type A", "Type B", "Type C"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# muscle area
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 <- pairwise.t.test(bodycomp_metadata$delta_MuscleArea, bodycomp_metadata$cluster_name, p.adjust.method = "none", pool.sd = FALSE)
bodycomp_metadata %>% cohens_d(delta_MuscleArea ~ cluster_name)

comparisons                  <- list(c("Type B", "Type C"), c("Type A", "Type B"), c("Type A", "Type C"))
p <- ggplot(bodycomp_metadata,  aes(x = cluster_name, y = pmin(delta_MuscleArea , 40))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ylab(expression(Delta*"(%)")) +
  theme(legend.position = "none",
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right =  element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  ggtitle("SKM") +
  #scale_y_break(breaks = c(41, 120), ticklabels = 112) +
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT") +
 # stat_compare_means(comparisons = list(c("Atrophy-Wasted", "Non-Wasted"), c("Inflammatory-Wasted", "Atrophy-Wasted"), c("Inflammatory-Wasted", "Non-Wasted")), linewidth = 0.1, size = 2, label = "p.signif", tip.length = 0)+
  scale_y_continuous(n.breaks = 6)
 
ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_muscle_clusters.pdf", width =1.75, height =2, units = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sat
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_SAT, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_SAT ~ cluster_name)

p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmin(delta_SAT, 100))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  labs(x = "") + 
  ggtitle("SAT") +
  ylab(expression(Delta*"(%)")) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_SAT_wastednonwasted.pdf", width =1.75, height =2, units = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# vat
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_VAT, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_VAT ~ cluster_name)

p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmin(delta_VAT, 100))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(%)")) +
  ggtitle("VAT") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) +
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_vat_wastednonwasted.pdf", width =1.74, height =2, units = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# skm density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_MuscleDensity, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_MuscleDensity ~ cluster_name)

p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmax(delta_MuscleDensity, -40))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(HU)")) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  ggtitle("SKM Density") +
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(method.args = list(var.equal = TRUE), comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")
ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_skmdensity_wastednonwasted.pdf", width =1.75, height =2, units = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  IMAT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_IMAT, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_IMAT ~ cluster_name)

p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmin(delta_IMAT, 100))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(%)")) +
  ggtitle("IMAT") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_imat_wastednonwasted.pdf", width = 1.75, height =2, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  liver
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_LiverArea, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_LiverArea ~ cluster_name)

p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmin(delta_LiverArea, 100))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(%)")) +
  ggtitle("Liver Volume") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_liver_wastednonwasted.pdf", width =1.75, height = 2, units = "in", dpi = 300)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  liver density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_LiverDensity,bodycomp_metadata$cluster_name,p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_LiverDensity~ cluster_name)

p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = delta_LiverDensity)) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(HU)")) +
  ggtitle("Liver Density") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_liverdensity_wastednonwasted.pdf", width =1.75, height = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# kidney
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_KidneyVolume, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_KidneyVolume~ cluster_name)
p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmin(delta_KidneyVolume, 75))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(%)")) +
  ggtitle("Kidney Volume") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  #scale_y_break(breaks = c(255, 355), ticklabels = 355) +
  #scale_y_break(breaks = c(380, 420), ticklabels = 420) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_kidneyvolume_wastednonwasted.pdf", width =1.75, height = 2, unit = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# kidney density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_KidneyDensity, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_KidneyDensity~ cluster_name)
p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = delta_KidneyDensity)) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(HU)")) +
  ggtitle("Kidney Density") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_kidneydensity_wastednonwasted.pdf", width = 1.75, height = 2, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# spleen 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_SpleenVolume, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_SpleenVolume~ cluster_name)
p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmin(delta_SpleenVolume, 100))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(%)")) +
  ggtitle("Spleen Volume") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_spleenvolume_wastednonwasted.pdf", width = 1.75, height = 2, unit = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# spleen density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_SpleenDensity, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_SpleenDensity~ cluster_name)
p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmax(pmin(delta_SpleenDensity, 100), -100))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(HU)")) +
  ggtitle("Spleen Density") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_spleendensity_wastednonwasted.pdf", width = 1.75, height = 2, unit = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pancreas volume
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_PancreasVolume, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_PancreasVolume~ cluster_name)
p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmin(delta_PancreasVolume, 100))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(%)")) +
  ggtitle("Pancreas Volume") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_pancreasvolume_wastednonwasted.pdf", width = 1.75, height = 2, unit = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pancreas density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_PancreasVolume, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_PancreasDensity~ cluster_name)
p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = pmin(pmax(delta_PancreasDensity, -100), 100))) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(HU)")) +
  ggtitle("Pancreas Density") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_pancreasdensity_wastednonwasted.pdf", width = 1.75, height = 2, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAT Density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_SATDensity, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_SATDensity~ cluster_name)
p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = delta_SATDensity)) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(HU)")) +
  ggtitle("SAT Density") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_satdensity_wastednonwasted.pdf", width =1.75, height = 2, units = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VAT Density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pairwise.t.test(bodycomp_metadata$delta_VATDensity, bodycomp_metadata$cluster_name, p.adjust.method = "BH")
bodycomp_metadata %>% cohens_d(delta_VATDensity~ cluster_name)
p <- ggplot(bodycomp_metadata, aes(x = cluster_name, y = delta_VATDensity)) + 
  geom_quasirandom(alpha = 0.1, aes(colour = cluster_name), width = 0.3, stroke = NA) + 
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, width = 0.5, fill = NA) + 
  theme_std() + 
  labs(x = "") + 
  ylab(expression(Delta*"(HU)")) +
  ggtitle("VAT Density") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.1) +
  theme(legend.position = "none",
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold", family = "ArialMT")) + 
  scale_colour_manual(values = c("#B07AA1FF","#499894FF", "#A0CBE8FF")) + 
  stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "BH", label = "p.format", step.increase = 0.1, bracket.size = 0.1, tip.length = 0, size = 2 , family = "ArialMT")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/boxplot_vatdensity_wastednonwasted.pdf", width =1.75, height = 2, units = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sex barplot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tbl_sex <- (as.data.frame(table(bodycomp_metadata$cluster_name, bodycomp_metadata$GENDER)))
colnames(tbl_sex) <- c("CLUSTER", "SEX", "FREQ")

tbl_sex <- tbl_sex %>% group_by(SEX) %>% mutate(n  = sum(FREQ), prop = FREQ/n) %>% ungroup()
tbl_sex$SEX <- str_to_title(tbl_sex$SEX)
chisq_test <- chisq.test(as.matrix(table(bodycomp_metadata$cluster_name, bodycomp_metadata$GENDER), ncol = 2))
chisq_test$p.value

p <- ggplot(tbl_sex, aes(x = SEX, y = prop, fill = factor(CLUSTER, levels = rev(c("Type A", "Type B", "Type C"))))) + 
  geom_bar(stat = "identity", linewidth = 0.1, colour = "white") + 
  theme_std()+
  coord_flip() +
  ggtitle(paste0("Chisq-Test P = ", round(chisq_test$p.value, 3))) +
  scale_fill_manual(values = rev(c("#B07AA1FF","#499894FF", "#A0CBE8FF"))) + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Proportion of patients", fill = "") +
  theme(axis.ticks.y = element_blank(),
        plot.title = element_text(family = "ArialMT", size = 6, hjust = 0.5, face = "italic"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "bottom",
        axis.line.y = element_blank())

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision//results//barplot_sex_wastednonwasted.pdf", width = 1.5, height = 1.25)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# bmi barplot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bodycomp_metadata                                   <- bodycomp_metadata %>%
                                                        mutate(CCX_START_WEIGHT_CATEGORY = case_when(
                                                          CCX_START_BMI < 18.5 ~ "Underweight",
                                                          CCX_START_BMI >= 18.5 & CCX_START_BMI < 25 ~ "Normal",
                                                          CCX_START_BMI >= 25 & CCX_START_BMI < 30 ~ "Overweight",
                                                          CCX_START_BMI >= CCX_START_BMI ~ "Obese"))

tbl_bmi <- (as.data.frame(table(bodycomp_metadata$cluster_name, bodycomp_metadata$CCX_START_WEIGHT_CATEGORY)))
colnames(tbl_bmi) <- c("CLUSTER", "BMI", "FREQ")
tbl_bmi <- tbl_bmi %>% group_by(BMI) %>% mutate(n  = sum(FREQ), prop = FREQ/n) %>% ungroup()
tbl_bmi$BMI <- str_to_title(tbl_bmi$BMI)
chisq_test <- chisq.test(as.matrix(table(bodycomp_metadata$cluster_name, bodycomp_metadata$weight_category), ncol = 2))
chisq_test$p.value

p <- ggplot(tbl_bmi, aes(x = factor(BMI, levels = c("Underweight", "Normal", "Overweight", "Obese")), y = prop, fill = factor(CLUSTER, levels = rev(c("Type A", "Type B", "Type C"))))) + 
  geom_bar(stat = "identity", linewidth = 0.1, colour = 'white', width = 0.95) + 
  theme_std() + 
  coord_flip() +
  scale_fill_manual(values = rev(c("#B07AA1FF","#499894FF", "#A0CBE8FF"))) + 
  ggtitle(paste0("Chisq-Test P = ",round(chisq_test$p.value, 3))) +
  labs(y = "Proportion of patients", x = "",fill = "") + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) + 
  theme(axis.ticks.y = element_blank(),
        plot.title = element_text(family = "ArialMT", size = 6, hjust = 0.5, face = "italic"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "bottom",
        axis.line.y = element_blank())
ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision/results//barplot_bmi_wastednonwasted.pdf", width = 1.75, height = 1.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cancertype barplot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#systematically test for enrichement in each cancer type 
cancertypes <- unique(bodycomp_metadata$CANCER_TYPE_DETAILED)

# 1 vs. the rest
results                             <- data.frame()
for(cancer in cancertypes){
  sub_df                            <- bodycomp_metadata %>% mutate(is_cancer = ifelse(CANCER_TYPE_DETAILED == cancer, "Y", "N"))
  tbl                               <- table(sub_df$cluster_name, sub_df$is_cancer)
  print(tbl)
  posthoc                           <- chisq.posthoc.test(tbl, method = "bonferroni")
  for(cluster in rownames(tbl)){
    pvalue_sub                      <- posthoc %>%
                                       filter(Value == "p values") %>%
                                       filter(Dimension == cluster)
    
    a <- tbl[cluster, "Y"]                
    b <- tbl[cluster, "N"]                 
    c <- sum(tbl[-which(rownames(tbl) == cluster), "Y"]) 
    d <- sum(tbl[-which(rownames(tbl) == cluster), "N"])  
    
    OR <- (a * d) / (b * c)
    
    
    results                         <- rbind(results, data.frame(cancer_type = cancer,
                                                                 cluster = cluster, 
                                                                 pval = pvalue_sub$Y,
                                                                 OR = OR))
  }
}
results$p.adj <- p.adjust(results$pval, method = "BH")

sign_results                         <- results %>%
                                        filter(p.adj < 0.05)
tbl_cancertype <- (as.data.frame(table(bodycomp_metadata$cluster_name, bodycomp_metadata$CANCER_TYPE_DETAILED)))
colnames(tbl_cancertype) <- c("CLUSTER", "CancerType", "FREQ")

tbl_cancertype <- tbl_cancertype %>% group_by(CancerType) %>% mutate(n  = sum(FREQ), prop = FREQ/n)  %>% ungroup()

cancerlist <-tbl_cancertype %>%
  filter(CLUSTER == "Type C") %>%
  arrange((prop)) %>% 
  pull(CancerType)

tbl_cancertype$CancerType <- factor(tbl_cancertype$CancerType, levels = unique(cancerlist))
p <- ggplot(tbl_cancertype, aes(x = CancerType, y = prop, fill = factor(CLUSTER, levels = rev(c("Type A", "Type B", "Type C"))))) + 
  geom_bar(stat = "identity", linewidth = 0.1, colour = "white") + 
  scale_fill_manual(values = rev(c("#B07AA1FF","#499894FF", "#A0CBE8FF"))) +
  coord_flip() +
  labs(x = "", y = "Proportion of patients", fill = "") +
  theme_std() +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.ticks.y = element_blank(),
        legend.key.size = unit(0.25, "cm"))
ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision//results//barplot_cancertype_wasted_vsnonwasted.pdf", width = 3.75, height = 2.25)

write.csv(results, file = "~/Desktop/reznik/bodycomp_main/revision//tables/supplementary_4_cluster_cancertypeenrichment.csv", row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# stage barplot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tbl_stage <- (as.data.frame(table(bodycomp_metadata$cluster_name, bodycomp_metadata$STAGE_CCX)))
colnames(tbl_stage) <- c("CLUSTER", "STAGE", "FREQ")

tbl_stage <- tbl_stage %>% group_by(STAGE) %>% mutate(n  = sum(FREQ), prop = FREQ/n) %>% ungroup()
tbl_stage$STAGE <- str_to_title(tbl_stage$STAGE)
chisq_test <- chisq.test(as.matrix(table(bodycomp_metadata$cluster_name, bodycomp_metadata$STAGE_CCX), ncol = 2))
chisq_test$p.value

p <- ggplot(tbl_stage, aes(x = STAGE, y = prop, fill = factor(CLUSTER, levels = rev(c("Type A", "Type B", "Type C"))))) + 
  geom_bar(stat = "identity", linewidth = 0.1, colour = "white") + 
  theme_std()+
  coord_flip() +
  ggtitle(paste0("Chisq-Test P = ", round(chisq_test$p.value, 3))) +
  scale_fill_manual(values = rev(c("#B07AA1FF","#499894FF", "#A0CBE8FF"))) + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Proportion of patients", fill = "") +
  theme(axis.ticks.y = element_blank(),
        plot.title = element_text(family = "ArialMT", size = 6, hjust = 0.5, face = "italic"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "bottom",
        axis.line.y = element_blank())

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision/results//barplot_stage.pdf", width =2.1, height = 1)


### pairwise comparisons df 
delta_vars <- c("delta_MuscleDensity", "delta_SATDensity", "delta_VATDensity",
                "delta_KidneyDensity", "delta_PancreasDensity", "delta_SpleenDensity",
                "delta_LiverDensity", "delta_LiverArea", "delta_KidneyVolume",
                "delta_PancreasVolume", "delta_SpleenVolume", "delta_BMDLStandard",
                "delta_MuscleArea", "delta_SAT", "delta_VAT", "delta_IMAT")

all_pairwise_results <- map_dfr(delta_vars, function(var) {
  bodycomp_metadata %>%
    pairwise_t_test(
      as.formula(paste(var, "~ cluster_name")),
      p.adjust.method = "none",
      pool.sd = FALSE
    ) %>%
    mutate(variable = var)
})


all_pairwise_results$p.adj.all <- p.adjust(all_pairwise_results$p, method = "BH")

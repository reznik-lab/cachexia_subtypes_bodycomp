# Author: Sonia Boscenco
# Manual vs. Automatic Segmentation Scatterplots
# Figure 1 

rm(list = ls())
gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
segmentations      <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/crc_matched_oscar.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p                  <- ggplot(segmentations, aes(x = Auto_SAT, Manual_SAT)) + 
                      geom_point(alpha = 0.7, size = 1, stroke = NA) + 
                      theme_std() + 
                      labs(x = "OSCAR SAT", y = "Manual SAT") + 
                      stat_cor(label.x.npc = 0.4, family = "ArialMT", size = 2)

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/sat_auto_v_manual.pdf", width = 1.5, height = 1.5)

p1                  <- ggplot(segmentations, aes(x = Auto_VAT, Manual_VAT)) + 
                       geom_point(alpha = 0.7, size = 1, stroke = NA) + 
                       theme_std() + 
                       labs(x = "OSCAR VAT", y = "Manual VAT") + 
                       stat_cor(label.x.npc = 0.5, family = "ArialMT", size = 2)

ggsave(p1, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/vat_auto_v_manual.pdf", width = 1.5, height = 1.5)

p2                  <- ggplot(segmentations, aes(x = Auto_Muscle, Manual_Muscle)) + 
                       geom_point(alpha = 0.7, size = 1, stroke = NA) + 
                       theme_std() + 
                       labs(x = "OSCAR SKM", y = "Manual SKM") + 
                       stat_cor(label.x.npc = 0.5, family = "ArialMT", size = 2)

ggsave(p2, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/muscle_auto_v_manual.pdf", width = 1.5, height = 1.5)
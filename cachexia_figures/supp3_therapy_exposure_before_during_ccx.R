
rm(list = ls())
gc()

source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load files 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bodycomp_metadata                   <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0828.csv")
therapy                             <- read.csv("~/Desktop/reznik/bodycomp_main/data/metadata/ddp_chemo.tsv", sep = "\t")

therapy                             <- therapy %>%
                                       filter(MRN %in% bodycomp_metadata$MRN)

therapy$CDG_THERAPEUTIC_CATEGORY    <- trimws(therapy$CDG_THERAPEUTIC_CATEGORY)
therapy$PRE_CCX_DATE                <- bodycomp_metadata$START_CCX_SCAN_DATE[match(therapy$MRN, bodycomp_metadata$MRN)]
therapy$POST_CCX_DATE               <- bodycomp_metadata$END_CCX_SCAN_DATE[match(therapy$MRN, bodycomp_metadata$MRN)]
therapy$CLUSTER                     <- bodycomp_metadata$cluster_name[match(therapy$MRN, bodycomp_metadata$MRN)]
therapy$during_ccx                  <- ifelse(as.Date(therapy$APR_START_DTE) >= (as.Date(therapy$PRE_CCX_DATE)) & as.Date(therapy$APR_END_DTE) <= as.Date(therapy$POST_CCX_DATE), TRUE, FALSE)
therapy$before_ccx                  <- ifelse(as.Date(therapy$APR_START_DTE) < (as.Date(therapy$PRE_CCX_DATE)), TRUE, FALSE)

therapy$chemo                       <- ifelse(therapy$CDG_CHEMO_FLAG == "Y", 1, 0)
therapy$hormone                     <- ifelse(therapy$CDG_HORMONE_FLAG == "Y", 1, 0)
therapy$targeted                    <- ifelse(therapy$CDG_TARGETED_FLAG == "Y", 1, 0)
therapy$immuno                      <- ifelse(therapy$CDG_IMMUNOTHERAPY_FLAG == "Y", 1, 0)

exposure_before                     <- therapy %>%
                                       filter(before_ccx == TRUE) 


chemo                               <- exposure_before %>%
  group_by(MRN) %>%
  summarise(n_chemo = sum(chemo))
chemo$n_chemo                       <- ifelse(chemo$n_chemo >= 1, 1, 0)
chemo                               <- chemo %>% 
  filter(n_chemo == 1)

hormone                             <- exposure_before %>%
  group_by(MRN) %>%
  summarise(n_hormone = sum(hormone))

hormone$n_hormone                   <- ifelse(hormone$n_hormone >= 1, 1, 0)
hormone                             <- hormone %>% 
  filter(n_hormone == 1)

targeted                            <- exposure_before %>%
  group_by(MRN) %>%
  summarise(n_targeted = sum(targeted))

targeted$n_targeted                 <- ifelse(targeted$n_targeted >= 1, 1, 0)
targeted                            <- targeted %>% 
  filter(n_targeted == 1)

immuno                              <- exposure_before %>%
  group_by(MRN) %>%
  summarise(n_immuno = sum(immuno))

immuno$n_immuno                    <- ifelse(immuno$n_immuno >= 1, 1, 0)
immuno                             <- immuno %>% 
  filter(n_immuno == 1)

all_prior                          <- merge(chemo, hormone, by = "MRN", all = TRUE)
all_prior                          <- merge(all_prior, targeted, by = "MRN", all = TRUE)
all_prior                          <- merge(all_prior, immuno, by = "MRN", all = TRUE)

all_prior$cluster                  <- bodycomp_metadata$cluster_name[match(all_prior$MRN, bodycomp_metadata$MRN)]
colnames(all_prior)               <- c("MRN", "chemo_before", "hormone_before", "targeted_before", "immuno_before")
write.csv(all_prior, file = "~/Desktop/reznik/bodycomp_main/revision/tables/therapy_before.csv", row.names = FALSE)
chemo_tbl                          <- as.data.frame(table(all_prior$n_chemo, all_prior$cluster)) %>%
  select(-c(Var1))
colnames(chemo_tbl)                <- c("CLUSTER", "Chemo")

hormone_tbl                        <- as.data.frame(table(all_prior$n_hormone, all_prior$cluster)) %>%
  select(-c(Var1))
colnames(hormone_tbl)              <- c("CLUSTER", "Hormone")

targeted_tbl                       <- as.data.frame(table(all_prior$n_targeted, all_prior$cluster)) %>%
  select(-c(Var1))
colnames(targeted_tbl)             <- c("CLUSTER", "Targeted")

immuno_tbl                         <- as.data.frame(table(all_prior$n_immuno, all_prior$cluster)) %>%
  select(-c(Var1))
colnames(immuno_tbl)               <- c("CLUSTER", "Immuno")

all_prior_tbl                      <- merge(chemo_tbl, hormone_tbl, by = "CLUSTER")
all_prior_tbl                      <- merge(all_prior_tbl, targeted_tbl, by = "CLUSTER")
all_prior_tbl                      <- merge(all_prior_tbl, immuno_tbl, by = "CLUSTER")

table_all                          <- as.data.frame(table(bodycomp_metadata$cluster_name))
colnames(table_all)                <- c("CLUSTER", "N")

all_prior_tbl                      <- merge(all_prior_tbl, table_all, by = "CLUSTER")

long_tbl                           <- all_prior_tbl %>% pivot_longer(cols = c("Chemo", "Hormone", "Targeted", "Immuno"))
long_tbl$prop                      <- long_tbl$value / long_tbl$N
#long_tbl$name                      <- factor(long_tbl$name , levels = rev(c("Chemo", "Immuno", "Targeted", "Hormone")))
p <- ggplot(long_tbl, aes(x = CLUSTER, y = prop)) + 
  geom_bar(stat = "identity", fill = "black") +
  theme_std() +
  labs(x = "", y = "Prop. patients exposed prior ccx") + 
  theme(legend.position = "bottom",
        axis.ticks.x = element_blank(),
        strip.background = element_blank()) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~name, nrow = 1, scales = "free_y")
  
ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision/figures/prop_patients_exposed_prior_treatment.pdf", width = 4.5, height = 1.5)


pairs <- combn(unique(long_tbl$CLUSTER), 2, simplify = FALSE)

results <- map_dfr(unique(long_tbl$name), function(drug) {
  map_dfr(pairs, function(p) {
    g1 <- long_tbl %>% filter(CLUSTER == p[1], name == drug)
    g2 <- long_tbl %>% filter(CLUSTER == p[2], name == drug)
    
    tab <- matrix(c(g1$value, g1$N - g1$value,
                    g2$value, g2$N - g2$value),
                  nrow = 2, byrow = TRUE)
    
    test   <- prop.test(tab)
    fisher <- fisher.test(tab)
    
    tibble(
      drug = drug,
      group1 = p[1],
      group2 = p[2],
      prop1 = round(g1$value / g1$N, 4),
      prop2 = round(g2$value / g2$N, 4),
      chi_sq = round(unname(test$statistic), 3),
      p_value = test$p.value,
      OR = round(unname(fisher$estimate), 3),
      OR_lower = round(fisher$conf.int[1], 3),
      OR_upper = round(fisher$conf.int[2], 3)
    )
  })
})

results <- results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         sig_adj = case_when(
           p_adj < 0.001 ~ "***",
           p_adj < 0.01  ~ "**",
           p_adj < 0.05  ~ "*",
           TRUE          ~ "ns"
         ))

exposure_during                     <- therapy %>%
  filter(during_ccx == TRUE) 


chemo                               <- exposure_during %>%
  group_by(MRN) %>%
  summarise(n_chemo = sum(chemo))
chemo$n_chemo                       <- ifelse(chemo$n_chemo >= 1, 1, 0)
chemo                               <- chemo %>% 
  filter(n_chemo == 1)

hormone                             <- exposure_during %>%
  group_by(MRN) %>%
  summarise(n_hormone = sum(hormone))

hormone$n_hormone                   <- ifelse(hormone$n_hormone >= 1, 1, 0)
hormone                             <- hormone %>% 
  filter(n_hormone == 1)

targeted                            <- exposure_during %>%
  group_by(MRN) %>%
  summarise(n_targeted = sum(targeted))

targeted$n_targeted                 <- ifelse(targeted$n_targeted >= 1, 1, 0)
targeted                            <- targeted %>% 
  filter(n_targeted == 1)

immuno                              <- exposure_during %>%
  group_by(MRN) %>%
  summarise(n_immuno = sum(immuno))

immuno$n_immuno                    <- ifelse(immuno$n_immuno >= 1, 1, 0)
immuno                             <- immuno %>% 
  filter(n_immuno == 1)

all_prior                          <- merge(chemo, hormone, by = "MRN", all = TRUE)
all_prior                          <- merge(all_prior, targeted, by = "MRN", all = TRUE)
all_prior                          <- merge(all_prior, immuno, by = "MRN", all = TRUE)
colnames(all_prior)               <- c("MRN", "chemo_before", "hormone_before", "targeted_before", "immuno_before")
write.csv(all_prior, file = "~/Desktop/reznik/bodycomp_main/revision/tables/therapy_during.csv", row.names = FALSE)
all_prior$cluster                  <- bodycomp_metadata$cluster_name[match(all_prior$MRN, bodycomp_metadata$MRN)]

chemo_tbl                          <- as.data.frame(table(all_prior$n_chemo, all_prior$cluster)) %>%
  select(-c(Var1))
colnames(chemo_tbl)                <- c("CLUSTER", "Chemo")

hormone_tbl                        <- as.data.frame(table(all_prior$n_hormone, all_prior$cluster)) %>%
  select(-c(Var1))
colnames(hormone_tbl)              <- c("CLUSTER", "Hormone")

targeted_tbl                       <- as.data.frame(table(all_prior$n_targeted, all_prior$cluster)) %>%
  select(-c(Var1))
colnames(targeted_tbl)             <- c("CLUSTER", "Targeted")

immuno_tbl                         <- as.data.frame(table(all_prior$n_immuno, all_prior$cluster)) %>%
  select(-c(Var1))
colnames(immuno_tbl)               <- c("CLUSTER", "Immuno")

all_prior_tbl                      <- merge(chemo_tbl, hormone_tbl, by = "CLUSTER")
all_prior_tbl                      <- merge(all_prior_tbl, targeted_tbl, by = "CLUSTER")
all_prior_tbl                      <- merge(all_prior_tbl, immuno_tbl, by = "CLUSTER")

table_all                          <- as.data.frame(table(bodycomp_metadata$cluster_name))
colnames(table_all)                <- c("CLUSTER", "N")

all_prior_tbl                      <- merge(all_prior_tbl, table_all, by = "CLUSTER")

long_tbl                           <- all_prior_tbl %>% pivot_longer(cols = c("Chemo", "Hormone", "Targeted", "Immuno"))
long_tbl$prop                      <- long_tbl$value / long_tbl$N
#long_tbl$name                      <- factor(long_tbl$name , levels = rev(c("Chemo", "Immuno", "Targeted", "Hormone")))
p <- ggplot(long_tbl, aes(x = CLUSTER, y = prop)) + 
  geom_bar(stat = "identity", fill = "black") +
  theme_std() +
  labs(x = "", y = "Prop. patients exposed during ccx") + 
  theme(legend.position = "bottom",
        axis.ticks.x = element_blank(),
        strip.background = element_blank()) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~name, nrow = 1, scales = "free_y")

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision/figures/prop_patients_exposed_during_treatment_new.pdf", width = 4.5, height = 1.5)


pairs <- combn(unique(long_tbl$CLUSTER), 2, simplify = FALSE)

results <- map_dfr(unique(long_tbl$name), function(drug) {
  map_dfr(pairs, function(p) {
    g1 <- long_tbl %>% filter(CLUSTER == p[1], name == drug)
    g2 <- long_tbl %>% filter(CLUSTER == p[2], name == drug)
    
    tab <- matrix(c(g1$value, g1$N - g1$value,
                    g2$value, g2$N - g2$value),
                  nrow = 2, byrow = TRUE)
    
    test   <- prop.test(tab)
    fisher <- fisher.test(tab)
    
    tibble(
      drug = drug,
      group1 = p[1],
      group2 = p[2],
      prop1 = round(g1$value / g1$N, 4),
      prop2 = round(g2$value / g2$N, 4),
      chi_sq = round(unname(test$statistic), 3),
      p_value = test$p.value,
      OR = round(unname(fisher$estimate), 3),
      OR_lower = round(fisher$conf.int[1], 3),
      OR_upper = round(fisher$conf.int[2], 3)
    )
  })
})

results <- results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         sig_adj = case_when(
           p_adj < 0.001 ~ "***",
           p_adj < 0.01  ~ "**",
           p_adj < 0.05  ~ "*",
           TRUE          ~ "ns"
         ))
results <- results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),  # or "bonferroni"
         sig_adj = case_when(
           p_adj < 0.001 ~ "***",
           p_adj < 0.01  ~ "**",
           p_adj < 0.05  ~ "*",
           TRUE          ~ "ns"
         ))


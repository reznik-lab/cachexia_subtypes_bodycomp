# Author: Sonia Boscenco
# Create clusters for body composition data

rm(list = ls())
gc()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data and scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)
source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
bodycomp_deltas               <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_0302.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# remove aortic calculation 
bodycomp_deltas                         <- bodycomp_deltas %>%
                                           dplyr::select(-c("delta_TotalAortic"))

delta_cols                              <- colnames(bodycomp_deltas)[grepl("delta_", colnames(bodycomp_deltas))]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# impute missing values and scale
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bodycomp_deltas_knn                     <- VIM::kNN(bodycomp_deltas,k = 10, imp_var = FALSE, variable = delta_cols)
rownames(bodycomp_deltas_knn)           <- bodycomp_deltas_knn$MRN
bodycomp_deltas_knn                     <- bodycomp_deltas_knn[grepl("delta_", colnames(bodycomp_deltas_knn))]
bodycomp_deltas_knn_scaled              <- scale(bodycomp_deltas_knn, center = TRUE, scale = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# kmeans clustering
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# find optimal number of clusters
wss                                     <- sapply(1:10, function(k) {kmeans(bodycomp_deltas_knn_scaled, nstart = 50, centers = k)$tot.withinss})
wss_df                                  <- data.frame(K = 1:10, wss = wss)

# elbow plot 
p <- ggplot(wss_df, aes(x = K, y = wss)) +
  geom_point(shape = 19, size = 0.5) +
  geom_line(linewidth = 0.1) +
  scale_x_continuous(breaks = wss_df$K) +
  labs(x = "Number of clusters K",
       y = "Total within-clusters\nsum of squares") +
  theme_std()
ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/kmeans_elbowplot.pdf", width = 2 , height = 1.5)

# now we cluster with k = 3!
km                                      <- kmeans(bodycomp_deltas_knn_scaled, centers = 3, nstart = 50)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add annotation to the data frame
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cluster_df                                <- as.data.frame(km$cluster)
colnames(cluster_df)                      <- c("cluster")
cluster_df$MRN                            <- as.numeric(rownames(cluster_df))
cluster_df$cluster_name                   <- ifelse(cluster_df$cluster == 1, "Type A",
                                                    ifelse(cluster_df$cluster == 2, "Type B", "Type C"))
bodycomp_deltas$cluster                   <- cluster_df$cluster[match(bodycomp_deltas$MRN, cluster_df$MRN)]
bodycomp_deltas$cluster_name              <- cluster_df$cluster_name[match(bodycomp_deltas$MRN, cluster_df$MRN)]


write.csv(cluster_df, file = paste0("~/Desktop/reznik/bodycomp_main/data/cachexia/cluster_assignments_mrn_", curr_date, ".csv"), row.names = FALSE)
write.csv(bodycomp_deltas, file = paste0("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_cluster_", curr_date, ".csv"), row.names = FALSE)

# RNAseq liver biopsies 
rm(list = ls())
gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data and scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
pancreatic_subset2                  <- read.csv("~/Desktop/reznik/bodycomp_main/data/metadata/pancreatic_liver_biopsies_clusters.csv")
bulk                                <- read.csv("~/Downloads/GSE245535_counts_scaled_DESeq_GEO.txt", sep = "\t", skip = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# gene expression
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pancreatic_subset2$metabo_id        <- substr(pancreatic_subset2$UniqueID, 1, nchar(pancreatic_subset2$UniqueID) -2)
rownames(bulk)              <- bulk$GeneID

# fix ensembl ids
numbers                     <- gsub("s_([0-9]+).*", "\\1", colnames(bulk))
colnames(bulk)              <- numbers
numbers_we_have             <- numbers[numbers %in% pancreatic_subset2$metabo_id]
bulk_sub                    <- bulk[colnames(bulk) %in% pancreatic_subset2$metabo_id]
entrez_ids                  <- rownames(bulk_sub)
gene_ids_clean              <- sub("\\..*$", "", entrez_ids)
rownames(bulk_sub)          <- gene_ids_clean
gene_ids_clean              <- gene_ids_clean[grepl("^ENS", gene_ids_clean)]

ensembl                     <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping                     <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = gene_ids_clean, mart = ensembl)
      
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# limma
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samples                     <- data.frame(condition = pancreatic_subset2$cluster_name,
                                          row = pancreatic_subset2$metabo_id,
                                          row.names = pancreatic_subset2$metabo_id)

samples_order               <- samples[colnames(bulk_sub), ]
counts_matrix               <- as.matrix(log2(bulk_sub + 1))

samples_order$condition     <- relevel(factor(samples_order$condition), ref = "Non-Wasted")
design                      <- model.matrix(~ condition, data = samples_order)

fit                         <- lmFit(counts_matrix, design)
fit                         <- eBayes(fit, robust = TRUE, trend = TRUE)

res_atrophy                 <- topTable(fit, coef = "conditionAtrophy-Wasted", number = Inf, adjust.method = "BH") %>% 
                               as.data.frame() %>% 
                               rownames_to_column("entrezgene_id") %>%
                               mutate(clean_id = gsub("\\..*$", "", entrezgene_id)) 

res_hepaton                 <- topTable(fit, coef = "conditionInflammatory-Wasted", number = Inf, adjust.method = "BH") %>% 
                               as.data.frame() %>% 
                               rownames_to_column("entrezgene_id") %>%
                               mutate(clean_id = gsub("\\..*$", "", entrezgene_id)) 

res_atrophy$gene             <- mapping$hgnc_symbol[match(res_atrophy$clean_id, mapping$ensembl_gene_id)]
res_hepaton$gene             <- mapping$hgnc_symbol[match(res_hepaton$clean_id, mapping$ensembl_gene_id)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# atrophy gsea
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gs                           <- read.gmt("~/Desktop/reznik/bodycomp_main/data/reference/h.all.v2025.1.Hs.symbols.gmt")

ranks_atrophy                <- res_atrophy[!is.na(res_atrophy$gene), ]
ranks_atrophy                <- ranks_atrophy[!is.na(ranks_atrophy$P.Value), ]
ranks_atrophy                <- ranks_atrophy[ranks_atrophy$gene != "", ]
ranks_atrophy_list           <- ranks_atrophy$logFC
names(ranks_atrophy_list)    <- ranks_atrophy$gene
ranks_atrophy_list           <- sort(ranks_atrophy_list, decreasing = TRUE)
resgsea                      <- fgsea(pathways = gs, stats = ranks_atrophy_list, minSize = 15, maxSize = 500)
resgsea$padj                 <- as.numeric(resgsea$padj)
resgsea                      <- resgsea[order(resgsea$padj, decreasing = FALSE), ]

resgsea$comparison           <- c("Type B - Type C")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# inflammatory gsea
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ranks_hepato               <- res_hepaton[!is.na(res_hepaton$gene), ]
res_hepaton                <- res_hepaton[!is.na(res_hepaton$P.Value), ]
res_hepaton                <- res_hepaton[res_hepaton$gene != "", ]
res_hepaton_list           <- res_hepaton$logFC
names(res_hepaton_list)    <- res_hepaton$gene

res_hepaton_list           <- sort(res_hepaton_list, decreasing = TRUE)
resgsea_hepato             <- fgsea(pathways = gs, stats = res_hepaton_list, minSize = 15, maxSize = 500)
resgsea_hepato$padj        <- as.numeric(resgsea_hepato$padj)
resgsea_hepato             <- resgsea_hepato[order(resgsea_hepato$padj, decreasing = FALSE), ]
resgsea_hepato$comparison  <- c("Type A - Type C")

all_gsea                   <- rbind(resgsea, resgsea_hepato) %>% as.data.frame()
all_gsea$leadingEdge       <- sapply(all_gsea$leadingEdge, paste, collapse=",")
write.csv(all_gsea, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/tables/supp_liver_gsea.csv", row.names = FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot gsea
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clean_hallmark_names <- function(df){
  df$cleaned_term    <- df$pathway |>
                        gsub("HALLMARK_", "", x = _) |>
                        gsub("_", " ", x = _) |>
                       str_to_title()
  
   df$cleaned_term   <- ifelse(is.na(hall_mark_equivs[df$cleaned_term]), df$cleaned_term, hall_mark_equivs[df$cleaned_term])
}

hall_mark_equivs                              <- c("E2f Targets" = "E2F Targets",
                                                   "P53 Pathway" = "p53 Pathway",
                                                   "Mtorc1 Signaling" = "mTORC1 Signaling",
                                                   "Il6 Jak Stat3 Signaling" = "IL-6/JAK/STAT3 Signaling",
                                                   "Uv Response Up" = "UV Response Up",
                                                   "Uv Response Dn" = "UV Response Dn",
                                                   "Tgf Beta Signaling" = "TGF-beta Signaling",
                                                   "Dna Repair" = "DNA Repair",
                                                   "Tgf Beta Signaling" = "TGF-beta Signaling",
                                                   "Kras Signaling Up" = "KRAS Signaling Up",
                                                   "Il2 Stat5 Signaling" = "IL-2/STAT5 Signaling",
                                                   "G2m Checkpoint" = "G2-M Checkpoint",
                                                   "Tnfa Signaling Via Nfkb" = "TNF-alpha Signaling via NF-kB",
                                                   "Pi3k Akt Mtor Signaling" = "PI3K/AKT/mTOR Signaling",
                                                   "Wnt Beta Catenin Signaling" = "Wnt-beta Catenin Signaling")


resgsea$cleaned_term                    <- clean_hallmark_names(resgsea)
resgsea_hepato$cleaned_term             <- clean_hallmark_names(resgsea_hepato)

resgsea$colour                          <- ifelse(resgsea$NES > 0, "red", "blue")
resgsea$psigcolour                      <- ifelse(resgsea$padj < 0.05 & resgsea$NES > 0, "blue",
                                                  ifelse(resgsea$padj < 0.05 & resgsea$NES < 0, "red", "black"))

resgsea_hepato$psigcolour                      <- ifelse(resgsea_hepato$padj < 0.05 & resgsea_hepato$NES > 0, "blue",
                                                  ifelse(resgsea_hepato$padj < 0.05 & resgsea_hepato$NES < 0, "red", "black"))

volcano_atrophy <- ggplot(resgsea, aes(x = NES, y = -log10(padj), colour = psigcolour)) + 
  geom_point(stroke = NA, alpha = 0.8) +
  theme_std() +
  theme(legend.position = "none") +
  ylab(expression(-log[10](q))) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 18.6)) +
  scale_colour_manual(values = c("black", "#4E79A7FF", "#E15759FF")) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.1) + 
  geom_text_repel(aes(label = cleaned_term), family = "ArialMT", size = 1, segment.size = 0.1, data = (filter(resgsea, cleaned_term %in% c("E2F Targets",  "TNF-alpha Signaling via NF-kB", "Epithelial Mesenchymal Transition", "Oxidative Phosphorylation"))))

ggsave(volcano_atrophy, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/gsea_atrophy_liver_volcano.pdf", width = 2, height = 1.5)

volcano_inflam <- ggplot(resgsea_hepato, aes(x = NES, y = -log10(padj), colour = psigcolour)) + 
  geom_point(stroke = NA, alpha = 0.8) +
  theme_std() +
  theme(legend.position = "none") +
  ylab(expression(-log[10](q))) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 18.6)) +
  scale_colour_manual(values = c("black", "#4E79A7FF", "#E15759FF")) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.1) + 
  geom_text_repel(aes(label = cleaned_term), family = "ArialMT", size = 1, segment.size = 0.1, data = (filter(resgsea_hepato, cleaned_term %in% c("E2F Targets",  "TNF-alpha Signaling via NF-kB", "Inflammatory Response"))))

ggsave(volcano_inflam, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/gsea_inflam_liver_volcano.pdf", width = 2, height = 1.5)

atrophy_gsea <- ggplot(filter(resgsea, padj < 0.05), aes(x = reorder(cleaned_term, NES), y = NES, fill = colour)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_std() +
  scale_fill_manual(values = c("#4E79A7FF", "#E15759FF")) +
  theme(legend.position = "none",
        axis.ticks.y = element_blank()) + 
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "") 

resgsea_hepato$colour                   <- ifelse(resgsea_hepato$NES > 0, "red", "blue")
inflam_gsea <- ggplot(filter(resgsea_hepato, padj < 0.05), aes(x = reorder(cleaned_term, NES), y = NES, fill = colour)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_std() +
  scale_fill_manual(values = c("#4E79A7FF", "#E15759FF")) +
  theme(legend.position = "none",
        axis.ticks.y = element_blank()) + 
  labs(x = "") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

ggsave(atrophy_gsea, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/gsea_atrophy_liver.pdf", width = 2.75, height = 2)
ggsave(inflam_gsea, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/gsea_inflam_liver.pdf", width = 2.25, height = 1.5)


p <- ggplot(filter(resgsea, padj < 0.05 & NES > 0), aes(x = reorder(cleaned_term, NES), y = NES, colour = padj, size = size)) + 
  geom_point() + 
  coord_flip() + 
  theme_std() +
  scale_size(range = c(1, 4)) +
  scale_colour_gradient(high = "#4E79A7FF", low = "#E15759FF") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 7, family = "ArialMT"),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.25, "cm")) + 
  ggtitle("Upregulated in Atrophy") +
  labs(x = "", size = "Gene count", colour = "q value") 

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/gsea_atrophy_liver.pdf", width = 3.5, height = 2)

p <- ggplot(filter(resgsea, padj < 0.05 & NES < 0), aes(x = reorder(cleaned_term, NES), y = NES, colour = padj, size = size)) + 
  geom_point() + 
  coord_flip() + 
  theme_std() +
  scale_colour_gradient(high = "#4E79A7FF", low = "#E15759FF") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 7, family = "ArialMT"),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.5, "cm")) + 
  ggtitle("Downregulated in Atrophy") +
  labs(x = "", size = "Gene count", colour = "q value") 

p <- ggplot(filter(resgsea_hepato, padj < 0.05), aes(x = reorder(cleaned_term, NES), y = NES, colour = padj, size = size)) + 
  geom_point() + 
  coord_flip() + 
  theme_std() +
  scale_size(range = c(1, 4)) +
  scale_colour_gradient(high = "#4E79A7FF", low = "#E15759FF") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 7, family = "ArialMT"),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.25, "cm")) + 
  ggtitle("Upregulated in Inflammatory") +
  labs(x = "", size = "Gene count", colour = "q value") 


ggsave(p, file = "~/Desktop/reznik/bodycomp_main/results/cachexia/gsea_inflam_liver.pdf", width = 3.5, height = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# heatmap of bodycomps for patients
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clusters                               <- as.data.frame(as.factor(bodycomp_deltas$cluster_name)) 
clusters$`as.factor(bodycomp_deltas$cluster_name)` <- factor(clusters$`as.factor(bodycomp_deltas$cluster_name)`, levels = c("Inflammatory-Wasted", "Atrophy-Wasted", "Non-Wasted"))
rownames(clusters)                     <- bodycomp_deltas$MRN

delta_values                           <- bodycomp_deltas %>% 
  dplyr::select(starts_with(("delta_"))) %>%
  as.matrix()
rownames(delta_values)                 <- bodycomp_deltas$MRN

f1                                     <- colorRamp2(seq(-50, 50, length = 3), c("#4E79A7FF","white","#E15759FF"))

clean_labels                           <- gsub("delta_", "", colnames(delta_values))
clean_labels                           <- gsub("(?<=\\w)(Density)", " \\1", clean_labels, perl = TRUE)
clean_labels                           <- gsub("(?<=\\w)(Volume)", " \\1", clean_labels, perl = TRUE)
clean_labels                           <- gsub("Area", " Volume", clean_labels, perl = TRUE)
clean_labels                           <- gsub("Muscle Volume", "SKM", clean_labels, perl = TRUE)
clean_labels                           <- gsub("BMDLStandard", " Bone Mineral Density", clean_labels, perl = TRUE)

mat                 <- t(delta_values)
row_labels          <- clean_labels  
column_labels       <- rownames(delta_values) 

ht <- Heatmap(
  (mat),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,     
  show_column_names = FALSE, 
  row_names_side = "left",
  cluster_column_slices = FALSE, 
  row_labels = row_labels,
  column_labels = column_labels,
  col = f1,
  row_dend_gp = gpar(lwd = 0.1),
  row_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  row_title_gp = gpar(fontsize = 6, fontfamily = "ArialMT", fontface = "bold"),
  column_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontfamily = "ArialMT"),
    labels_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
    legend_height = unit(5, "mm"), 
    legend_width = unit(20, "mm"),     
    grid_height = unit(2, "mm"),      
    grid_width  = unit(2, "mm"),
    direction = 'horizontal',
    title = expression(Delta~Value),
    title_position = "topcenter"
  ),
  show_column_dend = FALSE,
  show_parent_dend_line = FALSE,
  show_row_dend = FALSE,
  # row split must now be indexed by patients (columns)
  column_split = factor(clusters$`as.factor(bodycomp_deltas$cluster_name)`, labels = c("Type A", "Type B", "Type C")),
  column_title_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  column_gap = unit(2, "mm"),
  height = unit(nrow(mat) * 2.5, "mm"),
  width = unit(30, "mm"),
  border = "black",
  border_gp = gpar(lwd = 0.1)
)

pdf(file = "~/Desktop/reznik/bodycomp_main/results/cachexia/liver_pdac_heatmap_bodycompchanges.pdf", width = 2, height = 3)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

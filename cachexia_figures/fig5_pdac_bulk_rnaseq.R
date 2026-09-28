# Author: Sonia Boscenco
# RNA-Seq for pancreatic patients

rm(list = ls())
gc()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data and scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("~/Desktop/reznik/bodycomp_main/analysis/prerequisites.R")
source("~/Desktop/reznik/bodycomp_main/analysis/cachexia/rnaseq/run_gsea_func.R")

clusters                     <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0828.csv")
key                          <- read.csv("~/Desktop/reznik/bodycomp_main/data/rnaseq/pancreatic/Cachexia_SP_6.30.25.csv")
counts                       <- read.csv("~/Desktop/reznik/bodycomp_main/data/rnaseq/pancreatic/counts_by_gene_cachexia_39_WP.txt", sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
bodycomp_deltas              <- read.csv("~/Desktop/reznik/bodycomp_main/data/cachexia/cachexia_deltas_w_metdata_0828.csv")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# some gene names are empty of have been saved incorrectly
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
counts                       <- counts %>% 
                                filter(Gene != "")
counts                       <- counts %>% filter(!grepl("-", Gene))
rownames(counts)             <- counts$Gene

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subset for overlap
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
key_sub                      <- key %>% 
                                filter(MRN %in% clusters$MRN)

key_sub                      <- merge(key_sub, clusters, by = "MRN")
key_sub$Path_date            <- as.Date(as.character(key_sub$Path_date), format = "%m/%d/%y")
key_sub$time_to_ccx          <- as.Date(key_sub$START_CCX_SCAN_DATE) - key_sub$Path_date
key_sub$time_to_end          <- as.Date(key_sub$END_CCX_SCAN_DATE) - key_sub$Path_date
key_sub                      <- key_sub %>% filter(abs(time_to_ccx) < 365)
counts_sub                   <- counts[(colnames(counts) %in% key_sub$PanstudyID)]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run gsea
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# gene volcano plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples                <- data.frame(condition = key_sub$cluster_name,
                                     row = key_sub$X,
                                     sex = key_sub$Sex,
                                     sample_site = key_sub$Organ,
                                     row.names = key_sub$X)

samples_order           <- samples[colnames(counts_sub), ]
dds                     <- DESeqDataSetFromMatrix(countData = counts_sub,
                                                  colData = samples_order,
                                                  design = ~ sex + condition)
dds$condition           <- relevel(dds$condition, ref = "Type C")
dds                     <- DESeq(dds)
vsd                     <- vst(dds, blind = TRUE)  
plotPCA(vsd, intgroup = "sample_site")

res_raw                 <- results(dds, alpha = 0.05, contrast = c("condition", "Type A", "Type C"))
res                     <- lfcShrink(dds, contrast = c("condition", "Type A", "Type C"),
                                     res = res_raw, type = "ashr")
res$stat                <- res_raw$stat
res <- as.data.frame(res)
res                          <- res %>%
                                filter(!is.na(padj))
res$gene                     <- rownames(res)
res$color                    <- ifelse(res$log2FoldChange >1 & res$padj < 0.05, "red", 
                                       ifelse(res$log2FoldChange < -1 & res$padj < 0.05, "blue", "grey"))
sum(res$padj < 0.05 & res$log2FoldChange > 0)
sum(res$padj < 0.05 & res$log2FoldChange < 0)

labels_up                    <- res %>%
                                dplyr::filter(padj < 0.05) %>%
                                slice_max(n = 10, order_by = log2FoldChange) %>%
                                pull(gene)
labels_down                  <- res %>%
                                dplyr::filter(padj < 0.05) %>%
                                slice_min(n = 10, order_by = log2FoldChange) %>%
                                pull(gene)

labels                        <- c(labels_up, labels_down)
labels = c("HBB", "HBA2", "HBA1", "HBD", "ROBO1", "IGFN1", "FABP4", "CCL21", "CD247", "MCAM")
p <- ggplot(res, aes(x = log2FoldChange, y =-log10(padj), colour = color)) + 
     geom_point(alpha = 0.5, stroke = NA) + 
    theme_std() + 
    #scale_y_continuous(expand = c(0,0), limits = c(0, 5.5)) + 
    scale_x_continuous(limits = c(-15,15)) + 
    geom_hline(yintercept = -log10(0.05), linewidth = 0.1, linetype = "dashed") + 
    geom_vline(xintercept = -1, linewidth = 0.1, linetype = "dashed") + 
    geom_vline(xintercept = 1, linewidth = 0.1, linetype = "dashed") +
    scale_colour_manual(values = c("#4E79A7", "#BAB0AC", "#E15759")) + 
    theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", family = "ArialMT", size = 6)) + 
    ylab(expression(-log[10](q))) +
    xlab(expression(log[2](Fold-Change))) + 
    ggtitle("PDAC Bulk RNA (n = 9)") +
    geom_text_repel(aes(label = gene), data = subset(res, gene %in% labels), size = 1.5, family = "ArialMT", segment.size = 0.1)

ggsave(p, file = "~/Desktop/reznik/bodycomp_main//revision/main_figures/pdac_bulk_volcano_rnaseq.pdf", width = 2.5, height = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# hallmark barplot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

res$gene                <- rownames(res)
sub_ranks               <- res[!is.na(res$gene), ]
sub_ranks               <- sub_ranks[!is.na(res$padj), ]
sub_ranks               <- sub_ranks[sub_ranks$gene != "", ]

sub_ranks_max           <- sub_ranks %>% as.data.frame() %>%
  group_by(gene) %>%
  slice_max(stat, n = 1) %>%  
  ungroup()

ranks                   <- sub_ranks_max$stat
names(ranks)            <- sub_ranks_max$gene

ranks                   <- sort(ranks, decreasing = TRUE)
gmtfile = "~/Desktop/reznik/bodycomp_main/data/reference/h.all.v2025.1.Hs.symbols.gmt"
gs                      <- read.gmt(gmtfile)
resgsea                 <- fgsea(pathways = gs, stats = ranks, minSize = 15, maxSize = 500)
resgsea$padj            <- as.numeric(resgsea$padj)
resgsea                 <- resgsea[order(resgsea$padj, decreasing = FALSE), ]
#resgsea$leadingEdge     <- as.character(resgsea$leadingEdge)
resgsea$leadingEdge     <- sapply(resgsea$leadingEdge, function(x) paste(unlist(x), collapse = ";"))
write.csv(resgsea, file = "~/Desktop/reznik/bodycomp_main/revision/tables/pdac_bulk_gsea_hallmark_typeatypec.csv", row.names = FALSE)
write.csv(res, file = "~/Desktop/reznik/bodycomp_main/revision/tables/pdac_bulk_deseq2_hallmark_typeatypec.csv", row.names = FALSE)
                                
resgsea$colour           <- ifelse(resgsea$NES > 0 & resgsea$padj < 0.05, "up",
                                  ifelse(resgsea$NES < 0 & resgsea$padj < 0.05, "down", "gray"))


clean_hallmark_names <- function(df){
  df$cleaned_term    <- df$pathway |>
    gsub("HALLMARK_", "", x = _) |>
    gsub("_", " ", x = _) |>
    str_to_title()
  
  df$cleaned_term   <- ifelse(is.na(hall_mark_equivs[df$cleaned_term]), df$cleaned_term, hall_mark_equivs[df$cleaned_term])
}

hall_mark_equivs                              <- c("E2f Targets" = "E2F Targets",
                                                   "P53 Pathway" = "p53 Pathway",
                                                   "Mtorc1 Signaling" = "mTORC1 Signalling",
                                                   "Il6 Jak Stat3 Signaling" = "IL-6/JAK/STAT3 Signalling",
                                                   "Uv Response Up" = "UV Response Up",
                                                   "Uv Response Dn" = "UV Response Dn",
                                                   "Tgf Beta Signaling" = "TGF-beta Signal'ing",
                                                   "Dna Repair" = "DNA Repair",
                                                   "Tgf Beta Signaling" = "TGF-beta Signaling",
                                                   "Kras Signaling Up" = "KRAS Signaling Up",
                                                   "Il2 Stat5 Signaling" = "IL-2/STAT5 Signalling",
                                                   "G2m Checkpoint" = "G2-M Checkpoint",
                                                   "Tnfa Signaling Via Nfkb" = "TNF-alpha Signalling via NF-kB",
                                                   "Pi3k Akt Mtor Signaling" = "PI3K/AKT/mTOR Signalling",
                                                   "Wnt Beta Catenin Signaling" = "Wnt-beta Catenin Signalling",
                                                   "Oxphos" = "OXPHOS",
                                                   "Oxidative Phosphorylation" = "OXPHOS",
                                                   "Epithelial Mesenchymal Transition" = "EMT",
                                                   "Emt" = "EMT",
                                                   "Inf-Γ Response" = "Interferon Gamma Response")


resgsea$cleaned_term                    <- clean_hallmark_names(resgsea)

p <- ggplot(resgsea, aes(x = NES, y = -log10(padj), colour = colour)) + 
  geom_point(alpha = 0.9, stroke = NA) + 
  theme_std() +
  scale_colour_manual(values = c("#4E79A7", "#BAB0AC", "#E15759")) + 
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.1, linetype = "dashed") + 
  geom_vline(xintercept = 0, linewidth = 0.1, linetype = "dashed") + 
  scale_y_continuous(expand = c(0,0)) + 
  ggtitle("PDAC Bulk RNA (n = 9)") +
  geom_text_repel(data = subset(resgsea, cleaned_term %in% c("Bile Acid Metabolism", "Inflammatory Response", "IL-6/JAK/STAT3 Signalling", "OXPHOS", "Adipogenesis", "EMT", "Interferon Gamma Response")), aes(label = cleaned_term), size = 1.5, family = "ArialMT", segment.size = 0.1) +
  labs(x = "Normalized enrichment score") +
  ylab(expression(-log[10](q))) + 
  theme(
        plot.title = element_text(hjust = 0.5, face = "bold", family = "ArialMT", size = 6))

ggsave(p, file = "~/Desktop/reznik/bodycomp_main/revision/main_figures/pdac_bulk_gsea_hallmark_barplot.pdf",width = 2.5, height = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# bodycomp heatmap 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bodycomp_deltas                        <- bodycomp_deltas %>% filter(MRN %in% key_sub$MRN)
clusters                               <- as.data.frame(as.factor(bodycomp_deltas$cluster_name)) 
clusters$`as.factor(bodycomp_deltas$cluster_name)` <- factor(clusters$`as.factor(bodycomp_deltas$cluster_name)`, levels = c("Type A", "Type C"))
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
clean_labels                           <- gsub("Muscle Density", "SKM Density", clean_labels, perl = TRUE)
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
  column_split = factor(clusters$`as.factor(bodycomp_deltas$cluster_name)`, labels = c("Type A", "Type C")),
  column_title_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
  column_gap = unit(2, "mm"),
  height = unit(nrow(mat) * 2.5, "mm"),
  width = unit(30, "mm"),
  border = "black",
  border_gp = gpar(lwd = 0.1)
)

pdf(file = "~/Desktop/reznik/bodycomp_main/revision//main_figures//bulk_pdac_heatmap_bodycompchanges.pdf", width = 2, height = 3)
draw(ht, heatmap_legend_side = "bottom")
dev.off()


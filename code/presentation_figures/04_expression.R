library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(ggrepel)

outdir <- "figures/presentation_results/expression"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

res <- read_csv("analyses/07_deseq2/deseq2_results_chr3.csv", show_col_types = FALSE)

gene_col <- names(res)[names(res) %in% c("gene", "Geneid", "gene_id", "id", "X1", "...1")][1]

if (is.na(gene_col)) {
  res <- res %>% mutate(gene_label = paste0("gene_", row_number()))
} else {
  res <- res %>% mutate(gene_label = .data[[gene_col]])
}

res_plot <- res %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    status = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

top_labels <- res_plot %>%
  filter(status != "Not Sig") %>%
  arrange(padj) %>%
  slice_head(n = 10)

p1 <- ggplot(res_plot, aes(log2FoldChange, neg_log10_padj)) +
  geom_point(aes(colour = status), alpha = 0.75, size = 2) +
  geom_text_repel(
    data = top_labels,
    aes(label = gene_label),
    colour = "black",
    size = 4,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.size = 0.3
  ) +
  scale_colour_manual(
    values = c(
      "Down" = "steelblue",
      "Not Sig" = "grey75",
      "Up" = "firebrick"
    ),
    breaks = c("Down", "Not Sig", "Up")
  ) +
  theme_classic(base_size = 16) +
  labs(
    title = "Differential expression under heat treatment",
    x = "log2FC",
    y = "-log10(Pvalue)",
    colour = ""
  )

ggsave(file.path(outdir, "volcano_heat_response.png"), p1, width = 8, height = 6, dpi = 300)

sig_summary <- res_plot %>%
  filter(status != "Not Sig") %>%
  count(status)

p2 <- ggplot(sig_summary, aes(status, n, fill = status)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("Down" = "steelblue", "Up" = "firebrick")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Significant differentially expressed genes",
    x = "",
    y = "Number of genes"
  )

ggsave(file.path(outdir, "significant_de_genes.png"), p2, width = 6, height = 5, dpi = 300)

message("Expression figures saved to: ", outdir)

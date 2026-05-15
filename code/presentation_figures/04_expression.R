library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggrepel)

outdir <- "figures/presentation_results/expression"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

res <- read_csv("analyses/07_deseq2/deseq2_results_chr3.csv", show_col_types = FALSE)

# Detect gene name/id column
gene_col <- names(res)[names(res) %in% c("gene", "Geneid", "gene_id", "id", "X1")][1]

if (is.na(gene_col)) {
  res <- res %>% mutate(gene_label = row_number())
} else {
  res <- res %>% mutate(gene_label = .data[[gene_col]])
}

res_plot <- res %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Label top 10 most significant genes
top_labels <- res_plot %>%
  filter(significance != "Not significant") %>%
  arrange(padj) %>%
  slice_head(n = 10)

p1 <- ggplot(res_plot, aes(log2FoldChange, -log10(padj), colour = significance)) +
  geom_point(alpha = 0.65, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_labels,
    aes(label = gene_label),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_colour_manual(
    values = c(
      "Not significant" = "grey70",
      "Upregulated" = "firebrick",
      "Downregulated" = "steelblue"
    )
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Differential expression under heat treatment",
    subtitle = "Top significant genes are labelled",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    colour = ""
  )

ggsave(file.path(outdir, "volcano_heat_response.png"), p1, width = 8, height = 5, dpi = 300)

sig_summary <- res_plot %>%
  filter(significance != "Not significant") %>%
  count(significance)

p2 <- ggplot(sig_summary, aes(significance, n, fill = significance)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_manual(
    values = c(
      "Upregulated" = "firebrick",
      "Downregulated" = "steelblue"
    )
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Significant differentially expressed genes",
    x = "",
    y = "Number of genes"
  )

ggsave(file.path(outdir, "significant_de_genes.png"), p2, width = 6, height = 5, dpi = 300)

p3 <- ggplot(res_plot, aes(baseMean, log2FoldChange, colour = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(
    values = c(
      "Not significant" = "grey70",
      "Upregulated" = "firebrick",
      "Downregulated" = "steelblue"
    )
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "MA plot for heat treatment response",
    x = "Mean normalized expression",
    y = "log2 fold change",
    colour = ""
  )

ggsave(file.path(outdir, "MAplot_heat_response.png"), p3, width = 8, height = 5, dpi = 300)

top_genes <- res_plot %>%
  filter(significance != "Not significant") %>%
  arrange(padj) %>%
  select(gene_label, any_of(c("baseMean", "log2FoldChange", "padj", "pvalue", "significance"))) %>%
  head(20)

write_csv(top_genes, file.path(outdir, "top_20_DE_genes.csv"))

message("Expression figures saved to: ", outdir)

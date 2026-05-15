library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

outdir <- "figures/presentation_results/expression"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

res <- read_csv("analyses/07_deseq2/deseq2_results_chr3.csv", show_col_types = FALSE)

res_plot <- res %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    significant = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

p1 <- ggplot(res_plot, aes(log2FoldChange, -log10(padj), colour = significant)) +
  geom_point(alpha = 0.65, size = 1.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw(base_size = 13) +
  labs(
    title = "Differential expression under heat treatment",
    subtitle = "Dashed lines show |log2FC| = 1 and adjusted p-value = 0.05",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    colour = ""
  )

ggsave(file.path(outdir, "volcano_heat_response.png"), p1, width = 8, height = 5, dpi = 300)

sig_summary <- res_plot %>%
  filter(padj < 0.05) %>%
  count(significant)

p2 <- ggplot(sig_summary, aes(significant, n, fill = significant)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  theme_bw(base_size = 13) +
  labs(
    title = "Significant differentially expressed genes",
    x = "",
    y = "Number of genes"
  )

ggsave(file.path(outdir, "significant_de_genes.png"), p2, width = 6, height = 5, dpi = 300)

p3 <- ggplot(res_plot, aes(baseMean, log2FoldChange, colour = significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 13) +
  labs(
    title = "MA plot for heat treatment response",
    x = "Mean normalized expression",
    y = "log2 fold change",
    colour = ""
  )

ggsave(file.path(outdir, "MAplot_heat_response.png"), p3, width = 8, height = 5, dpi = 300)

top_genes <- res_plot %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  select(any_of(c("gene", "Geneid", "id", "baseMean", "log2FoldChange", "padj"))) %>%
  head(20)

write_csv(top_genes, file.path(outdir, "top_20_DE_genes.csv"))

copy_if_exists <- function(from, to_dir) {
  if (file.exists(from)) file.copy(from, to_dir, overwrite = TRUE)
}

copy_if_exists("figures/PCA_chr3.pdf", outdir)
copy_if_exists("figures/MAplot_chr3.pdf", outdir)
copy_if_exists("figures/deseq2/PCA_chr3.pdf", outdir)
copy_if_exists("figures/deseq2/MAplot_chr3.pdf", outdir)

message("Expression presentation figures saved to: ", outdir)

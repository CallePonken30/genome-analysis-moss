library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

outdir <- "figures/presentation_results/chloroplast"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cp_quast <- "analyses/09_chloroplast/quast/report.tsv"

if (file.exists(cp_quast)) {
  cp <- read_tsv(cp_quast, show_col_types = FALSE)
} else {
  cp <- tibble(
    Assembly = "Chloroplast",
    `Total length` = 123717,
    N50 = 123717,
    `GC (%)` = 28.27,
    `# N's per 100 kbp` = 0
  )
}

cp_long <- cp %>%
  pivot_longer(-1, names_to = "metric", values_to = "value") %>%
  filter(metric %in% c("Total length", "N50", "GC (%)", "# N's per 100 kbp")) %>%
  mutate(value = parse_number(as.character(value)))

p1 <- ggplot(cp_long, aes(metric, value, fill = metric)) +
  geom_col(show.legend = FALSE) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "Chloroplast assembly summary",
    x = "",
    y = "Value"
  )

ggsave(file.path(outdir, "chloroplast_assembly_summary.png"), p1, width = 7, height = 5, dpi = 300)

density <- tibble(
  genome = c("Chloroplast", "Chromosome 3"),
  length_mb = c(123717 / 1e6, 16.93),
  annotated_features = c(355, 3447)
) %>%
  mutate(features_per_mb = annotated_features / length_mb)

p2 <- ggplot(density, aes(genome, features_per_mb, fill = genome)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(features_per_mb, 1)), vjust = -0.4, size = 4) +
  theme_bw(base_size = 13) +
  labs(
    title = "Annotated feature density comparison",
    subtitle = "Use as feature density unless these numbers are confirmed as unique genes",
    x = "",
    y = "Annotated features per Mb"
  )

ggsave(file.path(outdir, "gene_density_comparison.png"), p2, width = 7, height = 5, dpi = 300)

copy_if_exists <- function(from, to_dir) {
  if (file.exists(from)) file.copy(from, to_dir, overwrite = TRUE)
}

copy_if_exists("figures/chloroplast/OGDRAW_chloroplast.jpg", outdir)
copy_if_exists("figures/OGDRAW_chloroplast.jpg", outdir)
copy_if_exists("figures/GenAnaCalle_scaffold_1--1321-,1331 ,1337-,1331-,1333 _OGDRAW.jpg", outdir)

message("Chloroplast presentation figures saved to: ", outdir)

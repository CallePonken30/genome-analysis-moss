library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

outdir <- "figures/presentation_results/assembly"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

quast_file <- "analyses/04_evaluation/quast_chr3_polish/transposed_report.tsv"
quast <- read_tsv(quast_file, show_col_types = FALSE)

wanted <- c("# contigs", "# contigs (>= 10000 bp)", "Total length", "N50", "Largest contig", "GC (%)")

quast_long <- quast %>%
  pivot_longer(-Assembly, names_to = "metric", values_to = "value") %>%
  filter(metric %in% wanted) %>%
  mutate(value = parse_number(as.character(value)))

p1 <- ggplot(quast_long, aes(Assembly, value, fill = Assembly)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(
    title = "Chromosome 3 assembly quality before and after polishing",
    x = "",
    y = "Value"
  )

ggsave(file.path(outdir, "quast_polishing_summary.png"), p1, width = 10, height = 6, dpi = 300)

copy_if_exists <- function(from, to_dir) {
  if (file.exists(from)) file.copy(from, to_dir, overwrite = TRUE)
}

copy_if_exists("analyses/04_evaluation/quast_chr3_polish/Nx_plot.pdf", outdir)
copy_if_exists("analyses/04_evaluation/quast_chr3_polish/cumulative_plot.pdf", outdir)
copy_if_exists("analyses/04_evaluation/quast_chr3_polish/GC_content_plot.pdf", outdir)

busco_file <- "analyses/04_evaluation/busco_chr3/busco_chr3/short_summary.specific.embryophyta_odb10.busco_chr3.txt"
busco_lines <- readLines(busco_file)

extract_busco <- function(pattern) {
  as.numeric(str_match(busco_lines[str_detect(busco_lines, pattern)][1], paste0(pattern, "([0-9.]+)%"))[,2])
}

busco_summary <- tibble(
  category = c("Complete", "Fragmented", "Missing"),
  percent = c(extract_busco("C:"), extract_busco("F:"), extract_busco("M:"))
)

p2 <- ggplot(busco_summary, aes(category, percent, fill = category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(percent, "%")), vjust = -0.4, size = 4) +
  ylim(0, 100) +
  theme_bw(base_size = 13) +
  labs(
    title = "BUSCO completeness for chromosome 3 assembly",
    subtitle = "Low completeness is expected because only chromosome 3 was assembled",
    x = "",
    y = "BUSCOs (%)"
  )

ggsave(file.path(outdir, "busco_chr3_summary.png"), p2, width = 7, height = 5, dpi = 300)

message("Assembly presentation figures saved to: ", outdir)

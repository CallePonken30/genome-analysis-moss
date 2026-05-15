library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

outdir <- "figures/presentation_results/fastqc"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

mqc_dir <- "analyses/01_preprocessing/multiqc/multiqc_data"

read_mqc <- function(file) {
  path <- file.path(mqc_dir, file)
  if (!file.exists(path)) stop("Missing file: ", path)
  read_tsv(path, comment = "#", col_types = cols(.default = col_character()))
}

clean_sample <- function(x) {
  x %>%
    str_replace("_fastqc$", "") %>%
    str_replace("\\.fastq\\.gz|\\.fq\\.gz|\\.gz", "")
}

get_read <- function(x) {
  case_when(
    str_detect(x, "R1|f1|_1") ~ "R1",
    str_detect(x, "R2|r2|_2") ~ "R2",
    TRUE ~ "Other"
  )
}

get_type <- function(x) {
  case_when(
    str_detect(x, "trim|paired|clean|fastp") ~ "Trimmed",
    TRUE ~ "Raw"
  )
}

make_long <- function(file, xname, yname) {
  df <- read_mqc(file)

  a <- df %>%
    rename(sample = 1) %>%
    pivot_longer(-sample, names_to = xname, values_to = yname)

  b <- df %>%
    rename(!!xname := 1) %>%
    pivot_longer(-all_of(xname), names_to = "sample", values_to = yname)

  score_a <- sum(!is.na(parse_number(a[[xname]])) & !is.na(parse_number(a[[yname]])))
  score_b <- sum(!is.na(parse_number(b[[xname]])) & !is.na(parse_number(b[[yname]])))

  if (score_b > score_a) b else a
}

quality <- make_long("fastqc_per_base_sequence_quality_plot.txt", "base", "quality") %>%
  mutate(
    sample = clean_sample(sample),
    base = parse_number(base),
    quality = parse_number(quality),
    read = get_read(sample),
    type = get_type(sample)
  ) %>%
  filter(!is.na(base), !is.na(quality), read %in% c("R1", "R2"))

p1 <- ggplot(quality, aes(base, quality, colour = type, group = sample)) +
  geom_hline(yintercept = 30, linetype = "dashed") +
  geom_line(linewidth = 0.9, alpha = 0.8) +
  facet_wrap(~read) +
  theme_bw(base_size = 13) +
  labs(
    title = "FastQC per-base quality before and after trimming",
    subtitle = "Dashed line marks Phred Q30",
    x = "Base position",
    y = "Mean Phred quality",
    colour = ""
  )

ggsave(file.path(outdir, "fastqc_per_base_quality.png"), p1, width = 8, height = 5, dpi = 300)

gc <- make_long("fastqc_per_sequence_gc_content_plot_Counts.txt", "gc", "count") %>%
  mutate(
    sample = clean_sample(sample),
    gc = parse_number(gc),
    count = parse_number(count),
    read = get_read(sample),
    type = get_type(sample)
  ) %>%
  filter(!is.na(gc), !is.na(count), read %in% c("R1", "R2"))

p2 <- ggplot(gc, aes(gc, count, colour = type, group = sample)) +
  geom_line(linewidth = 0.9, alpha = 0.8) +
  facet_wrap(~read, scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(
    title = "FastQC GC content distribution",
    x = "GC content (%)",
    y = "Read count",
    colour = ""
  )

ggsave(file.path(outdir, "fastqc_gc_content.png"), p2, width = 8, height = 5, dpi = 300)

summary_file <- file.path(mqc_dir, "multiqc_fastqc.txt")
if (file.exists(summary_file)) {
  summary <- read_tsv(summary_file, show_col_types = FALSE)
  write_csv(summary, file.path(outdir, "fastqc_summary_table.csv"))
}

message("FastQC presentation figures saved to: ", outdir)

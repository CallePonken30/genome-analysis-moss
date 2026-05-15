library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

outdir <- "figures/presentation_results/fastqc"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

mqc_dir <- "analyses/01_preprocessing/multiqc/multiqc_data"

read_mqc <- function(file) {
  read_tsv(
    file.path(mqc_dir, file),
    comment = "#",
    col_types = cols(.default = col_character())
  )
}

clean_sample <- function(x) {
  x %>%
    str_replace("_fastqc$", "") %>%
    str_replace("\\.fastq\\.gz|\\.fq\\.gz|\\.gz", "")
}

sample_type <- function(x) {
  case_when(
    x %in% c("chr3_illumina_R1", "chr3_illumina_R2") ~ "Raw Illumina",
    x %in% c("R1_paired", "R2_paired") ~ "Trimmed Illumina",
    x == "chr3_clean_nanopore" ~ "Nanopore",
    TRUE ~ "Other"
  )
}

read_pair <- function(x) {
  case_when(
    x %in% c("chr3_illumina_R1", "R1_paired") ~ "R1",
    x %in% c("chr3_illumina_R2", "R2_paired") ~ "R2",
    TRUE ~ "Single read"
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

# -------------------------
# Per-base quality
# -------------------------

quality_all <- make_long("fastqc_per_base_sequence_quality_plot.txt", "base", "quality") %>%
  mutate(
    sample = clean_sample(sample),
    base = parse_number(base),
    quality = parse_number(quality),
    type = sample_type(sample),
    read = read_pair(sample)
  ) %>%
  filter(!is.na(base), !is.na(quality))

quality_illumina <- quality_all %>%
  filter(type %in% c("Raw Illumina", "Trimmed Illumina"))

p1 <- ggplot(quality_illumina, aes(base, quality, colour = type, group = sample)) +
  geom_hline(yintercept = 30, linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~read) +
  theme_bw(base_size = 13) +
  labs(
    title = "Illumina read quality before and after trimming",
    subtitle = "Chromosome 3 Illumina reads. Dashed line marks Phred Q30",
    x = "Base position",
    y = "Mean Phred quality",
    colour = ""
  )

ggsave(file.path(outdir, "illumina_raw_vs_trimmed_quality.png"), p1, width = 8, height = 5, dpi = 300)

quality_nanopore <- quality_all %>%
  filter(type == "Nanopore")

p2 <- ggplot(quality_nanopore, aes(base, quality)) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  theme_bw(base_size = 13) +
  labs(
    title = "Nanopore per-base sequence quality",
    subtitle = "Chromosome 3 Nanopore reads. Dashed line marks Phred Q20",
    x = "Base position",
    y = "Mean Phred quality"
  )

ggsave(file.path(outdir, "nanopore_per_base_quality.png"), p2, width = 8, height = 5, dpi = 300)

# -------------------------
# GC content
# -------------------------

gc_raw <- read_mqc("fastqc_per_sequence_gc_content_plot_Counts.txt")

gc_df <- gc_raw %>%
  pivot_longer(-Sample, names_to = "gc_header", values_to = "pair") %>%
  mutate(
    sample = clean_sample(Sample),
    gc_percent = parse_number(gc_header),
    count = as.numeric(str_match(pair, "\\([0-9.]+, *([0-9.]+)\\)")[,2]),
    type = sample_type(sample),
    read = read_pair(sample)
  ) %>%
  filter(!is.na(gc_percent), !is.na(count))

gc_illumina <- gc_df %>%
  filter(type %in% c("Raw Illumina", "Trimmed Illumina")) %>%
  group_by(sample) %>%
  mutate(relative_abundance = count / sum(count, na.rm = TRUE)) %>%
  ungroup()

p3 <- ggplot(gc_illumina, aes(gc_percent, relative_abundance, colour = type, group = sample)) +
  geom_line(linewidth = 1) +
  facet_wrap(~read) +
  theme_bw(base_size = 13) +
  labs(
    title = "Illumina GC content before and after trimming",
    subtitle = "Chromosome 3 Illumina reads",
    x = "GC content (%)",
    y = "Relative abundance",
    colour = ""
  )

ggsave(file.path(outdir, "illumina_raw_vs_trimmed_gc_content.png"), p3, width = 8, height = 5, dpi = 300)

gc_nanopore <- gc_df %>%
  filter(type == "Nanopore") %>%
  mutate(relative_abundance = count / sum(count, na.rm = TRUE))

p4 <- ggplot(gc_nanopore, aes(gc_percent, relative_abundance)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 13) +
  labs(
    title = "Nanopore GC content distribution",
    subtitle = "Chromosome 3 Nanopore reads",
    x = "GC content (%)",
    y = "Relative abundance"
  )

ggsave(file.path(outdir, "nanopore_gc_content.png"), p4, width = 8, height = 5, dpi = 300)

message("FastQC figures saved to: ", outdir)

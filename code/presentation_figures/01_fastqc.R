library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

outdir <- "figures/presentation_results/fastqc"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

mqc_dir <- "analyses/01_preprocessing/multiqc/multiqc_data"

clean_sample <- function(x) {
  x %>%
    str_replace("_fastqc$", "") %>%
    str_replace("\\.fastq\\.gz|\\.fq\\.gz|\\.gz", "")
}

get_read <- function(x) {
  case_when(
    str_detect(x, "R1|f1|_1") ~ "R1",
    str_detect(x, "R2|r2|_2") ~ "R2",
    TRUE ~ "Single read"
  )
}

read_mqc <- function(file) {
  read_tsv(
    file.path(mqc_dir, file),
    comment = "#",
    col_types = cols(.default = col_character())
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
    read = get_read(sample)
  ) %>%
  filter(!is.na(base), !is.na(quality))

quality_illumina <- quality_all %>%
  filter(sample %in% c("chr3_illumina_R1", "chr3_illumina_R2"))

p1 <- ggplot(quality_illumina, aes(base, quality, colour = sample, group = sample)) +
  geom_hline(yintercept = 30, linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~read) +
  theme_bw(base_size = 13) +
  labs(
    title = "Illumina per-base sequence quality",
    subtitle = "Chromosome 3 Illumina reads. Dashed line marks Phred Q30",
    x = "Base position",
    y = "Mean Phred quality",
    colour = ""
  )

ggsave(file.path(outdir, "illumina_per_base_quality.png"), p1, width = 8, height = 5, dpi = 300)

quality_nanopore <- quality_all %>%
  filter(sample == "chr3_clean_nanopore")

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
    gc_percent = parse_number(gc_header),
    count = as.numeric(str_match(pair, "\\([0-9.]+, *([0-9.]+)\\)")[,2]),
    sample = clean_sample(Sample),
    read = get_read(sample)
  ) %>%
  filter(!is.na(gc_percent), !is.na(count))

gc_illumina <- gc_df %>%
  filter(sample %in% c("chr3_illumina_R1", "chr3_illumina_R2")) %>%
  group_by(sample) %>%
  mutate(relative_abundance = count / sum(count, na.rm = TRUE)) %>%
  ungroup()

p3 <- ggplot(gc_illumina, aes(gc_percent, relative_abundance, colour = sample, group = sample)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 13) +
  labs(
    title = "Illumina GC content distribution",
    subtitle = "Chromosome 3 Illumina reads",
    x = "GC content (%)",
    y = "Relative abundance",
    colour = ""
  )

ggsave(file.path(outdir, "illumina_gc_content.png"), p3, width = 8, height = 5, dpi = 300)

gc_nanopore <- gc_df %>%
  filter(sample == "chr3_clean_nanopore") %>%
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

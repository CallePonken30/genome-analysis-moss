library(ggplot2)
library(dplyr)
library(readr)

outdir <- "figures/presentation_results/synteny"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

paf_file <- "analyses/10_synteny/chr3_vs_cpurpureus.paf"

if (!file.exists(paf_file)) {
  stop("PAF file not found: ", paf_file)
}

paf <- read_tsv(
  paf_file,
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
)

colnames(paf)[1:12] <- c(
  "query", "query_length", "query_start", "query_end",
  "strand", "target", "target_length", "target_start", "target_end",
  "matches", "alignment_length", "mapq"
)

paf <- paf %>%
  mutate(
    query_mid = (query_start + query_end) / 2,
    target_mid = (target_start + target_end) / 2,
    identity = matches / alignment_length
  )

# Keep only more reliable/larger alignments for clearer synteny
paf_filtered <- paf %>%
  filter(
    alignment_length >= 5000,
    mapq >= 20,
    identity >= 0.75
  )

# If filtering is too strict, fall back to top 25% longest alignments
if (nrow(paf_filtered) < 10) {
  paf_filtered <- paf %>%
    filter(alignment_length >= quantile(alignment_length, 0.75, na.rm = TRUE))
}

# -------------------------
# 1. Raw overview dotplot
# -------------------------

p1 <- ggplot(paf, aes(target_mid / 1e6, query_mid / 1e6)) +
  geom_point(alpha = 0.25, size = 0.8) +
  theme_bw(base_size = 13) +
  labs(
    title = "Synteny overview: N. japonicum chr3 vs C. purpureus",
    subtitle = "All minimap2 alignments",
    x = "C. purpureus genome position (Mb)",
    y = "N. japonicum chr3 assembly position (Mb)"
  )

ggsave(file.path(outdir, "synteny_overview_all_alignments.png"),
       p1, width = 8, height = 6, dpi = 300)

# -------------------------
# 2. Filtered dotplot
# -------------------------

p2 <- ggplot(paf_filtered, aes(target_mid / 1e6, query_mid / 1e6, colour = strand)) +
  geom_point(aes(size = alignment_length / 1000, alpha = identity)) +
  theme_bw(base_size = 13) +
  labs(
    title = "Filtered synteny: N. japonicum chr3 vs C. purpureus",
    subtitle = "Longer, higher-confidence alignments only",
    x = "C. purpureus genome position (Mb)",
    y = "N. japonicum chr3 assembly position (Mb)",
    colour = "Orientation",
    size = "Alignment length (kb)",
    alpha = "Identity"
  )

ggsave(file.path(outdir, "synteny_filtered_dotplot.png"),
       p2, width = 8, height = 6, dpi = 300)

# -------------------------
# 3. Segment-style synteny plot
# -------------------------

p3 <- ggplot(paf_filtered) +
  geom_segment(
    aes(
      x = target_start / 1e6,
      xend = target_end / 1e6,
      y = query_start / 1e6,
      yend = query_end / 1e6,
      colour = strand,
      linewidth = alignment_length / 1000
    ),
    alpha = 0.65
  ) +
  scale_linewidth_continuous(range = c(0.2, 1.4)) +
  theme_bw(base_size = 13) +
  labs(
    title = "Syntenic blocks: N. japonicum chr3 vs C. purpureus",
    subtitle = "Segments show aligned genomic blocks after filtering",
    x = "C. purpureus genome position (Mb)",
    y = "N. japonicum chr3 assembly position (Mb)",
    colour = "Orientation",
    linewidth = "Alignment length (kb)"
  )

ggsave(file.path(outdir, "synteny_filtered_segments.png"),
       p3, width = 8, height = 6, dpi = 300)

# -------------------------
# 4. Alignment quality summary
# -------------------------

summary_df <- tibble(
  category = c("All alignments", "Filtered alignments"),
  count = c(nrow(paf), nrow(paf_filtered)),
  total_aligned_mb = c(
    sum(paf$alignment_length, na.rm = TRUE) / 1e6,
    sum(paf_filtered$alignment_length, na.rm = TRUE) / 1e6
  )
)

write_csv(summary_df, file.path(outdir, "synteny_alignment_summary.csv"))

p4 <- ggplot(summary_df, aes(category, total_aligned_mb, fill = category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(total_aligned_mb, 2)), vjust = -0.4, size = 4) +
  theme_bw(base_size = 13) +
  labs(
    title = "Total aligned sequence in synteny comparison",
    x = "",
    y = "Total aligned length (Mb)"
  )

ggsave(file.path(outdir, "synteny_alignment_summary.png"),
       p4, width = 6, height = 5, dpi = 300)

message("Synteny figures saved to: ", outdir)
message("All alignments: ", nrow(paf))
message("Filtered alignments used: ", nrow(paf_filtered))

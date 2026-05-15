library(ggplot2)
library(dplyr)
library(readr)

outdir <- "figures/presentation_results/synteny"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

paf_file <- "analyses/10_synteny/chr3_vs_cpurpureus.paf"

if (file.exists(paf_file)) {
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

  p1 <- ggplot(paf, aes(target_mid / 1e6, query_mid / 1e6, colour = strand)) +
    geom_point(aes(size = alignment_length, alpha = identity)) +
    theme_bw(base_size = 13) +
    labs(
      title = "Synteny dotplot: N. japonicum chr3 vs C. purpureus",
      x = "C. purpureus genome position (Mb)",
      y = "N. japonicum chr3 position (Mb)",
      colour = "Strand",
      size = "Alignment length",
      alpha = "Identity"
    )

  ggsave(file.path(outdir, "chr3_vs_cpurpureus_dotplot.png"), p1, width = 8, height = 6, dpi = 300)

  zoom <- paf %>%
    filter(alignment_length >= quantile(alignment_length, 0.75, na.rm = TRUE))

  p2 <- ggplot(zoom, aes(target_mid / 1e6, query_mid / 1e6, colour = strand)) +
    geom_point(aes(size = alignment_length), alpha = 0.75) +
    theme_bw(base_size = 13) +
    labs(
      title = "Synteny dotplot, longest alignments only",
      x = "C. purpureus genome position (Mb)",
      y = "N. japonicum chr3 position (Mb)",
      colour = "Strand",
      size = "Alignment length"
    )

  ggsave(file.path(outdir, "zoom_synteny_long_alignments.png"), p2, width = 8, height = 6, dpi = 300)
}

copy_if_exists <- function(from, to_dir) {
  if (file.exists(from)) file.copy(from, to_dir, overwrite = TRUE)
}

copy_if_exists("analyses/10_synteny/chr3_vs_cpurpureus.pdf", outdir)
copy_if_exists("analyses/10_synteny/zoom_synteny.pdf", outdir)

message("Synteny presentation figures saved to: ", outdir)

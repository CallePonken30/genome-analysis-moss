library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

outdir <- "figures/presentation_results/annotation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

repeat_tbl <- "analyses/03_annotation/repeatmasker_chr3/polished_chr3.fasta.tbl"
tbl <- readLines(repeat_tbl)

repeat_lines <- tbl[str_detect(tbl, "SINEs|LINEs|LTR elements|DNA elements|Small RNA|Satellites|Simple repeats|Low complexity|Unclassified")]
repeat_df <- tibble(raw = repeat_lines) %>%
  mutate(
    category = str_trim(str_extract(raw, "^[A-Za-z ]+")),
    bp = parse_number(str_extract(raw, "[0-9]+ bp"))
  ) %>%
  filter(!is.na(bp), bp > 0) %>%
  arrange(desc(bp))

p1 <- ggplot(repeat_df, aes(reorder(category, bp), bp / 1e6)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 13) +
  labs(
    title = "RepeatMasker repeat composition on chromosome 3",
    x = "",
    y = "Masked sequence (Mb)"
  )

ggsave(file.path(outdir, "repeatmasker_composition.png"), p1, width = 8, height = 5, dpi = 300)

total_line <- tbl[str_detect(tbl, "total interspersed repeats|Total interspersed repeats|bases masked|total masked")]
writeLines(total_line, file.path(outdir, "repeatmasker_key_lines.txt"))

eggnog_file <- "analyses/03_annotation/eggnog_braker_chr3_masked/chr3_braker_masked_functional.emapper.annotations"
eggnog <- read_tsv(eggnog_file, comment = "##", show_col_types = FALSE)
names(eggnog)[1] <- "query"

cog_col <- names(eggnog)[str_detect(names(eggnog), "COG_category|COG")][1]

if (!is.na(cog_col)) {
  cog_df <- eggnog %>%
    filter(!is.na(.data[[cog_col]])) %>%
    separate_rows(all_of(cog_col), sep = "") %>%
    rename(COG = all_of(cog_col)) %>%
    filter(COG != "-", COG != "") %>%
    count(COG, sort = TRUE)

  p2 <- ggplot(cog_df, aes(reorder(COG, n), n)) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 13) +
    labs(
      title = "eggNOG functional COG categories",
      x = "COG category",
      y = "Number of annotations"
    )

  ggsave(file.path(outdir, "eggnog_cog_categories.png"), p2, width = 8, height = 5, dpi = 300)
}

desc_cols <- names(eggnog)[str_detect(names(eggnog), "Description|Preferred_name|GOs|KEGG|PFAM|eggNOG")]
annotation_col <- desc_cols[1]

if (!is.na(annotation_col)) {
  ann_summary <- eggnog %>%
    mutate(status = ifelse(is.na(.data[[annotation_col]]) | .data[[annotation_col]] == "-" | .data[[annotation_col]] == "", "No clear annotation", "Annotated")) %>%
    count(status)

  p3 <- ggplot(ann_summary, aes(status, n, fill = status)) +
    geom_col(show.legend = FALSE) +
    theme_bw(base_size = 13) +
    labs(
      title = "Annotated vs unclear gene functions",
      x = "",
      y = "Number of predicted proteins"
    )

  ggsave(file.path(outdir, "annotation_status_summary.png"), p3, width = 6, height = 5, dpi = 300)
}

message("Annotation presentation figures saved to: ", outdir)

library(DESeq2)

counts <- read.table(
  "analyses/06_counts/featurecounts_chr3.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
)

sample_info <- read.csv(
  "analyses/07_deseq2/sample_info.csv",
  row.names = 1
)

countdata <- counts[, c(
  "Geneid",
  "analyses/05_expression/mapping/control/Control_1.sorted.bam",
  "analyses/05_expression/mapping/control/Control_2.sorted.bam",
  "analyses/05_expression/mapping/control/Control_3.sorted.bam",
  "analyses/05_expression/mapping/heat_treated/Heat_1.sorted.bam",
  "analyses/05_expression/mapping/heat_treated/Heat_2.sorted.bam",
  "analyses/05_expression/mapping/heat_treated/Heat_3.sorted.bam"
)]

colnames(countdata) <- c(
  "Geneid", "Control_1", "Control_2", "Control_3",
  "Heat_1", "Heat_2", "Heat_3"
)

rownames(countdata) <- countdata$Geneid
countdata <- countdata[, -1]
countdata <- round(as.matrix(countdata))

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = sample_info,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Heat", "Control"))
res <- res[order(res$padj), ]

write.csv(as.data.frame(res), file = "analyses/07_deseq2/deseq2_results_chr3.csv")

sig <- subset(as.data.frame(res), padj < 0.05)
write.csv(sig, file = "analyses/07_deseq2/deseq2_significant_chr3.csv")

sink("analyses/07_deseq2/deseq2_summary.txt")
cat("Total genes after filtering:", nrow(dds), "\n")
cat("Significant genes (padj < 0.05):", nrow(sig), "\n")
sink()

pdf("analyses/07_deseq2/MAplot_chr3.pdf")
plotMA(res, main = "DESeq2 MA plot: Heat vs Control", ylim = c(-5, 5))
dev.off()

vsd <- vst(dds, blind = FALSE)
pdf("analyses/07_deseq2/PCA_chr3.pdf")
plotPCA(vsd, intgroup = "condition")
dev.off()

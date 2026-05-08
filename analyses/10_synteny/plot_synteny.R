paf <- read.delim(
  "chr3_vs_cpurpureus.paf",
  header=FALSE,
  sep="\t",
  fill=TRUE,
  quote=""
)

pdf(
  "chr3_vs_cpurpureus.pdf",
  width=8,
  height=8
)

plot(
  paf$V8,
  paf$V3,
  pch=16,
  cex=0.3,
  xlab="Ceratodon purpureus genome position",
  ylab="N. japonicum chr3 assembly position",
  main="Synteny analysis: chr3 vs C. purpureus"
)

dev.off()

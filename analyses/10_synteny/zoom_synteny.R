paf <- read.delim(
  "chr3_vs_cpurpureus.paf",
  header=FALSE,
  sep="\t",
  fill=TRUE,
  quote=""
)

max_chr3 <- max(paf$V3, na.rm=TRUE)

pdf(
  "zoom_synteny.pdf",
  width=8,
  height=8
)

plot(
  paf$V8,
  paf$V3,
  pch=16,
  cex=0.3,
  xlim=c(0, max_chr3),
  ylim=c(0, max_chr3),
  xlab="Ceratodon purpureus genome position (zoomed)",
  ylab="N. japonicum chr3 assembly position",
  main="Zoomed synteny analysis"
)

dev.off()

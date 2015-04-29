library(RColorBrewer)
library(gplots)

taxlevelMatrix <- function(taxsummary, level, minfraction, frequency = FALSE) {
  df <- subset(taxsummary, taxlevel==level)
  row.names(df) <- df$taxon
  df <- df[6:length(df)]
  df$X <- NULL
  df <- t(df)
  if (frequency) {
    df = prop.table(df, 1)
  }
  # remove the genera with less than minfraction as their maximum relative abundance
  maxcount <- apply(df, 2, max)
  if (frequency) {
    threshold = minfraction 
  }
  else {
    threshold <- minfraction*max(maxcount)
  }
  n1 <- names(which(maxcount >= threshold))
  df <- df[, which(colnames(df) %in% n1)]
  return(as.matrix(df))
}

taxheatmap <- function(taxsummary, level, minfraction, margins, fonts, main, file="") {
  data <- taxlevelMatrix(taxsummary, level, minfraction)
  data.dist <- vegdist(data, method = "bray")
  scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
  row.clus <- hclust(data.dist, "aver")
  data.dist.g <- vegdist(t(data), method = "bray")
  col.clus <- hclust(data.dist.g, "aver")
  print(file)
  pdf(file=file, width=8.5, height=11)
  heatmap.2(data, Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus),
  col = scaleyellowred, trace = "none", density.info = "none", cexCol=fonts[1], cexRow=fonts[2], margins = margins, main=main)
  dev.off()
}

makeTaxHeatMaps <- function(taxsummary) {
  levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (level in c(2,3,4,5,6,7)) {
    filename =  cat(levels[level], ".pdf", sep="")
    taxheatmap(taxsummary, level, 0.05, c(13,4), c(1,0.08), levels[level], filename)
  }
}

findCohousedMetadata <- function(samples, metadata) {
  subset(metadata, Mouse==Cohoused & Date==Date)
}


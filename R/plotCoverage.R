plotCoverage <- function(x, chrom, start=1, end=length(x[[chrom]]), col="blue", xlab="Index", ylab="Coverage", main=chrom) {
  xWindow <- as.vector(window(x[[chrom]], start, end))
  x <- start:end
  xlim <- c(start, end)
  ylim <- c(0, max(xWindow))
  plot(x = start, y = 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, type = "n")
  polygon(c(start, x, end), c(0, xWindow, 0), col = col)
  }
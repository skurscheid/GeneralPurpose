
plotCoverageStrands <- function(pos, neg, chrom, start=1, end=length(pos[[chrom]]), pos.col="blue", neg.col="red", xlab="Index", ylab="Coverage", main=chrom){
  posWindow <- as.vector(window(pos[[chrom]], start, end))
  negWindow <- as.vector(window(neg[[chrom]], start, end))
  x <- start:end
  xlim <- c(start, end)
  ylim <- c(-1, 1) * min(max(posWindow), max(negWindow))
  plot(x = start, y = 0, xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, main = main, type = "n")
  polygon(c(start, x, end), c(0, posWindow, 0), col = pos.col)
  polygon(c(start, x, end), c(0, - negWindow, 0), col = neg.col)
}

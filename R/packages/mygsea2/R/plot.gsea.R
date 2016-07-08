plot.gsea <- function(object, values=NULL,
                      ranklist="Ordered list", ... ) {

  if (! any( class(object)=="gsea"))
    stop("Error: object is not a gsea object.")

  small.list <- object$small.list
  big.list <- object$big.list

  if (is.null(values)) {
    layout( matrix( 1:2, ncol=1), heights=c( 2/3, 1/3))
  } else {
    layout( matrix( 1:3, ncol=1), heights=c( 1/2, 1/6, 2/6))
  }
  #  par(mfrow=c(3,1))

  weightup <- 1/length(small.list)
  weightdo <- 1/(length(big.list)-length(small.list))

  weights <- -rep(weightdo, length(big.list))
  weights[match(small.list, big.list)] <- weightup

  sums <- cumsum(weights)

  mar.orig <- par("mar")

  par(mar = c(0, mar.orig[2:4]))
  
  plot(sums, type="l", col="blue", axes=F, xlab="", ylab="Enrichment score", yaxs="i", ...)
  #  abline(h=0, col="lightgray")
  axis(2)
  box(col="gray")

  par(mar = c(0, mar.orig[2], 0, mar.orig[4]))

  #  plot(sums, type="n", axes=F, ylab="")
  #  xcoord <- match(small.list, big.list)
  #  height <- par("usr")[4]-par("usr")[3]
  #  segments(xcoord, rep(par("usr")[3]+0.1*height, length(xcoord)),
  #           xcoord, rep(par("usr")[4]-0.1*height, length(xcoord)))
  #  abline(v= match(small.list, big.list))
  xcoord <- match(small.list, big.list)
  plot(xcoord, rep(1,length(xcoord)), type="h", axes=F, ylab="",
       xlim=c(1,length(big.list)), ylim=c(0,1))
  box(col="gray")

  if (! is.null(values)) {
    par(mar = c( mar.orig[1:2], 0, mar.orig[4]))
    plot(1:length(values), values, type="h", col="gray", xlab="Rank", ylab=ranklist, axes=F, yaxs="i")

    polygon(c(1, 1:length(values), length(values)),
            c(0, values, 0), border="black", col="gray")
    abline(h=0, col="lightgray")
    axis(1)
    axis(2)
    box(col="gray")

    crosszero <- which(diff(sign(values)) != 0)
    if (length(crosszero)==1)
      abline(v=crosszero, lty=3)
  }

  par("mar" = mar.orig)
}

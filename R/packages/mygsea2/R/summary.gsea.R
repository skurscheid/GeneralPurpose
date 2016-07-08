summary.gsea <- function(object) {

  if (! any( class(object)=="gsea"))
    stop("Error: object is not a gsea object.")

}	

print.gsea <- function(object) {

  if (! any( class(object)=="gsea"))
    stop("Error: object is not a gsea object.")

  if (! is.null(object$weights)) {
    cat("Weighted ")
  }

  cat("GSEA analysis (", object$nperms," permutations)\n\n", sep="")
  cat("Small list: ", length(object$small.list),"\n",
      "  Big list: ", length(object$big.list),"\n\n", sep="")

  coefs <- cbind( c(object$ks.pos, object$ks.neg),
		  c(object$p.pos, object$p.neg) )

  colnames(coefs) <- c("Ks stat","P-value")
  rownames(coefs) <- c("+", "-")

  printCoefmat(coefs, P.values=TRUE, has.Pvalue=TRUE)

#  cat("Ks pos: ", object$ks.pos, object$p.pos, "\n")
#  cat("Ks neg: ", object$ks.neg, object$p.neg, "\n")
  
}


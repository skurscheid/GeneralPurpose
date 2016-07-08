print.genesets <- function(genesets) {

  if (class(genesets)!="genesets")
    stop("Error: object should be of class genesets.")

  cat("Name:", genesets$name)

  if (!is.null(genesets$version))
    cat(" (version: ", genesets$version,")",sep="")

  cat("\nLength:", length(genesets$genesets))

  cat("\n\nSource: ", genesets$source, "\n", sep="")

  if (!is.null(genesets$collections)) {
    cat("\nCollections:\n")
    print(summary( genesets$collections))
  }
}

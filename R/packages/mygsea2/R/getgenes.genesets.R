getgenes.genesets <- function(object, geneset) {

  index <- which(object$names == geneset)

  if (length(index)==0)
    stop("Unknown gene set.")

  # We assume there is only one geneset for a given name; is should
  # be constrained somewhere else.

  return( object$genesets[[index]])
}

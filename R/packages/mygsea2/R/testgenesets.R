testgenesets <- function(lists, big.list, weights=NULL, n.perm=10000,
                         min.size=15, max.size=500,
                         verbose=FALSE) {

  if (class(lists) != "character" && class(lists) != "genesets")
    stop("Error: lists should a list of gene symbols or a geneset object.")

  if (class(lists) == "genesets") {
    names <- lists$names
    lists <- lists$genesets
  }

  # Remove duplicates from the big.list (code taken from mygsea2,
  # so that it does not have to be executed at every iteration

  dups=duplicated(big.list) | is.na(big.list)
  if (length(big.list[dups])>0) {
    if (verbose)
      cat("WARNING: Duplicated symbols found in big list, occurrences after the first removed.\n",paste(unique(big.list[dups]),coll=" "),"\n",sep="")
    big.list=big.list[!dups]
    if (!is.null(weights)) {
      cat("WARNING: Weights vector adjusted to big list length.\n")
      weights=weights[!dups]
    }
  }

  geneset <- NULL
  ks.pos <- NULL
  p.pos <- NULL
  ks.neg <- NULL
  p.neg <- NULL
  
# Use lapply ?
  for (i in 1:length(lists)) {    
    small.list <- lists[[i]]
    small.list <- small.list[ small.list != "" ]

    # @@TODO: Test if small.list \in big.list

    if (verbose) {
      cat("Geneset #", i,": ", names[i],"/", length(lists),
          " (length ",length(small.list),")\n", sep="")
    }

    if (length(small.list)<min.size || length(small.list)>max.size) {
      if (verbose)
        cat(" length of list outside of ", min.size, "-", max.size,".\n")
      
      next
    }
    
    result <- mygsea2(small.list, big.list, weights=weights, n.perm=n.perm,
                      verbose=verbose)

    geneset <- c(geneset, names[i])
    ks.pos <- c(ks.pos,result$ks.pos)
    p.pos <- c(p.pos, result$p.pos)
    ks.neg <- c(ks.neg, result$ks.neg)
    p.neg <- c(p.neg, result$p.neg) 
    
  }

  data.frame(geneset, ks.pos, p.pos, ks.neg, p.neg)
}

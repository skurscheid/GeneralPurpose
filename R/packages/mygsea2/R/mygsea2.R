# @@TODO: use the R warning() system instead of cat.

mygsea2 = function(small.list, big.list, weights=NULL, n.perm=10000, seed=NULL, verbose=FALSE) {

  if (!is.null(weights) & length(big.list)!=length(weights)) {
    cat("ERROR: Number of weights not equal to number of gene in big list.\n")
    return(as.list(NA,NA,NA,NA))
  }
  
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

  dups=duplicated(small.list)
  if (length(small.list[dups])>0) {
    if (verbose)
      cat("WARNING: Duplicated symbols found in small list, occurrences after the first removed.\n",paste(unique(big.list[dups]),coll=" "),"\n",sep="")
    small.list=small.list[!dups]
  }

  g    =match(small.list,big.list)
  nas  =is.na(g)
  if (length(g[nas])!=0) {
    if (verbose)
      cat("WARNING: Small list symbols not found in big list, symbols removed:\n",paste(small.list[nas],coll=" "),"\n",sep="")
  }
  small.list=small.list[!nas]
  if (length(small.list)==0) {
    cat("WARNING: Small list empty, no calculation performed.\n",paste(g[nas]),"\n")
    return(list(ks.pos=NA, ks.neg=NA, p.pos=NA, p.neg=NA))
  }
  g    =sort(match(small.list,big.list))-1  # -1 because g array is 0-based in C code
  m    =as.integer(length(small.list))
  ngene=as.integer(length(big.list))
  r    =as.double(1:ngene)

  if (verbose) {
#    cat("Some diagnostic here...\n");
  }

  if (is.null(seed)) {
    seed=runif(1)*(2^31)
  }

  resks=as.double(c(0.0,0.0))
  resperm=as.double(c(0.0,0.0))

  if (is.null(weights)) {
    w=rep(1,length(r))
  } else {
    w=weights
  }

  res=
    .C("Rpermute_wKS",
       m      =as.integer(m),
       g      =as.integer(g),
       r      =as.double(r),
       ngene  =as.integer(ngene),
       weights=w,
       nperm  =as.integer(n.perm),
       seed   =as.integer(seed),
       resks  =as.double(resks),
       resperm=as.double(resperm),
       iwrk   =as.integer(rep(0,ngene)));
    
#  ks.pos=res$resks[1]
#  ks.neg=res$res[2]
#  p.pos=res[3]
#  p.neg=res[4]

#  res

  z <- list(ks.pos=res$resks[1], ks.neg=res$resks[2], p.pos=res$resperm[1], p.neg=res$resperm[2])

  z$nperms <- n.perm
  z$weights <- weights
  z$small.list <- small.list
  z$big.list <- big.list

  class(z) <- "gsea"

  return(z)
}

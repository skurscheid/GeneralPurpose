\name{mygsea2}
\alias{mygsea2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gene Set Enrichment Analysis }
\description{
  Type ?gs.mygsea2 to get help for mygsea2 gene set manipulation functions.
  
  mygsea2 computes p-values associated with the position of
  small-list list of genes in a ranked (and possibly weighted)
  list of genes in a big list.
  An empirical KS test is performed by shuffling the position
  of the small-list genes in the big list.

  For the description of gene-set related helper functions:
  ?gs.mygsea2
}
\usage{
mygsea2(small.list, big.list, weights=NULL, n.perm=10000,
seed = NULL, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{small.list}{\code{small.list} must be a non-empty vector of
    unique strings which exist in big.list}
  \item{big.list}{\code{big.list} must be a non-empty vector of unique
    strings representing ranked genes. Ranking is typically from large
  positive to near-zero to large-negative t-values.}
  \item{weights}{\code{weights} is the weighting associated to the genes
  it is typically the absolute value of t in a t-test. If weights=NULL,
  all genes are given a weight of 1.0.}
  \item{n.perm}{\code{n.perm} is the number of permutations to perform
    to compute the empirical p-value }
  \item{seed}{\code{seed} is the value to initialize the RNG. If
    seed=NULL are a random seed is picked. }
  \item{verbose}{\code{verbose} is used for debugging purpose. The
    output if verbose=TRUE is subject to change.}
}
\details{
  mygsea2 performs a permutation on the genes, thus not taking into
  account the structure that can exist in the expression profile.

  WARNING: In its current implementation, the computed p-value may be
  subject to numerical instabilities if 1. The big list contains a small
  number of genes (say less than about 100) AND 2. all genes have the
  same weight. 
}
\value{
  mygsea2 returns a list which components are:
  \item{ks.pos }{Maximum of the running sum for the actual realization
    of the score, with genes in the small list (hits) given the weight
    picked from the "weight" vector, and genes not in the small list
    given a weight of -1/(number of genes in big list not in small
    list). If t-basef ranking, ks.pos picks genes coregulated with the
    factor used to make the ranking.}
  \item{ks.neg }{Minimum of the running sum. Under t-test-based ranking,
  ks.neg picks genes anticorrelated with the ranking factor.}
  \item{p.pos }{p-value associated to ks.pos}
  \item{p.neg }{p-value associated to ks.neg}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ T. Sengstag with contributions of Asa Wirapati and Vlad Popovici}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ gs.mygsea2 }
\examples{
   set.seed(0)

   #
   #  MSigDB gene lists
   #  data(msigdb)
   #  names(msigdb)
   #

   big.list=paste("G",1:20000,sep="")
   #
   #  Genes distributed randomly 
   #
   small.list=paste("G",sample(1:20000,size=20),sep="")
   mygsea2(small.list,big.list,n.perm=10000)

   #
   #  Genes with strong enrichment in positive t-values
   #
   small.list=paste("G",sample(1:10000,size=20),sep="")
   mygsea2(small.list,big.list,n.perm=10000)

}
\keyword{ misc }% at least one, from doc/KEYWORDS

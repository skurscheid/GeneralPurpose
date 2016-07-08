\name{testgenesets}
\alias{testgenesets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gene Set Enrichment Analysis of a set of genesets }
\description{
  testgenesets is a wrapper around the mygsea2 function for testing
  several genesets one after the other.
}
\usage{
testgenesets(lists, big.list, weights=NULL, n.perm=10000,
seed = NULL, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lists}{\code{small.list}: a list of a genesets object containing
    genesets}
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

}
\value{
  testgenesets returns a matrix, with each row containing results for a
  given geneset, in the following form:
  \item{geneset }{The geneset name}
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
\author{Frederic Schutz}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
   # Assume big.list contains a ranked list of genes

   data(msigdb)

   results <- testgenesets( msigdb, big.list )

}
\keyword{ misc }% at least one, from doc/KEYWORDS

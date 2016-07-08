\name{gs.mygsea2}
\alias{gs.mygsea2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gene Set Enrichment Analysis, Gene Set Helper Functions }
\description{
  gs.mygsea2 provides functions to retrieve gene sets from a database.
  (Currently, development MySQL database on sib-pc17.)
}
\usage{

Gene Set functions:
-------------------

  gs.getGeneSetList(class=NULL, taxid=9606)

  gs.getGeneSet(geneset=NULL, taxid=9606)

  gs.findSymbolInGeneSets(symbol=NULL, class=NULL, taxid=9606)

  gs.findGeneIDInGeneSets(geneid=NULL, taxid=9606)

  gs.annotateGeneSets(geneset=NULL)

Helper functions:
-----------------

  gs.getSymbolFromGeneID(geneid=NULL)

  gs.getGeneIDFromSymbol(symbol=NULL, taxid=9606)

  gs.getChrLocFromGeneID(geneid=NULL)

  gs.getChrLocFromSymbol(symbol=NULL, taxid=9606)

  gs.getTaxIDFromGeneID(geneid=NULL)

  gs.findAlternativeSymbols(symbol = NULL, taxid = 9606, strict=FALSE,
                            fullname=FALSE, verbose=FALSE)

# DO NOT USE THE FOLLOWING FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE DOING!!!

  gs.SQL(query=NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  Arguments should be clear from the definition of the functions and of
  the examples below. Less intuitive arguments are:

  \item{taxid}{
  "taxid" is NCBI's taxonomy ID for the organism to search symbols for
  (e.g. human=9606, mouse=10090, rat=10116, which are the currently
  supported organisms ask the package maintainers to extend the list
  if required).}
  
  \item{class}{The general category of a gene set, e.g. "GO" or
    "KEGG". This parameter can be used in some functions to restrict
    the scope of search.}
  
  Currently only atomic parameters can be passed to the functions,
  i.e. no lists are accepted (exception: gs.annotateGeneSets
  can accept a vector of gene set names.)
  
%  \item{small.list}{\code{small.list} must be a non-empty vector of
%    unique strings which exist in big.list}
}
\details{
  Primary information are retrieved from NCBI's and KEGG's FTP site.
  Gene symbols are as provided by NCBI. In most cases these are HUGO
  approved symbols, yet minor differences may occur (e.g. if a
  gene has not yet received an official symbol, a temporary 
  one can be provided; in some releases of NCBI's gene_info file
  HUGO-approved mitochondrial gene symbols were discarded to use
  the common symbols).

  gs.findAlternativeSymbols(symbol, taxid, strict, fullname, verbose):

  The string provided as "symbol" is case insensitive, both upper and
  lower case will match.

  If "strict" is TRUE then only symbols matching exactly the query are
  returned (though case remains irrelevant).
  
  If "fullname", searches and print hits in the full name of the genes
  as well.
  
  If "verbose", prints the MySQL query.
}

\value{
  gs.mygsea2 return values should be obvious.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ T. Sengstag and F. Schutz }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
#
#  Retrieving gene sets from the database
#  --------------------------------------

  list = gs.getGeneSetList("KEGG")
  list[1:10,]                               # you get the number of
                                            # genes for free...
  gs.annotateGeneSets(list[1:10,1])
  
  gs.getGeneSet("KEGG:05212")               # default organism Human...
  gs.getGeneSet("KEGG:05212", taxid=10090)  # ... the mouse version

  gs.findSymbolInGeneSets(c("A1BG","NAT2"))

  gs.findGeneIDInGeneSets(c(1,10))

  gs = gs.findSymbolInGeneSets("BRAF", class="KEGG")
  gs.annotateGeneSets(gs[,"geneset"])

#
#  Helper functions for mapping GeneIDs with symbols
#  -------------------------------------------------

  gs.getSymbolFromGeneID(2099)

  my.genes <-  c("esr1","A1BG","NAT2")
  gs.getGeneIDFromSymbol(my.genes)               # ... search is case insensitive...
  gs.getGeneIDFromSymbol(my.genes,taxid=10090)   # the mouse homologs... 

# ... but here is a big warning:
#
# WARNING: HOMOLOGY INFERRED FROM GENE SYMBOLS CAN BE WRONG !!!
# Consider using HomoloGene to find actual gene homologs.

#
#  Finding alternative symbols in annotation
#  -----------------------------------------
# ... in human
  gs.findAlternativeSymbols("erbb")
  gs.findAlternativeSymbols("erbb2")
  gs.findAlternativeSymbols("tcf4",strict=TRUE)

# ... in mouse
  gs.findAlternativeSymbols("tcf4", strict=TRUE, taxid=10090)

# ... searching the full names as well
  gs.findAlternativeSymbols("epidermal", fullname=TRUE)


#
#  Getting chromosomal locations
#  -----------------------------

  gs.getChrLocFromSymbol(my.genes)                # The location of our genes...
  gs.getChrLocFromSymbol(my.genes, taxid=10090)   # ... in mouse...
  gs.getChrLocFromSymbol(my.genes, taxid=10116)   # ... and in rat.

  my.geneids <- c(2099, 13982, 24890)  #  Mixing organisms is even possible...
  gs.getTaxIDFromGeneID(my.geneids)    #  ... but here you are on your own...
  gs.getChrLocFromGeneID(my.geneids)

#
#  Advanced usage using the gs.SQL functions
#  -----------------------------------------
#  WARNING: Using this function there is a potential to destroy
#  the database. Do not use the gs.SQL function unless you are
#  truly sure you will not harm the resource.
#
gs.SQL("show tables;")
gs.SQL("describe gene_info;")

#
#  Find all genes on Human Chromosome 9 and sort them accroding to position
#
  chrom9 = gs.SQL("
select distinct
  gene2accession.GeneID,Symbol,
  start_position_on_the_genomic_accession,
  end_position_on_the_genomic_accession
from gene2accession,gene_info
where gene2accession.GeneID = gene_info.GeneID
  and gene_info.tax_id = 9606
  and chromosome='9'
  and genomic_nucleotide_accession like 'NC_\%'
  order by start_position_on_the_genomic_accession;")
  chrom9[1:10,]

#
# Statistics on the number of genes in each chromosome
#
gs.SQL("
select chromosome,count(*)
from gene_info where tax_id=9606
group by chromosome order by count(*);")

#
#  Statistics on annotated type of genes in mouse
#
gs.SQL("
select type_of_gene,count(*)
from gene_info where tax_id=10090
group by type_of_gene order by count(*);");

}
\keyword{ misc }% at least one, from doc/KEYWORDS

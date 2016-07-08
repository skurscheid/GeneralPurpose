.RGS.host     = "sib-pc17.unil.ch"
.RGS.user     = "bcf"
.RGS.password = NULL
.RGS.dbname    = "genesets"

gs.getGeneSet = function(geneset=NULL, taxid=9606) {
  if (is.null(geneset)) {
    warning("Must provide a GeneSet ID.\n");
    return(NULL)
  }
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  query = paste('select geneset2geneid.geneid,symbol,weight from geneset2geneid,gene_info where ',
          sprintf('geneset="%s"',geneset),
          ' and tax_id = ',sprintf("%d",taxid),
          ' and geneset2geneid.geneid=gene_info.GeneID order by symbol;',
          sep='')
  rs  = dbSendQuery(con,query)
  genes = fetch(rs, n=-1)
  dbDisconnect(con)
  if (nrow(genes)!=0) {
    return(list(geneids=genes[,1],symbols=genes[,2],weights=genes[,3]))
  } else {
    return(list(geneids=NULL,symbols=NULL,weights=NULL))
  }
}

#-------------------------------------------------------------------------------

gs.getGeneSetList = function(class=NULL, taxid=9606) {
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  query = paste('select geneset,tax_id,count(*) ',
    'from geneset2geneid,gene_info where ',
    'geneset2geneid.geneid = gene_info.GeneID ',
    ifelse(!is.null(class),sprintf('and geneset like "%s:%%" ',class),''),
    'and tax_id = ',sprintf("%d ",taxid),
    'group by geneset ',
    'order by geneset;',sep='')
#  cat(query,"\n")
  rs   = dbSendQuery(con,query)
  list = fetch(rs, n=-1)
  if (nrow(list)>0) {
    list = list[,c(1,3)]
    colnames(list) = c("geneset","number_of_genes")
  }
  dbDisconnect(con)
  return(list)
}

#-------------------------------------------------------------------------------

gs.findSymbolInGeneSets = function(symbol=NULL, class=NULL, taxid=9606) {
  if (is.null(symbol)) {
    warning("Must provide a gene symbol.\n");
    return(NULL)
  }
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  symbol = unlist(symbol)
  symbol = paste(symbol,sep="",collapse="\",\"")
  query = paste("
select distinct gene_info.GeneID,Symbol,geneset
from geneset2geneid,gene_info
where geneset2geneid.geneid=gene_info.GeneID
      and Symbol in  (\"",symbol,"\")",
      ifelse(!is.null(class),sprintf(' and geneset like "%s:%%" ',class),''),
      " and tax_id=",sprintf('"%d"',taxid),";",sep='')
  rs       = dbSendQuery(con,query)
  genesets = fetch(rs, n=-1)
  dbDisconnect(con)
  return(genesets)
}

#-------------------------------------------------------------------------------

gs.getTaxIDFromGeneID = function(geneid=NULL) {
  if (is.null(geneid)) {
    warning("Must provide a gene id.\n");
    return(NULL)
  }
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  geneid = unlist(geneid)
  geneid = paste(geneid,sep="",collapse="\",\"")
  query = paste("
select distinct GeneID,tax_id
from gene_info 
where GeneID in (\"",geneid,"\");",sep='')
  rs       = dbSendQuery(con,query)
  genesets = fetch(rs, n=-1)
  dbDisconnect(con)
  return(genesets)
}

#-------------------------------------------------------------------------------

gs.findGeneIDInGeneSets = function(geneid=NULL) {
  if (is.null(geneid)) {
    warning("Must provide a Gene ID.\n");
    return(NULL)
  }
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  geneid=unlist(geneid)
  geneid=paste(geneid,sep="",collapse="\",\"")
  query = paste("select distinct tax_id 
from gene_info where GeneID in (\"",geneid,"\");
",sep="");
  rs       = dbSendQuery(con,query)
  res      = fetch(rs, n=-1)
  if (length(unlist(res[,1]))>1) {
    warning("More than one organism implicitly specified by Gene IDs; tax_id:",
            paste(unlist(res),collapse=" "),collapse=" ")
  }
  query = paste('
select distinct geneset2geneid.geneid,Symbol,geneset
from geneset2geneid,gene_info
where geneset2geneid.geneid = gene_info.GeneID and geneset2geneid.geneid in ("',geneid,'");',sep='')
  rs       = dbSendQuery(con,query)
  genesets = fetch(rs, n=-1)
  dbDisconnect(con)
  return(genesets)
}

#-------------------------------------------------------------------------------

gs.getSymbolFromGeneID = function(geneid=NULL) {
  if (is.null(geneid)) {
    warning("Must provide a Gene ID.\n");
    return(NULL)
  }
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  geneid=unlist(geneid)
  geneid=paste(geneid,sep="",collapse="\",\"")
  query = paste("select distinct tax_id 
from gene_info where GeneID in (\"",geneid,"\");
",sep="");
  rs       = dbSendQuery(con,query)
  res      = fetch(rs, n=-1)
  if (length(unlist(res))>1) {
    warning("More than one organism implicitly specified by Gene IDs; tax_id:",
            paste(unlist(res),collapse=" "),collapse=" ")
  }
  query = paste('select distinct GeneID,Symbol from gene_info',
    ' where GeneID in ("',geneid,'");',sep='')
  rs       = dbSendQuery(con,query)
  symbols = fetch(rs, n=-1)
  dbDisconnect(con)
  return(symbols)
}

#-------------------------------------------------------------------------------

gs.getChrLocFromGeneID = function(geneid=NULL) {
  if (is.null(geneid)) {
    warning("Must provide a (list of) Gene ID(s).\n");
    return(NULL)
  }
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  geneid=unlist(geneid)
  geneid=paste(geneid,sep="",collapse="\",\"")
  query = paste("select distinct tax_id 
from gene_info where GeneID in (\"",geneid,"\");
",sep="");
  rs       = dbSendQuery(con,query)
  res      = fetch(rs, n=-1)
  if (length(unlist(res))>1) {
    warning("More than one organism implicitly specified by Gene IDs; tax_id:",
            paste(unlist(res),collapse=" "),collapse=" ")
  }
  query = paste("
select distinct
  gene2accession.GeneID,
  Symbol,
  chromosome,
  genomic_nucleotide_accession,
  start_position_on_the_genomic_accession,
  end_position_on_the_genomic_accession
from gene2accession,gene_info
where gene2accession.GeneID = gene_info.GeneID
  and gene2accession.GeneID in (\"",geneid,"\")
  and genomic_nucleotide_accession like 'NC_%';",sep="")
#  order by GeneID;",sep="")
  rs       = dbSendQuery(con,query)
  res      = fetch(rs, n=-1)
  dbDisconnect(con)
  if (nrow(res)>1) {
    if (length(which(duplicated(res[,1])))>0) {
      warning(paste("Duplicate entries found for Gene IDs: ",
                    paste(unique(res[which(duplicated(res[,1])),1]),sep=" ",collapse=""),
                    sep="",collapse=""))
    }
  }
  return(res)
}
# gids=c(2,1,3,4,5,6,7,8,9,10)
# gs.getChrLocFromGeneID(gids)

#-------------------------------------------------------------------------------

gs.getChrLocFromSymbol = function(symbol=NULL, taxid=9606) {
  if (is.null(symbol)) {
    warning("Must provide a (list of) Gene ID(s).\n");
    return(NULL)
  }
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  symbol=unlist(symbol)
  symbol=paste(symbol,sep="",collapse="\",\"")
  query = paste("
select distinct
  gene2accession.GeneID,
  Symbol,
  chromosome,
  genomic_nucleotide_accession,
  start_position_on_the_genomic_accession,
  end_position_on_the_genomic_accession
from gene2accession,gene_info
where gene2accession.GeneID = gene_info.GeneID
  and gene_info.tax_id = ",taxid,"
  and Symbol in (\"",symbol,"\")
  and genomic_nucleotide_accession like 'NC_%';",sep="")
#  order by GeneID;",sep="")
  rs       = dbSendQuery(con,query)
  res      = fetch(rs, n=-1)
  dbDisconnect(con)
  if (nrow(res)>1) {
    if (length(which(duplicated(res[,1])))>0) {
      warning(paste("Duplicate entries found for Gene IDs: ",
                    paste(unique(res[which(duplicated(res[,1])),1]),sep=" ",collapse=""),
                    sep="",collapse=""))
    }
  }
  return(res)
}
#syms=c("ERBB2","NAT2")
#gs.getChrLocFromSymbol(syms)

#-------------------------------------------------------------------------------

gs.getGeneIDFromSymbol = function(symbol=NULL, taxid=9606) {
  if (is.null(symbol)) {
    warning("Must provide a Gene ID.\n");
    return(NULL)
  }
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  symbol = unlist(symbol)
  symbol = paste(symbol,sep="",collapse="\",\"")
  query = paste("
select distinct GeneID,Symbol
from gene_info
where Symbol in (\"",symbol,"\")
and tax_id=",sprintf("%d",taxid),
    ";",sep='')
  rs       = dbSendQuery(con,query)
  geneids = fetch(rs, n=-1)
  dbDisconnect(con)
  return(geneids)
}

#-------------------------------------------------------------------------------

gs.annotateGeneSets = function(geneset=NULL) {
  annot=NULL
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  for (i in 1:length(geneset)) {
    query = sprintf('select geneset,description,comment from geneset_annot where geneset="%s";',geneset[i]);
    rs    = dbSendQuery(con,query)
    res   = fetch(rs, n=-1)
    if (nrow(res)==0) {
      warning(sprintf("No geneset %s\n",geneset[i]))
      annot[i]=""
    } else {
      if(nrow(res)>1) {
        warning(sprintf("Multiple hits for %s onyly first one returned\n",geneset[i]))
      }
      annot[i]=paste(res[1,2],ifelse( !is.na(res[1,3]), sprintf(' ( %s )',res[1,3]), '' ) ,sep="")
    }
  }
  dbDisconnect(con)
  return(as.data.frame(cbind(geneset,annotation=annot)))
}

#-------------------------------------------------------------------------------

gs.findAlternativeSymbols = function(symbol = NULL, taxid = 9606, strict=FALSE, fullname=FALSE, verbose=FALSE) {
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  if (is.null(symbol)) {
    warning("Must provide a symbol.\nCommon taxid for organisms are:\n 9606 human\n10090 mouse\n10116 rat\n");
    return(NULL)
  }

  query = paste("select Symbol,Synonyms,Symbol_from_nomenclature_authority,GeneID",
    ifelse(fullname,",Full_name_from_nomenclature_authority",""),
    " from gene_info ",
    sep="")

  if (strict) {

    query = paste(query,
      sprintf("where tax_id=%d and (Symbol like \"%s\" or Synonyms like \"%s\" or Synonyms like \"%s|%%\" or Synonyms like \"%%|%s\" or Synonyms like \"%%|%s|%%\" or Symbol_from_nomenclature_authority like \"%s\"",taxid,symbol,symbol,symbol,symbol,symbol,symbol),
      ifelse(fullname,sprintf(" or Full_name_from_nomenclature_authority like \"%s\"",symbol),""),
      sep="");

  } else {
    
    query = paste(query,
      sprintf("where tax_id=%d and (Symbol like \"%%%s%%\" or Synonyms like \"%%%s%%\" or Symbol_from_nomenclature_authority like '%s'",
              taxid,symbol,symbol,symbol,symbol),
      ifelse(fullname,sprintf(" or Full_name_from_nomenclature_authority like \"%%%s%%\"",symbol),""),
      sep="")
  }
  
  query = paste (query, ");",sep="")

  if (verbose) {
    cat("MySQL query:\n",query,"\n")
  }
  
  rs    = dbSendQuery(con,query)
  res   = fetch(rs, n=-1)
  
  dbDisconnect(con)
  if (nrow(res)==0) {return(NULL)}
  return(res)
}

#-------------------------------------------------------------------------------

gs.SQL = function(sql=NULL) {
  if (is.null(sql)) {
    return()
  }
  library(RMySQL)
  con   = dbConnect(MySQL(),
    user      =.RGS.user,
    password  =.RGS.password,
    host      =.RGS.host,
    dbname    =.RGS.dbname)
  rs  = dbSendQuery(con,sql)
  res = fetch(rs, n=-1)
  dbDisconnect(con)
  return(res)
}

#-------------------------------------------------------------------------------

if(FALSE) {
# Examples:
  list = gs.getGeneSetList("KEGG")
  list[1:10,]
  
  gs.annotateGeneSets(list[1:10,])
  
  gs.getGeneSet("05212")
  
  gs=gs.findSymbolInGeneSets("KRAS")
  gs=gs[which(gs[,1]=="KEGG"),]
  gs.annotateGeneSets(gs)
  
  gs=gs.findSymbolInGeneSets("BRAF")
  gs=gs[which(gs[,1]=="KEGG"),]
  gs.annotateGeneSets(gs)

  gids=c(2,1,3,4,5,6,7,8,9,10)
  gs.getChrLocFromGeneID(gids)
  
  syms=c("ERBB2","NAT2")
  gs.getChrLocFromSymbol(syms)

  chrom9 = gs.SQL("select distinct gene2accession.GeneID,Symbol,start_position_on_the_genomic_accession,end_position_on_the_genomic_accession from gene2accession,gene_info where gene2accession.GeneID = gene_info.GeneID and gene_info.tax_id = 9606 and chromosome='9' and genomic_nucleotide_accession like 'NC_%' order by start_position_on_the_genomic_accession;")

  
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

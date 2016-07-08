amILocal <- function(machinename = NULL){
  if(is.null(machinename)) stop("Machinename is missing")
  m <- Sys.info()["nodename"]
  mn <- unlist(lapply(strsplit(m, "\\."), function(x) x[1]))
  if (mn == machinename) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


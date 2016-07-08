amILocal <- function(){
  m <- Sys.info()["nodename"]
  mn <- unlist(lapply(strsplit(m, "\\."), function(x) x[1]))
  if (mn == "JCSMR027564ML") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


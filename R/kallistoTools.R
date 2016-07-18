# creates matrix of kallisto output files
makeCountMatrix <- function(files = NULL, column = NULL){
  if (is.null(files)) {stop("Filenames missing")}
  if (is.null(column)) {stop("Column name missing")}
  n <- length(files)
  if (!is.null(column)){
    header = T
  } else {
    header = F
  }
  for (i in 1:n){
    s <- unlist(lapply(strsplit(names(files)[i], "\\."), function(x) x[1]))
    if (i == 1) {
      mat0 <- read.table(files[i], header = header, as.is = T, sep = ("\t"))
      rn <- mat0[,i]
      mat0 <- as.matrix(mat0[, column])
      colnames(mat0) <- s
      rownames(mat0) <- rn
    } else {
      mat1 <- read.table(files[i], header = header, as.is = T, sep = ("\t"))
      mat1 <- as.matrix(mat1[, column])
      colnames(mat1) <- s
      mat0 <- cbind(mat0, mat1)
    }
  }
  return(mat0)
}

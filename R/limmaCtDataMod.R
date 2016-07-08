limmaCtDataMod <- function (q, design = NULL, contrasts, sort = TRUE, stringent = TRUE, 
          ndups = 1, spacing = NULL, dupcor, ret.fit = FALSE, topTableOut = FALSE, ...) 
{
  data <- exprs(q)
  featPos <- featurePos(q)
  if (missing(dupcor)) {
    if (ndups > 1) {
      dup.cor <- duplicateCorrelation(data, ndups = ndups, 
                                      spacing = spacing, design = design)
      temp <- unwrapdups(featPos, ndups = ndups, spacing = spacing)
      featPos <- apply(temp, 1, paste, collapse = ";")
    }
    else {
      dup.cor <- NULL
    }
  }
  fit <- lmFit(data, design = design, ndups = ndups, spacing = spacing, 
               correlation = dup.cor$consensus, ...)
  if (!missing(contrasts)) 
    fit <- contrasts.fit(fit, contrasts = contrasts)
  fit2 <- eBayes(fit)
  out <- list()
  if (!missing(contrasts)) {
    coefs <- colnames(contrasts)
    cont <- design %*% contrasts
  }
  else {
    coefs <- colnames(design)
    cont <- design
  }
  for (coef in coefs) {
    res <- topTable(fit2, coef = coef, number = nrow(fit2), 
                    sort.by = "none", confint = TRUE, ...)
    if (topTableOut == TRUE){
      tt1 <- res
    }
    both.means <- both.cats <- array(0, c(nrow(res), 2), 
                                     list(rownames(res), c("Test", "Reference")))
    for (i in c(-1, 1)) {
      sample <- ifelse(i == 1, "Test", "Reference")
      index <- cont[, coef] == i
      mean <- rowMeans(unwrapdups(data[, index], ndups = ndups, 
                                  spacing = spacing))
      new.cat <- rep("OK", length(mean))
      old.cat <- unwrapdups(featureCategory(q)[, index], 
                            ndups = ndups, spacing = spacing)
      count.cat <- apply(old.cat, 1, function(x) sum(x %in% 
                                                       c("Undetermined", "Unreliable")))
      cutoff <- ifelse(stringent, 1, ceiling(sum(index)/2))
      new.cat[count.cat >= cutoff] <- "Undetermined"
      both.means[, sample] <- mean
      both.cats[, sample] <- new.cat
    }
    res.out <- cbind(rownames(res), featPos, res[, c("t", 
                                                     "P.Value", "adj.P.Val", "logFC")], 2^(-res$logFC), 
                     both.means, both.cats)
    colnames(res.out) <- c("genes", "feature.pos", "t.test", 
                           "p.value", "adj.p.value", "ddCt", "FC", "meanTarget", 
                           "meanCalibrator", "categoryTarget", "categoryCalibrator")
    if (sort) 
      res.out <- res.out[order(res.out$adj.p.value), ]
    out[[coef]] <- res.out
  }
  res <- decideTests(fit2, ...)
  rownames(res) <- rownames(topTable(fit2, sort = "none", n = nrow(fit2)))
  if (ret.fit == TRUE){
    out[["Fit"]] <- fit2   
  }
  if (topTableOut == TRUE){
    out[["topTable"]] <- tt1
  }
  out[["Summary"]] <- res
  out
}

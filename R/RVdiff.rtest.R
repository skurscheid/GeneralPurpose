# ==============================================================================
# Title: RVdiff.rtest.R
# Description:selection of optimal number of variables based on differences of RV coefficient 
# Author:  Pierre Bady <pierre.bady@unil.ch>
# Date : Dec 11, 2012 12:03:40 PM
# Version: 0.1
# Revision: Dec 11, 2012 12:03:40 PM
# Comments: RAS
# License: GPL version 2 or newer
# Copyright (C) 2012  Pierre Bady
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# ==============================================================================

RVdiff.rtest <- function (df1, df2,df3,df4=NULL, nrepet = 99,center=TRUE){
	if (!is.data.frame(df1)) 
		stop("data.frame expected")
	if (!is.data.frame(df2)) 
		stop("data.frame expected")
	if (!is.data.frame(df3)) 
		stop("data.frame expected")
	l1 <- nrow(df1)
	if (nrow(df2) != l1) 
		stop("Row numbers are different")
	if (nrow(df3) != l1) 
		stop("Row numbers are different")
	if (any(row.names(df2) != row.names(df1))) 
		stop("row names are different")
	if(is.null(df4))
		df4 <- df3
	if (!is.data.frame(df4)) 
		stop("data.frame expected")	
	if (nrow(df4) != l1) 
		stop("Row numbers are different")	
	if(center){    
		X <- scale(df1, scale = FALSE)
		Y <- scale(df2, scale = FALSE)
		Z <- scale(df3, scale = FALSE)
		W <- scale(df4, scale = FALSE)
	}else{
		X <- df1
		Y <- df2
		Z <- df3
		W <- df4
	}
	X <- X/(sum(svd(X)$d^4)^0.25)
	Y <- Y/(sum(svd(Y)$d^4)^0.25)
	Z <- Z/(sum(svd(Z)$d^4)^0.25)
	W <- W/(sum(svd(W)$d^4)^0.25)	
	X <- as.matrix(X)
	Y <- as.matrix(Y)
	Z <- as.matrix(Z)
	W <- as.matrix(W)	
	# absolute difference between the two observed RVs
	obs <- abs(sum(svd(t(X) %*% Z)$d^2)-sum(svd(t(Y) %*% W)$d^2))
	if (nrepet == 0) 
		return(obs)
	perm <- matrix(0, nrow = nrepet, ncol = 1)
	# empricial distribution of the absolute difference between two 'random' RVs
	perm <- apply(perm, 1, function(x) abs(sum(svd(t(X) %*% Z[sample(l1),])$d^2)-sum(svd(t(Y) %*% W[sample(l1),])$d^2)))
	w <- as.randtest(obs = obs, sim = perm, call = match.call(),alter="greater")
	return(w)
}

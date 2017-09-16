# Copyright (C) 2014 - 2017
#
# This file is part of pdc R Methods Package
#
# The `pdc` R package is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# The `pdc` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



#' @title Compute Ordering Sequence
#'
#' @description Compute variable ordered sequence based the greedy algorithm described
#' in  Avella Medina et al., (2017).
#' @param y A \code{vector} corresponding to the response vector of dimension \eqn{n}.
#' @param X A \code{matrix} corresponding to the design matrix of dimension \eqn{n x p}.
#' @param intercept A \code{logical} value used to indicate if an intercept should be included
#' in the model.
#' @return A \code{dataframe} with the following structure:
#' \describe{
#' \item{seqence}{Sequence of ordered variables}
#' \item{yhat}{Predictions matrix. The \eqn{j}-th column of this matrix
#' corresponds to the predictions made from the model including the first \eqn{j}
#' variables (according the sequence of ordered variables).}
#' }
#' @author Marco Andres Avella Medina and Stéphane Guerrier
#' @references "\emph{A Prediction Divergence Criterion for Model Selection}", Avella
#' Medina, M., Guerrier, S. & Victoria-Feser, M.-P. ArXiv link:
#' \href{https://arxiv.org/abs/1511.04485}{https://arxiv.org/abs/1511.04485}
#'
order_variables <- function(y, X, intercept = FALSE, adaptive = TRUE){

  # Compute sample size and number of varibales
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  # Center covariates
  if (intercept == TRUE){
    a0 <- mean(y)
    y0 <- y-rep(a0,n)
    mu <- t(rep(1/n,n))%*%X
    sX <- X-rep(1,n)%*%mu
    D <- sqrt(diag(t(sX)%*%sX/n) )
    sX <- sX%*%diag(1/D)
  }else{
    sX <- X
    y0 <- y
    a0 <- 0
  }

  if (adaptive == TRUE){
    w <- as.vector(abs(solve(t(sX)%*%sX,t(sX)%*%y0)))
    norm.col <- as.vector(sqrt(t(rep(1,n))%*%sX^2))
    wX <- X%*%diag(w/norm.col)
  }else{
    wX <- X
  }

  # Start Greedy search
  seq <- NULL
  ind0 <- ind <- 1:p
  beta <- rep(0, p)
  yhat <- matrix(NA, n, p)

  for(i in 1:p){
    r <- y0-wX%*%beta
    ind.i <- which(abs(t(wX[, ind0])%*%r) == max(abs(t(wX[, ind0])%*%r)))
    seq[i] <- ind0[ind.i]
    ind0 <- ind[-seq]
    beta[seq] <- solve(t(wX[, seq])%*%wX[ ,seq], t(wX[, seq])%*%y0)
    yhat[ ,i] <- sX[ ,seq]%*%solve(t(sX[ ,seq])%*%sX[, seq], t(sX[, seq])%*%y0) + rep(a0, n)
  }

  # Return sequence
  out <- NULL
  out <- list(yhat = yhat, sequence = seq)
  out
}

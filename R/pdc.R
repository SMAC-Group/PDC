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


#' @title Model Selection based on PDC algorithm
#'
#' @description Model selection based on PDC algorithm for the case \eqn{n > p}
#' see Avella Medina et al., (2017) for more details.
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
#' @references "A Prediction Divergence Criterion for Model Selection", Avella Medina,
#' M., Guerrier, S. & Victoria-Feser, M.-P. ArXiv link:
#' \href{https://arxiv.org/abs/1511.04485}{https://arxiv.org/abs/1511.04485}

pdc <- function(y, X, intercept = FALSE){

  # Intial checks
  n   <- length(y)
  n.x <- dim(X)
  p   <- n.x[2]

  # Check dimensions
  if (n != n.x[1]){ stop("Incompatible dimensions") }
  if (p >= n){ stop("p >= n") }

  # Compute sequence
  seq.ord <- order_variables(y, X, intercept = intercept)
  yhat <- seq.ord$yhat

  # Compute estimated residual variance
  sig2 <- 1/(n - p)*L2.2(y - yhat[, p])

  # Compute PDC criterion
  pdc.crit <- rep(NA, (p-1))
  for (i in 1:(p-1)){
    pdc.crit[i] <- L2.2(yhat[, i] - yhat[, i+1]) + 2*sig2*i
  }

  # Selected model
  seq.select <- sort(seq.ord$sequence[1:which.min(pdc.crit)])

  # Output
  beta.pdc <- rep(0, p)
  if (intercept == FALSE){
    inter <- as.vector(lm(y ~ X[, seq.select] - 1)$coef)
    beta.pdc[seq.select] <- inter
    intercept <- 0
  }else{
    inter <- as.vector(lm(y ~ X[, seq.select])$coef)
    beta.pdc[seq.select] <- inter[-1]
    intercept <- inter[1]
  }

  out <- list(beta.hat = beta.pdc, y.hat = yhat[,length(seq.select)],
             variables = seq.select, intercept = intercept)
  out
}

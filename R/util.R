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


#' Squared L2 norm
#'
#' @param y vector
#' @return squared L2 norm of \code{y}
#' @examples
#' x = 1:10
#' L2.2(x)
L2.2 <- function(x){sum(x^2)}

#' Get ggplot2-like colors
#'
#' @param n number of colors.
#' @param alpha transparency.
#' @return list of colors
#' @export
ggplot_like_colors <- function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}


#' Estimator of the mean squared error
#'
#' @description The function compute the empirical version of \eqn{E[ (x - \mu)^2 ]}
#' @param x vector of random variables
#' @param mu true value
#' @return Empirical mean squared error: \eqn{1/n \sum (xi - \mu)^2}
#' @examples
#' set.seed(1982)
#' x = rnorm(1000)
#' MSE(x = x, mu = 0)
MSE <- function(x,mu){
  (mean(x)-mu)^2 + var(x)
}

# find selected variables
find.select <- function(x){
  return((1:length(x))*(x != 0))
}


is.correct <- function(x){
  if( sum((true.beta != 0) == (x != 0)) == length(x)){
    return(1)
  }else{
    return(0)
  }
}


hit <- function(x){
  sum((x != 0)[true.beta != 0])
}


is.pf = function(x){
  sum((x != 0)[true.beta == 0])
}

nb.regress = function(x){
  return(sum(x != 0))
}

construct.cor = function(cor.in = 0.5, nb.reg = 8){
  out = matrix(NA,nb.reg,nb.reg)
  for (i in 1:nb.reg){
    for (j in 1:nb.reg){
      out[i,j] = cor.in^(abs(i-j))
    }
  }
  return(out)
}


get.X <- function(n = 100, seed = 1, cor.in = 0.5, nb.reg = 8, cor = "abs"){

  # Set seed
  set.seed(seed)

  # Sigma
  Sigma = construct.cor(cor.in, nb.reg = nb.reg)

  # Simulate RMN
  X = rmvnorm(n = n, mean = rep(0,nb.reg) , sigma = Sigma)

  # Return X
  return(X)
}

is.included <- function(x){
  if( sum((true.beta != 0)[true.beta != 0] == (x != 0)[true.beta != 0]) == length((true.beta != 0)[true.beta != 0])){
    return(1)
  }else{
    return(0)
  }
}


boot.med <- function(x, B = 500){
  med = rep(NA,B)
  n = length(x)
  for (i in 1:B){
    med[i] = median(x[sample(1:n,n, replace = TRUE)])
  }
  mean.med = mean(med)
  return(sqrt(sum((med - mean.med)^2)/(B-1)))
}






#' @export
STEP.AIC = function(y, X, method = "AIC", intercept = FALSE, adaptive = TRUE){
  # y = response vector, dim = n x 1
  # X = design matrix, dim = n x p

  # Intial checks
  n = length(y)
  n.x = dim(X)
  p = n.x[2]
  if (n != n.x[1]){ stop("Incompatible dimensions") }
  if (p >= n){ stop("p >= n") }
  if (sum(method == c("AIC","BIC","HQ")) == 0){ stop("Method not supported")}

  # Compute sequence
  seq.ord = order_variables(y, X, intercept = intercept, adaptive = adaptive)
  yhat = seq.ord$yhat

  # Compute penality
  if (method == "AIC"){pen = 2*(1:p + 1)/n}

  if (method == "BIC"){pen = (1:p)*log(n)/n}

  if (method == "HQ"){pen = 2*(1:p)*log(log(n))/n}

  # Compute criterion
  crit = rep(NA,p)
  for (i in 1:p){
    crit[i] = log(1/(n - i)*L2.2(y - yhat[,i])) + pen[i]
  }

  # Select model
  seq.select = sort(seq.ord$sequence[1:which.min(crit)])

  # Output
  beta.meth = rep(0, p)
  if (intercept == FALSE){
    inter = as.vector(lm(y ~ X[,seq.select] - 1)$coef)
    beta.meth[seq.select] = inter
    intercept = 0
  }else{
    inter = as.vector(lm(y ~ X[,seq.select])$coef)
    beta.meth[seq.select] = inter[-1]
    intercept = inter[1]
  }

  out = list(beta.hat = beta.meth,
             y.hat = yhat[,length(seq.select)],
             variables = seq.select,
             intercept = intercept)
  out
}

#' @export
enet.get.lambda = function(y,X){
  # Initial grid
  grid.val = c(0,0.01,0.1,1,10,100,1000)
  cv.val = rep(1,length(grid.val))

  # CV
  for (i in 1:length(grid.val)){
    cv = cv.enet(X,y,lambda = grid.val[i], s = seq(0,1,length=100), mode = "fraction", plot.it = FALSE, se = FALSE,intercept = FALSE, normalize = FALSE)
    cv.val[i] = min(cv$cv)
  }

  min.index = which.min(cv.val)
  ord.index = order(cv.val)
  if (min.index == 1){
    dw = 0
    up = grid.val[1]
  }else{
    if (min.index == length(grid.val)){
      dw = grid.val[(length(grid.val) - 1)]
      up = grid.val[(length(grid.val))]
    }
    else{
      dw = grid.val[(min.index - 1)]
      up = grid.val[(min.index + 1)]
    }
  }

  new.grid = s = seq(dw, up, length = 4)
  cv.val = rep(1,length(new.grid ))

  # CV
  for (i in 1:length(grid.val)){
    cv = cv.enet(X,y,lambda = new.grid[i], s = seq(0,1,length=10), mode = "fraction", plot.it = FALSE, se = FALSE)
    cv.val[i] = min(cv$cv)
  }

  # Best lambda
  return(new.grid[which.min(cv.val)])
}




#' @export
enet.get.lambda.fast = function(y,X){
  # Initial grid
  grid.val = c(0.1,1,10)
  cv.val = rep(1,length(grid.val))

  # CV
  for (i in 1:length(grid.val)){
    cv = cv.enet(X,y,lambda = grid.val[i], s = seq(0,1,length=100), mode = "fraction", plot.it = FALSE, se = FALSE,intercept = FALSE, normalize = FALSE)
    cv.val[i] = min(cv$cv)
  }

  grid.val[which.min(cv.val)]
}



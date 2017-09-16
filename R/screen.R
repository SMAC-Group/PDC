#' @export
sure_screen = function(y, X, gamma = 0.8){
  n = nrow(X)
  p = ncol(X)
  k = floor(n^gamma)
  abs.beta = rep(NA,p)
  for (i in 1:p){
    abs.beta[i] = abs(lm(y ~ X[,i])$coefficients[2])
  }
  ranking=order(abs.beta,decreasing=T)
  ranking[1:k]
}

#' @export
sure_screen_perturbed = function(y, X, B = 100, gamma = 0.8){
  n = nrow(X)
  p = ncol(X)
  k = floor(n^gamma)
  beta = matrix(NA,B, p)
  for (j in 1:B){
    for (i in 1:p){
      beta[j,i] = as.numeric(lm(y ~ X[,i], weights = rexp(n,1))$coefficients[2])
    }
  }
  ranking = order(abs(apply(beta,2,mean)),decreasing=T)
  ranking[1:k]
}

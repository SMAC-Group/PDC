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



#' @title Simulation y and X
#'
#' @description This function allows to simulate \eqn{y = X \beta + \epsilon}, where
#' each row of the design matrix is simulated multivariate normal with unit variance
#' and AR1-like correlation structure (i.e. rho^{|i - j|}) and \eqn{\epsilon} is
#' simulated from a zero-mean normal distribution with variance \eqn{\sigma^2} which
#' is user specific.
#'
#' This function was used in the simulation study of Avella Medina et al. (2017).
#' @param beta A \code{vector} denoting the true parameter vector of dimension \eqn{p}.
#' @param n An \code{integer} denoting the sample size (default value = 100).
#' @param sigma2 A \code{double} representing the residual variance \eqn{\epsilon}.
#' @param cor_design A \code{double} used for the correlation structure of the \code{X}
#' matrix (default value = 0.5).
#' @param seed A \code{scalar} using as "seed" for random number generation
#' (default value = 1982).
#' @return \code{$y}: A \code{vector} corresponding to the simulated response vector \eqn{y}.
#' @return \code{$X}: A \code{matrix} corresponding to the simulated design matrix \eqn{X}.
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
#' @references "\emph{A Prediction Divergence Criterion for Model Selection}", Avella Medina,
#' M., Guerrier, S. & Victoria-Feser, M.-P. ArXiv link:
#' \href{https://arxiv.org/abs/1511.04485}{https://arxiv.org/abs/1511.04485}
simulate_data <- function(beta, n = 100, sigma2 = 1,
                          cor_design = 0.5, seed = 1982, intercept = 0){

  if (abs(cor_design) > 0.999999){stop("Correlation in cor_design is outside of the
                                       range (-1,1).")}
  set.seed(seed = seed)
  p <- length(beta)
  X <- get_X(cor_design = cor_design, n = n, seed = seed, p = p)
  y <- X%*%beta + rnorm(n, sd = sigma2) + intercept

  out <- list(y = y, X = X)
  out
}


#' @title Covariance matrix used for X
#'
#' @description Construction of the covariance matrix used for \eqn{X}.
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
construct_cov_matrix <- function(cor_design, p){
  out <- matrix(NA, p, p)

  for (i in 1:nb.reg){
    for (j in 1:nb.reg){
      out[i,j] <- cor_design^(abs(i-j))
    }
  }
  out
}

#' @title Simulate design matrix X
#'
#' @description Simulation of the design matrix. Each row of this matrix is simulated
#' multivariate normal with unit variance and AR1-like correlation (i.e. rho^{|i - j|})
#' structure.
#' @param n A \code{scalar} denoting the sample size (default value = 100).
#' @param cor_design A \code{scalar} denoting the correlation used to construct \code{X}.
#' @param p An \code{integer} denoting the number of parameters (i.e. length of the vector \eqn{\beta}).
#' @return A \code{matrix} corresponding to the simulated design matrix \eqn{X}.
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
get_X = function(n, seed, cor_design, p){
  # Set seed
  set.seed(seed)

  # Sigma
  Sigma = construct_cov_matrix(cor_design, nb.reg = p)

  # Simulate RMN
  X = rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)
  X
}

#' @title Num
#'
#' @description somehting
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
hit = function(x, beta){
  sum((x != 0)[beta != 0])
}


#' @title Find index of selected variables
#'
#' @description This function is used in the simulation study
#' of Avella Medina et al. (2017) and allow to find non-zero indices.
#' @param x A \code{vector} of estimated parameters.
#' @return A \code{vector} containing the indices of the selected variables.
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
find_select <- function(x){ (1:length(x))*(x != 0) }

#' @title Verify if the correct model is selected
#'
#' @description This function is used in the simulation study
#' of Avella Medina et al. (2017) and allows to check if the "true"
#' model is selected.
#' @param x A \code{vector} of estimated parameters.
#' @param beta A \code{vector} containing the "true" parameter vector.
#' @return A \code{logical} value (0 = true model NOT selected;
#' 1 = true model selected).
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
is_correct <- function(x, beta){
  if( sum((beta != 0) == (x != 0)) == length(x)){
    return(1)
  }else{
    return(0)
  }
}

#' @title Number of correct variables selected
#'
#' @description This function is used in the simulation study
#' of Avella Medina et al. (2017) and allows to compute the number
#' "true" variables contained in a candidate model.
#' @param x A \code{vector} of estimated parameters.
#' @param beta A \code{vector} containing the "true" parameter vector.
#' @return An \code{integer} corresponding to the number of correct
#' variables selected.
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
nb_correct <- function(x, beta){ sum((x != 0)[beta != 0]) }


#' @title Number of non-significant variables selected
#'
#' @description This function is used in the simulation study
#' of Avella Medina et al. (2017) and allows to compute the number
#' non-significant variables contained in a candidate model.
#' @param x A \code{vector} of estimated parameters.
#' @param beta A \code{vector} containing the "true" parameter vector.
#' @return An \code{integer} corresponding to the number of
#' non-significant variables selected.
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
nb_false_pos <- function(x, beta){ sum((x != 0)[beta == 0]) }


#' @title Number of false positive
#'
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
nb_selected <- function(x){ sum(x != 0) }

#' @title Number of false positive
#'
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
nb_signif <- function(x, beta){ sum((x != 0)[beta != 0]) }



#' @title Verify if the selected model is included in the correct model
#'
#' @description This function is used in the simulation study
#' of Avella Medina et al. (2017) and allows to check if the correct model
#' is nested within a candidate model.
#' @param x A \code{vector} of estimated parameters.
#' @param beta A \code{vector} containing the "true" parameter vector.
#' @return A \code{logical} value (0 = true model is NOT included;
#' 1 = true model included).
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
is_included <- function(x, beta){
  if( sum((beta != 0)[beta != 0] == (x != 0)[beta != 0]) == length((beta != 0)[beta != 0])){
    return(1)
  }else{
    return(0)
  }
}


#' @title Bootstrap estimation of the standard error of the median
#'
#' @description This function is used in the simulation study
#' of Avella Medina et al. (2017) and allows to compute various standard
#' errors.
#' @param x A \code{vector} of values on which bootstrap procedure will be
#' applied.
#' @param B An \code{integer}  corresponding to then umber of bootstrap replications
#' (default value = 500).
#' @return A \code{double}  corresponding to the estimated standard error.
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
median_se_boostrap <- function(x, B = 500){
  med <- rep(NA, B)
  n <- length(x)

  for (i in 1:B){
    med[i] <- median(x[sample(1:n, n, replace = TRUE)])
  }

  mean.med <- mean(med)
  sqrt(sum((med - mean.med)^2)/(B-1))
}


#' @title Bootstrap estimation of the covariance matrix of median vector
#'
#' @description This function is used in the simulation study
#' of Avella Medina et al. (2017) to compute "confidence regions".
#' @param X A \code{matrix} of dimension \eqn{n x p}, where \eqn{n > p}.
#' @param B A \code{integer} corresponding to the number of
#' bootstrap replications (default value = 500).
#' @return A \code{matrix} corresponding to the estimated covariance matrix.
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
median_multivariate_cov_bootstrap <- function(X, B = 500){
  res <- matrix(NA, B, 2)
  n <- dim(X)[1]
  for (i in 1:B){
    X.star <- X[sample(1:n, n, replace = TRUE), ]
    res[i, ] <- apply(X.star, 2, median)
  }
  cov(res)
}



#' @title Generic function used to run each simulation
#'
#' @description This function is used to perform the simulation study
#' of Avella Medina et al. (2017).
#' @param nb_simu An \code{integer} used to denote the number of Monte Carlo replications.
#' @param true_beta An \code{vector} corresponding to the true parameter vector.
#' @param n An \code{integer} corresponding to the sample size.
#' @param n_star An \code{integer} corresponding to the sample size of the testing set.
#' @param sigma A \code{double} used for the standard deviation of the residuals.
#' @param cor_X A \code{double} used to defined the correlation structure of the design matrix
#' \eqn{X}.
#' @return A \code{matrix} corresponding to the estimated covariance matrix.
#' @export
#' @author Marco Andres Avella Medina and Stéphane Guerrier
run_simulation <- function(nb_simu, true_beta, n, n_star, sigma, cor_X,
                           adaptive = TRUE,
                           omitted_variables = NULL,
                           H1 = 10, H2 = 500,
                           x_inter = NULL){

  # True parameters
  true_beta <- true_beta
  p <- length(true_beta) - length(omitted_variables)
  m <- n + n_star

  # Initialisation
  # Estimated beta
  beta.LS <- beta.lasso <- beta.lasso <- beta.enet <- beta.alasso <- beta.mcp <-
    beta.scad <- beta.pdc <- beta.apdc <- beta.step.AIC <- beta.step.BIC <-
    beta.step.HQ <- matrix(0, nb_simu, p)

  # Prediciton error
  PE.pdc <- PE.apdc <- PE.lasso <- PE.enet <- PE.alasso <- PE.mcp <-
    PE.scad <- PE.step.AIC <- PE.step.BIC <- PE.step.HQ <-
    PE.true <- PE.LS <- rep(NA, nb_simu)

  # Beta MSE
  MSE.pdc <- MSE.apdc <- MSE.lasso <- MSE.enet <- MSE.alasso <-
    MSE.mcp <- MSE.scad <- MSE.step.AIC <- MSE.step.BIC <-
    MSE.step.HQ <- MSE.true <- MSE.LS <- rep(NA, nb_simu)

  # Redefine "true" beta in the case of missing variables
  if (!is.null(omitted_variables)){
    true_beta2 <- 0
    for (i in 1:H1){
      X <- get.X(cor.in = cor_X, n = H2*m, seed = i, nb.reg = length(true_beta))

      if (!is.null(x_inter)){
        X[,omitted_variables] = X[,omitted_variables]*x_inter[1:(H2*m)]
      }

      y <- X%*%true_beta+rnorm(H2*m,sd=sigma)
      true_beta2 <- true_beta2 + lm(y ~ X[,-omitted_variables] - 1)$coef
    }
    true_beta2 <- true_beta2/H1
  }else{
    true_beta2 <- true_beta
  }

  # Start bootstrap
  for (i in 1:nb_simu){

    # Simulate data set
    X <- get.X(cor.in = cor_X, n = m, seed = i, nb.reg = (p+length(omitted_variables)))

    if (!is.null(x_inter)){
      X[,omitted_variables] = X[,omitted_variables]*x_inter[1:m]
    }

    y <- X%*%true_beta+rnorm(m,sd=sigma)

    # Omitted variables?
    if (is.null(omitted_variables)){
      # Training and testing sets
      y_train <- y[1:n]
      y_test  <- y[(n + 1):m]
      X_train <- X[1:n, ]
      X_test  <- X[(n + 1):m, ]
    }else{
      # Training and testing sets
      y_train <- y[1:n]
      y_test  <- y[(n + 1):m]
      X_train <- X[1:n, -omitted_variables]
      X_test  <- X[(n + 1):m, -omitted_variables]
    }


    # Sure Screening if p > n
    if (p > (n - 1)){
      K <- sure_screen(y_train, X_train)
    }else{
      K <- 1:p
    }

    # PDC
    res.pdc        <- pdc(y_train, X_train[, K], intercept = FALSE, adaptive = adaptive)
    beta.pdc[i, K] <- res.pdc$beta.hat
    PE.pdc[i]      <- (mean((X_test%*%beta.pdc[i, ] - y_test)^2))
    MSE.pdc[i]     <- sqrt(sum((beta.pdc[i, ] - true_beta2)^2))

    #res.apdc      <- adaptive.pdc(y_train, X_train)
    #beta.apdc[i, ] <- res.apdc$beta.hat
    #PE.apdc[i]    <- (mean((X_test%*%beta.apdc[i, ]-y_test)^2))
    #MSE.apdc[i]   <- sqrt(sum((beta.apdc[i, ]-true_beta)^2))

    # Step BIC
    beta.step.BIC[i, K] <- STEP.AIC(y_train, X_train[, K], method = "BIC",
                                    adaptive = adaptive)$beta.hat
    PE.step.BIC[i]      <- (mean((X_test%*%beta.step.BIC[i, ] - y_test)^2))
    MSE.step.BIC[i]     <- sqrt(sum((beta.step.BIC[i, ] - true_beta2)^2))

    # Step AIC
    beta.step.AIC[i, K] <- STEP.AIC(y_train, X_train[, K], method = "AIC",
                                   adaptive = adaptive)$beta.hat
    PE.step.AIC[i]      <- (mean((X_test%*%beta.step.AIC[i, ] - y_test)^2))
    MSE.step.AIC[i]     <- sqrt(sum((beta.step.AIC[i, ] - true_beta2)^2))

    # Step HQ
    beta.step.HQ[i, K]  <- STEP.AIC(y_train, X_train[, K], method = "HQ",
                                   adaptive = adaptive)$beta.hat
    PE.step.HQ[i]       <- (mean((X_test%*%beta.step.HQ[i, ] - y_test)^2))
    MSE.step.HQ[i]      <- sqrt(sum((beta.step.HQ[i, ] - true_beta2)^2))

    # Lasso
    fitnet  <- cv.glmnet(y = y_train, x = X_train, intercept = FALSE, standardize = FALSE)
    ind     <- which(fitnet$lambda == fitnet$lambda.min)
    fitnet0 <- glmnet(y = y_train, x = X_train, intercept = FALSE, standardize = FALSE, lambda = fitnet$lamdba)

    beta.lasso[i, ] <- fitnet0$beta[ ,ind]
    PE.lasso[i]     <- (mean((X_test%*%beta.lasso[i, ] - y_test)^2))
    MSE.lasso[i]    <- sqrt(sum((beta.lasso[i, ] - true_beta2)^2))

    # Fit elastic net
    fitnet1 <- cv.glmnet(y = y_train,x = X_train, alpha = 0.25, intercept = FALSE, standardize = FALSE)
    fitnet2 <- cv.glmnet(y = y_train,x = X_train, alpha = 0.5, intercept = FALSE, standardize = FALSE)
    fitnet3 <- cv.glmnet(y = y_train,x = X_train, alpha = 0.75, intercept = FALSE, standardize = FALSE)

    inter <- which.min(c(min(fitnet1$cvm),min(fitnet2$cvm),min(fitnet3$cvm)))

    if (inter == 1){
      ind <- which(fitnet1$lambda == fitnet1$lambda.min)
      fitnet0 <- glmnet(y = y_train, x = X_train, alpha = 0.25,
                        intercept = FALSE, standardize = FALSE,
                        lambda = fitnet1$lamdba)
    }

    if (inter == 2){
      ind <- which(fitnet2$lambda == fitnet2$lambda.min)
      fitnet0 <- glmnet(y = y_train, x = X_train, alpha = 0.5,
                        intercept = FALSE, standardize = FALSE,
                        lambda = fitnet2$lamdba)
    }

    if (inter == 3){
      ind <- which(fitnet3$lambda == fitnet3$lambda.min)
      fitnet0 <- glmnet(y = y_train, x = X_train, alpha = 0.75,
                        intercept=FALSE, standardize = FALSE,
                        lambda = fitnet3$lamdba)
    }

    beta.enet[i, ] <- fitnet0$beta[ ,ind]
    PE.enet[i]     <- (mean((X_test%*%beta.enet[i, ] - y_test)^2))
    MSE.enet[i]    <- sqrt(sum((beta.enet[i, ] - true_beta2)^2))


    # MCP
    cvfit <- cv.ncvreg(X_train,y_train, familiy = "gaussian", intercept = FALSE)
    fit <- cvfit$fit
    beta.mcp[i, ] <- as.vector(fit$beta[2:(p+1), cvfit$min])
    PE.mcp[i]   <- (mean((fit$beta[1] + X_test%*%beta.mcp[i, ] - y_test)^2))
    MSE.mcp[i]  <- sqrt(sum((beta.mcp[i, ] - true_beta2)^2))

    # SCAD
    cvfit <- cv.ncvreg(X_train, y_train, familiy = "gaussian",
                       penalty = "SCAD", intercept = FALSE)
    fit <- cvfit$fit
    beta <- fit$beta[ ,cvfit$min]
    beta.scad[i, ] <- fit$beta[2:(p+1), cvfit$min]
    PE.scad[i]  <- (mean((fit$beta[1] + X_test%*%beta.scad[i, ] - y_test)^2))
    MSE.scad[i] <- sqrt(sum((beta.scad[i, ] - true_beta2)^2))

    # Adaptive lasso
    fit <- adalasso(X_train, y_train, k = 10)
    beta.alasso[i, ] <- fit$coefficients.adalasso
    PE.alasso[i]  <- (mean((X_test%*%beta.alasso[i, ] - y_test)^2))
    MSE.alasso[i]  <- sqrt(sum((beta.alasso[i, ] - true_beta2)^2))

    # Full LS
    mod.full      <- lm(y_train ~ X_train[ ,K] - 1)
    beta.LS[i, K] <- coef(mod.full)
    PE.LS[i]      <- (mean((X_test%*%beta.LS[i, ] - y_test)^2))
    MSE.LS[i]     <- sqrt(sum((beta.LS[i, ] - true_beta2)^2))
  }


  # Output results
  simu_out <- matrix(NA, 10, 10)
  dimnames(simu_out)[[1]] <- c("LS", "AIC", "BIC", "HQ",
                                  "lasso", "PDC", "enet",
                                  "alasso", "MCP", "SCAD")
  dimnames(simu_out)[[2]] <- c("Median MSEy", "SD Med MSEy",
                                  "Median MSEb", "SD Med MSEb",
                                  "correct", "included", "true+",
                                  "false+", "# sel.", "SD # sel.")

  if (!is.null(omitted_variables)){
    true_beta <- true_beta2
  }
  simu_out[1, ]  <- fill_result_matrix(PE.LS, MSE.LS, beta.LS, true_beta)
  simu_out[2, ]  <- fill_result_matrix(PE.step.AIC, MSE.step.AIC, beta.step.AIC, true_beta)
  simu_out[3, ]  <- fill_result_matrix(PE.step.BIC, MSE.step.BIC, beta.step.BIC, true_beta)
  simu_out[4, ]  <- fill_result_matrix(PE.step.HQ, MSE.step.HQ, beta.step.HQ, true_beta)
  simu_out[5, ]  <- fill_result_matrix(PE.lasso, MSE.lasso, beta.lasso, true_beta)
  simu_out[6, ]  <- fill_result_matrix(PE.pdc, MSE.pdc, beta.pdc, true_beta)
  simu_out[7, ]  <- fill_result_matrix(PE.enet, MSE.enet, beta.enet, true_beta)
  simu_out[8, ]  <- fill_result_matrix(PE.alasso, MSE.alasso, beta.alasso, true_beta)
  simu_out[9, ]  <- fill_result_matrix(PE.mcp, MSE.mcp, beta.mcp, true_beta)
  simu_out[10, ] <- fill_result_matrix(PE.scad, MSE.scad, beta.scad, true_beta)


  output = list(summary_table = simu_out,
                PE.scad = PE.scad,
                PE.enet = PE.enet,
                PE.alasso = PE.alasso,
                PE.lasso = PE.lasso,
                PE.mcp = PE.mcp,
                PE.pdc = PE.pdc,
                PE.step.AIC = PE.step.AIC,
                PE.step.BIC = PE.step.BIC,
                PE.step.HQ = PE.step.HQ,
                PE.LS = PE.LS,
                MSE.scad = MSE.scad,
                MSE.enet = MSE.enet,
                MSE.alasso = MSE.alasso,
                MSE.lasso = MSE.lasso,
                MSE.mcp = MSE.mcp,
                MSE.pdc = MSE.pdc,
                MSE.step.AIC = MSE.step.AIC,
                MSE.step.BIC = MSE.step.BIC,
                MSE.step.HQ = MSE.step.HQ,
                MSE.LS = MSE.LS,
                beta.LS = beta.LS,
                beta.lasso = beta.lasso,
                beta.enet = beta.enet,
                beta.alasso = beta.alasso,
                beta.mcp = beta.mcp,
                beta.scad = beta.scad,
                beta.pdc = beta.pdc,
                beta.step.AIC = beta.step.AIC,
                beta.step.BIC = beta.step.BIC,
                beta.step.HQ = beta.step.HQ,
                beta0 = true_beta)
  class(output) = "simulation"
  output
}

#' @export
add_point <- function(PE, myMSE, couleur, couleur_trans, point_pch){
  point.meth <- c(median(PE), median(myMSE))
  mat.meth <- as.matrix(cbind(PE, myMSE))
  cov.meth <- boot.conf.region(mat.meth, B = 1000)
  el.meth <- ellipse(cov.meth, centre = point.meth)
  points(point.meth[1],point.meth[2], pch = point_pch,
         col = couleur, cex = 1.5)
  lines(el.meth, col = couleur)
  polygon(el.meth, col = couleur_trans, border = NA)
}

#' @export
fill_result_matrix <- function(PE.meth, MSE.meth, beta.meth, true_beta){
  output <- rep(NA, 10)
  output[1] <- median(PE.meth)
  output[2] <- boot.med(PE.meth)
  output[3] <- median(MSE.meth^2)
  output[4] <- boot.med(MSE.meth^2)
  output[5] <- 100*mean(apply(beta.meth, 1, is_correct, beta = true_beta))
  output[6] <- 100*mean(apply(beta.meth, 1, is_included, beta = true_beta))
  output[7] <- mean(apply(beta.meth, 1, nb_signif, beta = true_beta))
  output[8] <- mean(apply(beta.meth, 1, nb_false_pos, beta = true_beta))
  output[9] <- mean(apply(beta.meth, 1, nb_selected))
  output[10] <- boot.mean(apply(beta.meth, 1, nb_selected))
  output
}

#' @export
boot.conf.region = function(X, B = 500){
  res = matrix(NA,B,2)
  n = dim(X)[1]
  for (i in 1:B){
    X.star = X[sample(1:n, n, replace = TRUE),]
    res[i,1] = median(X.star[,1])
    res[i,2] = mean(X.star[,2])
  }
  cov(res)
}

#' @export
add_legend = function(x, y, xlim, ylim){
  (x > xlim[1] & x < xlim[2]) & (y > ylim[1] & y < ylim[2])
}

#' @export
boot.mean <- function(x, B = 500){
  med = rep(NA,B)
  n = length(x)
  for (i in 1:B){
    med[i] = mean(x[sample(1:n, n, replace = TRUE)])
  }
  mean.med = mean(med)
  return(sqrt(sum((med - mean.med)^2)/(B-1)))
}

#' @export
make_correct_graph = function(simu, main = "Frequency of True Model Selection",
                              ymax = NULL, add_max = NULL, coef = 1.05,
                              cex.names = 1, cex.axis = 1, cex = 1){
  correct = simu$summary_table[,5]
  coleur = ggplot_like_colors(10, alpha = 0.7)
  ord = c(10, 8, 5, 9, 6, 7, 2, 3, 4)
  if (is.null(ymax)){
    ymax = max(correct)*coef
  }
  a = barplot(correct[ord], col = coleur[1:9], main = main, ylim = c(0,ymax),
          cex.names = cex.names, cex.axis = cex.axis,
          cex.main = cex)
  if (!is.null(add_max)){
    points(a[which.max(correct[ord])], max(correct[ord]), col = "red", pch = 16, cex = 3)
  }
}

#' @export
make_included_graph = function(simu, main = "Frequency Model Including True Model", ymax = NULL){
  correct = simu$summary_table[,6]
  coleur = ggplot_like_colors(10, alpha = 0.7)
  ord = c(10, 8, 5, 9, 6, 7, 2, 3, 4)
  if (is.null(ymax)){
    ymax = max(correct[ord])
  }
  barplot(correct[ord], col = coleur[1:9], main = main, ylim = c(0,ymax))
}

#' @export
make_selected_graph = function(simu, main = "Average Number of Variables Selected",
                               ymax = NULL, add_min = NULL,
                               cex.names = 1, cex.axis = 1, cex = 1){
  correct = simu$summary_table[,9]
  coleur = ggplot_like_colors(10, alpha = 0.7)
  ord = c(10, 8, 5, 9, 6, 7, 2, 3, 4)
  if (is.null(ymax)){
    ymax = max(correct[ord])
  }
    a = barplot(correct[ord], col = coleur[1:9], main = main, ylim = c(0,ymax),
          cex.names = cex.names, cex.axis = cex.axis,
          cex.main = cex)
  if (!is.null(add_min)){
    points(a[which.min(correct[ord])], min(correct[ord]), col = "red", pch = 16, cex = 3)
  }
}

#' @export
make_main_graph = function(simu, xlim = NULL, ylim = NULL, xlab = NULL,
                           ylab = NULL, title = NULL, leg_pos = NULL,
                           cex.legend = 1.2, cex.axis = 1,
                           cex.lab = 1,
                           cex.main = 1){
  # Make graph
  coleur = ggplot_like_colors(10)
  coleurTrans = ggplot_like_colors(10, alpha = 0.08)
  point.pch = c(15:18,21:25,7)

  if (is.null(xlim)){
    xlim <- range(c(median(simu$PE.scad), median(simu$PE.pdc),
                    median(simu$PE.lasso), median(simu$PE.enet),
                    median(simu$PE.alasso), median(simu$PE.mcp),
                    median(simu$PE.step.AIC),
                    median(simu$PE.step.BIC), median(simu$PE.step.HQ),
                    median(simu$PE.LS)))
  }

  if (is.null(ylim)){
    ylim <- range(c(median(simu$MSE.scad), median(simu$MSE.pdc),
                    median(simu$MSE.lasso), median(simu$MSE.enet),
                    median(simu$MSE.alasso), median(simu$MSE.mcp),
                    median(simu$MSE.step.AIC),
                    median(simu$MSE.step.BIC), median(simu$MSE.step.HQ),
                    median(simu$MSE.LS)))
  }

  if (is.null(xlab)){ xlab = "Med(PE)" }

  if (is.null(ylab)){ ylab = "Med(MSE)" }

  if (is.null(title)){ title= "Simulation"}

  #plot(NA, xlim = xlim, ylim = ylim,
  #     xlab = xlab, ylab = " ",
  #     main = title, cex.lab = cex.lab,
  #     cex.main = cex.main,
  #     cex.axis = cex.axis)


  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab,
       ylab = ylab, main = title,
       cex.axis = cex.axis,
       cex.lab = cex.lab,
       cex.main = cex.main)

  #mtext(ylab, side = 2, line = 2.5, cex = cex.lab)

  grid()

  # Add points
  add_point(simu$PE.scad, simu$MSE.scad, coleur[1], coleurTrans[1], point.pch[1])
  add_point(simu$PE.alasso, simu$MSE.alasso, coleur[2], coleurTrans[2], point.pch[2])
  add_point(simu$PE.lasso, simu$MSE.lasso, coleur[3], coleurTrans[3], point.pch[3])
  add_point(simu$PE.mcp, simu$MSE.mcp, coleur[4], coleurTrans[4], point.pch[4])
  add_point(simu$PE.pdc, simu$MSE.pdc, coleur[5], coleurTrans[5], point.pch[5])
  add_point(simu$PE.enet, simu$MSE.enet, coleur[6], coleurTrans[6], point.pch[6])
  add_point(simu$PE.step.AIC, simu$MSE.step.AIC, coleur[7], coleurTrans[7], point.pch[8])
  add_point(simu$PE.step.BIC, simu$MSE.step.BIC, coleur[8], coleurTrans[8], point.pch[8])
  add_point(simu$PE.step.HQ, simu$MSE.step.HQ, coleur[9], coleurTrans[9], point.pch[9])
  add_point(simu$PE.LS, simu$MSE.LS, coleur[10], coleurTrans[10], point.pch[10])

  # Check which lab to display in legend
  leg_indic = rep(NA, 10)
  leg_indic[1]  <- add_legend(median(simu$PE.scad), median(simu$MSE.scad), xlim, ylim)
  leg_indic[2]  <- add_legend(median(simu$PE.alasso), median(simu$MSE.alasso), xlim, ylim)
  leg_indic[3]  <- add_legend(median(simu$PE.lasso), median(simu$MSE.lasso), xlim, ylim)
  leg_indic[4]  <- add_legend(median(simu$PE.mcp), median(simu$MSE.mcp), xlim, ylim)
  leg_indic[5]  <- add_legend(median(simu$PE.pdc), median(simu$MSE.pdc), xlim, ylim)
  leg_indic[6]  <- add_legend(median(simu$PE.enet), median(simu$MSE.enet), xlim, ylim)
  leg_indic[7]  <- add_legend(median(simu$PE.step.AIC), median(simu$MSE.step.AIC), xlim, ylim)
  leg_indic[8]  <- add_legend(median(simu$PE.step.BIC), median(simu$MSE.step.BIC), xlim, ylim)
  leg_indic[9]  <- add_legend(median(simu$PE.step.HQ), median(simu$MSE.step.HQ), xlim, ylim)
  leg_indic[10] <- add_legend(median(simu$PE.LS), median(simu$MSE.LS), xlim, ylim)

  leg_lab = c("SCAD", "alasso", "lasso", "MCP", "PDC",
              "enet", "AIC", "BIC", "HQ", "LS")
  leg_pt = c(1.5, 1.5, 1.4, 1.5, 1.8, rep(1.4, 5))

  if (sum(leg_indic) > 0){

    if (is.null(leg_pos)){
      leg_pos = "topleft"
    }

    legend(leg_pos, leg_lab[leg_indic],
           pch = point.pch[leg_indic], col = coleur[leg_indic],
           pt.cex = leg_pt[leg_indic],
           bty = "n", cex = cex.legend)
  }

}

#' @export
plot.simulation = function(simu,
                           main.correct = "Frequency of True Model Selection",
                           ymax.correct = NULL,
                           add_max = TRUE,
                           coef.correct = 1.05,
                           cex.correct = 1,
                           cex.axis.correct = 1,
                           cex.names.correct = 1,
                           cex.selected = 1,
                           cex.axis.selected = 1,
                           cex.names.selected = 1,
                           main.selected = "Average Number of Variables Selected",
                           ymax.selected = NULL,
                           add_min = TRUE,
                           xlim = NULL,
                           ylim = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           title = NULL,
                           leg_pos = NULL,
                           add_main = TRUE,
                           add_correct = TRUE,
                           add_selected = TRUE,
                           cex.legend = 1.2,
                           cex.axis = 1,
                           cex.lab = 1,
                           cex.main = 1){

  if (add_correct || add_selected){
    layout(matrix(c(1,1,1,1,2,3,2,3), 2, 4, byrow = FALSE))
  }


  if (add_main){
    make_main_graph(simu, xlim = xlim,
                    ylim = ylim,
                    xlab = xlab,
                    ylab = ylab,
                    title = title,
                    leg_pos = leg_pos,
                    cex.legend = cex.legend,
                    cex.axis = cex.axis,
                    cex.lab = cex.lab,
                    cex.main = cex.lab)
  }

  if (add_correct){
    make_correct_graph(simu, main = main.correct,
                       coef = coef.correct,
                       add_max = add_max,
                       cex.names = cex.names.correct,
                       cex.axis = cex.axis.correct,
                       cex = cex.correct)
  }

  if (add_selected){
    make_selected_graph(simu, add_min = add_min,
                        main = main.selected,
                        ymax = ymax.selected,
                        cex.names = cex.names.selected,
                        cex.axis = cex.axis.selected,
                        cex = cex.selected)
  }

}




#' Elastic Net Coordinate Descent
#'
#' This function estimates coefficients for elastic net and also has function to visualize solution path.
#'
#' @param Y = response vector
#' @param X = covariate matrix
#' @param alpha = parameter between [0,1]
#' @param B =  initial vector
#' @return a matrix with lambda values in the first column, and corresponding coefficient estimates
#'
#' @keywords elastic net, coordinate descent
#' @export
#' @examples
#' library(MASS)
#'  set.seed(123)
#'  n <- 100
#'  mu <- rep(0, 20)
#'  Sigma <- diag(20)
#' Sigma[1,2] <- Sigma[2,1] <- Sigma[5,6] <- Sigma[6,5] <- 0.8
#' TrueBeta <- c(2, 0, -2, 0, 1, 0, -1, 0, rep(0, times=12))
#' X <- mvrnorm(n, mu, Sigma)
#' Y <- X%*%TrueBeta + rnorm(n)
#'
#'require(ggplot2)
#' require(reshape2)
#' output <- elnet_coord(Y, X, alpha=0, B=rep(0, times=ncol(X)))
#' visualize_coord(output)
#' output <- elnet_coord(Y, X, alpha=0.5, B=rep(0, times=ncol(X)))
#' visualize_coord(output)
#' output <- elnet_coord(Y, X, alpha=1, B=rep(0, times=ncol(X)))
#' visualize_coord(output)
#'
#'
#'
#'

elnet_coord <- function(Y, X, alpha, B){

soft_threshold <- function(z, gamma){
  if (z >0 & gamma < abs(z)){
    output <- z-gamma
  }
  else if (z <= 0 & gamma < abs(z)){
    output <- z+ gamma
  }
  else{ output <- 0}
  return(output)
}

coord_desc <- function(B, X, Y, alpha, lambda, n){
  for (j in 1:20){
    oldB <- B[j]
    Xpartial <- X[, -j]
    Bpartial <- B[-j]
    rj <- Y - (Xpartial%*%Bpartial)
    num <- t(X[, j])%*%rj
    Bstar <- num/n
    newB <- soft_threshold(Bstar, lambda*alpha)/(1 + lambda*(1-alpha))
    if (abs(newB - oldB) > 1e-3){
      B[j] <- newB
    }
  }
  return(B)
}

  p <- ncol(X)
  n <- nrow(X)
  Y <- scale(Y) #need to standardize variables first
  X <- scale(X)

  if (alpha == 0){
    lmax <- max(abs(t(Y - mean(Y)*(1-mean(Y))) %*% X ))/(0.0005 * n)
  } else {lmax <- max(abs(t(Y - mean(Y)*(1-mean(Y))) %*% X ))/(alpha * n)}
  lgrid <- seq(0.01, lmax, length.out=100) #grid search

  solpath <- matrix(0, length(lgrid), p)

  for (l in 1:length(lgrid)){
    lambda <- lgrid[l]
    niter = 1

    while (niter < 501) {
      newB <- c(coord_desc(B, X, Y, alpha, lambda, n))
      if (niter == 500) {
        warning("Reached 500 iterations")
        break
      }
      niter <- niter + 1
      B <- newB
    }
    solpath[l, ] <- B
  }
  return(data.frame(lgrid, solpath))
}


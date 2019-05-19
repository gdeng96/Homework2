#' Solving Linear Systems
#'
#' This function solves a linear system using Gauss-Seidel, Jacobi method, or Jacobi (parallel) method.
#'
#' @param A = a symmetric matrix
#' @param b = a vector
#' @param x0 = initial vector for iterative process
#' @param method = one of three methods: "Gauss-Seidel", "Jacobi", or "Jacobi-Parallel"
#' @return a list with x, the final estimated vector, and total_iter, the number of iterations
#' @keywords linear system
#' @export
#' @examples
#' A3 = tridiagonal(3, -1, 100)
#' v = rep(c(1, 0), times=50)
#' b3 = A3%*%v
#'
#' solve_ols(A=A3, b=b3, x0=rep(0, nrow(A3)), method="Gauss-Seidel")
#' solve_ols(A=A3, b=b3, x0=rep(0, nrow(A3)), method="Jacobi")
#' solve_ols(A=A3, b=b3, x0=rep(0, nrow(A3)), method="Jacobi-Parallel")
#'

solve_ols<- function(A, b, x0, method, ncores){

library(foreach)
library(doParallel)

GaussSeidel <- function(A, b, x) {
  a <- diag(A)
  diag(A) <- 0
  for (i in 1:length(x)){
    x[i] <- (b[i] - crossprod(A[i, ], x))/a[i]
  }
  return(x)
}

GS_itersolve <- function(A, b, x0, niter = 1000) {
  n = 1
  while (n < niter) {
    x <- c(GaussSeidel(A, b, x0))
    if (n == niter) {
      warning("Reached 1000 iterations")
      break
    }
    n <- n + 1
    x0 <- x
  }
  return(list(x = x, total_iter = (n-1)))
}


Jacobi_itersolve <- function(A, b, x0, eps=10e-7, niter = 1000) {
  nall = nrow(A)
  D  <- diag(diag(A), nall, nall)
  L <- U <- diag(0, nall, nall)
  L[(row(L) - col(L)) == 1] <- -1
  U[(col(U) - row(U)) == 1] <- -1

  n = 1
  Dinv <- solve(D)
  LU <- L+U

  while (n < niter) {
    x <- crossprod(Dinv, (b - crossprod(LU,x0))) #c(Jacobi(L, D, U, b, x0))
    if (n == niter) {
      warning("Reached 1000 iterations")
      break
    }
    n <- n + 1
    x0 <- x
  }
  return(list(x = x, total_iter=(n-1)))
}

Jacobi_paritersolve <- function(A, b, x0, niter = 1000) {
  n = 1
  a <- diag(A)
  diag(A) <- 0

  while (n < niter) {
    x <- foreach(i=1:length(x0), .combine=c)%dopar%{
      (b[i] - crossprod(A[i, ], x0))/a[i]}
    if (n == niter) {
      warning("Reached 1000 iterations")
      break
    }
    n <- n + 1
    x0 <- x

  return(list(x = x, total_iter=(n-1)))
  }
}

  if (method ==  "Gauss-Seidel"){
    output <- GS_itersolve(A, b, x0)
  }
  else if (method  == "Jacobi"){
    output <- Jacobi_itersolve(A, b, x0)
  }

  else if (method  ==  "Jacobi-Parallel"){
    cores=detectCores()
    cl <- makeCluster(cores[1]-2) #not to overload your computer
    registerDoParallel(cl)
    output <- Jacobi_paritersolve(A, b, x0)
    stopCluster(cl)
  }
  else{
    print("Error: Must Specify Method")
  }
  return(output)
}


#' Tridiagonal matrix
#'
#' Generates a symmetric tridiagonal matrix with "offd" on off diagonal, "d" on the diagonal
#' @param d the value on the diagonal
#' @param offd the value along the off diagonal
#' @param n the dimension of the matrix
#' @return A, an n x n tridiagonal matrix
#' @export
#'
tridiagonal <- function(d, offd, n){
  x = diag(d, n, n)
  x[abs(row(x) - col(x)) == 1] <- offd
  return(x)
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ols_solve.R
\name{solve_ols}
\alias{solve_ols}
\title{Solving Linear Systems}
\usage{
solve_ols(A, b, x0, method, ncores)
}
\arguments{
\item{A}{= a symmetric matrix}

\item{b}{= a vector}

\item{x0}{= initial vector for iterative process}

\item{method}{= one of three methods: "Gauss-Seidel", "Jacobi", or "Jacobi-Parallel"}
}
\value{
a list with x, the final estimated vector, and total_iter, the number of iterations
}
\description{
This function solves a linear system using Gauss-Seidel, Jacobi method, or Jacobi (parallel) method.
}
\examples{
A3 = tridiagonal(3, -1, 100)
v = rep(c(1, 0), times=50)
b3 = A3\%*\%v

solve_ols(A=A3, b=b3, x0=rep(0, nrow(A3)), method="Gauss-Seidel")
solve_ols(A=A3, b=b3, x0=rep(0, nrow(A3)), method="Jacobi")
solve_ols(A=A3, b=b3, x0=rep(0, nrow(A3)), method="Jacobi-Parallel")

}
\keyword{linear}
\keyword{system}

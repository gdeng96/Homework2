# Homework2

This is an example R package for 3 functions: ols_solve(), algo_leverage(), and elnet_coord().

## Installing the R package
```
library(devtools)
install_github("gdeng96/Homework2") #install the package
library(Homework2) #loads package
```

## Make sure you also have these other R packages in your local folder

```
library(ggplot2)
library(reshape2)
library(doParallel)
library(foreach)

```
# Tutorial Code

## ols_solve()
This function solves a linear system using Gauss-Seidel, Jacobi method, or Jacobi (parallel) method, for a given tri-diagonal matrix A.

```
A3 = tridiagonal(3, -1, 100)
v = rep(c(1, 0), times=50)
b3 = A3%*%v

solve_ols(A=A3, b=b3, x0=rep(0, nrow(A3)), method="Gauss-Seidel")
solve_ols(A=A3, b=b3, x0=rep(0, nrow(A3)), method="Jacobi")
solve_ols(A=A3, b=b3, x0=rep(0, nrow(A3)), method="Jacobi-Parallel")
```

## algo_leverage()
This function estimates regression coefficients using a randomly selected subsample of an original dataset by two methods: "Unweighted" or "Leveraging". Original paper by Ma and Sun (2014), "Leveraging for big data regression".

```
set.seed(123)
fullX <- rt(500, df=6)
fullY <- -1*fullX + rnorm(100)
fullDF <- data.frame(fullY, fullX)
plot(fullX, fullY, main="Full Sample")

summary(lm(fullY ~ fullX))$coef #True coefficients, compare with uniform/leverage output

algo_leverage(fullY, data.frame(fullX), rsize = 200, method="Uniform")
algo_leverage(fullY, data.frame(fullX), rsize = 200, method="Leverage")
```

## elnet_coord()
This function estimates coefficients for elastic net based on coordinate gradient descent algorithm. 

```
library(MASS)
set.seed(123)
n <- 100
mu <- rep(0, 20)
Sigma <- diag(20)
Sigma[1,2] <- Sigma[2,1] <- Sigma[5,6] <- Sigma[6,5] <- 0.8 
TrueBeta <- c(2, 0, -2, 0, 1, 0, -1, 0, rep(0, times=12))
X <- mvrnorm(n, mu, Sigma)
Y <- X%*%TrueBeta + rnorm(n)

output <- elnet_coord(Y, X, alpha=0, B=rep(0, times=ncol(X)))
visualize_coord(output)
output <- elnet_coord(Y, X, alpha=0.5, B=rep(0, times=ncol(X)))
visualize_coord(output)
output <- elnet_coord(Y, X, alpha=1, B=rep(0, times=ncol(X)))
visualize_coord(output)
```

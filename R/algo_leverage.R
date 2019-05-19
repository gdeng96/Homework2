#' Leveraging
#'
#' This function estimates regression coefficients using a randomly selected subset of the full dataset.
#'
#'
#' @param Y = response vector
#' @param X = covariate data frame
#' @param rsize = the subset sample size
#' @param method = either "Uniform" or "Leverage" for weighting
#' @return output of coefficient estimates, with t-statistic and p-value
#' @keywords leveraging
#' @export
#' @examples
#' set.seed(123)
#' fullX <- rt(500, df=6)
#' fullY <- -1*fullX + rnorm(100)
#' fullDF <- data.frame(fullY, fullX)
#' plot(fullX, fullY, main="Full Sample")
#'
#' summary(lm(fullY ~ fullX))$coef #True coefficients, compare with uniform/leverage output
#'
#' algo_leverage(fullY, data.frame(fullX), rsize = 200, method="Uniform")
#' algo_leverage(fullY, data.frame(fullX), rsize = 200, method="Leverage")
#'

algo_leverage <- function(Y, X, rsize, method){
  if (method == "Uniform"){
    index <- sample(1:length(Y), size=rsize, replace=TRUE, prob=rep(1/length(Y), length(Y)))
    newY <- Y[index]
    newX <- X[index, ]
    subweights <- rep(1/length(Y), length(index))
    linearmodel <- summary(lm(newY ~ newX, weights=1/subweights))
    estimate <- linearmodel$coef
  }

  else if (method == "Leverage"){

    newX <- as.matrix(X)
    H <- newX%*% solve(t(newX)%*%newX) %*% t(newX)
    h <- diag(H)
    normalizedh <- h/sum(h)

    index <- sample(1:length(Y), size=rsize, replace=TRUE, prob=normalizedh)
    newY <- Y[index]
    newX <- X[index, ]
    subweights <- normalizedh[index]
    linearmodel <- summary(lm(newY ~ newX, weights=1/subweights))
    estimate <- linearmodel$coef
  }
  else{
    print ("Error: Method not specified.")
    }
  return(estimate)
}

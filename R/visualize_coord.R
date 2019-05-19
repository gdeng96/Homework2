#' Visualizing Elastic Net
#'
#' This function visualizes the solution path, given by the output of elnet_coord()
#'
#' @param solpath = the output of elnet_coord()
#' @return a plot of the solution path similar to enet package
#' @keywords visualizing solution path
#' @export
#' @examples
#' require(ggplot2)
#' require(reshape2)
#' output <- elnet_coord(Y, X, alpha=0, B=rep(0, times=ncol(X)))
#' visualize_coord(output)
#' output <- elnet_coord(Y, X, alpha=0.5, B=rep(0, times=ncol(X)))
#' visualize_coord(output)
#' output <- elnet_coord(Y, X, alpha=1, B=rep(0, times=ncol(X)))
#' visualize_coord(output)

visualize_coord <- function(solpath){
  library(ggplot2)
  library(reshape2)
  solpath$lgrid <- log(solpath$lgrid)
  df <- melt(solpath, id.vars="lgrid", variable.name="X")
  ggplot(df, aes(lgrid,value)) + geom_line(aes(colour=as.factor(X))) + ggtitle("Solution Path") + xlab("log(lambda)") + ylab("Coefficients")
}

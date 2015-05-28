#' Dimension Reduction Methods for Multivariate Time Series.
#'
#' BigVAR contains a series of functions that allow for the estimation of Penalized Vector Autoregressive models.  This package originated as a 2014 Google "Summer of Code" Project.
#'
#' @author Will Nicholson \email{wbn8@@cornell.edu},
#'
#' @docType package
#' @useDynLib BigVAR
#' @name BigVAR
#' @details To use the facilities of this package, starting with an \eqn{k\times T} multivariate time series and run \code{\link{constructModel}} to create an object of class \code{\link{BigVAR}}.  \code{\link{cv.BigVAR}} creates an object of class \code{\link{BigVAR.results}}, which chooses an optimal penalty parameter based on minimizing h-step ahead forecasts on a specified cross-validation period over a grid of values as well as comparisons against AIC, unconditional mean, and a random walk.  There are plot functions for both BigVAR (\code{\link{plot.BigVAR}}) and BigVAR.results (\code{\link{plot}}) as well as a predict function for BigVAR.results (\code{\link{predict}}).
#' @seealso \code{\link{constructModel}}, \code{\link{cv.BigVAR}}, \code{\link{BigVAR.results}}, \code{\link{plot}}, \code{\link{predict}}
#' @examples 
#' data(Y)
#' head(Y)
#' T1=floor(nrow(Y)/3)
#' T2=floor(2*nrow(Y)/3)
#' m1=constructModel(Y,p=4,struct="None",gran=c(50,10),verbose=FALSE,T1=T1,T2=T2)
#' plot(m1)
#' results=cv.BigVAR(m1)
#' plot(results)
#' predict(results,n.ahead=1)
#' @references Lutkepohl "New Introduction to Multivariate Time Series",
#' William B Nicholson, Jacob Bien, and David S Matteson. "Hierarchical vector autoregression." arXiv preprint arXiv:1412.5250, 2014.
#' William B Nicholson, David S. Matteson, and Jacob Bien (2015), "Structured regularization for large vector
#' autoregressions with exogenous variables," \url{http://www.wbnicholson.com/Nicholsonetal2015.pdf}.

#' @import Rcpp
#' @import RcppArmadillo
#' @import expm
#' @import methods


NULL

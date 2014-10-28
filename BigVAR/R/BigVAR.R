#' Dimension Reduction Methods for Multivariate Time Series.
#'
#' BigVAR contains a series of functions that allow for the estimation of Penalized Vector Autoregressive models.
#'
#' @author Will Nicholson \email{wbn8@@cornell.edu},
#'
#' @docType package
#' @useDynLib BigVAR
#' @name BigVAR
#' @details To use the facilities of this package, starting with an \eqn{k\times T} multivariate time series and run \code{\link{constructModel}} to create an object of class \code{\link{BigVAR}}.  \code{\link{cv.BigVAR}} creates an object of class \code{\link{BigVAR.results}}, which chooses an optimal penalty parameter based on minimizing h-step ahead forecasts on a specified cross-validation period over a grid of values as well as comparisons against AIC, unconditional mean, and a random walk.  There are plot functions for both BigVAR (\code{\link{plot.BigVAR}}) and BigVAR.results (\code{\link{plot}}) as well as a predict function for BigVAR.results (\code{\link{predict}}).
#' @seealso \code{\link{constructModel}}, \code{\link{cv.BigVAR}}, \code{\link{BigVAR.results}}, \code{\link{plot}}), (\code{\link{predict}}
#' @examples 
#' data(Y)
#' Y=Y[1:100,]
#' m1=constructModel(Y,p=4,struct="None",gran=c(50,10),
#' RVAR=FALSE,MN=FALSE,verbose=FALSE,h=1,cv="Rolling")
#' plot(m1)
#' results=cv.BigVAR(m1)
#' plot(results)
#' predict(results,n.ahead=1)
#' @references Lutkepohl "New Introduction to Multivariate Time Series", Nicholson et al (2014)
#' @import Rcpp
#' @import RcppArmadillo
#' @import zoo
#' @import expm
#' @import methods
#' @import lattice


NULL

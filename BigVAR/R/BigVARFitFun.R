# BigVAR fit used in rolling cv, out of sample forecast evaluation and
# BigVAR.est
.BigVAR.fit <- function(group, beta, trainZ, trainY, lambda, tol, p, m = 0, k1, k,
    s = 0, s1 = 0, MN = FALSE, C, intercept = TRUE, separate_lambdas, dual, activeset = NULL,
    starting_eigvals = NULL, groups = NULL, compgroups = NULL, VARX = FALSE, alpha = NULL,
    palpha=NULL, gamma=3) {
    if (is.null(s)) {
        s <- 0
    }
    if (is.null(s1)) {
        s1 <- 0
    }
    if (is.null(m)) {
        m <- 0
    }
    pre_proc <- pre_process(trainY, trainZ, C, MN, intercept)

    if (separate_lambdas) {
        if (is.vector(lambda)) {
            lambda <- matrix(lambda, nrow = 1)
        }
        gran2 <- nrow(lambda)
    }


    trainY <- pre_proc$Y
    trainZ <- pre_proc$Z
    YMean <- pre_proc$YMean
    ZMean <- pre_proc$ZMean
    C <- pre_proc$C

    if (group == "Basic") {
                    
         beta <- .lassoVARFistX(beta, trainZ, trainY, lambda, tol, p, MN, k,
            k1, s + s1, m, C, YMean, ZMean, separate_lambdas)

            
        }


    if(group=="MCP"|group=="SCAD")
    {
        
        beta <- .MCPFit(beta,trainZ,trainY,lambda,tol,p,MN,k,k1,s,m,C,group,gamma, YMean, ZMean)

    }

    if (group == "BasicEN") {
        
        beta <- .lassoVARFistXEN(beta, trainZ, trainY, lambda, alpha, tol, p, MN, k,
            k1, s + s1, m, C, YMean, ZMean, separate_lambdas)

    }


    if (group == "Lag") {

        GG <- .GroupLassoVAR1(beta, groups, compgroups, trainY, trainZ, lambda, activeset,
            tol, p, MN, k, k1, s + s1, C, YMean, ZMean)

        beta <- GG$beta

        activeset <- GG$active

    }

    if (group == "SparseLag") {

        if (VARX) {

            if (!dual) {
                GG <- .SparseGroupLassoVARX(beta, groups, compgroups, trainY, trainZ,
                  lambda, alpha, INIactive = activeset, tol, starting_eigvals, p, MN,
                  k, s + s1, k1, C, YMean, ZMean)

            } else {

                GG <- .SparseGroupLassoVARXDual(beta, groups, compgroups, trainY,
                  trainZ, lambda, alpha, INIactive = activeset, tol, starting_eigvals,
                  p, MN, k, s + s1, k1, C, YMean, ZMean)

            }
        } else {

            if (!dual) {
                GG <- .SparseGroupLassoVAR(beta, groups, compgroups, trainY, trainZ,
                  lambda, alpha, INIactive = activeset, tol, starting_eigvals, p, MN,
                  C, YMean, ZMean)

            } else {
                GG <- .SparseGroupLassoVARDual(beta, groups, compgroups, trainY,
                  trainZ, lambda, alpha, INIactive = activeset, tol, starting_eigvals,
                  p, MN, C, YMean, ZMean)


            }
        }

        beta <- GG$beta

        activeset <- GG$active

        starting_eigvals <- GG$q1

    }

    if (group == "OwnOther") {

        if (VARX) {

            GG <- .GroupLassoOOX(beta, groups, compgroups, trainY, trainZ, lambda,
                activeset, tol, p, MN, k, k1, s + s1, C, YMean, ZMean)

        } else {


            GG <- .GroupLassoOO(beta, groups, compgroups, trainY, trainZ, lambda, activeset,
                tol, p, MN, C, YMean, ZMean)
        }

        beta <- GG$beta

        activeset <- GG$active

    }

    if (group == "SparseOO") {
        if (VARX) {

            GG <- .SparseGroupLassoVAROOX(beta, groups, compgroups, trainY, trainZ,
                lambda, alpha, INIactive = activeset, tol, p, MN, k1, s + s1, k, dual,
                C, YMean, ZMean)

        } else {

            GG <- .SparseGroupLassoVAROO(beta, groups, compgroups, trainY, trainZ,
                lambda, alpha, INIactive = activeset, tol, starting_eigvals, p, MN,
                dual, C, YMean, ZMean)

            starting_eigvals <- GG$q1

        }

        beta <- GG$beta

        activeset <- GG$active

    }

    if (group == "Tapered") {

        beta <- .lassoVARTL(beta, trainZ, trainY, lambda, tol, p, MN, palpha, C, YMean,
            ZMean)
    }

    if (group == "EFX") {

        beta <- .EFVARX(beta, trainY, trainZ, lambda, tol, MN, k1, s, m, p, C, YMean,
            ZMean)

    }

    if (group == "HLAGC") {
        beta <- .HLAGCAlg(beta, trainY, trainZ, lambda, tol, p, MN, C, YMean, ZMean,
            separate_lambdas)

    }

    if (group == "HLAGOO") {
        beta <- .HLAGOOAlg(beta, trainY, trainZ, lambda, tol, p, MN, C, YMean, ZMean,
            separate_lambdas)
    }

    if (group == "HLAGELEM") {

        beta <- .HLAGElemAlg(beta, trainY, trainZ, lambda, tol, p, MN, C, YMean, ZMean,
            separate_lambdas)

    }

    if (group == "BGR") {
        trainZ <- rbind(1, trainZ)
        beta <- BGRGridSearch(trainY, trainZ, p, lambda, as.numeric(MN))
    }

    ## if (group == "MCP" | group == "SCAD") {

    ##     beta <- .MCPFit(beta, trainZ, trainY, lambda, tol, p, MN, k, k1, s, m, C, YMean,
    ##         ZMean, group, lambda)

    ## }



    if (!exists("activeset")) {
        activeset <- NULL
    }
    if (!exists("starting_eigvals")) {
        starting_eigvals <- NULL
    }
    return(list(beta = beta, activeset = activeset, starting_eigvals = starting_eigvals))

}



#' Simple function to fit BigVAR model with fixed penalty parameter
#' @param Y \eqn{T \times k} multivariate time series or Y \eqn{T \times (k+m)} endogenous and exogenous series, respectively 
#' @param p Predetermined maximal lag order (for modeled series)
#' @param struct The choice of penalty structure (see details).
#' @param lambda vector or matrix of penalty parameters.
#' @param intercept True or False: option to fit an intercept
#' @param RVAR True or False: option to refit based upon the support selected using the Relaxed-VAR procedure
#' @param refit_fraction fraction of least squares refit to incorporate (default 1)
#' @param MN Minnesota Prior Indicator
#' @param VARX List containing VARX model specifications. 
#' @param alpha grid of candidate parameters for the alpha in the Sparse Lag and Sparse Own/Other VARX-L 
#' @param C vector of coefficients to shrink toward a random walk (if \code{MN} is \code{TRUE})
#' @param tf transfer function indicator (i.e. VARX in which p=0 & s>0) (default false)
#' @param tol optimization tolerance (default 1e-4)
#' @param separate_lambdas indicator for separate penalty parameters for each time series (default \code{FALSE})
#' @param beta optional \eqn{k\times (k\times p + m\times s +1)} coefficient matrix to use as a 'warm start' (default \code{NULL})
#' @param gamma additional parameter for SCAD/MCP penalty (default 3)
#' 
#'  @details The choices for 'struct' are as follows
#' \itemize{
#' \item{  'Basic' (Basic VARX-L)}
#' \item{  'BasicEN' (Basic Elastic Net VARX-L)}
#' \item{  'Lag' (Lag Group VARX-L)} 
#' \item{  'SparseLag' (Lag Sparse Group VARX-L)} 
#' \item{  'OwnOther' (Own/Other Group VARX-L) }
#' \item{  'SparseOO' (Own/Other Sparse Group VARX-L) }
#' \item{  'EFX' (Endogenous First VARX-L)}
#' \item{  'HLAGC' (Componentwise HLAG) }
#' \item{  'HLAGOO' (Own/Other HLAG) }
#' \item{  'HLAGELEM' (Elementwise HLAG)}
#' \item{  'Tapered' (Lag weighted Lasso VAR)}
#' \item{  'BGR' (Bayesian Ridge Regression (cf. Banbura et al))}
#' \item{  'MCP' (Minimax Concave Penalty (cf. Breheny and Huang))}
#' \item{  'SCAD' (Smoothly Clipped Absolute Deviation (cf. Breheny and Huang))}
#' }
#'
#' VARX specifications consist of a list with entry k denoting the series that are to be modeled and entry s to denote the maximal lag order for exogenous series.
#'
#' The argument alpha is ignored unless the structure choice is 'SparseLag' or 'Lag.'  By default 'alpha' is set to \code{NULL} and will be initialized as 1/(k+1) in \code{cv.BigVAR} and \code{BigVAR.est}.  Any user supplied values must be between 0 and 1.  

#' @note The specifications 'Basic', 'Lag,' 'SparseLag,' 'SparseOO,' and 'OwnOther' can accommodate both VAR and VARX models.  EFX only applies to VARX models.  'HLAGC,' 'HLAGOO,' 'HLAGELEM,' and 'Tapered' can only be used with VAR models.  Our implementation of the SCAD and MCP penalties is heavily influenced by the implementation in \code{ncvreg}.
#'
#' @seealso \code{\link{cv.BigVAR}},\code{\link{BigVAR.est}},\code{\link{constructModel}}
#' 
#' @references  
#' Banbura, Marta, Domenico Giannone, and Lucrezia Reichlin. 'Large Bayesian vector auto regressions.' Journal of Applied Econometrics 25.1 (2010): 71-92.
#' Breheny P, Huang J (2011). “Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection.” Annals of Applied Statistics, 5(1), 232–253.
#' William B Nicholson, Jacob Bien, and David S Matteson. 'High Dimensional Forecasting via Interpretable Vector Autoregression.' arXiv preprint arXiv:1412.5250, 2016.
#' William B. Nicholson, David S. Matteson, Jacob Bien,VARX-L: Structured regularization for large vector autoregressions with exogenous variables, International Journal of Forecasting, Volume 33, Issue 3, 2017, Pages 627-651,
#' William B Nicholson, David S. Matteson, and Jacob Bien (2016), 'BigVAR: Tools for Modeling Sparse High-Dimensional Multivariate Time Series' arxiv:1702.07094
#' @examples
#' # VARX Example
#' # Fit a Basic VARX-L with k=2, m=1, s=2, p=4, lambda=1e-2
#' VARX=list()
#' VARX$k=2 # indicates that the first two series are modeled
#' VARX$s=2 # sets 2 as the maximal lag order for exogenous series
#' data(Y)
#' BigVAR.fit(Y,p=4,'Basic',lambda=1e-2,VARX=VARX)
#' @export
BigVAR.fit <- function(Y, p, struct, lambda, alpha = NULL, VARX = list(), separate_lambdas = F,
    MN = F, C = as.double(NULL), intercept = TRUE, tf = F, tol = 1e-04, RVAR = F, refit_fraction=1,
    beta = NULL, gamma=3) {
    if (!is.matrix(Y)) {
        stop("Y needs to be a matrix")
    }
    if (is.null(lambda)) {
        stop("Must include penalty parameter")
    }
    if (is.null(alpha)) {
        alpha <- 1/(ncol(Y) + 1)
    }
    dual <- FALSE
    k <- ncol(Y)
    # runs constructModel just to check for errors
    temp.bv <- constructModel(Y, p, struct = struct, gran = c(lambda), ownlambdas = TRUE,
        VARX = VARX, cv = "LOO", model.controls=list(MN=MN, RVAR=RVAR, C=C, intercept = intercept,gamma=gamma))
    gran2 <- length(lambda)
    ## structures=c('Basic','Lag','SparseLag','OwnOther','SparseOO','HLAGC','HLAGOO','HLAGELEM','Tapered','EFX','BGR')

    group <- struct

    groups <- 0
    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }
    s <- ifelse(VARX, VARX$s, 0)
    T <- nrow(Y) - max(p, s)
    if (length(VARX) != 0) {


        k1 <- VARX$k
        s <- VARX$s

        if (!is.null(VARX$contemp)) {

            contemp <- TRUE
            s1 <- 1

        } else {

            contemp <- FALSE
            s1 <- 0
        }
        VARX <- TRUE
        m <- k - k1
        Y1 <- matrix(Y[, 1:k1], ncol = k1)
        X <- matrix(Y[, (ncol(Y) - m + 1):ncol(Y)], ncol = m)

        if (!tf) {

            trainZ <- VARXCons(Y1, X, k1, p, m, s, contemp = contemp)

        } else {

            trainZ <- VARXCons(matrix(0, ncol = 1, nrow = nrow(X)), matrix(X, ncol = m),
                k = 0, p = 0, m = m, s = s, contemp = contemp, oos = FALSE)

        }
        trainZ <- trainZ[2:nrow(trainZ), ]

        trainY <- matrix(Y[(max(c(p, s)) + 1):nrow(Y), 1:k1], ncol = k1)

        grps <- create_group_indexes(group, p, k, gran2 * length(alpha), VARX, k1,
            s + s1)
        groups <- grps$groups
        compgroups <- grps$compgroups
        activeset <- grps$activeset
        starting_eigvals <- grps$starting_eigvals

        if (group == "BGR") {
            Grid <- seq(1, 5, by = 0.025)
            grid <- Grid * sqrt(k * p)
            MSFE <- matrix(0, nrow = 1, ncol = length(grid))
        }


        # Initial Coefficient Matrix
        beta1 <- array(0, dim = c(k1, k1 * p + (k - k1) * (s + s1) + 1, gran2 * length(alpha)))
        if (!is.null(beta) & all(dim(beta) == dim(beta1))) {
            beta <- beta
        } else {
            beta <- beta1
        }

    } else {
        # No VARX
        VARX <- FALSE
        s <- p
        s1 <- 0
        m <- 0
        
        Z1 <- VARXCons(Y, matrix(0, nrow = nrow(Y)), k, p, 0, 0)

        trainZ <- Z1[2:nrow(Z1), , drop = F]

        trainY <- matrix(Y[(p + 1):nrow(Y), ], ncol = k)

        k1 <- k
        s <- 0
        grps <- create_group_indexes(group, p, k, gran2)
        groups <- grps$group
        compgroups <- grps$compgroups
        activeset <- grps$activeset
        starting_eigvals <- grps$starting_eigvals
        beta <- array(0, dim = c(k, k * p + 1, gran2 * length(alpha)))

    }

    if (group == "Tapered")

    {
        palpha <- seq(0, 1, length = 10)
        palpha <- rev(palpha)
        gran2 <- length(lambda) * length(palpha)
        beta <- array(0, dim = c(k, k * p + 1, gran2))
    }


    if (group == "BGR") {
        trainZ <- rbind(1, trainZ)
        beta <- BGRGridSearch(trainY, trainZ, p, lambda, as.numeric(MN))
    } else {
        temp <- .BigVAR.fit(group, beta, trainZ, trainY, lambda, tol, p, m, k1, k,
            s, s1, MN, C, intercept, separate_lambdas, dual, activeset, starting_eigvals,
            groups, compgroups, VARX, alpha, palpha)
        beta <- temp$beta
        
    }
    
    # refit if varx
    if (RVAR) {
        for(i in dim(beta)[3]) {
            beta_rls<- RelaxedLS(cbind(t(trainZ), trainY), beta[,,i])
            beta[,,i] <- (1-refit_fraction)*beta[,,i] +refit_fraction*beta_rls
        }
    }

    return(beta)
}

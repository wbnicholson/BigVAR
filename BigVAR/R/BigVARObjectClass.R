                                        # Ensures that the created BigVAR object is valid
check.BigVAR <- function(object) {

    errors <- character()

    VARX <- object@VARX
    Y <- object@Data

    if (any(is.na(Y))) {
        msg <- c("Remove NA values before running ConstructModel")

        errors <- c(errors, msg)
    }
    if(!is.matrix(Y)){
        msg('Y must be coercible to a matrix')
        errors <- c(errors,msg)
    }
    if (dim(Y)[2] > dim(Y)[1] & length(VARX) == 0) {
        msg <- paste("k is", ncol(Y), "which is greater than T, is Y formatted correctly (k x T)?")
    }
    if (object@lagmax < 0) {
        msg <- c("Maximal lag order must be at least 0")
        errors <- c(errors, msg)
    }
    if (object@lagmax == 0 & !object@Structure %in% c("Basic", "BasicEN")) {
        msg <- c("Only Basic VARX-L supports a transfer function")
        errors <- c(errors, msg)
    }
    structures <- c("Basic", "Lag", "SparseLag", "OwnOther", "SparseOO", "HLAGC", "HLAGOO", "HLAGELEM", "Tapered", "EFX",
                    "BGR", "BasicEN", "MCP", "SCAD")
    cond1 <- object@Structure %in% structures
    if (cond1 == FALSE) {
        msg <- paste("struct must be one of", structures)
        errors <- c(errors, msg)
    }
    if (object@horizon < 1) {
        msg <- paste("Forecast Horizon is ", object@horizon, " must be at least 1")

    }


    if (object@delta < 0) {
        msg <- paste("Huber Delta is ", object@delta, " must be positive")
    }


    if (object@gamma < 0) {
        msg <- paste("Gamma is ", object@gamma, " must be positive")
    }


    if (object@refit_fraction < 0 | object@refit_fraction > 1) {
        msg <- paste("Refit fraction is ", object@refit_fraction, " must be between 0 and 1")
    }

    if (!(object@crossval %in% c("Rolling", "LOO", "None"))) {
        msg <- c("Cross-Validation type must be one of Rolling, LOO or None")
        errors <- c(errors, msg)
    }
    if (length(object@Granularity) != 2 & object@ownlambdas == FALSE) {
        msg("Granularity must have two parameters")
        errors <- c(errors, msg)

    }
    if (any(object@Granularity <= 0)) {
        msg <- c("Granularity parameters must be positive")
        errors <- c(errors, msg)
    }
    structure2 <- c("Basic", "Lag", "HLAGC", "BasicEN", "MCP", "SCAD")
    cond2 <- object@Structure %in% structure2
    k1 <- 0
    if (length(VARX) != 0) {

        k1 <- VARX$k
        if (k1 > ncol(Y)) {
            msg <- c("k is greater than the number of columns in Y")

            errors <- c(errors, msg)
        }
    } else {
        k <- 0
    }
    m <- ncol(Y) - k1

    if (object@tf & object@lagmax != 0) {
        msg <- c("p must be 0 if fitting a transfer function")
        errors <- c(errors, msg)
    }
    nseries <- ncol(Y) - ifelse(m < ncol(Y), m, 0)
    if (nseries == 1 & cond2 == FALSE) {
        msg <- c("Univariate support is only available for Basic VARX-L, Lag Group VARX-L, SCAD, MCP and Componentwise HLAG")

        errors <- c(errors, msg)

    }
    if (length(VARX) == 0 & object@Structure == "EFX") {

        msg <- c("EFX is only supported in the VARX framework")

        errors <- c(errors, msg)

    }
    if (!object@ownlambdas & object@Granularity[2] == 1 & object@separate_lambdas) {
        stop("separate lambda estimation requires more than one candidate penalty parameter")
    }

    if (is.list(VARX) & length(VARX) > 0 & !(exists("k", where = VARX) & exists("s", where = VARX))) {

        msg <- c("VARX Specifications entered incorrectly")

        errors <- c(errors, msg)
    }


    if (object@Structure == "EFX" & !is.null(VARX$contemp)) {
        if (VARX$contemp) {
            msg <- c("EFX does not support contemporaneous dependence")
            errors <- c(errors, msg)
        }

    }
    structs <- c("HLAGC", "HLAGOO", "HLAGELEM")
    if (length(VARX) != 0 & object@Structure %in% structs) {
        msg <- c("EFX is the only nested model supported in the VARX framework")

        errors <- c(errors, msg)

    }
    if (object@T1 > nrow(Y) | object@T2 > nrow(Y) | object@T2 < object@T1) {
        msg <- c("Training dates exceed series length")
        errors <- c(errors, msg)

    }

    if (any(object@alpha < 0) || any(object@alpha > 1)) {
        msg <- c("alpha must be between zero and 1")
        errors <- c(errors, msg)
    }
    if (object@recursive & length(VARX) > 0) {
        msg <- c("recursive forecasts can only be used with VAR models")
        errors <- c(errors, msg)
    }

    if (length(errors) == 0)
        TRUE else errors

}

#' BigVAR Object Class
#'
#' An object class to be used with cv.BigVAR
#' 
#' @slot Data a \eqn{T \times k} multivariate time series
#' @slot model_data processed time series and lag matrix
#' @slot lagmax Maximal lag order for modeled series
#' @slot intercept Indicator as to whether an intercept should be included 
#' @slot Structure Penalty Structure
#' @slot Relaxed Indicator for relaxed VAR
#' @slot Granularity Granularity of penalty grid
#' @slot horizon Desired Forecast Horizon
#' @slot crossval Cross-Validation Procedure
#' @slot Minnesota Minnesota Prior Indicator
#' @slot verbose Indicator for Verbose output
#' @slot dates dates extracted from an xts object 
#' @slot ic Indicator for including AIC and BIC benchmarks
#' @slot VARX VARX Model Specifications
#' @slot VARXI VARX Indicator 
#' @slot T1 Index of time series in which to start cross validation
#' @slot T2  Index of times series in which to start forecast evaluation
#' @slot ONESE Indicator for 'One Standard Error Heuristic'
#' @slot ownlambdas Indicator for user-supplied lambdas
#' @slot tf Indicator for transfer function
#' @slot alpha Grid of candidate alpha values (applies only to Sparse VARX-L and Elastic Net models)
#' @slot recursive Indicator as to whether recursive multi-step forecasts are used (applies only to multiple horizon VAR models)
#' @slot constvec vector indicating variables to shrink toward a random walk instead of toward zero (valid only if Minnesota is \code{TRUE})
#' @slot tol optimization tolerance
#' @slot window.size size of rolling window.  If set to NULL an expanding window will be used.
#' @slot separate_lambdas indicator to use separate penalty parameter for each time series (default \code{FALSE})
#' @slot loss Loss function to select penalty parameter (one of 'L1','L2','Huber').
#' @slot delta delta parameter for Huber loss (default 2.5)
#' @slot gamma gamma parameter for SCAD or MCP penalty (default 3)
#' @slot rolling_oos True or False: indicator to update the penalty parameter over the evaluation period (default \code{False})
#' @slot linear indicator for linearly decrementing penalty grid (FALSE is log-linear).
#' @slot refit_fraction fraction of least squares refit to incorporate (default is 1).
#' @details To construct an object of class BigVAR, use the function \code{\link{constructModel}}
#' @seealso \code{\link{constructModel}}
#' @export
setClass(Class = "BigVAR", representation(Data = "matrix", model_data = "list", lagmax = "numeric", Structure = "character",
                                          Relaxed = "logical", Granularity = "numeric", intercept = "logical", Minnesota = "logical", horizon = "numeric", verbose = "logical",
                                          crossval = "character", ic = "logical", VARX = "list", T1 = "numeric", T2 = "numeric", ONESE = "logical", ownlambdas = "logical",
                                          tf = "logical", alpha = "numeric", recursive = "logical", dates = "character", constvec = "numeric", tol = "numeric",
                                          window.size = "numeric", separate_lambdas = "logical", loss = "character", delta = "numeric", gamma = "numeric", rolling_oos = "logical",
                                          VARXI = "logical", linear = "logical", refit_fraction = "numeric"), validity = check.BigVAR)


#' Construct an object of class BigVAR
#' @param Y \eqn{T \times k} multivariate time series or Y \eqn{T \times (k+m)} endogenous and exogenous series, respectively. 
#' @param p Predetermined maximal lag order (for modeled series).
#' @param struct The choice of penalty structure (see details).
#' @param gran vector of penalty parameter specifications.
#' @param h Desired forecast horizon.
#' @param cv Cross-validation approach, either 'Rolling' for rolling cross-validation or 'LOO' for leave-one-out cross-validation. 'None' for use with BigVAR.fit.
#' @param verbose Verbose output while estimating.
#' @param IC True or False: whether to include AIC and BIC benchmarks.
#' @param VARX List containing VARX model specifications. 
#' @param T1 Index of time series in which to start cross validation.
#' @param T2  Index of times series in which to start forecast evaluation.
#' @param ONESE True or False: whether to use the 'One Standard Error Heuristic.'
#' @param ownlambdas True or False: Indicator for user-supplied penalty parameters.
#' @param recursive True or False: Indicator as to whether iterative multi-step predictions are desired in the VAR context if the forecast horizon is greater than 1.
#' @param loss Loss function to select penalty parameter (one of 'L1','L2','Huber')
#' @param dates optional vector of dates corresponding to \eqn{Y}.
#' @param separate_lambdas indicator for separate penalty parameters for each time series (default \code{FALSE}).
#' @param window.size size of rolling window.  If set to 0 an expanding window will be used. 
#' @param linear indicator for linearly decrementing penalty grid (FALSE is log-linear; default \code{FALSE}).
#' @param rolling_oos True or False: indicator to update the penalty parameter over the evaluation period (default \code{False})
#' @param model.controls named list of control parameters for BigVAR model estimation (see details).
#' @details The choices for 'struct' are as follows
#' \itemize{
#' \item{  'Basic' (Basic VARX-L)}
#' \item{  'BasicEN' (Elastic Net VARX-L)}
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
#' \item{  'MCP' (Minimax Concave Penalty  (cf. Breheny and Huang))}
#' \item{  'SCAD' (Smoothly Clipped Absolute Deviation Penalty (cf. Breheny and Huang))}
#' }
#'
#' The first number in the vector 'gran' specifies how deep to construct the penalty grid and the second specifies how many penalty parameters to use  If ownlambas is set to TRUE, gran should contain the user-supplied penalty parameters.
#' 
#' VARX specifications consist of a named list with entry k denoting the series that are to be modeled and entry s to denote the maximal lag order for exogenous series.
#'
#' As the capabilities of BigVAR have expanded, we have decided to consolidate parameters in the list model.controls.  These parameters include: 
#' \itemize{
#' \item{ 'alpha:' grid of candidate parameters for the alpha in the Basic Elastic Net, Sparse Lag, Sparse Own/Other VARX-L.}
#' \item{ 'C:' vector of coefficients to shrink toward a random walk (if \code{MN} is \code{TRUE}).}
#' \item{ 'delta:' parameter for Huber loss (default 2.5)} 
#' \item{ 'intercept:' option to fit an intercept, default \code{TRUE}}
#' \item{ 'loss:' Loss function to select penalty parameter (one of 'L1','L2','Huber')}
#' \item{ 'MN:' Minnesota Prior Indicator, default \code{FALSE}}
#' \item{ 'RVAR:'  option to refit based upon the support selected using the Relaxed-VAR procedure (default FALSE).}
#' \item{ 'refit_fraction:'  If RVAR is \code{TRUE}, proportional tradeoff between least squares fit and penalized fit (default 1).}
#' \item{ 'tol:' optimization tolerance (default 1e-4)}
#' }
#' 
#' The argument alpha is ignored unless the structure choice is 'SparseLag' or 'Lag.'  By default 'alpha' is set to \code{NULL} and will be initialized as 1/(k+1) in \code{cv.BigVAR} and \code{BigVAR.est}.  Any user supplied values must be between 0 and 1.  

#' @note The specifications 'Basic','BasicEN', 'Lag,' 'SparseLag,' 'SparseOO','OwnOther', 'MCP', and 'SCAD.' can accommodate both VAR and VARX models.  EFX only applies to VARX models.  'HLAGC,' 'HLAGOO,' 'HLAGELEM,' and 'Tapered' can only be used with VAR models.  Our implementation of the SCAD and MCP penalties is heavily influenced by the package \code{ncvreg}. 
#'
#' @seealso \code{\link{cv.BigVAR}},\code{\link{BigVAR.est}}
#' 
#' @references
#' Banbura, Marta, Domenico Giannone, and Lucrezia Reichlin. 'Large Bayesian vector auto regressions.' Journal of Applied Econometrics 25.1 (2010): 71-92.
#' Breheny P, Huang J (2011). “Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection.” Annals of Applied Statistics, 5(1), 232–253.
#' Nicholson, William, I. Wilms, J. Bien, and D. S. Matteson. High dimensional forecasting via interpretable vector autoregression. Journal of Machine Learning Research, 21(166):1–52, 2020.
#' William B. Nicholson, David S. Matteson, Jacob Bien,VARX-L: Structured regularization for large vector autoregressions with exogenous variables, International Journal of Forecasting, Volume 33, Issue 3, 2017, Pages 627-651,
#' William B Nicholson, David S. Matteson, and Jacob Bien (2016), 'BigVAR: Tools for Modeling Sparse High-Dimensional Multivariate Time Series' arxiv:1702.07094
#' @examples
#' # VARX Example
#' # Create a Basic VARX-L with k=2, m=1, s=2, p=4
#' VARX=list()
#' VARX$k=2 # indicates that the first two series are modeled
#' VARX$s=2 # sets 2 as the maximal lag order for exogenous series
#' data(Y)
#' T1=floor(nrow(Y)/3)
#' T2=floor(2*nrow(Y)/3)
#' Model1=constructModel(Y,p=4,struct='Basic',gran=c(50,10),verbose=FALSE,VARX=VARX,T1=T1,T2=T2)
#' @export
constructModel <- function(Y, p, struct, gran, h = 1, cv = "Rolling", verbose = TRUE, IC = TRUE, VARX = list(), T1 = floor(nrow(Y)/3),
                           T2 = floor(2 * nrow(Y)/3), ONESE = FALSE, ownlambdas = FALSE, recursive = FALSE, dates = as.character(NULL), window.size = 0,
                           separate_lambdas = FALSE, linear = FALSE, loss = "L2", rolling_oos = FALSE, model.controls = list()) {


    if (exists("RVAR", where = model.controls)) {
        RVAR <- model.controls$RVAR
    } else {
        RVAR <- FALSE
    }

    if (exists("refit_fraction", where = model.controls)) {
        refit_fraction <- model.controls$refit_fraction
    } else {
        refit_fraction <- 1
    }


    ## if (exists("linear")) {
    ##     linear <- linear
    ## } else {
    ##     linear <- FALSE
    ## }
    linear <- linear
    
    if (exists("alpha", where = model.controls)) {
        alpha = model.controls$alpha
    } else {
        alpha <- as.double(NULL)
    }

    if (exists("tol", where = model.controls)) {
        tol <- model.controls$tol
    } else {
        tol <- 1e-04
    }

    if (exists("MN", where = model.controls)) {
        MN <- model.controls$MN
    } else {
        MN <- FALSE
    }
    if (exists("C", where = model.controls)) {
        C <- model.controls$C
    } else {
        C <- as.double(NULL)
    }
    if (exists("delta", where = model.controls)) {
        delta <- model.controls$delta
    } else {
        delta <- 2.5
    }
    if (exists("gamma", where = model.controls)) {
        gamma <- model.controls$gamma
    } else {
        gamma <- 3
    }
    if (exists("intercept", where = model.controls)) {
        intercept <- model.controls$intercept
    } else {
        intercept <- TRUE
    }
                                        # remove alpha if it's not used
    if (length(alpha) != 0 & !struct %in% c("BasicEN", "SparseLag", "SparseOO")) {
        alpha <- as.double(NULL)
    }
    if (any(is.na(Y))) {
        stop("Remove NA values before running constructModel")
    }
    if (dim(Y)[2] > dim(Y)[1] & length(VARX) == 0) {
        warning("k is greater than T, is Y formatted correctly (k x T)?")
    }
    if (p < 0) {
        stop("Maximal lag order must be at least 0")
    }
    if (p == 0 & !struct %in% c("Basic", "BasicEN")) {
        stop("Only Basic VARX-L supports a transfer function")
    }
    oldnames <- c("None", "Diag", "SparseDiag")
    if (struct %in% oldnames)
        stop("Naming Convention for these structures has changed. Use Basic, OwnOther, and SparseOO.")
    structures <- c("Basic", "Lag", "SparseLag", "OwnOther", "SparseOO", "HLAGC", "HLAGOO", "HLAGELEM", "Tapered", "EFX",
                    "BGR", "BasicEN", "MCP", "SCAD")
    if (struct == "BasicEN" & length(alpha) > 1 & separate_lambdas) {
        stop("Cannot perform separate lambdas per series and range of alphas simultaneously")
    }
    cond1 <- struct %in% structures


    if (!loss %in% c("L1", "L2", "Huber")) {
        stop("loss must be one of L1,L2,Huber")
    }
    if (!cond1) {
        stop(paste("struct must be one of", structures))
    }
    if (h < 1) {
        stop("Forecast Horizon must be at least 1")
    }
    if (cv != "None" & T1 - 2 * h - p < 0) {
        stop("Forecast Horizon too long; increase T1 or decrease h ")
    }
    if (!(cv %in% c("Rolling", "LOO", "None"))) {
        stop("Cross-Validation type must be one of Rolling, LOO or None")
    }
    if (length(gran) != 2 & ownlambdas == FALSE) {
        stop("Granularity must have two parameters")
    }
    if (any(gran <= 0)) {
        stop("Granularity parameters must be positive")
    }
    if (tol < 0 | tol > 0.1) {
        stop("Tolerance must be positive")
    }
    if (window.size > nrow(Y) | window.size < 0) {
        stop("window size must be shorter than the series length")
    }
    ## ws_init<- max(c(T1 - window.size - h-p, 1))
    ## start_index<- (ws1 + h):(T1 - h)   
    ## if (length(start_index)<10) {
    ##     stop("window.size and forecast horizon")
    ## }
    
    if (delta < 0) {
        stop("huber delta must be positive")
    }
    bss <- c("Basic", "HLAGC", "HLAGOO", "HLAGELEM", "BasicEN", "SCAD", "MCP")
    if (separate_lambdas & !struct %in% c("Basic", "HLAGC", "HLAGOO", "HLAGELEM", "BasicEN", "SCAD", "MCP")) {
        stop(paste("separate lambda estimation only available for ", toString(bss)))
    }
    start_ind <- (T1 - p - h + 1)
    if (cv == "Rolling" & start_ind < 5) {
        stop("too few rows for rolling validation, try running BigVAR.fit")
    }
    if (MN & intercept) {
        intercept <- FALSE
    }
    structure2 <- c("Basic", "Lag", "HLAGC", "BasicEN", "MCP", "SCAD")
    cond2 <- struct %in% structure2
    if (length(VARX) != 0) {
        k <- VARX$k
        if (k > ncol(Y)) {
            stop("k is greater than the number of columns in Y")
        }
    } else {
        k <- ncol(Y)
    }
    m <- ncol(Y) - k
    nseries <- ncol(Y) - ifelse(m < ncol(Y), m, 0)
    if (p == 0) {
        tf <- TRUE
    } else {
        tf <- FALSE
    }
    if (nseries == 1 & cond2 == FALSE) {
        stop("Univariate support is only available for Basic, Elastic Net, Lag Group, and Componentwise HLAG")
    }
    if (length(VARX) == 0 & struct == "EFX") {
        stop("EFX is only supported in the VARX framework")
    }
    if (struct == "EFX" & !is.null(VARX$contemp)) {
        if (VARX$contemp) {
            stop("EFX does not support contemporaneous dependence")
        }
    }
    structs <- c("HLAGC", "HLAGOO", "HLAGELEM")
    if (length(VARX) != 0 & struct %in% structs) {
        stop("EFX is the only nested model supported in the VARX framework")
    }
    if (length(VARX) != 0 & struct == "BGR") {
        stop("BGR is only available in the VAR framework")
    }
    if (length(VARX) != 0 & struct == "Tapered") {
        stop("Lag-Weighted Lasso only available in VAR framework")
    }
    if (T1 > nrow(Y) | T2 > nrow(Y) | T2 < T1) {
        stop("Training dates exceed series length")
    }
    if (is.list(VARX) & length(VARX) > 0 & !(exists("k", where = VARX) & exists("s", where = VARX))) {
        stop("VARX Specifications entered incorrectly")
    }
    if (!is.null(alpha)) {
        if (any(alpha < 0) || any(alpha > 1)) {
            stop("alpha must be between 0 and 1")
        }
    }
    if (length(C) != 0) {
        if (length(C) != k) {
            stop("C must have length k")
        }
        if (!all(C %in% c(0, 1))) {
            stop("Values of C must be either 0 or 1")
        }
    } else {
        C <- rep(1, k)
    }
    if (length(dates) != 0) {
        ind <- dates
    } else {
        ind <- as.character(NULL)
    }
    model_data <- VARXConsModel(Y, p, VARX, tf)
    var_data <- list(trainY = model_data$trainY, trainZ = model_data$trainZ)
    VARXI <- model_data$VARXI
    if (VARXI) {
        VARX$s1 <- model_data$s1
        VARX$contemp <- model_data$contemp
    }
    if (VARXI & window.size > 0) {
        window.size <- window.size + 1
    }
    if (ncol(model_data$trainY) == 1 & separate_lambdas) {
        separate_lambdas <- FALSE
    }
                               (BV1 <- new("BigVAR", Data = Y, model_data = var_data, lagmax = p, Structure = struct, Relaxed = RVAR, Granularity = gran,
                                           Minnesota = MN, verbose = verbose, horizon = h, crossval = cv, ic = IC, VARX = VARX, VARXI = VARXI, T1 = T1, T2 = T2,
                                           ONESE = ONESE, ownlambdas = ownlambdas, tf = tf, alpha = alpha, recursive = recursive, dates = ind, constvec = C,
                                           intercept = intercept, tol = tol, window.size = window.size, separate_lambdas = separate_lambdas, loss = loss, delta = delta,
                                           gamma = gamma, rolling_oos = rolling_oos, linear = linear, refit_fraction = refit_fraction))
    return(BV1)
}


                                        # show-default method to show an object when its name is printed in the console.
#' Default show method for an object of class BigVAR
#'
#' @param object \code{BigVAR} object created from \code{ConstructModel}
#' @return Displays the following information about the BigVAR object:
#' \itemize{
#' \item{Prints the first 5 rows of \code{Y}}
#' \item{ Penalty Structure}
#' \item{ Relaxed Least Squares Indicator}
#' \item{Maximum lag order} 
#' \item{ VARX Specifications (if applicable)}
#' \item{Start, end of cross validation period}
#' }
#' @seealso \code{\link{constructModel}} 
#' @name show.BigVAR
#' @aliases show,BigVAR-method
#' @docType methods
#' @rdname show-methods
#' @export
setMethod("show", "BigVAR", function(object) {
    T1P <- ifelse(length(object@dates) != 0, object@dates[object@T1], object@T1)
    T2P <- ifelse(length(object@dates) != 0, object@dates[object@T2], object@T2)
    nrowShow <- min(5, nrow(object@Data))
    cat("*** BIGVAR MODEL *** \n")
    cat("Data (First 5 Observations):\n")
    print(formatC(object@Data[1:nrowShow, ]), quote = FALSE)
    cat("Structure\n")
    print(object@Structure)
    cat("Forecast Horizon \n")
    print(object@horizon)
    cat("Relaxed VAR \n")
    print(object@Relaxed)
    cat("Minnesota Prior \n")
    print(object@Minnesota)
    cat("Maximum Lag Order \n")
    print(object@lagmax)
    if (length(object@VARX) != 0) {
        cat("VARX Specs \n")
        print(object@VARX)
    }
    cat("Start of Cross Validation Period \n")
    print(T1P)
    cat("End of Cross Validation Period \n")
    print(T2P)
})


#' Plot a BigVAR object
#' 
#' @param x BigVAR object created from \code{ConstructModel}
#' @param y NULL
#' @param ... additional plot arguments
#' @return NA, side effect is graph
#' @details Uses plot.zoo to plot each indivdual series of \code{Y} on a single plot
#' @name plot.BigVAR
#' @import methods
#' @seealso \code{\link{constructModel}}
#' @aliases plot,BigVAR-method
#' @docType methods
#' @method plot method
#' @rdname plot.BigVAR-methods
#' @export
#' @importFrom zoo plot.zoo
#' @importFrom zoo as.zoo
#' @importFrom zoo zoo
#' @importFrom zoo as.yearqtr
#' @importFrom zoo index
#' @importFrom graphics abline
#' @importFrom graphics legend
setMethod(f = "plot", signature = "BigVAR", definition = function(x, y = NULL, ...) {
    T1P <- ifelse(length(x@dates) != 0, x@dates[x@T1], x@T1)
    T2P <- ifelse(length(x@dates) != 0, x@dates[x@T2], x@T2)
    g <- ncol(x@Data)
    names <- ifelse(rep(!is.null(colnames(x@Data)), ncol(x@Data)), colnames(x@Data), as.character(1:g))
    if (length(x@dates) != 0) {
        dates <- as.yearqtr(x@dates)
    } else {
        dates <- seq_len(nrow(x@Data))
    }
    Yzoo <- zoo(as.matrix(x@Data), order.by = dates)
    plot.zoo(Yzoo, plot.type = "single", col = 1:g)
    legend("topright", names, lty = 1, col = 1:g)
    abline(v = index(Yzoo[as.yearqtr(T1P)]))
    abline(v = index(Yzoo[as.yearqtr(T2P)]))
})

#' Cross Validation for BigVAR
#' 
#' @usage cv.BigVAR(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details The main function of the BigVAR package. Performs cross validation to select penalty parameters over a training sample (as the minimizer of in-sample MSFE), then evaluates them over a test set.  Compares against sample mean, random walk, AIC, and BIC benchmarks.  Creates an object of class \code{BigVAR.results}
#' @return An object of class \code{BigVAR.results}.
#' @seealso \code{\link{constructModel}}, \code{\link{BigVAR.results}},\code{\link{BigVAR.est}} 
#' @name cv.BigVAR
#' @aliases cv.BigVAR,BigVAR-method
#' @docType methods
#' @rdname cv.BigVAR-methods
#' @examples
#' data(Y)
#' # Fit a Basic VARX-L with rolling cross validation 
#' Model1=constructModel(Y,p=4,struct='Basic',gran=c(50,10), verbose=FALSE)
#' results=cv.BigVAR(Model1)
#' @importFrom abind adrop
#' @importFrom abind abind
#' @export
setGeneric(name = "cv.BigVAR", def = function(object) {
    standardGeneric("cv.BigVAR")
})
                                        # Cross-validation and evaluation function
setMethod(f = "cv.BigVAR", signature = "BigVAR", definition = function(object) {
    p <- object@lagmax
    Y <- object@Data
    k <- ncol(Y)
    alpha <- object@alpha
    gamma <- object@gamma
    RVAR <- object@Relaxed
    refit_fraction <- object@refit_fraction
    group <- object@Structure
    cvtype <- object@crossval
    if(cvtype=='None'){
        stop("set cv to rolling or LOO to run cv.BigVAR")
    }
    intercept <- object@intercept
    recursive <- object@recursive
    VARX <- object@VARX
    VARXI <- object@VARXI
    tol <- object@tol
    window.size <- object@window.size
    verbose <- object@verbose
    loss <- object@loss
    delta <- object@delta
    linear <- object@linear
    rolling_oos <- object@rolling_oos
    if (length(alpha) == 0) {
        if (length(VARX) > 0) {
            alpha <- 1/(VARX$k + 1)
        } else {
            alpha <- 1/(k + 1)
        }
    }
    C <- object@constvec
    if (length(alpha) > 1 & group %in% c("SparseLag", "SparseOO", "BasicEN")) {
        dual <- TRUE
    } else {
        dual <- FALSE
    }
    MN <- object@Minnesota
    h <- object@horizon
    jj <- 0
    separate_lambdas <- object@separate_lambdas
    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }
    if (object@crossval == "Rolling") {
        T1 <- object@T1
        T1 <- floor(T1)
    } else {
        T1 <- p + h + 1
    }
    T2 <- object@T2
    T2 <- floor(T2)
    ONESE <- object@ONESE
    if (object@ownlambdas) {
        lambda <- object@Granularity
        gran2 <- length(lambda)
        if (gran2 == 1) {
            ONESE <- FALSE
        }
    } else if (object@Granularity[2] == 1) {
        warning("only one penalty parameter; more efficient to run BigVAR.est instead of cv.BigVAR")
    }
    if (!object@ownlambdas) {
        gran2 <- object@Granularity[2]
        gran1 <- object@Granularity[1]
    }
    trainZ <- object@model_data$trainZ
    trainY <- object@model_data$trainY
    if (VARXI) {
        k1 <- object@VARX$k
        s <- object@VARX$s
        contemp <- object@VARX$contemp
        s1 <- object@VARX$contemp
        m <- ncol(Y) - ncol(trainY)
        beta <- array(0, dim = c(k1, k1 * p + (k - k1) * (s + s1) + 1, gran2 * length(alpha)))
    } else {
        k1 <- 0
        s <- 0
        contemp <- FALSE
        s1 <- 0
        m <- 0
        beta <- array(0, dim = c(k, k * p + 1, gran2 * length(alpha)))
    }
    if (object@crossval == "Rolling") {
        T1 <- T1 - max(p, s)
        T2 <- T2 - max(p, s)
    }
    grps <- create_group_indexes(group, p, k, gran2 * length(alpha), VARXI, k1, s + s1)
    groups <- grps$groups
    compgroups <- grps$compgroups
    activeset <- grps$activeset
    starting_eigvals <- grps$starting_eigvals
    if (object@ownlambdas) {
        lambda <- object@Granularity
    } else {
        gran2 <- object@Granularity[2]
        gran1 <- object@Granularity[1]
        
        lambda <- create_lambda_grid(trainY[1:T2, , drop = FALSE], trainZ[, 1:T2], lapply(groups, function(x) {
            x + 1
        }), gran1, gran2, group, p, k1, s + s1, m, k, MN, alpha, C, intercept, tol, VARXI, separate_lambdas, dual, gamma,
        linear, verbose)
    }
    h <- object@horizon
    ZFull <- list()
    if (!is.matrix(trainZ)) {
        trainZ <- matrix(trainZ, ncol = 1)
    }
    if (!is.matrix(trainY)) {
        trainY <- matrix(trainY, ncol = 1)
    }
    ZFull$Z <- trainZ
    ZFull$Y <- trainY
    T3 <- nrow(trainY)
    if (object@ownlambdas) {
        lambda <- object@Granularity
    }
    if (group == "Tapered") {
        palpha <- seq(0, 1, length = 10)
        palpha <- rev(palpha)
        gran2 <- length(lambda) * length(palpha)
        beta <- array(0, dim = c(k, k * p + 1, gran2))
    } else {
        palpha <- NULL
    }
    if (!is.matrix(ZFull$Y)) {
        ZFull$Y <- matrix(ZFull$Y, ncol = 1)
    }
    if (!dual) {
        if (separate_lambdas) {
            if (!VARXI) {
                MSFE <- array(0, dim = c(T2 - T1 + 1, gran2, k))
            } else {
                MSFE <- array(0, dim = c(T2 - T1 + 1, gran2, k1))
            }
        } else {
            MSFE <- matrix(0, nrow = T2 - T1 + 1, ncol = gran2)
            lambda <- as.matrix(lambda)
        }
    } else {
        nalpha <- length(alpha)
        MSFE <- matrix(0, nrow = T2 - T1 + 1, ncol = gran2 * nalpha)
    }
    if (verbose) {
        pb <- txtProgressBar(min = T1, max = T2, style = 3)
        cat("Validation Stage:", group)
    }
    YT <- Y[1:T2, , drop = FALSE]
    betaWS <- beta
    for (v in (T1 - h + 1):T2) {
        if (v + h - 1 > T2) {
            break
        }
        if (cvtype == "Rolling") {
            if (h > 1 & !recursive) {
                if (window.size != 0) {
                    ws1 <- max(c(v - window.size - h, 1))
                    index_y <- (ws1 + h-1):(v - 1)
                    index_z <- (ws1 ):(v - h)
                    trainY <- ZFull$Y[index_y, , drop = FALSE]
                    trainZ <- ZFull$Z[, index_z,drop=FALSE]
                } else {
                    index_y <- (h):(v - 1)
                    index_z <- 1:(v - h)
                    trainY <- ZFull$Y[index_y, , drop = FALSE]
                    trainZ <- ZFull$Z[, index_z, drop = F]
                }
            } else {
                if (window.size != 0) {
                    ws1 <- max(c(v - window.size, 1))
                    trainY <- ZFull$Y[(ws1):(v - 1), , drop = FALSE]
                    trainZ <- ZFull$Z[, (ws1):(v - 1), drop = FALSE]
                } else {
                    trainY <- ZFull$Y[(1):(v - 1), , drop = FALSE]
                    trainZ <- ZFull$Z[, (1):(v - 1), drop = FALSE]
                }
            }
        } else {
            if (VARXI) {
                YT2 <- YT[-v, , drop = FALSE]
                Y1 <- YT2[, 1:k1, drop = FALSE]
                X <- YT2[, (ncol(YT2) - m + 1):ncol(YT2), drop = FALSE]
                trainZ <- VARXCons(Y1, X, k1, p, m, s, contemp = contemp)
                trainZ <- trainZ[2:nrow(trainZ), , drop = FALSE]
                trainY <- YT2[(max(c(p, s)) + 1):nrow(YT2), 1:k1, drop = FALSE]
            } else {
                YT2 <- YT[-v, , drop = FALSE]
                Z1 <- VARXCons(YT2, matrix(0, nrow = nrow(YT2)), k, p, 0, 0)
                trainZ <- Z1[2:nrow(Z1), ]
                trainY <- YT2[(p + 1):nrow(YT2), , drop = FALSE]                
            }
        }
        temp <- .BigVAR.fit(group, betaWS, trainZ, trainY, lambda, tol, p, m, k1, k, s, s1, MN, C, intercept, separate_lambdas,
                            dual, activeset, starting_eigvals, groups, compgroups, VARXI, alpha, palpha, gamma)
        beta <- temp$beta
        betaWS <- temp$beta
        if (MN) {
            for (i in 1:dim(betaWS)[3]) {
                submat <- adrop(betaWS[1:dim(betaWS)[1], 1:dim(betaWS)[1], i, drop = F], drop = 3)
                diag(submat) <- diag(submat) - C
                betaWS[1:dim(betaWS)[1], 1:dim(betaWS)[1], i] <- submat
            }
        }
        activeset <- temp$activeset
        q1a <- temp$q1a
        eZ <- c(1, ZFull$Z[, v])
        msfe_index <- v - (T1 - h)

        temp_results <- refine_and_forecast(beta, as.matrix(eZ), trainZ, trainY, ZFull$Y[v + h - 1, , drop = F], lambda = lambda,
                                            h = h, recursive = recursive, MN = MN, RVAR = RVAR, refit_fraction = refit_fraction, separate_lambdas = separate_lambdas,
                                            inds = NULL, loss = loss, delta = delta, k = k, p = p, k1 = k1, s = s, oos = FALSE)
        if (separate_lambdas) {
            MSFE[msfe_index, , ] <- temp_results$MSFE
        } else {
            MSFE[msfe_index, ] <- temp_results$MSFE
        }
        if (verbose) {
            setTxtProgressBar(pb, v)
        }
    }
    tapered <- group == "Tapered"

    optimal_indices <- find_optimal_lambda(MSFE, lambda, palpha, alpha, ONESE, dual, separate_lambdas, tapered)
    optind <- optimal_indices$optind
    lambdaopt <- optimal_indices$lambdaopt
    alphaopt <- optimal_indices$alphaopt
    int_results <- new("BigVAR.intermediate", ZFull = ZFull, InSampMSFE = MSFE, InSampSD = apply(MSFE, 2, sd), index = optind,
                       OptimalLambda = lambdaopt, dual = dual, contemp = contemp, LambdaGrid = as.matrix(lambda), object, T1 = T2, T2 = T3,
                       alpha = alphaopt)
    OOSEval <- BigVAR.Eval(int_results)
    MSFEOOSAgg <- na.omit(OOSEval$MSFE)
    betaPred <- OOSEval$betaPred
    betaArray <- OOSEval$betaArray
    lambda_evolve <- OOSEval$lambda_evolve
    preds <- OOSEval$predictions
    Y <- object@Data
    if (VARXI) {
        if (contemp) {
            OOS <- FALSE
        } else {
            OOS <- TRUE
        }
        if (!object@tf) {
            Zvals <- VARXCons(matrix(Y[, 1:k1], ncol = k1), matrix(Y[, (ncol(Y) - m + 1):ncol(Y)], ncol = m), k1, p, m, s,
                              oos = OOS, contemp = contemp)
        } else {
            Zvals <- VARXCons(matrix(0, ncol = 1, nrow = nrow(Y)), matrix(Y[, (ncol(Y) - m + 1):ncol(Y)], ncol = m), 0, 0,
                              m, s, oos = FALSE, contemp = contemp)
        }
    } else {
        m <- 0
        s <- 0
        k1 <- k
        Zvals <- VARXCons(matrix(Y[, 1:k1], ncol = k1), matrix(0, nrow = nrow(Y)), k1, p, m, s, oos = TRUE)
    }
    Zvals <- matrix(Zvals[, ncol(Zvals)], ncol = 1)
    if (ncol(Y) == 1 | k1 == 1) {
        betaPred <- matrix(betaPred, nrow = 1)
    }
    lagmatrix <- rbind(rep(1, ncol(ZFull$Z)), ZFull$Z)
    fitted <- t(betaPred %*% lagmatrix)
    resids <- ((ZFull$Y) - fitted)
    MSFEOOS <- mean(na.omit(MSFEOOSAgg))
    seoos <- sd(na.omit(MSFEOOSAgg))/sqrt(length(na.omit(MSFEOOSAgg)))
    if (!VARXI) {
        k1 <- k
    }
    meanbench <- .evalMean(ZFull$Y[, 1:k1], T2, T3, h = h, loss = loss, delta = delta)
    RWbench <- .evalRW(ZFull$Y[, 1:k1], T2, T3, h = h, loss = loss, delta = delta)
    if (object@ic == FALSE | object@tf) {
        AICbench <- list()
        AICbench$Mean <- as.double(NA)
        AICbench$SD <- as.double(NA)
        AICbench$preds <- as.matrix(NA)
        AICbench$pvec <- as.double(NA)
        AICbench$svec <- as.double(NA)
        BICbench <- list()
        BICbench$Mean <- as.double(NA)
        BICbench$SD <- as.double(NA)
        BICbench$preds <- as.matrix(NA)
        BICbench$pvec <- as.double(NA)
        BICbench$svec <- as.double(NA)
    } else {
        if (!VARXI) {
            X <- matrix(0, nrow = nrow(Y), ncol = k)
            AICbench1 <- VARXForecastEval(matrix(ZFull$Y, ncol = k), X, p, 0, T2, T3, "AIC", h, loss = loss, delta = delta)
            AICbench <- list()
            AICbench$Mean <- mean(AICbench1$MSFE)
            AICbench$SD <- sd(AICbench1$MSFE)/sqrt(length(AICbench1$MSFE))
            AICbench$preds <- AICbench1$pred
            AICbench$pvec <- AICbench1$p
            AICbench$svec <- AICbench1$s
            BICbench1 <- VARXForecastEval(matrix(ZFull$Y, ncol = k), X, p, 0, T2, T3, "BIC", h, loss = loss, delta = delta)
            BICbench <- list()
            BICbench$Mean <- mean(BICbench1$MSFE)
            BICbench$SD <- sd(BICbench1$MSFE)/sqrt(length(BICbench1$MSFE))
            BICbench$preds <- BICbench1$pred
            BICbench$pvec <- BICbench1$p
            BICbench$svec <- BICbench1$s
        } else {
            offset <- max(c(p, s))
            X <- matrix(Y[(offset + 1):nrow(Y), (k1 + 1):ncol(Y)], ncol = m)
            AICbench1 <- VARXForecastEval(matrix(ZFull$Y[, 1:k1], ncol = k1), as.matrix(X), p, s, T2, T3, "AIC", h = h, loss = loss,
                                          delta = delta)
            AICbench <- list()
            AICbench$Mean <- mean(AICbench1$MSFE)
            AICbench$SD <- sd(AICbench1$MSFE)/sqrt(length(AICbench1$MSFE))
            AICbench$preds <- AICbench1$pred
            AICbench$pvec <- AICbench1$p
            AICbench$svec <- AICbench1$s
            BICbench1 <- VARXForecastEval(matrix(ZFull$Y[, 1:k1], ncol = k1), X, p, s, T2, T3, "BIC", h = h, loss = loss,
                                          delta = delta)
            BICbench <- list()
            BICbench$Mean <- mean(BICbench1$MSFE)
            BICbench$SD <- sd(BICbench1$MSFE)/sqrt(length(BICbench1$MSFE))
            BICbench$preds <- BICbench1$pred
            BICbench$pvec <- BICbench1$p
            BICbench$svec <- BICbench1$s
        }
    }
    if (!VARXI) {
        contemp <- FALSE
    }
    if (separate_lambdas) {
        tmean <- t(apply(MSFE, 3, colMeans))
    }
    if (separate_lambdas) {
        isMSFE <- MSFE
    } else {
        isMSFE <- array(MSFE, dim = c(nrow(MSFE), ncol(MSFE), 1))
    }
    sparse_count <- function(x) {
        x_ss <- x[, 2:ncol(x)]
        sc <- length(x_ss[x_ss != 0])/length(x_ss)
        sc
    }
    sc <- mean(apply(betaArray, 3, sparse_count))
    results <- new("BigVAR.results", InSampMSFE = isMSFE, InSampSD = apply(MSFE, 2, sd)/sqrt(nrow(MSFE)), LambdaGrid = as.matrix(lambda),
                   index = optind, OptimalLambda = lambdaopt, OOSMSFE = as.matrix(MSFEOOSAgg), seoosmsfe = seoos, MeanMSFE = meanbench$Mean,
                   AICMSFE = AICbench$Mean, AICpvec = AICbench$pvec, AICsvec = AICbench$svec, AICPreds = AICbench$preds, BICpvec = BICbench$pvec,
                   BICsvec = BICbench$svec, BICPreds = BICbench$preds, RWMSFE = RWbench$Mean, RWPreds = RWbench$preds, MeanSD = meanbench$SD,
                   MeanPreds = meanbench$preds, AICSD = AICbench$SD, BICMSFE = BICbench$Mean, BICSD = BICbench$SD, RWSD = RWbench$SD,
                   sparse_count = sc, betaPred = betaPred, Zvals = Zvals, resids = resids, VARXI = VARXI, preds = preds, alpha = alphaopt,
                   fitted = fitted, lagmatrix = lagmatrix, betaArray = betaArray, dual = dual, contemp = contemp, lambda_evolve_path = lambda_evolve,
                   object)
    return(results)
})



#' BigVAR.intermediate
#' This class contains the in-sample results for cv.BigVAR
#'
#' It inherits the class BigVAR, and contains the results from rolling validation 
#' @field ZFull List containing full lag matrix and time series
#' @field InSampMSFE In-sample MSFE from optimal value of lambda
#' @field LambdaGrid Grid of candidate lambda values
#' @field index Index order of optimal lambda value 
#' @field OptimalLambda Value of lambda that minimizes MSFE
#' @field Data a \eqn{T \times k} or \eqn{T\times k + m} multivariate time Series
#' @field lagmax Maximal lag order
#' @field Structure Penalty structure
#' @field Relaxed Indicator for relaxed VAR
#' @field Granularity Granularity of penalty grid
#' @field horizon Desired forecast horizon
#' @field crossval Cross-Validation procedure
#' @field alpha additional penalty parameter for Sparse Lag Group or Sparse Own/Other methods. Will contain either the heuristic choice of \eqn{1/(k+1)} or the value selected by cross validation if the argument \code{dual} is set to \code{TRUE}
#' @field Minnesota Minnesota Prior Indicator
#' @field verbose  verbose indicator
#' @field dual indicator as to whether dual cross validation was conducted
#' @field contemp indicator if contemporaneous exogenous predictors are used
#'
#' @note One can also access any object of class BigVAR from BigVAR.intermediate
#' @name BigVAR.intermediate
#' @rdname BigVAR.intermediate
#' @aliases BigVAR.intermediate-class
#' @exportClass BigVAR.intermediate
#' @author Will Nicholson
#' @export
setClass("BigVAR.intermediate", representation(ZFull = "list", InSampMSFE = "array", InSampSD = "numeric", LambdaGrid = "matrix",
                                               index = "numeric", OptimalLambda = "numeric", dual = "logical", contemp = "logical"), contains = "BigVAR")
setGeneric(name = "BigVAR.Eval", def = function(object) {
    standardGeneric("BigVAR.Eval")
})


setGeneric(name = "BigVAR.Eval", def = function(object) {
    standardGeneric("BigVAR.Eval")
})

setMethod(f = "BigVAR.Eval", signature = "BigVAR.intermediate", definition = function(object) {
    VARXI <- object@VARXI
    VARX <- object@VARX
    p <- object@lagmax

    gamma <- object@gamma
    tol <- object@tol
    T2 <- object@T1
    T3 <- object@T2
    MN <- object@Minnesota
    h <- object@horizon
    loss <- object@loss
    delta <- object@delta
    rolling_oos <- object@rolling_oos
    if (rolling_oos) {
        lambda <- object@LambdaGrid
        optind <- object@index
    } else {
        lambda <- object@OptimalLambda
        optind <- 1
    }
    gran2 <- nrow(as.matrix(lambda))
    lambdaopt <- object@OptimalLambda
    delta <- object@delta
    verbose <- object@verbose
    window.size <- object@window.size
    recursive <- object@recursive
    group <- object@Structure
    relaxed <- object@Relaxed
    alpha <- object@alpha
    ZFull <- object@ZFull
    RVAR <- object@Relaxed
    refit_fraction <- object@refit_fraction
    intercept <- object@intercept
    Y <- ZFull$Y
    k <- ncol(object@Data)
    n_oos_obs <- length((T2 + 1):T3)
    separate_lambdas <- object@separate_lambdas
    if (rolling_oos) {
        MSFE <- object@InSampMSFE
        if (separate_lambdas) {
            MSFE_oos <- matrix(NA, nrow = n_oos_obs, ncol = ncol(Y))
            lambda_evolve <- matrix(NA, nrow = n_oos_obs, ncol = ncol(Y))
        } else {
            MSFE_oos <- matrix(NA, nrow = n_oos_obs, ncol = 1)
            lambda_evolve <- matrix(NA, nrow = n_oos_obs, ncol = 1)
        }
    } else {
        MSFE_oos <- rep(NA, n_oos_obs)
        lambda_evolve <- matrix(lambda, nrow = n_oos_obs, ncol = ncol(object@LambdaGrid))
        gran2 <- 1
    }
    ONESE <- object@ONESE
    if (VARXI) {
        k1 <- object@VARX$k
        s <- object@VARX$s
        contemp <- object@VARX$contemp
        s1 <- object@VARX$contemp
        m <- k - k1
        beta <- array(0, dim = c(k1, k1 * p + (k - k1) * (s + s1) + 1, gran2 * length(alpha)))
    } else {
        k1 <- 0
        s <- 0
        contemp <- FALSE
        s1 <- 0
        m <- 0
        beta <- array(0, dim = c(k, k * p + 1, gran2 * length(alpha)))
    }
    C <- object@constvec
    s1 <- 0
    if (exists("contemp", where = VARX)) {
        if (VARX$contemp) {
            contemp <- TRUE
            s1 <- 1
        } else {
            contemp <- FALSE
            s1 <- 0
        }
    } else {
        contemp <- FALSE
        s1 <- 0
    }
    if (group == "Tapered") {
        palpha <- object@alpha
        gran2 <- length(lambda) * length(palpha)
        beta <- array(0, dim = c(k, k * p + 1, gran2))
        tapered <- TRUE
    } else {
        palpha <- NULL
        tapered <- FALSE
    }
    preds <- matrix(NA, nrow = n_oos_obs, ncol = ncol(Y))
    if (verbose) {
        print("Evaluation Stage")
        pb <- txtProgressBar(min = T2 - h + 2, max = T3 - h, style = 3)
    }
    betaArray <- array(0, dim = c(dim(beta)[1], dim(beta)[2], n_oos_obs))
    grps <- create_group_indexes(group, p, k, gran2 * length(alpha), VARXI, k1, s + s1)
    groups <- grps$groups
    compgroups <- grps$compgroups
    activeset <- grps$activeset
    starting_eigvals <- grps$starting_eigvals
    betaWS <- beta
    for (v in (T2 - h + 2):T3) {
        ## if (h > 1 & !recursive) {
        ##     if (window.size != 0) {
        ##         ws1 <- max(c(v - window.size - h, 1))
        ##         trainY <- ZFull$Y[(ws1 + h):(v - 1), , drop = F]
        ##         trainZ <- ZFull$Z[, (ws1 + h):(v - h),drop=F]
        ##     } else {
        ##         trainY <- ZFull$Y[(h):(v - 1), , drop = F]
        ##         trainZ <- ZFull$Z[, 1:(v - h),drop=F]
        ##     }
        ## } else {
        ##     if (window.size != 0) {
        ##         ws1 <- max(c(v - window.size, 1))
        ##         trainY <- ZFull$Y[(ws1):(v - 1), , drop = F]
        ##         trainZ <- ZFull$Z[, (ws1):(v - 1),drop=F]
        ##     } else {
        ##         trainY <- ZFull$Y[(1):(v - 1), , drop = F]
        ##         trainZ <- ZFull$Z[, (1):(v - 1),drop=F]
        ##     }
        ## }
        if (h > 1 & !recursive) {
            if (window.size != 0) {
                ws1 <- max(c(v - window.size - h, 1))
                index_y <- (ws1 + h-1):(v - 1)
                index_z <- (ws1 ):(v - h)
                trainY <- ZFull$Y[index_y, , drop = FALSE]
                trainZ <- ZFull$Z[, index_z,drop=FALSE]
            } else {
                index_y <- (h):(v - 1)
                index_z <- 1:(v - h)
                trainY <- ZFull$Y[index_y, , drop = FALSE]
                trainZ <- ZFull$Z[, index_z, drop = F]
            }
        } else {
            if (window.size != 0) {
                ws1 <- max(c(v - window.size, 1))
                trainY <- ZFull$Y[(ws1):(v - 1), , drop = FALSE]
                trainZ <- ZFull$Z[, (ws1):(v - 1), drop = FALSE]
            } else {
                trainY <- ZFull$Y[(1):(v - 1), , drop = FALSE]
                trainZ <- ZFull$Z[, (1):(v - 1), drop = FALSE]
            }
        }

        if (v + h - 1 > T3) {
            break
        }
        dual <- FALSE
        temp <- .BigVAR.fit(group, betaWS, trainZ, trainY, lambda, tol, p, m, k1, k, s, s1, MN, C, intercept, separate_lambdas,
                            dual, activeset, starting_eigvals, groups, compgroups, VARXI, alpha, palpha, gamma)
        eZ <- c(1, ZFull$Z[, v])
        beta <- temp$beta
        betaWS <- temp$beta
        if (MN) {
            for (i in 1:dim(betaWS)[3]) {
                submat <- adrop(betaWS[1:dim(betaWS)[1], 1:dim(betaWS)[1], i, drop = F], drop = 3)
                diag(submat) <- diag(submat) - C
                betaWS[1:dim(betaWS)[1], 1:dim(betaWS)[1], i] <- submat
            }
        }
        activeset <- temp$activeset
        msfe_index <- v - (T2 - h) - 1
        temp_results <- refine_and_forecast(beta, as.matrix(eZ), trainZ, trainY, ZFull$Y[v + h - 1, , drop = FALSE], lambda = lambda,
                                            h = h, recursive = recursive, MN = MN, RVAR = RVAR, refit_fraction = refit_fraction, separate_lambdas = separate_lambdas,
                                            C = C, inds = NULL, loss = loss, delta = delta, k = k, p = p, k1 = k1, s = s, oos = !rolling_oos)
        if (!rolling_oos) {
            if (separate_lambdas) {
                MSFE_oos[msfe_index] <- temp_results$MSFE[optind]
                for (col in seq_len(ncol(temp_results$MSFE))) {
                    preds[msfe_index, col] <- temp_results$preds[col, optind[col]]
                    betaArray[col, , msfe_index] <- temp$beta[col, , optind[col]]
                }
            } else {
                MSFE_oos[msfe_index] <- temp_results$MSFE
                preds[msfe_index, ] <- temp_results$preds[, optind]
                betaArray[, , msfe_index] <- temp_results$betaArray
            }
        } else {
            if (!separate_lambdas) {
                MSFE_temp <- temp_results$MSFE[optind]
                MSFE_oos[msfe_index, ] <- MSFE_temp
                MSFE <- rbind(MSFE, temp_results$MSFE)
                preds[msfe_index, ] <- temp_results$preds[, optind]
                betaArray[, , msfe_index] <- temp_results$betaArray
            } else {
                opts <- c()
                for (col in seq_len(ncol(temp_results$MSFE))) {
                    opts[col] <- temp_results$MSFE[optind[col], col]
                    betaArray[col, , msfe_index] <- temp$beta[col, , optind[col]]
                    preds[msfe_index, col] <- temp_results$preds[col, optind[col]]
                }
                MSFE_oos[msfe_index, ] <- opts
                MSFE <- abind(MSFE, temp_results$MSFE, along = 1)
            }
            opt_lambdas <- find_optimal_lambda(MSFE, lambda, palpha, alpha, ONESE, dual, separate_lambdas, tapered)
            optind <- opt_lambdas$optind
            lambdaopt <- opt_lambdas$lambdaopt
            lambda_evolve[msfe_index, ] <- lambdaopt
        }
    }
    temp <- .BigVAR.fit(group, beta, ZFull$Z, ZFull$Y, lambda, tol, p, m, k1, k, s, s1, MN, C, intercept, separate_lambdas,
                        dual, activeset, starting_eigvals, groups, compgroups, VARXI, alpha, palpha)
    betas_full <- temp$beta
    if (separate_lambdas) {
        betaPred <- matrix(0, nrow = ncol(Y), ncol = dim(betas_full)[2])
        for (col in seq_len(ncol(Y))) betaPred[col, ] <- betas_full[col, , optind[col]]
    } else {
        betaPred <- as.matrix(betas_full[, , optind])
    }
    return(list(MSFE = MSFE_oos, betaPred = betaPred, predictions = preds, betaArray = betaArray, lambda_evolve = lambda_evolve))
})

#' BigVAR Estimation
#' @description
#' Fit a BigVAR object with a structured penalty (VARX-L or HLAG).
#' @usage BigVAR.est(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details Fits HLAG or VARX-L model on a BigVAR object.  Does not perform cross-validation.  This method allows the user to construct their own penalty parameter selection procedure.
#' @return An array of \eqn{k \times kp \times n} or \eqn{k\times kp+ms \times n} coefficient matrices; one for each of the n values of lambda.
#' @seealso \code{\link{constructModel}}, \code{\link{BigVAR.results}},\code{\link{cv.BigVAR}} 
#' @name BigVAR.est
#' @aliases BigVAR.est,BigVAR-method
#' @docType methods
#' @rdname BigVAR.est-methods
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' #construct a Basic VAR-L
#' Model1=constructModel(Y,p=4,struct='Basic',gran=c(50,10))
#' BigVAR.est(Model1)
#' @export
setGeneric(name = "BigVAR.est", def = function(object) {
    standardGeneric("BigVAR.est")
})
setMethod(f = "BigVAR.est", signature = "BigVAR", definition = function(object) {
    p <- object@lagmax
    s1 <- 0
    Y <- object@Data
    k <- ncol(Y)
    alpha <- object@alpha
    gamma <- object@gamma
    RVAR <- object@Relaxed
    refit_fraction <- object@refit_fraction
    group <- object@Structure
    intercept <- object@intercept
    VARX <- object@VARX
    tol <- object@tol
    loss <- object@loss
    linear <- object@linear
    delta <- object@delta
    if (length(alpha) == 0) {
        if (length(VARX) > 0) {
            alpha <- 1/(VARX$k + 1)
        } else {
            alpha <- 1/(k + 1)
        }
    }
    C <- object@constvec
    if (length(alpha) > 1 & group %in% c("SparseLag", "SparseOO", "BasicEN")) {
        dual <- TRUE
    } else {
        dual <- FALSE
    }
    MN <- object@Minnesota
    h <- object@horizon
    separate_lambdas <- object@separate_lambdas
    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }
    s <- ifelse(length(object@VARX) != 0, object@VARX$s, 0)
    if (object@ownlambdas) {
        lambda <- object@Granularity
        gran2 <- length(lambda)
        if (gran2 == 1) {
            ONESE <- FALSE
        }
    } else {
        gran2 <- object@Granularity[2]
        gran1 <- object@Granularity[1]
    }
    trainZ <- object@model_data$trainZ
    trainY <- object@model_data$trainY
    VARXI <- object@VARXI
    if (VARXI) {
        k1 <- object@VARX$k
        s <- object@VARX$s
        contemp <- object@VARX$contemp
        s1 <- object@VARX$contemp
        m <- ncol(Y) - ncol(trainY)
        beta <- array(0, dim = c(k1, k1 * p + (k - k1) * (s + s1) + 1, gran2 * length(alpha)))
    } else {
        k1 <- 0
        s <- 0
        contemp <- FALSE
        s1 <- 0
        m <- 0
        beta <- array(0, dim = c(k, k * p + 1, gran2 * length(alpha)))
    }
    grps <- create_group_indexes(group, p, k, gran2 * length(alpha), VARXI, k1, s + s1)
    groups <- grps$groups
    compgroups <- grps$compgroups
    activeset <- grps$activeset
    starting_eigvals <- grps$starting_eigvals
    if (object@ownlambdas) {
        lambda <- object@Granularity
    } else {
        gran2 <- object@Granularity[2]
        gran1 <- object@Granularity[1]
        lambda <- create_lambda_grid(trainY, trainZ, lapply(groups, function(x) {
            x + 1
        }), gran1, gran2, group, p, k1, s + s1, m, k, MN, alpha, C, intercept, tol, VARXI, separate_lambdas, dual, gamma,
        linear, verbose = FALSE)
    }
    h <- object@horizon
    ZFull <- list()
    if (!is.matrix(trainZ)) {
        trainZ <- matrix(trainZ, ncol = 1)
    }
    if (!is.matrix(trainY)) {
        trainY <- matrix(trainY, ncol = 1)
    }
    ZFull$Z <- trainZ
    ZFull$Y <- trainY
    T <- nrow(trainY)
    if (group == "Tapered") {
        palpha <- seq(0, 1, length = 10)
        palpha <- rev(palpha)
        gran2 <- length(lambda) * length(palpha)
        beta <- array(0, dim = c(k, k * p + 1, gran2))
    } else {
        palpha <- NULL
    }
    ## temp <- .BigVAR.fit(group, beta, trainZ, trainY, lambda, tol, p, m, k1, k, s, s1, MN, C, intercept,
    ## separate_lambdas, dual, activeset, q1a, jj, jjcomp, VARX, alpha, kk, palpha)
    temp <- .BigVAR.fit(group, beta, trainZ, trainY, lambda, tol, p, m, k1, k, s, s1, MN, C, intercept, separate_lambdas,
                        dual, activeset, starting_eigvals, groups, compgroups, VARXI, alpha, palpha, gamma)
    beta <- temp$beta
    if (RVAR) {
        beta_rls <- RelaxedLS(cbind(t(trainZ), trainY), beta)
        beta <- (1 - refit_fraction) * beta + refit_fraction * beta_rls

    }
    activeset <- temp$activeset
    q1a <- temp$q1a
    return(list(B = beta, lambdas = lambda))
})



## bigVAR results, inherits class bigVAR, prints results from cv.bigVAR


#' BigVAR.results
#' This class contains the results from cv.BigVAR.
#'
#' It inherits the class BigVAR, but contains substantially more information. 
#' 
#' @field InSampMSFE In-sample MSFE from optimal value of lambda
#' @field LambdaGrid Grid of candidate lambda values
#' @field index Rank of optimal lambda value 
#' @field OptimalLambda Value of lambda that minimizes MSFE
#' @field OOSMSFE Average Out of sample MSFE of BigVAR model with optimal lambda
#' @field seoosfmsfe Standard error of out of sample MSFE of BigVAR model with optimal lambda
#' @field MeanMSFE Average out of sample MSFE of unconditional mean forecast
#' @field MeanSD Standard error of out of sample MSFE of unconditional mean forecast
#' @field MeanPreds predictions from conditional mean model
#' @field RWMSFE Average out of sample MSFE of random walk forecast
#' @field RWPreds Predictions from random walk model
#' @field RWSD Standard error of out of sample MSFE of random walk forecast
#' @field AICMSFE Average out of sample MSFE of AIC forecast
#' @field AICSD Standard error of out of sample MSFE of AIC forecast
#' @field AICPreds Predictions from AIC VAR/VARX model
#' @field AICpvec Lag orders selected from AIC VAR model
#' @field AICpvec Lag orders selected from AIC VARX model
#' @field BICMSFE Average out of sample MSFE of BIC forecast
#' @field BICSD Standard error of out of sample MSFE of BIC forecast
#' @field BICPreds Predictions from BIC VAR/VARX model
#' @field BICpvec Lag orders selected from BIC VAR model
#' @field BICpvec Lag orders selected from BIC VARX model
#' @field betaPred The final estimated \eqn{k\times kp+ms+1} coefficient matrix, to be used for prediction
#' @field Zvals The final lagged values of \code{Y}, to be used for prediction
#' @field fitted fitted values obtained from betaPred
#' @field resids residuals obtained from betaPred
#' @field Data a \eqn{T \times k} or \eqn{T\times k + m} multivariate time Series
#' @field lagmax Maximal lag order
#' @field Structure Penalty structure
#' @field Relaxed Indicator for relaxed VAR
#' @field Granularity Granularity of penalty grid
#' @field horizon Desired forecast horizon
#' @field crossval Cross-Validation procedure
#' @field alpha additional penalty parameter for Sparse Lag Group or Sparse Own/Other methods. Will contain either the heuristic choice of \eqn{1/(k+1)} or the value selected by cross validation if the argument \code{dual} is set to \code{TRUE}
#' @field VARXI VARX Indicator 
#' @field Minnesota Minnesota Prior Indicator
#' @field verbose  verbose indicator
#' @field dual indicator as to whether dual cross validation was conducted
#' @field contemp indicator if contemporaneous exogenous predictors are used
#' @field lagmatrix matrix of lagged values used to compute residuals (of which Zvals is the final column)
#' @field betaArray array of VAR/VARX coefficients from out of sample forecasts
#' @field sparse_count average fraction of active coefficients in validation period
#' @field lambda_evolve_path evolution of lambda over evaluation period
#'
#' @note One can also access any object of class BigVAR from BigVAR.results
#' @name BigVAR.results 
#' @rdname BigVAR.results
#' @aliases BigVAR.results-class
#' @exportClass BigVAR.results
#' @author Will Nicholson
#' @export
setClass("BigVAR.results", representation(InSampMSFE = "array", InSampSD = "numeric", LambdaGrid = "matrix", index = "numeric",
                                          OptimalLambda = "numeric", OOSMSFE = "matrix", seoosmsfe = "numeric", MeanMSFE = "numeric", AICMSFE = "numeric", AICPreds = "matrix",
                                          BICMSFE = "numeric", BICpvec = "numeric", BICsvec = "numeric", AICpvec = "numeric", AICsvec = "numeric", BICSD = "numeric",
                                          BICPreds = "matrix", RWMSFE = "numeric", RWPreds = "matrix", MeanSD = "numeric", MeanPreds = "matrix", AICSD = "numeric",
                                          RWSD = "numeric", betaPred = "matrix", Zvals = "matrix", VARXI = "logical", resids = "matrix", preds = "matrix", dual = "logical",
                                          contemp = "logical", fitted = "matrix", lagmatrix = "matrix", betaArray = "array", sparse_count = "numeric", lambda_evolve_path = "matrix"),
         contains = "BigVAR")


#' Plot an object of class BigVAR.results
#' 
#' @param x BigVAR.results object created from \code{cv.BigVAR}
#' @param y NULL
#' @param ... additional arguments
#' @details Plots the in sample MSFE of all values of lambda with the optimal value highlighted.
#' @name plot
#' @import methods
#' @aliases plot,BigVAR.results-method
#' @aliases plot-methods
#' @docType methods
#' @method plot method
#' @rdname BigVAR.results-plot-methods
#' @importFrom graphics abline
#' @importFrom graphics par
#' @export
setMethod(f = "plot", signature = "BigVAR.results", definition = function(x, y = NULL, ...) {
    if (!x@separate_lambdas) {
        plot(x@LambdaGrid, colMeans(x@InSampMSFE[, , 1]), type = "o", xlab = "Value of Lambda", ylab = "MSFE", log = "x")
        abline(v = x@OptimalLambda, col = "green")
    } else {
        k <- ncol(x@Data)
        par(mfrow = c(k, 1))
        for (i in 1:k) {
            plot(x@LambdaGrid[, i], colMeans(x@InSampMSFE[, , i]), type = "o", xlab = "Value of Lambda", ylab = "MSFE", log = "x")
            abline(v = x@OptimalLambda[i], col = "green")
        }
    }
})

#' Default show method for an object of class BigVAR.results
#' 
#' @param object BigVAR.results object created from \code{cv.BigVAR}
#' @details prints forecast results and additional diagnostic information as well as comparisons with mean, random walk, and AIC, and BIC benchmarks
#' @seealso \code{\link{cv.BigVAR}},\code{\link{BigVAR.results}} 
#' @name show
#' @aliases show,BigVAR.results-method
#' @docType methods
#' @method show method
#' @rdname show-methods-BigVAR.results
#' @export
setMethod("show", "BigVAR.results", function(object) {
    cat("*** BIGVAR MODEL Results *** \n")
    cat("Structure\n")
    print(object@Structure)
    if (object@Relaxed == TRUE) {
        cat("Relaxed VAR \n")
        print(object@Relaxed)
    }
    cat("Loss \n")
    print(object@loss)
    cat("Forecast Horizon \n")
    print(object@horizon)
    cat("Minnesota VAR\n")
    print(object@Minnesota)
    if (object@VARXI) {
        cat("VARX Specs \n")
        print(object@VARX)
    }
    cat("Maximum Lag Order \n")
    print(object@lagmax)
    cat("Optimal Lambda \n")
    print(signif(object@OptimalLambda, digits = 4))
    if (object@dual) {
        cat("Optimal Alpha \n")
        print(signif(object@alpha, digits = 2))
    }
    cat("Grid Depth \n")
    print(object@Granularity[1])
    cat("Index of Optimal Lambda \n")
    print(object@index)
    cat("Fraction of active coefficients \n")
    print(signif(object@sparse_count, digits = 4))
    if (!object@separate_lambdas) {
        cat("In-Sample Loss\n")
        print(signif(mean(object@InSampMSFE[, object@index, ]), digits = 3))
    } else {
        cat("In-Sample Loss\n")
        print(signif(apply(object@InSampMSFE[, object@index, ], 2, mean), digits = 3))
    }
    cat("BigVAR Out of Sample Loss\n")
    print(signif(mean(object@OOSMSFE), digits = 3))
    cat("*** Benchmark Results *** \n")
    cat("Conditional Mean Out of Sample Loss\n")
    print(signif(object@MeanMSFE, digits = 3))
    cat("AIC Out of Sample Loss\n")
    print(signif(object@AICMSFE, digits = 3))
    cat("BIC Out of Sample Loss\n")
    print(signif(object@BICMSFE, digits = 3))
    cat("RW Out of Sample Loss\n")
    print(signif(object@RWMSFE, digits = 3))
})



#' Forecast using a BigVAR.results object
#' 
#' @usage predict(object,...)
#' @param object BigVAR.results object from \code{cv.BigVAR}
#' @param ... additional arguments affecting the predictions produced (e.g. \code{n.ahead}, \code{confint})
#' @details Provides \code{n.ahead} step forecasts using the model produced by cv.BigVAR.  If \code{confint} is set to \code{TRUE}, a 95 percent confidence interval will also be returned. 
#' @seealso \code{\link{cv.BigVAR}} 
#' @name predict
#' @aliases predict,BigVAR.results-method
#' @docType methods
#' @method predict method
#' @rdname predict-methods-BigVAR.results
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' Model1=constructModel(Y,p=4,struct='Basic',gran=c(50,10),verbose=FALSE)
#' results=cv.BigVAR(Model1)
#' predict(results,n.ahead=1)
#' @export
setMethod("predict", "BigVAR.results", function(object, n.ahead = 1, newxreg = NULL, predict_all = FALSE, confint = FALSE,
                                                ...) {
    MN <- object@Minnesota
    eZ <- object@Zvals
    betaPred <- object@betaPred
    Y <- object@Data
    k <- object@VARX$k
    m <- ncol(object@Data) - k
    p <- object@lagmax
    s <- object@VARX$s
    VARX <- object@VARXI
    contemp <- object@contemp
    s1 <- 0
    if (confint) {
        YZ <- object@model_data
        ci = create_sigma_u(YZ$trainY, YZ$trainZ, betaPred, n.ahead)
    }
    fcst <- matrix(betaPred %*% eZ, ncol = 1)
    fcst_full <- fcst
    if (n.ahead == 1) {
        if (confint) {
            lower <- fcst + ci[, 1]
            upper <- fcst + ci[, 2]
            fcst <- as.data.frame(cbind(fcst, lower, upper))
            names(fcst) <- c("forecast", "lower", "upper")
        }
        return(fcst)
    } else {
        if (!VARX) {
            fcst <- predictMS(matrix(fcst, nrow = 1), Y[(nrow(Y) - p + 1):nrow(Y), ], n.ahead - 1, betaPred, p, MN, predict_all = predict_all)
            if (predict_all) {
                row.names(fcst) <- paste0("T+", 1:n.ahead)
            } else {
                fcst <- t(fcst)
            }
        } else {
            if (is.null(newxreg)) {
                stop("Need new data for multi-step VARX forecasts.  Re-run with new data in newxreg")
            } else {
                if (nrow(newxreg) < n.ahead - 1) {
                    stop(paste("Need at least ", n.ahead - 1, "rows of new data"))
                }
                C <- max(p, s)
                if (contemp) {
                    C <- C + 3
                }
                fcst <- matrix(predictMSX(matrix(fcst, nrow = 1), as.matrix(Y[(nrow(Y) - C + 1):nrow(Y), 1:(k)]), n.ahead -
                                                                                                                  1, betaPred, p, newxreg, matrix(Y[(nrow(Y) - C + 1):nrow(Y), (ncol(Y) - m + 1):ncol(Y)], ncol = m), m,
                                          s, 1, MN, contemp), ncol = 1)
            }
        }
    }
    if (confint) {
        lower <- fcst + ci[, 1]
        upper <- fcst + ci[, 2]
        fcst <- as.data.frame(cbind(fcst, lower, upper))
        names(fcst) <- c("forecast", "lower", "upper")
    }
    return(fcst)
})

#' Sparsity Plot of a BigVAR.results object 
#'
#' @param object BigVAR.results object
#' @return NA, side effect is graph
#' @details Uses \code{levelplot} from the \code{lattice} package to plot the magnitude of each coefficient in the last coefficient estimated by \code{cv.BigVAR}.
#' @name SparsityPlot.BigVAR.results
#' @aliases SparsityPlot.BigVAR.results,BigVAR.results-method
#' @seealso \code{\link{cv.BigVAR}}, \code{\link{BigVAR.results}}
#' @docType methods
#' @rdname SparsityPlot.BigVAR.results-methods
#' @examples
#' data(Y)
#' Y <- Y[1:100,]
#' Model1 <- constructModel(Y,p=4,struct='Basic',gran=c(50,10),verbose=FALSE)
#' SparsityPlot.BigVAR.results(cv.BigVAR(Model1))
#' @export
#' @importFrom lattice levelplot
#' @importFrom lattice panel.abline
#' @importFrom lattice panel.levelplot
#' @importFrom grDevices colorRampPalette
setGeneric(name = "SparsityPlot.BigVAR.results", def = function(object) {
    standardGeneric("SparsityPlot.BigVAR.results")
})
setMethod(f = "SparsityPlot.BigVAR.results", signature = "BigVAR.results", definition = function(object) {
    B <- object@betaPred
    if (nrow(B) == 1) {
        B <- matrix(B[, 2:ncol(B)], nrow = 1)
    } else {
        B <- B[, 2:ncol(B)]
    }
    k <- nrow(B)
    p <- object@lagmax
    s1 <- 0
    if (length(object@VARX != 0)) {
        s <- object@VARX$s
        m <- ncol(object@Data) - object@VARX$k
        contemp <- object@VARX$contemp
        if (!is.null(contemp)) {
            if (contemp) {
                s1 <- 1
            }
        } else {
            s1 <- 0
        }
    } else {
        m <- 0
        s <- 0
    }
    s <- s + s1
    text <- c()
    for (i in 1:p) {
        text1 <- as.expression(bquote(bold(Phi)^(.(i))))
        text <- append(text, text1)
    }
    if (m > 0) {
        for (i in (p + 1):(p + s + 1)) {
            text1 <- as.expression(bquote(bold(beta)^(.(i - p - s1))))
            text <- append(text, text1)
        }
    }
    f <- function(m) t(m)[, nrow(m):1]
    rgb.palette <- colorRampPalette(c("white", "blue"), space = "Lab")
    at <- seq(k/2 + 0.5, p * (k) + 0.5, by = k)
    if (m > 0) {
        at2 <- seq(p * k + m/2 + 0.5, p * k + s * m + 0.5, by = m)
    } else {
        at2 <- c()
    }
    at <- c(at, at2)
    se2 <- seq(1.75, by = k, length = k)
    L2 <- levelplot(as.matrix(f(abs(B))), col.regions = rgb.palette, colorkey = NULL, xlab = NULL, ylab = NULL, main = list(label = "Sparsity Pattern Generated by BigVAR",
                                                                                                                            cex = 1), panel = function(...) {
                                                                                                                                panel.levelplot(...)
                                                                                                                                panel.abline(a = NULL, b = 1, h = seq(1.5, m * s + p * k + 0.5, by = 1), v = seq(1.5, by = 1, length = p * k + m *
                                                                                                                                                                                                                                  s), lwd = 0.5)
                                                                                                                                bl1 <- seq(k + 0.5, p * k + 0.5, by = k)
                                                                                                                                b23 <- seq(p * k + 0.5, p * k + 0.5 + s * m, by = m)
                                                                                                                                b1 <- c(bl1, b23)
                                                                                                                                panel.abline(a = NULL, b = 1, v = p * k + 0.5, lwd = 3)
                                                                                                                                panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
                                                                                                                            }, scales = list(x = list(alternating = 1, labels = text, cex = 1, at = at, tck = c(0, 0)), y = list(alternating = 0,
                                                                                                                                                                                                                                 tck = c(0, 0))))
    return(L2)
})


                                        # show-default method to show an object when its name is printed in the console.
#' Default coef method BigVAR-results, returns the last coefficient matrix from the evaluation period
#'
#' @param object BigVAR.results object created from \code{cv.BigVAR}
#' @details displays formatted coefficient matrix
#' @name coef
#' @import methods
#' @aliases coef,BigVAR.results-method
#' @aliases coef-methods
#' @docType methods
#' @method coef method
#' @rdname BigVAR.results-coef-methods
#' @export
setMethod(f = "coef", signature = "BigVAR.results", definition = function(object) {
    B <- data.frame(object@betaPred)
    k <- nrow(B)
    p <- object@lagmax
    Y <- object@Data
    intercept <- object@intercept
    if (!intercept) {
        B <- B[, 2:ncol(B)]
    }
    row.names(B) <- paste0("Y", 1:k)
    if (!object@VARXI) {
        bnames <- c(outer(X = paste0("Y", 1:k), Y = paste0("L", 1:p), paste0))
        if (intercept) {
            bnames <- c("intercept", bnames)
        }
        names(B) <- bnames
    } else {
        if (p > 0) {
            bnames <- c(outer(X = paste0("Y", 1:k), Y = paste0("L", 1:p), paste0))
        } else {
            bnames <- NULL
        }
        m <- ncol(Y) - k
        s <- object@VARX$s
        if (!is.null(object@VARX$contemp)) {
            bnamesX <- c(outer(X = paste0("X", 1:m), Y = paste0("L", 1:s), paste0))
        } else {
            bnamesX <- c(outer(X = paste0("X", 1:m), Y = paste0("L", 0:s), paste0))
        }
        bnames <- c(bnames, bnamesX)
        if (intercept) {
            bnames <- c("intercept", bnames)
        }
        names(B) <- bnames
    }
    return(B)
})

### Support functions for BigVAR Package These are mostly utility functions that will not be seen by the user




#' Construct a VAR or VARX lag matrix
#' 
#' @param Y a \eqn{T \times k} matrix of endogenous (modeled) series
#' @param X a \eqn{T \times m} matrix of exogenous (unmodeled) series (default NULL)
#' @param p Endogenous Lag order
#' @param s exogenous lag order (default zero)
#' @param oos indicator as to whether the data should be constructed for out of sample prediction (i.e. last available entries of Y as final lags default FALSE)
#' @param contemp indicator as to whether to use contemporaneous exogenous predictors (for example, if exogenous series become available before exogenous default FALSE).
#' @return list with two entries:
#' \itemize{
#' \item{'Z'}{\eqn{kp+ms+1\times T-max(p,s)} VARX lag matrix}
#' \item{'Y'}{adjusted \eqn{k\times T-max(p,s)} endogenous series}
#' }
#' @details This function is not required unless you which to design your own cross validation routine. 
#' @references
#' See page 15 of Lutkepohl, 'A New Introduction to Multiple Time Series Analysis
#' @seealso \code{\link{MultVarSim}}
#' @examples
#' data(Y)
#' # construct VAR lag matrix with p=4
#' ZZ<-VARXLagCons(Y,X=NULL,p=4,s=0)
#' @export
VARXLagCons <- function(Y, X = NULL, p, s = 0, oos = FALSE, contemp = FALSE) {
    if (is.null(X)) {
        X <- matrix(0, nrow = nrow(Y))
    }
    if (nrow(Y) != nrow(X)) {
        stop("Y and X must have same dimensions")
    }
    if (s == 0 & !contemp) {
        m <- 0
    } else {
        m <- ncol(X)
    }

    if (p < 0 | m < 0) {
        stop("lag orders must be positive")
    }
    k <- ifelse(is.null(Y), 0, ncol(Y))
    XX <- VARXCons(Y, X, k, p, m, s, oos, contemp)
    Y <- t(Y[(max(c(p, s)) + 1):nrow(Y), ])
    return(list(Z = XX, Y = Y))
}

# VARX Construction for constructModel
VARXConsModel <- function(Y, p, VARX, tf) {
    k <- ncol(Y)
    if (length(VARX) != 0) {

        VARXI <- TRUE
        k1 <- VARX$k
        s <- VARX$s

        if (!is.null(VARX$contemp)) {

            contemp <- TRUE
            s1 <- 1

        } else {

            contemp <- FALSE
            s1 <- 0
        }

        m <- k - k1
        Y1 <- matrix(Y[, 1:k1, drop = F], ncol = k1)
        X <- matrix(Y[, (ncol(Y) - m + 1):ncol(Y), drop = F], ncol = m)

        if (!tf) {
            trainZ <- VARXCons(Y1, X, k1, p, m, s, contemp = contemp)

        } else {
            trainZ <- VARXCons(matrix(0, ncol = 1, nrow = nrow(X)), matrix(X, ncol = m), k = 0, p = 0, m = m, s = s, contemp = contemp,
                oos = FALSE)

        }

        trainZ <- trainZ[2:nrow(trainZ), , drop = F]

        trainY <- matrix(Y[(max(c(p, s)) + 1):nrow(Y), 1:k1, drop = F], ncol = k1)

    } else {
        # VAR setting
        VARXI <- FALSE
        contemp <- FALSE
        k1 <- 0
        s <- 0
        m <- 0
        s1 <- 0

        Z1 <- VARXCons(Y, matrix(0, nrow = nrow(Y)), k, p, 0, 0)

        trainZ <- Z1[2:nrow(Z1), ]

        trainY <- matrix(Y[(p + 1):nrow(Y), , drop = FALSE], ncol = k)

    }
    return(list(trainY = trainY, trainZ = trainZ, s1 = s1, contemp = contemp, VARXI = VARXI))
}


.check_is_matrix <- function(x) {
    !is.null(attr(x, "dim"))
}

.huber_loss <- function(r, delta) {
    l <- ifelse(abs(r) < delta, 1/2 * abs(r)^2, delta * (abs(r) - 1/2 * delta))
    return(l)
}

.calc.loss <- function(x, univ = FALSE, loss, delta) {
    if (loss == "L1") {
        l <- sum(abs(x))
    } else if (loss == "L2") {
        if (univ) {
            l <- x^2
        } else {
            l <- norm2(x)^2
        }

    } else if (loss == "Huber") {
        l <- .huber_loss(x, delta)

    }
    return(sum(l))
}
# mean benchmark
.evalMean <- function(Y, T1, T2, h = 1, loss = "L2", delta = 2.5) {

    ypredF <- NULL
    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }

    MSFE <- c()

    k <- ncol(Y)

    for (u in (T1 - h + 2):T2) {

        if (h + u - 1 > T2) {
            break
        }

        trainY1 <- Y[1:(u - 1), ]

        if (k > 1) {

            ypred <- colMeans(trainY1)
            ypredF <- rbind(ypredF, ypred)
        } else {
            ypred <- mean(trainY1)
            ypredF <- c(ypredF, ypred)
        }

        uhat <- matrix(Y[u + h - 1, ] - ypred, ncol = k)

        MSFE <- c(MSFE, .calc.loss(uhat, univ = FALSE, loss, delta))

    }
    ypredF <- unname(ypredF)
    return(list(Mean = mean(na.omit(MSFE)), SD = sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE))), preds = as.matrix(ypredF)))
}

# random walk benchmark
.evalRW <- function(Y, T1, T2, h = 1, loss = "L2", delta = 2.5) {

    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }
    ypredF <- NULL
    MSFE <- c()

    k <- ncol(Y)

    for (u in (T1 - h + 2):T2) {

        if (h + u - 1 > T2) {
            break
        }

        trainY1 <- Y[u - 1, ]
        ypredF <- rbind(ypredF, trainY1)
        uhat <- matrix(Y[u + h - 1, ] - trainY1, ncol = k)

        MSFE <- c(MSFE, .calc.loss(uhat, univ = FALSE, loss, delta))

    }
    ypredF <- unname(ypredF)
    return(list(Mean = mean(na.omit(MSFE)), SD = sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE))), preds = as.matrix(ypredF)))
}

# Construct Lambda Grid:
.LambdaGrid <- function(gran1, gran2, groups, Y, Z, group, p, k1, s, m, k, MN, alpha, C, intercept, tol, VARX = FALSE, separate_lambdas = FALSE,
    verbose = FALSE, gamma = 3, linear) {

    nseries <- ifelse(VARX == TRUE, k1, k)

    if (group == "Lag") {

        mat <- list()

        for (i in seq_len(length(groups))) {

            if (k > 1) {

                mat[[i]] <- norm2(Z[groups[[i]], ] %*% Y)

            } else {
                mat[[i]] <- norm2(t(Y) %*% Z[groups[[i]], ])
            }

        }

        lambdastart <- max(unlist(mat))

    }

    if (group == "Basic" | group == "BasicEN" | group == "Tapered") {
        if (!separate_lambdas) {
            if (group == "Basic" | group == "Tapered") {
                lambdastart <- max(abs(t(Y) %*% t(Z)))
            } else {
                lambdastart <- max(abs(t(Y) %*% t(Z)))/max(c(alpha, 0.01))
            }
        } else {
            lambdastart <- c()
            for (i in 1:nseries) {

                if (group == "Basic" | group == "Tapered") {
                  lambdastart[i] <- max(abs(t(Y[, i, drop = F]) %*% t(Z)))
                } else {
                  lambdastart[i] <- max(abs(t(Y[, i, drop = F]) %*% t(Z)))/max(c(alpha, 0.01))
                }
            }

        }
    }


    if (group == "MCP" | group == "SCAD") {
        if (!separate_lambdas) {
            lambdastart <- max(abs(crossprod(t(Z), Y)))/nrow(Y)
        } else {
            lambdastart <- c()
            for (i in 1:nseries) {
                lambdastart[i] <- max(abs(crossprod(t(Z), Y[, i])))/nrow(Y)
            }
        }
    }


    if (group == "SparseLag") {

        mat <- list()

        if (alpha > 0) {

            for (i in seq_len(length(groups))) {

                if (k > 1) {

                  mat[[i]] <- norm2(Z[groups[[i]], ] %*% Y * (1/(alpha)))

                } else {

                  mat[[i]] <- norm2(t(Y) %*% Z[groups[[i]], ])

                }

            }

            lambdastart <- max(unlist(mat))

        } else {

            lambdastart <- max(t(Y) %*% t(Z))

        }


    }

    if (group == "OwnOther") {

        mat <- list()

        ZZ <- kronecker(t(Z), diag(nseries))

        for (i in seq_len(length(groups))) {

            mat[[i]] <- norm(as.vector(t(Y)) %*% ZZ[, groups[[i]]]/sqrt(length(groups[[i]])), "F")

        }

        lambdastart <- max(unlist(mat))

    }

    if (group == "SparseOO") {

        mat <- list()

        ZZ <- kronecker(t(Z), diag(nseries))

        if (alpha > 0) {
            for (i in seq_len(length(groups))) {

                mat[[i]] <- norm(1/(alpha) * as.vector(t(Y)) %*% ZZ[, groups[[i]]]/sqrt(length(groups[[i]])), "F")
            }

            lambdastart <- max(unlist(mat))

        } else {

            lambdastart <- max(t(Y) %*% t(Z))

        }

    }

    if (group == "EFX") {

        gmax <- c()

        for (i in 1:nseries) {

            gmax[i] <- norm2(Z %*% Y[, i])/sqrt(k * p)

        }

        lambdastart <- max(gmax)

    }

    if (group == "HLAGC" | group == "HLAGOO" | group == "HLAGELEM") {

        gmax <- c()

        for (i in 1:nseries) {

            gmax[i] <- norm2(Z %*% Y[, i])

        }

        if (!separate_lambdas) {
            lambdastart <- max(gmax)
        } else {
            lambdastart <- gmax
        }
    }

    if (VARX) {
        beta <- array(0, dim = c(k1, k1 * p + s * m + 1, 1))
    } else if (group == "Tapered") {

        beta <- array(0, dim = c(k, k * p + 1, 1))

    } else {
        beta <- array(0, dim = c(k, k * p + 1, 1))
    }

    if (!separate_lambdas) {
        lambdastart <- LGSearch(lambdastart, Y, Z, beta, group, k1, p, s, m, groups, k, MN, alpha, C, intercept, tol, VARX,
            gamma)

        if (!linear) {
            lambda <- exp(seq(from = log(lambdastart), to = log(lambdastart/gran1), length = gran2))
        } else {
            lambda <- exp(seq(from = lambdastart, to = lambdastart/gran1, length = gran2))

        }

    } else {

        lambda <- matrix(NA, nrow = c(gran2), ncol = ncol(Y))

        for (i in seq_len(ncol(lambda))) {
            lambdastart[i] <- LGSearch(lambdastart[i], Y, Z, beta, group, k1, p, s, m, groups, k, MN, alpha, C, intercept,
                tol, VARX, gamma)
            lambdastart[i] <- ifelse(lambdastart[i] == 0, 1e-04, lambdastart[i])
            if (verbose & i%%20 == 0) {
                print(sprintf("determined lambda grid for series %s", i))
            }

            if (!linear) {
                lambda[, i] <- exp(seq(from = log(lambdastart[i]), to = log(lambdastart[i]/gran1), length = gran2))
            } else {
                lambda[, i] <- exp(seq(from = lambdastart[i], to = lambdastart[i]/gran1, length = gran2))
            }

        }


    }


    return(lambda)

}

#' Converts a VAR coefficient matrix of order p to multiple companion form
#' 
#' @param B a \eqn{k \times kp} coefficient matrix
#' @param p Lag order
#' @param k Number of Series
#' @return Returns a \eqn{kp \times kp} coefficient matrix representing all coefficient matrices contained in Ai as a VAR(1).
#' @references See page 15 of Lutkepohl, 'A New Introduction to Multiple Time Series Analysis'
#' @seealso \code{\link{MultVarSim}}
#' @examples
#' k=3;p=6
#' B=matrix(0,nrow=k,ncol=p*k)
#' A1<- matrix(c(.4,-.02,.01,-.02,.3,.02,.01,.04,.3),ncol=3,nrow=3)
#' A2 <- matrix(c(.2,0,0,0,.3,0,0,0,.13),ncol=3,nrow=3)
#' B[,1:k]=A1
#' B[,(4*k+1):(5*k)]=A2
#' A <- VarptoVar1MC(B,p,k)
#' @export
VarptoVar1MC <- function(B, p, k) {

    Fp <- matrix(0, nrow = k * p, ncol = k * p)

    Fp[1:k, ] <- B

    Fp[-(1:k), 1:(k * (p - 1))] <- diag(k * (p - 1))
    # We require that the coefficient matrix generates a stationary VAR
    if (max(Mod(eigen(Fp)$values)) > 1) {
        warning("Coefficient Matrix is not stationary")
    }

    return(Fp)

}


#' Simulate a VAR
#' 
#' @param k Number of Series
#' @param A1 Either a \eqn{k \times k} coefficient matrix or a \eqn{kp \times kp} matrix created using \code{\link{VarptoVar1MC}}. 
#' @param p Maximum Lag Order
#' @param Sigma Residual Covariance Matrix of dimension \eqn{k\times k}
#' @param T Number of simulations
#' @return Returns a \eqn{T \times k} of realizations from a VAR.
#' @references Lutkepohl, 'A New Introduction to Multiple Time Series Analysis'
#' @seealso \code{\link{VarptoVar1MC}}
#' @examples
#' k=3;p=6
#' B=matrix(0,nrow=k,ncol=p*k)
#' A1<- matrix(c(.4,-.02,.01,-.02,.3,.02,.01,.04,.3),ncol=3,nrow=3)
#' A2 <- matrix(c(.2,0,0,0,.3,0,0,0,.13),ncol=3,nrow=3)
#' B[,1:k]=A1
#' B[,(4*k+1):(5*k)]=A2
#' A <- VarptoVar1MC(B,p,k)
#' Y <-MultVarSim(k,A,p,.1*diag(k),100)
#' @export
#' @importFrom MASS mvrnorm
MultVarSim <- function(k, A1, p, Sigma, T) {

    if (max(Mod(eigen(A1)$values)) > 1) {
        stop("Error: Generator Matrix is not stationary")
    }

    # add 500 observations for initialization purposes

    Y <- matrix(0, nrow = T + 500 + p, ncol = k)

    YY <- as.vector(Y)

    for (i in seq(from = (k * p + 1), to = (nrow(Y) * k - 1), by = k)) {

        u <- as.vector(c(mvrnorm(1, rep(0, k), Sigma), rep(0, k * p - k)))

        YY[(i + k):(i - k * p + 1 + k)] <- A1 %*% YY[(i):(i - k * p + 1)] + as.matrix(u, ncol = 1)

    }

    YY <- YY[YY != 0]

    Y <- matrix(YY, ncol = k, byrow = TRUE)

    Y <- Y[, c(ncol(Y):1)]

    Y <- Y[501:nrow(Y), ]

    return(Y)
}

# function to create subsets for lag group VARX-L
.groupfun <- function(p, k) {

    jjj <- list()
    jjj[[1]] <- 1:k

    if (p > 1) {

        for (i in 2:p) {
            jjj[[i]] <- jjj[[i - 1]] + k

        }

    }

    return(jjj)

}

# C++ groupings to account for 0 indexing
.groupfuncpp <- function(p, k) {
    jjj <- list()

    jjj[[1]] <- 0:(k - 1)

    if (p > 1) {
        for (i in 2:p) {

            jjj[[i]] <- jjj[[i - 1]] + k
        }
    }
    return(jjj)
}

# subsetting complement of groups in rcpp
.groupfuncomp <- function(p, k) {


    ownoth <- .groupfuncpp(p, k)

    kk2 <- list()

    pmax <- max(unlist(ownoth))

    to <- 0:(pmax)

    for (i in seq_len(length(ownoth))) {

        kk2[[i]] <- to[is.na(pmatch(to, ownoth[[i]]))]
    }

    return(kk2)

}

# Group indexing for own/other VARX-L
.lfunction2 <- function(p, k) {

    kk <- list()

    kk[[1]] <- 1:(k^2)

    if (p > 1) {
        for (i in 2:p) {
            kk[[i]] <- 1:(k^2) + tail(kk[[i - 1]], 1)

        }
    }
    return(kk)
}

.lfunction2cpp <- function(p, k) {

    kk <- list()

    kk[[1]] <- 0:(k^2 - 1)

    if (p > 1)

    {

        for (i in 2:p) for (i in 2:p) {

            kk[[i]] <- 0:(k^2 - 1) + tail(kk[[i - 1]], 1) + 1


        }

    }

    return(kk)


}

.lfunction3 <- function(p, k) {

    kk <- .lfunction2(p, k)

    oo <- list()

    pp <- list()

    for (i in seq_len(length(kk))) {
        j <- 0
        oo[[i]] <- kk[[i]][(seq(1, length(kk[[1]]), k + 1) + (j * k^2))]
        pp[[i]] <- kk[[i]][-(seq(1, length(kk[[1]]), k + 1) + (j * k^2))]
        j <- j + 1

    }

    ownoth <- c(oo, pp)
    return(ownoth)
}

.lfunction3cpp <- function(p, k) {

    kk <- .lfunction2cpp(p, k)
    oo <- list()
    pp <- list()

    for (i in seq_len(length(kk))) {
        j <- 0
        oo[[i]] <- kk[[i]][(seq(1, length(kk[[1]]), k + 1) + (j * k^2))]
        pp[[i]] <- kk[[i]][-(seq(1, length(kk[[1]]), k + 1) + (j * k^2))]
        j <- j + 1
    }

    ownoth <- c(oo, pp)

    return(ownoth)

}


.lfunctioncomp <- function(p, k) {

    ownoth <- .lfunction3cpp(p, k)
    kk2 <- list()
    pmax <- max(unlist(ownoth))
    to <- 0:(pmax)
    for (i in seq_len(length(ownoth))) {
        kk2[[i]] <- to[is.na(pmatch(to, ownoth[[i]]))]
    }

    return(kk2)
}


# This function should work for arbitrary groups
.lfunction <- function(groups, p) {

    H <- as.vector(do.call("cbind", groups))
    kk <- list()
    kk[[1]] <- H
    if (p > 1) {
        for (i in 2:p) {
            kk[[i]] <- as.vector(do.call("cbind", groups)) + tail(kk[[i - 1]], 1)
        }
    }
    return(kk)

}

# Indexing for HLAG
.vsubs <- function(p, k) {
    vi <- list()

    for (i in p:1) {
        g <- max(k * i - k, 0)
        vi[[i]] <- g:(k * (p) - 1)
    }
    return(vi)

}

# indexing for HLAG OO
.oofun <- function(p, k) {

    kk <- .lfunction2(p, k)
    oo <- list()
    pp <- list()

    for (i in seq_len(length(kk))) {
        j <- 0
        oo[[i]] <- kk[[i]][(seq(1, length(kk[[1]]), k + 1) + (j * k^2))]
        pp[[i]] <- kk[[i]][-(seq(1, length(kk[[1]]), k + 1) + (j * k^2))]
        j <- j + 1
    }
    ownoth <- list()
    jj <- .lfunction3(p, k)

    for (i in seq_len(length(jj))) {

        if (i == 1)

        {

            ownoth[[i]] <- jj[[i]]
            oo[[1]] <- NULL

        }

        if (i == length(jj)) {
            ownoth[[i]] <- tail(jj, 1)
            pp[[1]] <- NULL
        }

        if (i != 1 & i%%2 != 0) {

            ownoth[[i]] <- head(oo, 1)
            oo[[1]] <- NULL
        }

        if (i != length(jj) & i%%2 == 0) {
            ownoth[[i]] <- head(pp, 1)
            pp[[1]] <- NULL
        }


    }

    return(rev(ownoth))

}

.oocumfun <- function(p, k) {

    kk <- rev(.oofun(p, k))
    oogroups <- list()
    oogroups[[1]] <- unlist(kk)
    for (i in 2:length(kk)) {
        oogroups[[i]] <- unlist(kk[-(1:(i - 1))])

    }

    return(oogroups)
}


# indexing function for Own/Other HLAG
.vecoovars <- function(p, k, k1) {

    vv <- list()

    vv[[1]] <- 1:(p * k)

    vv[[2]] <- vv[[1]][-k1]

    q1 <- 1

    if (p > 1) {
        for (i in 3:(2 * p)) {
            if (i%%2 != 0) {
                vv[[i]] <- (q1 * k + 1):(k * p)
                q1 <- q1 + 1
            } else {
                vv[[i]] <- vv[[i - 1]][-k1]
            }
        }
    }
    return(vv)

}

# indexing to start at zero for use within rcpp
.vecoovarscpp <- function(p, k, k1) {

    vv <- list()
    vv[[1]] <- 0:(p * k - 1)
    vv[[2]] <- vv[[1]][-(k1)]
    q1 <- 1

    if (p > 1) {
        for (i in 3:(2 * p)) {
            if (i%%2 != 0) {
                vv[[i]] <- (q1 * k):(k * p - 1)
                q1 <- q1 + 1
            } else {

                vv[[i]] <- vv[[i - 1]][-(k1)]

            }

        }

    }

    return(vv)

}


# VARX Lag Group function
groupfunVARX <- function(p, k, k1, s) {

    jj <- list()
    m <- k - k1
    jj <- .groupfuncpp(p, k1)
    kp <- k1 * p + m * s - 1
    jj2 <- list()
    startjj <- max(unlist(jj)) + 1
    for (i in seq(startjj, kp, by = 1)) {
        jj[[i]] <- i
    }
    jj[sapply(jj, is.null)] <- NULL

    return(jj)
}

groupfunVARXcomp <- function(p, k, k1, s) {
    ownoth <- groupfunVARX(p, k, k1, s)
    kk2 <- list()
    pmax <- max(unlist(ownoth))
    to <- 0:(pmax)
    for (i in seq_len(length(ownoth))) {
        kk2[[i]] <- to[is.na(pmatch(to, ownoth[[i]]))]
    }
    return(kk2)

}


diaggroupfunVARX <- function(p, k, k1, s) {
    m <- k - k1
    jj <- list()
    jj <- .lfunction3cpp(p, k1)
    kp <- k1 * (p * k1 + s * m) - 1
    jj2 <- list()
    startjj <- max(unlist(jj)) + 1

    for (i in seq(startjj, kp, by = k1)) {
        jj[[i]] <- i:(i + k1 - 1)

    }

    jj[sapply(jj, is.null)] <- NULL


    return(jj)

}

diaggroupfunVARXcomp <- function(p, k, k1, s) {

    ownoth <- diaggroupfunVARX(p, k, k1, s)
    kk2 <- list()
    pmax <- max(unlist(ownoth))
    to <- 0:(pmax)

    for (i in seq_len(length(ownoth))) {
        kk2[[i]] <- to[is.na(pmatch(to, ownoth[[i]]))]
    }

    return(kk2)
}


diaggroupfunVARXL <- function(p, k, k1) {
    jj <- list()
    jj <- .lfunction3cpp(p, k1)
    kp <- k1 * p * k - 1
    jj2 <- list()
    startjj <- max(unlist(jj)) + 1

    for (i in seq(startjj, kp, by = 1)) {
        jj[[i]] <- i

    }

    jj[sapply(jj, is.null)] <- NULL
    return(jj)

}

diaggroupfunVARXcompL <- function(p, k, k1) {

    ownoth <- diaggroupfunVARXL(p, k, k1)
    kk2 <- list()
    pmax <- max(unlist(ownoth))

    to <- 0:(pmax)

    for (i in seq_len(length(ownoth))) {

        kk2[[i]] <- to[is.na(pmatch(to, ownoth[[i]]))]

    }

    return(kk2)

}


# iterative procedure to find a tighter bound for lambda starting value via binary search
LGSearch <- function(gstart, Y, Z, BOLD, group, k1, p, s, m, gs, k, MN, alpha, C, intercept, tol, VARX, gamma) {
    s1 <- 0
    palpha <- NULL
    tk <- 1/max(Mod(eigen(Z %*% t(Z))$values))
    lambdah <- gstart
    lambdal <- 0
    activeset <- list(rep(rep(list(0), length(gs))))

    gran2 <- 1
    grps <- create_group_indexes(group, p, k, gran2 * length(alpha), VARX, k1, s)
    groups <- grps$groups
    compgroups <- grps$compgroups
    activeset <- grps$activeset
    starting_eigvals <- grps$starting_eigvals
    nseries <- nrow(BOLD)
    while (max(abs(lambdah - lambdal)) > 10 * tol) {
        lambda <- (lambdah + lambdal)/2
        dual <- FALSE
        separate_lambdas <- FALSE
        temp <- .BigVAR.fit(group, BOLD, Z, Y, lambda, tol, p, m, k1, k, s, s1, MN, C, intercept, separate_lambdas, dual,
            activeset, starting_eigvals, groups, compgroups, VARX, alpha, palpha, gamma)
        # remove intercept from consideration
        if (group == "Tapered") {
            param <- adrop(BOLD[, -1, , drop = F], drop = 3)

        } else {

            BOLD <- temp$beta
            param <- adrop(BOLD[, -1, , drop = F], drop = 3)
        }

        activeset <- temp$activeset
        q1a <- temp$q1a

        if (MN) {
            submat <- param[1:nseries, 1:nseries, drop = F]
            diag(submat) <- ifelse(C == 0, diag(param[1:nseries, 1:nseries, drop = F]), diag(param[1:nseries, 1:nseries,
                drop = F]) - C)
            param[1:nseries, 1:nseries] <- submat
            BOLD[, -1, ] <- param

        }
        if (max(abs(param)) < tol) {
            lambdah <- lambda
        } else {
            lambdal <- lambda
        }

    }

    lambdah
}

#' Evaluate forecasts from a VAR or VARX with lag orders selected by AIC/BIC
#' 
#' @param Y a \eqn{T \times k} multivariate time series 
#' @param X a \eqn{T \times m} multivariate time series of unmodeled exogenous variables
#' @param p maximum lag order for endogenous series
#' @param s maximum lag order for exogenous series
#' @param T1 start of forecast evaluation period.
#' @param T2 end of forecast evaluation period
#' @param IC specifies whether to select lag order according to 'AIC' or 'BIC'
#' @param h desired forecast horizon
#' @param loss loss function (default 'L2', one of 'L1','L2','Huber')
#' @param delta delta for Huber loss function (default 2.5)
#' @param iterated indicator as to whether to use iterated or direct multistep forecasts (if applicable, VAR context only)
#' @return Returns the one-step ahead MSFE as well as the forecasts over the evaluation period and lag order selected.
#' @details This function evaluates the one-step ahead forecasts of a VAR or VARX fit by least squares over an evaluation period.  At every point in time, lag orders for the endogenous and exogenous series are selected according to AIC or BIC.  This function is run automatically when \code{\link{cv.BigVAR}} is called unless \code{ic} is set to \code{FALSE} in \code{\link{constructModel}}.      
#' @references Neumaier, Arnold, and Tapio Schneider. 'Estimation of parameters and eigenmodes of multivariate autoregressive models.' ACM Transactions on Mathematical Software (TOMS) 27.1 (2001): 27-57.
#' @seealso \code{\link{VARXFit}},\code{\link{constructModel}}, \code{\link{cv.BigVAR}}
#' @examples
#' data(Y)
#'
#' # Evaluate the performance of a VAR with lags selected by BIC.
#' p <- 4
#' T1 <- floor(nrow(Y))/3
#' T2 <- floor(2*nrow(Y))/3
#' # Matrix of zeros for X
#' X <- matrix(0,nrow=nrow(Y),ncol=ncol(Y))
#' BICMSFE <- VARXForecastEval(Y,X,p,0,T1,T2,'BIC',1)
#' 
#' @export
VARXForecastEval <- function(Y, X, p, s, T1, T2, IC, h, iterated = FALSE, loss = "L2", delta = 2.5) {


    if (T1 > nrow(Y) | T2 > nrow(Y) | T2 < T1) {
        stop("Training dates exceed series length")
    }

    if (!IC %in% c("AIC", "BIC")) {

        stop("IC must either be AIC or BIC")

    }

    MSFE <- c()
    predF <- NULL
    pvec <- NULL
    svec <- NULL
    k <- ncol(Y)
    m <- ifelse(s != 0, ncol(X), 0)
    for (i in (T1 - h + 2):T2) {

        if (h + i - 1 > T2) {
            break
        }

        testY <- as.matrix(Y[1:(i - 1), ])
        testX <- as.matrix(X[1:(i - 1), ])

        if (!iterated) {
            hd <- h

        } else {

            hd <- 1

        }
        if (IC == "BIC") {
            popt <- ICX(testY, testX, k, p, s, m, "BIC", h = hd)
        }
        if (IC == "AIC") {

            popt <- ICX(testY, testX, k, p, s, m, "AIC", h = hd)
        }
        B1 <- popt$B

        if (popt$p == 0 & popt$s == 0) {

            eZ <- matrix(rep(1, 1), ncol = 1)

            pred <- B1 %*% eZ

        } else {


            C <- max(popt$p, popt$s)

            ## # possibly memory leak in VARX lag matrix construction in Eigen if maxlag is 1.  # to be on the safe
            ## side, we will perform it in R
            eZ <- VARXCons(as.matrix(Y[(i - C):(i), ]), as.matrix(X[(i - C):(i), ]), k, popt$p, m, popt$s)

            ## }

            pred <- B1 %*% eZ

            # iterated multistep forecasts (if VAR and horizon greater than 1)
            if (h > 1 & s == 0 & iterated) {

                pred <- predictMS(matrix(pred, nrow = 1), Y, h - 1, B1, C, FALSE)

            }

        }
        predF <- rbind(predF, t(pred))
        MSFEi <- .calc.loss(Y[i + h - 1, ] - pred, univ = FALSE, loss, delta)

        MSFE <- c(MSFE, MSFEi)
        svec <- c(svec, popt$s)
        pvec <- c(pvec, popt$p)
    }

    return(list(MSFE = MSFE, pred = as.matrix(predF), p = pvec, s = svec))


}

#' Fit a VAR or VARX model by least squares
#' 
#' @param Y a \eqn{t \times k} multivariate time series
#' @param p maximum lag order
#' @param IC Information criterion indicator, if set to \code{NULL}, it will fit a least squares VAR(X) of orders p and s.  Otherwise, if set to 'AIC' or 'BIC' it return the model with lag orders that minimize the given IC. 
#' @param VARX a list of VARX specifications (as in \code{\link{constructModel}} (or NULL )
#' @return Returns a list with four entries:
#' \itemize{
#' \item{'Bhat'}{Estimated \eqn{k\times kp+ms} coefficient matrix}
#' \item{'SigmaU}{Estimated \eqn{k\times k} residual covariance matrix}
#' \item{'phat'}{Selected lag order for VAR component}
#' \item{'shat'}{Selected lag order for VARX component}
#' \item{'Y'}{multivariate time series retained for prediction purposes}
#' \item{'Y'}{number of endogenous (modeled) time series}
#' }
#' @details This function uses a modified form of the least squares technique proposed by Neumaier and Schneider (2001).  It fits a least squares VAR or VARX via a QR decomposition that does not require explicit matrix inversion.  This results in improved computational performance as well as numerical stability over the conventional least squares approach. 
#' @references Neumaier, Arnold, and Tapio Schneider. 'Estimation of parameters and eigenmodes of multivariate autoregressive models.' ACM Transactions on Mathematical Software (TOMS) 27.1 (2001): 27-57.
#' @seealso \code{\link{constructModel}}, \code{\link{cv.BigVAR}},\code{\link{BigVAR.fit}}
#' @examples
#' data(Y)
#' # fit a VAR_3(3)
#' mod <- VARXFit(Y,3,NULL,NULL)
#' # fit a VAR_3 with p= 6 and lag selected according to AIC
#' modAIC <- VARXFit(Y,6,'AIC',NULL)
#' # Fit a VARX_{2,1} with p=6, s=4 and lags selected by BIC
#' modXBIC <- VARXFit(Y,6,'BIC',list(k=1,s=4))
#' 
#' @export
VARXFit <- function(Y, p, IC, VARX = NULL) {

    if (!is.null(VARX)) {

        if (is.list(VARX) & !(exists("k", where = VARX) & exists("s", where = VARX))) {

            stop("VARX Specifications entered incorrectly")

        }

    }
    if (is.list(VARX) & (length(VARX) != 0)) {

        k1 <- VARX$k
        s <- VARX$s
        Y1 <- matrix(Y[, 1:k1], ncol = k1)
        m <- ncol(Y) - k1
        X <- matrix(Y[, (k1 + 1):ncol(Y)], ncol = m)

        if (exists("contemp", where = VARX)) {
            contemp <- VARX$contemp
        } else {
            contemp <- FALSE
        }
        Z <- VARXCons(Y1, X, k1, p, m, s, contemp = contemp)
        offset <- max(p, s) + 1
        YT <- matrix(Y1[offset:nrow(Y), ], ncol = k1)
        X <- matrix(X[offset:nrow(X), ], ncol = m)
    } else {

        k <- ncol(Y)
        k1 <- k
        s <- 0
        m <- 0
        offset <- p + 1
        X <- matrix(0, nrow = nrow(Y))

        Z <- VARXCons(Y, X, k, p, m, s)
        YT <- matrix(Y[(offset):nrow(Y), ], ncol = ncol(Y))

    }
    if (is.null(IC)) {

        Res <- ARFitVARXR(cbind(t(Z), YT), k1, p, m, s)

        shat <- s
        phat <- p

    } else {
        if (!IC %in% c("AIC", "BIC")) {

            stop("IC must either be AIC,BIC, or set to NULL")

        }

        Res <- ICX(YT, X, k1, p, s, m, IC)

        shat <- Res$s
        phat <- Res$p

    }
    if (is.null(VARX)) {
        k <- ncol(Y)
    } else {
        k <- VARX$k
    }
    list(Bhat = Res$B, SigmaU = Res$SigmaU, phat = phat, shat = shat, Y = Y, k = k)

}

#' One-step ahead predictions for VARX models
#' 
#' @param VARXRes the results from \code{\link{VARXFit}}
#' @return Returns a vector consisting of the out-of-sample forecasts for the provided \code{\link{VARXFit}} model.
#' @seealso \code{\link{VARXFit}}
#' @examples
#' data(Y)
#' # fit a VAR_3(3)
#' mod <- VARXFit(Y,3,NULL,NULL)
#' pred <-PredictVARX(mod)
#' 
#' @export
PredictVARX <- function(VARXRes) {

    B <- VARXRes$Bhat
    Y <- VARXRes$Y
    k <- VARXRes$k
    m <- ncol(Y) - k

    if (k < ncol(Y)) {
        Z <- VARXCons(Y[, 1:k, drop = FALSE], Y[, (k + 1):ncol(Y), drop = FALSE], k, VARXRes$phat, m, VARXRes$shat, oos = TRUE)
    } else {
        Z <- VARXCons(Y[, 1:k, drop = FALSE], matrix(0, nrow = nrow(Y)), k, VARXRes$phat, m, VARXRes$shat, oos = TRUE)
    }

    return(as.numeric(tail(t(B %*% Z), 1)))


}

# Recursive multi-step predictions
predictMS <- function(pred, Y, n.ahead, B, p, MN = FALSE, predict_all = FALSE, n.ahead_full = n.ahead) {

    # Augment Y with predictions, create lag matrix (no intercept if MN)
    Y <- rbind(Y, pred)

    Z <- VARXCons(Y, matrix(0, nrow = nrow(Y), ncol = 1), ncol(Y), p, 0, 0, oos = TRUE)

    if (MN) {
        Z <- Z[2:nrow(Z), , drop = F]
    }
    Z <- Z[, ncol(Z), drop = F]
    pred <- matrix(B %*% Z, ncol = ncol(Y), nrow = 1)

    if (n.ahead == 1) {
        if (predict_all) {
            return(rbind(Y[((nrow(Y) - n.ahead_full) + 1):nrow(Y), ], pred))
        } else {
            return(pred)
        }
    }

    predictMS(pred, Y, n.ahead - 1, B, p, MN, predict_all, n.ahead_full = n.ahead)

}

# Multi-step VARX with new data.
predictMSX <- function(pred, Y, n.ahead, B, p, newxreg, X, m, s, cumulative, MN, contemp = FALSE) {

    Y <- rbind(Y, pred)
    X <- rbind(X, matrix(newxreg[cumulative, ], ncol = m))

    if (nrow(Y) != nrow(X)) {
        stop("error, dimension issue")
    }
    if (!contemp) {
        Z <- VARXCons(as.matrix(Y), X, ncol(Y), p, m, s, oos = TRUE)
    } else {
        Z <- VARXCons(as.matrix(Y), as.matrix(X), ncol(Y), p, m, s, oos = FALSE, contemp = TRUE)
    }
    Z <- Z[, ncol(Z), drop = F]
    if (MN) {

        Z <- as.matrix(Z[2:nrow(Z), drop = F])
        pred <- matrix(B[, 2:ncol(B), drop = F] %*% Z, ncol = ncol(Y), nrow = 1)

    } else {

        pred <- matrix(B %*% Z, ncol = ncol(Y), nrow = 1)
    }
    if (n.ahead == 1) {
        return(pred)
    }

    predictMSX(pred, Y, n.ahead - 1, B, p, newxreg, X, m, s, cumulative + 1, MN)

}




# Find optimal values in 2-d gridsearch
findind <- function(opt, lambda1, lambda2) {
    if (opt < length(lambda2)) {
        lambda1ind <- 1
    } else {
        lambda1ind <- ceiling(opt/length(lambda2))
    }
    if (lambda1ind == 1) {
        jind <- opt
    } else {
        jind <- opt - (length(lambda2)) * (lambda1ind - 1)
    }
    return(c(lambda1ind, jind))
}



# Bayesian VAR with MN prior
BVARLitterman <- function(Y, Z, p, tau, mu, H, iRW) {
    T <- nrow(Y)
    k <- ncol(Y)

    # prior covariance based on univariate AR models
    sigmas <- c()
    for (i in 1:k) {
        Z1 <- VARXCons(Y[, i, drop = F], matrix(0, nrow = nrow(Y), ncol = 1), 1, p, 0, 0)
        K <- cbind(t(Z1), Y[(p + 1):nrow(Y), i])
        sigmas[i] <- sqrt(ARFitVARXR(K, 1, p, 0, 0)$SigmaU)

    }
    MMO <- colMeans(Y)

    # create prior random walk dummy
    Yrw1 <- diag(sigmas * iRW)
    Yrw2 <- matrix(0, nrow = k * (p - 1), ncol = k)
    Yrw <- tau * (rbind(Yrw1, Yrw2))
    Zrw <- tau * cbind(kronecker(diag(1:p), diag(sigmas)), matrix(0, nrow = k * p, ncol = 1))


    # create dummy for intercept
    epsilon <- 1e-05
    Ycs <- 1e-05 * matrix(0, nrow = 1, ncol = k)
    Zcs <- epsilon * cbind(matrix(0, ncol = k * p, nrow = 1), 1)

    # dummy on the sums of coefficients
    Ylr <- mu * diag(MMO * iRW)
    Zlr1 <- kronecker(matrix(1, nrow = 1, ncol = p), diag(MMO) * iRW)
    Zlr <- mu * (cbind(Zlr1, matrix(0, nrow = k, ncol = 1)))


    # Dummy for residual covariance matrix
    Ycv <- diag(sigmas)
    Zcv <- matrix(0, nrow = k, ncol = k * p + 1)

    Yprior <- rbind(Yrw, Ylr, Ycv, Ycs)
    Zprior <- rbind(Zrw, Zlr, Zcv, Zcs)

    Tstar <- nrow(Yprior)
    Z <- t(Z)
    Z <- cbind(Z[, 2:ncol(Z)], 1)

    ZZinv <- solve(t(Zprior) %*% Zprior + t(Z) %*% Z)
    ZY <- t(Zprior) %*% Yprior + t(Z) %*% Y
    beta <- ZZinv %*% ZY

    return(t(beta))

}

# grid search for BGR
BGRGridSearch <- function(Y, Z, p, grid, RWIND) {
    preds <- list()
    for (i in seq_len(length(grid))) {
        pi <- grid[i]
        mu <- pi * 0.1  # used in BGR paper
        preds[[i]] <- BVARLitterman(Y, Z, p, pi, mu, -1, RWIND)

    }
    preds <- array(unlist(preds), dim = c(nrow(preds[[1]]), ncol(preds[[1]]), length(preds)))

    return(preds)

}

# process MN prior
MN_prior <- function(Y, Z, C1) {
    C <- matrix(0, nrow = ncol(Y), ncol = nrow(Z))
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y - C %*% Z
    Y <- t(Y)
    return(list(Y = Y, C = C))
}

# pre-process data (subtract intercept, adjust for shrinking toward constants)
pre_process <- function(Y, Z, C1, MN, intercept) {
    k <- ncol(Y)
    if (MN) {
        YC <- MN_prior(Y, Z, C1)
        Y <- YC$Y
        C <- YC$C
    }

    Y <- t(Y)

    if (intercept) {
        YMean <- c(apply(Y, 1, mean))
        ZMean <- c(apply(Z, 1, mean))
        if (k > 1) {
            Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
        } else {
            Y <- Y - mean(Y)
        }
        Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))

    } else {
        YMean <- rep(0, nrow(Y))
        ZMean <- rep(0, nrow(Z))
    }
    Y <- t(Y)
    return(list(Y = Y, Z = Z, YMean = YMean, ZMean = ZMean, C = C))

}


# s needs to be s+ s1 create indices for group structures
create_group_indexes <- function(group, p, k, gran2, VARX = FALSE, k1 = NULL, s = NULL) {
    starting_eigvals <- NULL
    groups <- NULL
    compgroups <- NULL
    activeset <- NULL
    if (VARX) {
        if (group == "Lag") {

            groups <- groupfunVARX(p, k, k1, s)
            compgroups <- groupfunVARXcomp(p, k, k1, s)
            activeset <- rep(list(rep(rep(list(0), length(groups)))), gran2)

        } else if (group == "SparseLag") {

            groups <- groupfunVARX(p, k, k1, s)
            compgroups <- groupfunVARXcomp(p, k, k1, s)
            starting_eigvals <- list()

            for (i in 1:(p + s)) {

                starting_eigvals[[i]] <- matrix(runif(k1, -1, 1), ncol = 1)

            }

            activeset <- rep(list(rep(rep(list(0), length(groups)))), gran2)


        } else if (group == "OwnOther") {
            groups <- diaggroupfunVARX(p, k, k1, s)
            compgroups <- diaggroupfunVARXcomp(p, k, k1, s)
            activeset <- rep(list(rep(rep(list(0), length(groups)))), gran2)


        } else if (group == "SparseOO") {
            groups <- diaggroupfunVARX(p, k, k1, s)
            compgroups <- diaggroupfunVARXcomp(p, k, k1, s)
            activeset <- rep(list(rep(rep(list(0), length(groups)))), gran2)

        }


    } else {
        if (group == "Lag") {

            groups <- .groupfuncpp(p, k)


            compgroups <- .groupfuncomp(p, k)


            activeset <- rep(list(rep(rep(list(0), length(groups)))), gran2)

        } else if (group == "SparseLag") {

            groups <- .groupfuncpp(p, k)
            compgroups <- .groupfuncomp(p, k)

            starting_eigvals <- list()

            for (i in 1:p) {

                starting_eigvals[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
            }


            activeset <- rep(list(rep(rep(list(0), length(groups)))), gran2)



        } else if (group == "OwnOther") {

            groups <- .lfunction3cpp(p, k)
            compgroups <- .lfunctioncomp(p, k)
            activeset <- rep(list(rep(rep(list(0), length(groups)))), gran2)

        } else if (group == "SparseOO") {

            groups <- .lfunction3cpp(p, k)
            compgroups <- .lfunctioncomp(p, k)
            activeset <- rep(list(rep(rep(list(0), length(groups)))), gran2)
            starting_eigvals <- list()

            for (i in 1:(2 * p)) {

                starting_eigvals[[i]] <- matrix(runif(length(groups[[i]]), -1, 1), ncol = 1)
            }

        }
    }
    return(list(groups = groups, compgroups = compgroups, starting_eigvals = starting_eigvals, activeset = activeset))
}


create_lambda_grid <- function(trainY, trainZ, groups, gran1, gran2, group, p, k1, s, m, k, MN, alpha, C, intercept, tol,
    VARX, separate_lambdas, dual, gamma, linear, verbose) {
    # Constructs penalty grid if both alpha and lambda are selected
    if (dual) {
        lambda <- matrix(0, nrow = gran2, ncol = length(alpha))

        for (i in 1:length(alpha)) {

            lambda[, i] <- .LambdaGrid(gran1, gran2, groups, trainY, trainZ, group, p, k1, s, m, k, MN, alpha[i], C, intercept,
                tol, VARX = VARX, linear = FALSE)

        }

    } else {
        # Penalty parameter grid for just lambda
        if (group != "BGR") {

            lambda <- .LambdaGrid(gran1, gran2, groups, trainY, trainZ, group, p, k1, s, m, k, MN, alpha, C, intercept, tol,
                VARX = VARX, separate_lambdas, verbose, linear = FALSE)
        } else {
            # special handling for BGR
            lambda <- seq(1, 5, length = gran2)
            lambda <- lambda * sqrt(k * p)

        }
    }


    return(lambda)

}

# runs at each iteration of penalty parameter selection and validation stage
refine_and_forecast <- function(betaArray, eZ, trainZ, trainY, testY, lambda, h, recursive, MN, RVAR, refit_fraction, separate_lambdas,
    C, inds, loss, delta, k, p, k1, s, oos = FALSE) {
    if (separate_lambdas) {
        MSFE_temp <- matrix(0, nrow = dim(betaArray)[3], ncol = ncol(trainY))
    } else {
        MSFE_temp <- c()
    }
    if (MN) {

        eZ <- eZ[2:nrow(eZ), , drop = FALSE]

    }

    nlambdas <- dim(betaArray)[3]
    preds <- matrix(0, nrow = ncol(trainY), ncol = nlambdas)
    for (i in 1:nlambdas) {

        beta <- adrop(betaArray[, , i, drop = F], drop = 3)
        if (RVAR) {

            beta_rls <- RelaxedLS(cbind(t(trainZ), trainY), beta)
            beta <- (1 - refit_fraction) * beta + refit_fraction * beta_rls

        }
        if (!is.null(inds)) {
            if (i == inds) {
                beta_return <- beta
            }
        } else {
            # if we don't need beta matrix (i.e validation stage, just return the last value)
            beta_return <- beta
        }

        if (MN) {
            beta <- beta[, 2:ncol(beta)]
        }

        if (h > 1 & recursive) {
            ptemp <- beta %*% eZ

            pred <- matrix(ptemp, nrow = 1)

            pred <- predictMS(pred, trainY, h - 1, beta, p, MN)
            pred <- matrix(pred, ncol = 1)
        } else {
            pred <- beta %*% eZ
        }
        if (separate_lambdas) {
            if (!oos) {
                for (j in seq_len(ncol(testY))) {
                  MSFE_temp[i, j] <- .calc.loss(testY[, j] - pred[j, ], univ = TRUE, loss, delta)
                }
            } else {
                MSFE_temp[i] <- .calc.loss(t(testY) - beta %*% eZ, univ = FALSE, loss, delta)
            }
        } else {
            MSFE_temp[i] <- .calc.loss(t(testY) - beta %*% eZ, univ = FALSE, loss, delta)
        }


        preds[, i] <- pred
    }

    ## # just return optimal prediction in validation stage
    return(list(MSFE = MSFE_temp, betaArray = beta_return, preds = preds))
}

# determines optimal lambda from a grid of values
find_optimal_lambda <- function(MSFE, lambda, palpha, alpha, ONESE, dual, separate_lambdas, tapered) {

    if (tapered) {

        indopt <- which.min(colMeans(MSFE))

        if (indopt < length(lambda)) {

            alphaind <- 1

            alphaopt <- palpha[1]

        } else {

            alphaopt <- palpha[floor(indopt/length(lambda))]

            alphaind <- floor(indopt/length(lambda))

        }

        if (alphaind == 1) {

            lambdaopt <- lambda[indopt]

        } else if (indopt%%length(lambda) == 0) {

            lambdaopt <- lambda[length(lambda)]

        } else {
            lambdaind <- indopt - length(lambda) * alphaind
            lambdaopt <- lambda[lambdaind]

        }

        palpha <- alphaopt
        optind <- indopt
    } else {
        # one standard error correction
        if (ONESE & !dual & !separate_lambdas) {

            MSFE2 <- MSFE
            G2 <- colMeans(na.omit(MSFE2))
            G3 <- sd(na.omit(MSFE2))/sqrt(nrow(na.omit(MSFE2)))
            optind <- min(which(G2 < (min(G2) + G3)))
            lambdaopt <- lambda[optind]
        } else {

            if (!tapered & !dual) {
                # in rare cases in which MSFE is equal, the smaller penalty parameter is chosen.  This prevents
                # extremely sparse solutions
                if (separate_lambdas) {
                  if (ONESE) {
                    MSFES <- t(apply(MSFE, 3, colMeans))
                    sds <- t(apply(MSFE, 3, function(x) sd(na.omit(x))/sqrt(nrow(na.omit(x)))))
                    lambdaopt <- c()
                    optinds <- c()
                    for (i in 1:nrow(MSFES)) {
                      optinds[i] <- min(which(MSFES[i, ] < sds[i] + min(MSFES[i, ])))
                      lambdaopt[i] <- lambda[optinds[i], i, drop = F]
                    }
                    optind = optinds
                  } else {
                    MSFES <- t(apply(MSFE, 3, colMeans))
                    optinds <- apply(MSFES, 1, which.min)
                    lambdaopt <- c()
                    for (i in 1:nrow(MSFES)) {
                      lambdaopt[i] <- lambda[optinds[i], i]
                    }
                    optind = optinds


                  }
                } else {

                  optind <- max(which(colMeans(na.omit(MSFE)) == min(colMeans(na.omit(MSFE)))))
                  lambdaopt <- lambda[optind]
                }
            } else if (dual) {
                if (!ONESE) {

                  optind <- max(which(colMeans(na.omit(MSFE)) == min(colMeans(na.omit(MSFE)))))
                  inds <- findind(optind, lambda[, 1], alpha)
                } else {
                  G2 <- colMeans(na.omit(MSFE))

                  G3 <- sd(na.omit(MSFE))/sqrt(nrow(na.omit(MSFE)))

                  optind <- min(which(G2 < (min(G2) + G3)))
                  inds <- findind(optind, lambda[, 1], alpha)

                }
                lambdaopt <- lambda[inds[1], inds[2]]
                lambda <- lambda[, inds[2]]
                alphaopt <- alpha[inds[2]]
                optind <- inds

            }
        }
        if (!dual) {
            alphaopt <- alpha
        }

    }
    if (!exists("alphaopt")) {
        alphaopt <- NULL
    }
    return(list(optind = optind, lambdaopt = lambdaopt, alphaopt = alphaopt))
}



create_sigma_u <- function(Y, Z, B, h) {
    trainY <- Y[(h):(nrow(Y) - 1), , drop = F]
    trainZ <- Z[, 1:(ncol(Z) - h)]
    Sigma = 1/nrow(trainY) * sqrt(diag(tcrossprod(t(Y) - B %*% rbind(1, Z))))
    cis <- c(qnorm(0.025), qnorm(0.975))
    lower <- cis[1] * Sigma
    upper <- cis[2] * Sigma
    return(cbind(lower, upper))

}

# add C back to coefficient matrix in MN setting
adjust_mn_var <- function(beta, C) {
    for (i in 1:(dim(beta)[3])) beta[, 2:dim(beta)[2], i] <- beta[, 2:dim(beta)[2], i] + C
    return(beta)
}

.MCPFit <- function(B, Z, Y, lambda, eps, p, MN, k, k1, s, m, C,  group, gamma = 3, YMean, ZMean) {


    tk <- 1/max(Mod(eigen(Z %*% t(Z), only.values = TRUE)$values))

    B1 <- abind::adrop(B[, 2:dim(B)[2], 1, drop = F], 3)
    nc <- apply(B, 3, ncol)[1]
    BINI <- B[, 2:nc, , drop = F]
    if (group == "MCP") {
        mcp = TRUE
    } else {
        mcp = FALSE
    }

    beta <- gamloopMCP(BINI, Y, Z, as.matrix(lambda), eps, as.matrix(YMean), as.matrix(ZMean), gamma = gamma, mcp)

    if (MN) {
        beta <- adjust_mn_var(beta, C)
    }

    return(beta)

}


                                        # Sparse Own/Other (VAR)
.SparseGroupLassoVAROO <- function(beta, groups, compgroups, Y, Z, lambda, alpha, INIactive, eps, q1a, p, MN, dual = FALSE,
                                   C, YMean, ZMean) {

    k <- ncol(Y)
    Y <- t(Y)
    m <- 0
    ZZ <- kronecker(t(Z), diag(k))
    M1f <- list()
    M2f <- list()
    jj <- .lfunction3(p, k)
    eigs <- c()
    q1 <- list()
                                        # get step size from inverse of max eigenvalue via power method
    for (j in seq_len(length(jj))) {
        M1f[[j]] <- ZZ[, jj[[j]]]
        M2f[[j]] <- crossprod(ZZ[, jj[[j]]])
        gg1 <- powermethod(M2f[[j]], q1a[[j]])
        eigs[j] <- gg1$lambda
        q1[[j]] <- gg1$q1
    }
    jj <- .lfunction3cpp(p, k)
    jjfull <- jj
    jjcomp <- .lfunctioncomp(p, k)
    dims <- dim(beta)
    beta <- array(beta[, 2:ncol(beta[, , 1]), ], dim = c(dims[1], dims[2] - 1, dims[3]))
    if (!dual) {

        BB <- GamLoopSGLOO(beta, INIactive, lambda, alpha, Y, ZZ, jj, jj, jjcomp, eps, YMean, ZMean, k, p * k, M2f, eigs,
                           m)


    } else {
        BB <- GamLoopSGLOODP(beta, INIactive, lambda, alpha, Y, ZZ, jj, jj, jjcomp, eps, YMean, ZMean, k, p * k, M2f, eigs,
                             m)

    }

    BB$q1 <- q1

    if (MN) {
        BB$beta <- adjust_mn_var(BB$beta, C)
    }
    return(BB)
}



                                        # Sparse Lag (VAR)
.SparseGroupLassoVAR <- function(beta, groups, compgroups, Y, Z, lambda, alpha, INIactive, eps, q1a, p, MN, C, YMean, ZMean) {
    k <- ncol(Y)

    Y <- t(Y)

    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()
    jj <- .groupfun(p, k)

                                        # get step size from inverse of max eigenvalue via power method
    for (j in seq_len(length(jj))) {
        M1f[[j]] <- Z[jj[[j]], ]
        M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])
        gg1 <- powermethod(M2f[[j]], q1a[[j]])
        eigs[j] <- gg1$lambda
        q1[[j]] <- gg1$q1
    }
    jj <- .groupfuncpp(p, k)
    jjfull <- jj
    jjcomp <- .groupfuncomp(p, k)
    dims <- dim(beta)
    beta <- array(beta[, 2:ncol(beta[, , 1]), ], dim = c(dims[1], dims[2] - 1, dims[3]))
    BB <- GamLoopSGL(beta, INIactive, lambda, alpha, Y, Z, jj, jjfull, jjcomp, eps, YMean, as.matrix(ZMean), k, p * k, M1f,
                     M2f, eigs)
    BB$q1 <- q1

    if (MN) {
        BB$beta <- adjust_mn_var(BB$beta, C)

    }

    return(BB)
}

                                        # Sparse Lag (VAR) Dual Search
.SparseGroupLassoVARDual <- function(beta, groups, compgroups, Y, Z, lambda, alpha, INIactive, eps, q1a, p, MN, C, YMean,
                                     ZMean) {
    k <- ncol(Y)
    Y <- t(Y)
    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()
    jj <- .groupfun(p, k)
                                        # get step size from inverse of max eigenvalue via power method
    for (j in seq_len(length(jj))) {
        M1f[[j]] <- Z[jj[[j]], ]
        M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])
        gg1 <- powermethod(M2f[[j]], q1a[[j]])
        eigs[j] <- gg1$lambda
        q1[[j]] <- gg1$q1
    }

    jj <- .groupfuncpp(p, k)
    jjfull <- jj
    jjcomp <- .groupfuncomp(p, k)
    ngp <- length(alpha) * length(lambda)
    dims <- dim(beta)
    beta <- array(beta[, 2:ncol(beta[, , 1]), ], dim = c(dims[1], dims[2] - 1, dims[3]))
    BB <- GamLoopSGLDP(beta, INIactive, lambda, alpha, Y, Z, jj, jjfull, jjcomp, eps, YMean, ZMean, k, p * k, M1f, M2f, eigs)
    BB$q1 <- q1

    if (MN) {
        BB$beta <- adjust_mn_var(BB$beta, C)

    }

    return(BB)
}


                                        # Lag Group (VAR/VARX-L)
.GroupLassoVAR1 <- function(beta, groups, compgroups, Y, Z, lambda, INIactive, eps, p, MN, k, k1, s, C, YMean, ZMean) {

    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }
    k1 <- ncol(Y)

    m <- k - k1

    Y <- t(Y)
    fullgroups <- groups
    Eigsys <- Eigencomp(Z, groups, length(groups), k1)
    M2 <- Eigsys$M3
    eigvals <- Eigsys$eigval
    eigvecs <- Eigsys$eigvec
    beta <- beta[, 2:dim(beta)[2], , drop = F]
    BB <- GamLoopGL2(beta, INIactive, lambda, Y, Z, groups, fullgroups, compgroups, eps, YMean, ZMean, k1, p * k1 + m * s,
                     M2, eigvals, eigvecs)

    if (MN) {
        BB$beta <- adjust_mn_var(BB$beta, C)

    }
    return(BB)
}


                                        # Do I really need a separate function for X vs not?  my guess is no.... should be the case for all of these Group
                                        # Lasso Own/Other (VARXL)
.GroupLassoOOX <- function(beta, groups, compgroups, Y, Z, lambda, INIactive, eps, p, MN, k, k1, s, C, YMean, ZMean) {
    m <- k - k1


    Y <- t(Y)

    ZZ <- kronecker(t(Z), diag(k1))

    Eigsys <- EigencompOO(ZZ, groups, length(groups), k1)
    M2 <- Eigsys$M3
    eigvals <- Eigsys$eigval
    eigvecs <- Eigsys$eigvec

    groups_full <- groups

    beta <- array(beta[, 2:ncol(as.matrix(beta[, , 1])), ], dim = c(k1, (k1) * p + m * s, length(lambda)))

    BB <- GamLoopGLOO(beta, INIactive, lambda, Y, ZZ, groups, groups_full, compgroups, eps, YMean, ZMean, k1, p * (k1) +
                                                                                                              m * s, M2, eigvals, eigvecs, k1)

    if (MN) {
        BB$beta <- adjust_mn_var(BB$beta, C)
    }

    return(BB)
}


                                        # Own/Other Group VAR-L
.GroupLassoOO <- function(beta, groups, compgroups, Y, Z, lambda, INIactive, eps, p, MN, C, YMean, ZMean) {

    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }

    k <- ncol(Y)

    Y <- t(Y)

    fullgroups <- groups

    ZZ <- kronecker(t(Z), diag(k))

    Eigsys <- EigencompOO(ZZ, groups, length(groups), k)
    M2 <- Eigsys$M3
    eigvals <- Eigsys$eigval
    eigvecs <- Eigsys$eigvec

    beta <- array(beta[, 2:ncol(as.matrix(beta[, , 1])), ], dim = c(k, k * p, length(lambda)))


    BB <- GamLoopGLOO(beta, INIactive, lambda, Y, ZZ, groups, fullgroups, compgroups, eps, YMean, ZMean, k, p * k, M2, eigvals,
                      eigvecs, k)

    if (MN) {
        BB$beta <- adjust_mn_var(BB$beta, C)

    }

    return(BB)
}

                                        # Sparse Lag VARX-L
.SparseGroupLassoVARX <- function(beta, groups, compgroups, Y, Z, lambda, alpha, INIactive, eps, starting_eigvals, p, MN,
                                  k, s, k1, C, YMean, ZMean) {

    m <- k - k1

    Y <- t(Y)
    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()
    grps_lg <- lapply(groups, function(x) {
        x + 1
    })
    for (j in seq_len(length(grps_lg))) {
        M2f[[j]] <- Z[grps_lg[[j]], ] %*% t(Z[grps_lg[[j]], ])

        if (j <= p) {

            gg1 <- powermethod(M2f[[j]], starting_eigvals[[j]])
            eigs[j] <- gg1$lambda
            q1[[j]] <- gg1$q1
        } else {

            M2f[[j]] <- as.vector(Z[grps_lg[[j]], ]) %*% as.vector(t(Z[grps_lg[[j]], ]))
            eigs[j] <- M2f[[j]]

        }

    }

    groups_full <- groups

    beta <- array(beta[, 2:ncol(as.matrix(beta[, , 1])), ], dim = c(k1, (k1) * p + (s * m), length(lambda)))
    BB <- GamLoopSGLX(beta, INIactive, lambda, alpha, Y, Z, groups, groups_full, compgroups, eps, YMean, ZMean, k1, (k1) *
                                                                                                                    p + (s * m), M2f, eigs, k1)

    BB$q1 <- q1

    if (MN) {

        BB$beta <- adjust_mn_var(BB$beta, C)

    }

    return(BB)
}



.SparseGroupLassoVARXDual <- function(beta, groups, compgroups, Y, Z, lambda, alpha, INIactive, eps, starting_eigvals, p,
                                      MN, k, s, k1, C, YMean, ZMean) {

    m <- k - k1

    Y <- t(Y)

    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()
    jj <- lapply(groups, function(x) {
        x + 1
    })
    for (j in seq_len(length(jj))) {

        M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])

        if (j <= p) {

            gg1 <- powermethod(M2f[[j]], starting_eigvals[[j]])
            eigs[j] <- gg1$lambda
            q1[[j]] <- gg1$q1
        } else {

            M2f[[j]] <- as.vector(Z[jj[[j]], ]) %*% as.vector(t(Z[jj[[j]], ]))
            eigs[j] <- M2f[[j]]

        }

    }



    beta <- array(beta[, 2:ncol(as.matrix(beta[, , 1])), ], dim = c(k1, (k1) * p + (s * m), nrow(lambda) * length(alpha)))

    BB <- GamLoopSGLXDP(beta, INIactive, lambda, alpha, Y, Z, groups, groups, compgroups, eps, YMean, ZMean, k1, (k1) * p +
                                                                                                                 (s * m), M2f, eigs, k1)

    BB$q1 <- q1

    if (MN) {

        BB$beta <- adjust_mn_var(BB$beta, C)

    }

    return(BB)
}


                                        # Sparse Own/Other (VARX)
.SparseGroupLassoVAROOX <- function(beta, groups, compgroups, Y, Z, lambda, alpha, INIactive, eps, p, MN, k1, s, k, dual = FALSE,
                                    C, YMean, ZMean) {

    m <- k - k1

    Y <- t(Y)
    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()

                                        # function for R calculations jj <- diaggroupfunVARXLG(p, k, k1, s)
    rgroups <- lapply(groups, function(x) {
        x + 1
    })
                                        # for c++ calculations kk <- diaggroupfunVARX(p, k, k1, s) jjcomp <- diaggroupfunVARXcomp(p, k, k1, s)

    ZZ <- kronecker(t(Z), diag(k1))


    for (j in seq_len(length(rgroups))) {

        M2f[[j]] <- crossprod(ZZ[, rgroups[[j]]])

        eigs[j] <- max(Mod(eigen(M2f[[j]], only.values = TRUE)$values))
    }

    fullgroups <- groups

    if (!dual) {

        gran2 <- length(lambda)

    } else {

        gran2 <- nrow(lambda) * ncol(lambda)

    }

    beta <- array(beta[, 2:ncol(as.matrix(beta[, , 1])), ], dim = c(k1, k1 * p + m * s, gran2))
    if (dual) {
        BB <- GamLoopSGLOODP(beta, INIactive, lambda, alpha, Y, ZZ, groups, groups, compgroups, eps, YMean, ZMean, k1, p *
                                                                                                                       k1 + m * s, M2f, eigs, m)
    } else {
        BB <- GamLoopSGLOO(beta, INIactive, lambda, alpha, Y, ZZ, groups, fullgroups, compgroups, eps, YMean, ZMean, k1,
                           p * k1 + m * s, M2f, eigs, m)
    }
    if (MN) {
        BB$beta <- adjust_mn_var(BB$beta, C)

    }


    return(BB)
}


                                        # Elementwise HLAG
.HLAGElemAlg <- function(beta, Y, Z, lambda, eps, p, MN, C, YMean, ZMean, separate_lambdas = FALSE) {

    k <- ncol(Y)

    betafin <- beta

    tk <- 1/max(Mod(eigen(Z %*% t(Z))$values))
    lambda <- as.matrix(lambda)
    betaini <- array(beta[, 2:ncol(beta[, , 1]), ], dim = c(k, k * p, nrow(lambda)))
    betafin <- gamloopElem(betaini, Y, Z, lambda, eps, YMean, ZMean, as.matrix(betaini[, , 1]), k, p, separate_lambdas)

    if (MN) {
        betafin <- adjust_mn_var(betafin, C)
    }

    return(betafin)
}

.lassoVARFistX <- function(B, Z, Y, lambda, eps, p, MN, k, k1, s, m, C, YMean, ZMean, separate_lambdas = FALSE) {
    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }

    tk <- 1/max(Mod(eigen(Z %*% t(Z), only.values = TRUE)$values))

    B1 <- abind::adrop(B[, 2:dim(B)[2], 1, drop = F], 3)


    nc <- apply(B, 3, ncol)[1]

    BINI <- B[, 2:nc, , drop = F]
    beta <- gamloopFista(BINI, Y, Z, as.matrix(lambda), eps, as.matrix(YMean), as.matrix(ZMean), B1, k, p, tk, k1, s, separate_lambdas)

    if (MN) {
        beta <- adjust_mn_var(beta, C)
    }
    return(beta)

}

                                        # general basic/elastic net
.lassoVARFistXEN <- function(B, Z, Y, lambda, alpha, eps, p, MN, k, k1, s, m, C, YMean, ZMean, separate_lambdas = FALSE) {

    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }


    tk <- 1/max(Mod(eigen(Z %*% t(Z), only.values = TRUE)$values))
    BFOO1 <- abind::adrop(B[, 2:dim(B)[2], 1, drop = F], 3)

    nc <- apply(B, 3, ncol)[1]
    BFOO <- B[, 2:nc, , drop = F]
    ## if (length(alpha) == 1) { alpha <- rep(alpha, dim(B)[3]) }
    beta <- gamloopFistaEN(BFOO, Y, Z, as.matrix(lambda), as.matrix(alpha), eps, as.matrix(YMean), as.matrix(ZMean), BFOO1,
                           k, p, tk, k1, s, separate_lambdas)

    if (MN) {
        beta <- adjust_mn_var(beta, C)
    }

    return(beta)

}


                                        # Componentwise HLAG
.HLAGCAlg <- function(beta, Y, Z, lambda, eps, p, MN, C, YMean, ZMean, separate_lambdas = FALSE) {
    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }

    k <- ncol(Y)


    betafin <- beta

    tk <- 1/max(Mod(eigen(Z %*% t(Z))$values))

    lambda <- as.matrix(lambda)
    if (separate_lambdas) {
        betaini <- array(beta[, 2:ncol(as.matrix(beta[, , 1])), ], dim = c(k, k * p, nrow(lambda)))
    } else {
        betaini <- array(beta[, 2:ncol(as.matrix(beta[, , 1])), ], dim = c(k, k * p, length(lambda)))
    }
    betafin <- gamloopHLAG(betaini, Y, Z, lambda, eps, YMean, ZMean, as.matrix(betaini[, , 1]), k, p, separate_lambdas)

    if (MN) {
        betafin <- adjust_mn_var(betafin, C)
    }
    return(betafin)
}

                                        # Endogenous First VARX-L
.EFVARX <- function(beta, Y, Z, lambda, eps, MN, k1, s, m, p, C, YMean, ZMean) {



    betafin <- beta

    tk <- 1/max(Mod(eigen(Z %*% t(Z))$values))


    if (k1 == 1) {

        betaini <- array(beta[, 2:(k1 * p + m * s + 1), ], dim = c(1, k1 * p + m * s, dim(beta)[3]))

    } else {

        betaini <- beta[, 2:ncol(beta[, , 1]), ]

    }

    for (i in seq_len(length(lambda))) {

        if (dim(beta)[3] > 1) {

            betaF <- fistaX(Y, Z, matrix(betaini[, , i], nrow = k1), p, k1, lambda[i], eps, tk, m, s)
        } else {

            betaF <- fistaX(Y, Z, matrix(betaini, nrow = k1), p, k1, lambda[i], eps, tk, m, s)

        }

        nu <- YMean - betaF %*% ZMean

        betafin[, , i] <- cbind(nu, betaF)

    }

    if (MN) {

        betafin <- adjust_mn_var(betafin, C)

    }

    return(betafin)
}

                                        # HLAG Own/Other
.HLAGOOAlg <- function(beta, Y, Z, lambda, eps, p, MN, C, YMean, ZMean, separate_lambdas = FALSE) {

    k <- ncol(Y)
    betafin <- beta
    YOLD <- Y
    ZOLD <- Z


    weights <- sqrt(c(rep(c(1, k - 1), length = 2 * p)))
    groups <- list()
    for (i in 1:k) {

        groups[[i]] <- .vecoovarscpp(p, k, i)

    }

    lambda <- as.matrix(lambda)
    betaini <- array(beta[, 2:ncol(as.matrix(beta[, , 1])), ], dim = c(k, k * p, nrow(lambda)))

    betafin <- gamloopOO(betaini, Y, Z, lambda, eps, YMean, ZMean, as.matrix(betaini[, , 1]), k, p, weights, groups, separate_lambdas)


    if (MN) {

        betafin <- adjust_mn_var(betafin, C)

    }


    return(betafin)

}


                                        # indexing for efx
vxsubs <- function(i, k, m, p, s) {


    vv <- c(((i - 1) * k + 1):((i - 1) * k + k), ((i - 1) * m + k * p + 1):((i - 1) * m + k * p + m))

    vv

}

prox2 <- function(v, lambda, k, p, m, s) {

    for (i in 1:p) {

        if (i <= s) {

            vv <- vxsubs(i, k, m, p, s)
            F1 <- 0

        }

        if (i > s) {
            vv <- ((i - 1) * k + 1):((i - 1) * k + k)
            F1 <- 1
        }

        v2 <- proxvx2(v[vv], p, lambda, m, k, F1)
        v[vv] <- v2


    }

    v

}


fistaX <- function(Y, Z, beta, p, k1, lambda, eps, tk, m, s) {

    for (i in 1:k1) {

        phiOLD <- beta[i, ]
        phiOLDOLD <- beta[i, ]
        j <- 1
        thresh <- 10 * eps

        while (thresh > eps) {

            v <- matrix(phiOLD + ((j - 2)/(j + 1)) * (phiOLD - phiOLDOLD), nrow = 1)
            phiR <- prox2(v + tk * as.vector((Y[, i] - v %*% Z) %*% t(Z)), tk * lambda, k1, p, m, s)
            thresh <- max(abs(phiR - v))
            phiOLDOLD <- phiOLD
            phiOLD <- phiR
            j <- j + 1
        }

        beta[i, ] <- phiR

    }

    beta
}



                                        # Lag weighted lasso: VAR only

.lassoVARTL <- function(B, Z, Y, lambda, eps, p, MN, alpha, C, YMean, ZMean) {

    if (!is.matrix(Y)) {
        Y <- matrix(Y, ncol = 1)
    }

    k <- ncol(Y)


    tk <- 1/max(Mod(eigen(Z %*% t(Z), only.values = TRUE)$values))
    BFOO1 <- as.matrix(B[, 2:ncol(B[, , 1]), 1])
    BFOO <- array(B[, 2:ncol(as.matrix(B[, , 1])), ], dim = c(k, k * p, length(lambda) * length(alpha)))
    p2 <- 1:p

    gran2 <- length(lambda)

    for (i in seq_len(length(alpha))) {

        W <- rep(p2^(-alpha[i]), each = k)
        ZADJ <- diag(W) %*% Z

        B[, , (1 + (i - 1) * gran2):(i * length(lambda))] <- gamloopFista(array(BFOO[, , (1 + (i - 1) * gran2):(i * length(lambda))],
                                                                                dim = c(k, k * p, length(lambda))), Y, ZADJ, as.matrix(lambda), eps, as.matrix(YMean), as.matrix(ZMean), BFOO1,
                                                                          k, p, tk, k, p)


        for (j in (1 + (i - 1) * gran2):(i * length(lambda))) {

            B[, 2:(k * p + 1), j] <- B[, 2:(k * p + 1), j] %*% diag(W)

        }


    }
    if (MN) {
        B <- adjust_mn_var(B, C)

    }

    return(B)
}


# Block Group Lasso

.GroupLassoVAROpt <- 
function (beta, groups, Y, Z, gamm, INIactive, eps,p,MN) 
{
        if(class(Y)!="matrix")
        {Y <- matrix(Y,ncol=1)}
    k <- ncol(Y)
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    if(k>1){    
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))}
    else{Y <- Y-mean(Y)}    
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    jj = groups
    jjfull = groups
   Eigsys <- Eigencomp(Z,jj,length(jj),k)
    M2 <- Eigsys$M3
    eigvals <- Eigsys$eigval
    eigvecs <- Eigsys$eigvec

    jjcomp = .groupfuncomp(p, k)
    beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k,k*p,length(gamm)))
        
       ## beta <- beta[,2:ncol(beta[,,1]),]
BB <- GamLoopGL(beta,INIactive,gamm,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k,p*k,M2,eigvals,eigvecs)
    return(BB)
}

# Group Lasso Own/Other

.GroupLassoOONew <- function (beta, groups, Y, Z, gamm, INIactive, eps,p,MN) 
{
            if(class(Y)!="matrix")
        {Y <- matrix(Y,ncol=1)}
    k <- ncol(Y)
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
if(k>1){
    YMEAN <- apply(Y, 1, mean)}
   else{YMEAN <- mean(Y)}         
    ZMEAN <- apply(Z,1,mean)

        if(k>1){    
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))}
    else{Y <- Y-mean(Y)}    
  Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    kk = groups
    kkfull = groups
    kkcomp = .lfunctioncomp(p, k)
    ZZ <- kronecker(t(Z), diag(k))
    Eigsys <- EigencompOO(ZZ,kk,length(kk),k)
    M2 <- Eigsys$M3
    eigvals <- Eigsys$eigval
    eigvecs <- Eigsys$eigvec
## beta <- beta[,2:ncol(beta[,,1]),]
beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k,k*p,length(gamm)))
            
BB <- GamLoopGLOO(beta,INIactive,gamm,Y,ZZ,kk,kkfull,kkcomp,eps,YMEAN,ZMEAN,k,p*k,M2,eigvals,eigvecs)

return(BB)
}


# Block Sparse Group Lasso
.SparseGroupLassoVAROpt<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN) 
{
    k <- ncol(Y)
   # p <- nrow(Z)/2
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    M1f = list()
    M2f = list()
    eigs <- c()
    q1 = list()
   jj=.groupfun(p,k)
    for (j in 1:length(jj)) {
        M1f[[j]] <- Z[jj[[j]], ]
        M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])
        gg1 <- powermethod(M2f[[j]], q1a[[j]])
        eigs[j] <- gg1$lambda
        q1[[j]] <- gg1$q1
    }
    
    jj=.groupfuncpp(p,k)
    jjfull=jj
    jjcomp = .groupfuncomp(p, k)
    beta <- beta[,2:ncol(beta[,,1]),]
BB <- GamLoopSGL(beta,INIactive,gamm,alpha,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k,p*k,M1f,M2f,eigs)
BB$q1 <- q1
    return(BB)
}



# Sparse Group Lasso/Own/Other
.SparseGroupLassoVAROptOO<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN) 
{
    k <- ncol(Y)
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))

    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
   ZZ <- kronecker(t(Z),diag(k))  
    M1f = list()
    M2f = list()
    jj=.lfunction3(p,k)
    eigs <- c()
    q1 = list()
    for (j in 1:length(jj)) {
        M1f[[j]] <- ZZ[,jj[[j]] ]
        M2f[[j]] <- crossprod(ZZ[,jj[[j]]])
        gg1 <- powermethod(M2f[[j]], q1a[[j]])
        eigs[j] <- gg1$lambda
        q1[[j]] <- gg1$q1
    }    
    jj=.lfunction3cpp(p,k)
    jjfull=jj
    jjcomp = .lfunctioncomp(p, k)
    beta <- beta[,2:ncol(beta[,,1]),]
    BB <- GamLoopSGLOO(beta,INIactive,gamm,alpha,Y,ZZ,jj,jj,jjcomp,eps,YMean,ZMean,k,p*k,M1f,M2f,eigs)
   BB$q1 <- q1
    return(BB)
}



# Lasso-VAR Fista Implementation
.lassoVARFist <- function (B, Z, Y, gamm, eps,p,MN) 
{
   if(class(Y)!="matrix")
   {Y <- matrix(Y,ncol=1)}
  
   k <- ncol(Y)
   if(MN==TRUE)
        {
    C=matrix(0,nrow=k,ncol=k*p)        
    diag(C) <- rep(1,k)
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
        }
    
    Y <- t(Y)
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
     Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
     Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
    Y <- t(Y)
    ## tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))

    BFOO1 <- as.matrix(B[, 2:ncol(B[, , 1]), 1])
    BFOO <- array(B[,2:ncol(as.matrix(B[,,1])),],dim=c(k,k*p,length(gamm)))

    ## beta <- gamloopFista(B[, 2:ncol(B[, , 1]), ], Y, Z, gamm, eps, 
    ##     as.matrix(YMean), as.matrix(ZMean), BFOO,k,p)
       beta <- gamloopFista(BFOO, Y, Z, gamm, eps, 
        as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p)

    if(MN==TRUE)
         {
     for(i in 1:(dim(beta)[3]))
             
         beta[,2:ncol(beta[,,i]),i] <- beta[,2:ncol(beta[,,i]),i]+C
          }
    return(beta)
}



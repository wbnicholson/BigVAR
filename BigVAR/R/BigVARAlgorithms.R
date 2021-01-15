# BigVAR Algorithms
# Most of the computationally expensive portions of the code have been exported to C++


## ST <- function(z,gam)
## {
##     if(z>0 & gam<abs(z))
##     {return(z-gam)}
    
##     if(z<0 & gam<abs(z))
##     {return(z+gam)}
##     if(gam >=abs(z) )
##     {return(0)}

## }

## MCP_pen <- function(z,lambda,gamma){
##     if(abs(z)>gamma*lambda){
##         return(z)
##     }
##     if(abs(z)<=gamma*lambda){
##         return(ST(z,lambda)/(1-1/gamma))
##     }
## }

## MCP <- function(Z, Y, lambda,gamma, maxiter=1e3, tol, B=NULL,znorm2) {
##     pk <- nrow(Z)
##     T <- ncol(Z)
##     k <- nrow(Y)
##     if (is.null(B)) B <- matrix(0, nrow=length(lambda), ncol=pk)
##     BOLD <- B[1,]

##                                         # need to update the partial residual or bring in entire grid
##     ## r=as.numeric(Y[1,]-B%*%Z)
##     r=Y[1,]
##     gloss <- sum(r^2)
##     sdy <- sqrt(gloss/T)
##     tot_iter=0
##     max_iter=10000
##     for (l in seq(nrow(B))) {
##         thresh=10
##         if(l>1){
##             BOLD <- B[l-1,]
##         }
##         while (tot_iter < max_iter) {
            
##             tot_iter=tot_iter+1
##             for (m in seq(ncol(B))) {

##                 zhat <- as.numeric(r%*%Z[m,]/T)+BOLD[m]

##                 B[l,m] <- MCP_pen(zhat,lambda[l],gamma)
##                 if(l==30){
##                     print(sprintf("zhat %f",zhat))
##                     print(sprintf("b %f",B[l,m]))
##                 }

##                 shift <- B[l,m]-BOLD[m]
##                 if(shift!=0){
##                     r <- as.numeric(r-(shift)%*%Z[m,])
##                 }
                
##             }
##             thresh <- max(abs(B[l,]-BOLD))
##             BOLD <- B[l,]
##             if (thresh < tol*sdy) {
##                 break
##             }


##         }

##     }
##     print(tot_iter)
##     return(B)
## }

## .mcp_fit <- function(trainY,trainZ,gamm,beta){
##     ## browser()
##     trainZ <- t(trainZ)
##     ## mod <- scad(trainZ,trainY,lambda=gamm,penalty="SCAD")
##     for(i in 1:length(gamm)){
##         if(i<=ncol(mod$beta)){
##             beta[,,i] <- mod$beta[,i]
##             }
##     }
##     ## browser()
##     return(beta)
## }

## .MCPFit <- function (B, Z, Y, gamm, eps,p,MN,k,k1,s,m,C1,intercept,group,gamma=3) 
## {
    
##     if(!"matrix"%in%class(Y))
##    {
##        Y <- matrix(Y,ncol=1)
##    }
  

##     if(MN)
##         {
##             C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
##             diag(C) <- C1
##             Y <- t(Y)
##             Y <- Y-C%*%Z
##             Y <- t(Y)
##         }
    
##     Y <- t(Y)
##     intercept=TRUE
##     if(intercept){
##     YMean <- c(apply(Y, 1, mean))
##     ZMean <- c(apply(Z, 1, mean))
##     Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
##     Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
##     }else{
##         YMean <- rep(0,nrow(Y))
##         ZMean <- rep(0,nrow(Z))
##         }
##     Y <- t(Y)

##     tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))

##     BFOO1 <- matrix(B[, 2:dim(B)[2], 1],ncol=k1)
##     ## BFOO <- array(B[,2:ncol(as.matrix(B[,,1])),],dim=c(k1,k1*p+(k-k1)*s,length(gamm)))
    
##     nc <- apply(B,3,ncol)[1]
##     ## BFOO1 <-as.matrix(B[, 2:nc,1,drop=F])
##     BFOO <- B[,2:nc,,drop=F]
##     ## browser()
##     if(group=="MCP"){
##        mcp=TRUE
##     }else{
##        mcp=FALSE
##     }
##     ## browser()
##     ## if(length(gamm)==1){browser()}
##     beta <- gamloopMCP(BFOO,Y,Z,as.matrix(gamm),eps,as.matrix(YMean),as.matrix(ZMean),gamma=3,mcp)
##     ## foo2=MCP(Z,as.matrix(t(Y[,1])),lambda=gamm,gamma=3,maxiter=1000,tol=1e-4,B=NULL,NULL)

##     ## dyn.load('/home/will/Downloads/ncvreg_control/ncvreg/src/aaa')
##     ## X <- t(Z)
##     ## res2 <- .Call("rawfit_gaussian", t(Z), t(Y), BFOO[,,1], "MCP", gamm, 1e-4, as.integer(10000), as.double(3),rep(1,nrow(Z)), alpha=1)

##     ## res <- .Call("cdfit_gaussian", t(Z), Y, "MCP", gamm, as.double(1e-4), as.integer(10000), as.double(3), rep(1,ncol(X)),1, as.integer(ncol(X)+1), as.integer(FALSE))

##     ## beta <- gamloopFista(BFOO, Y, Z, as.matrix(gamm), eps, 
##     ##                      as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k1,s,separate_lambdas)

##     if(MN)
##         {
##             for(i in 1:(dim(beta)[3]))             
##                 beta[,2:dim(beta[,,i,drop=F])[2],i] <- beta[,2:dim(beta[,,i,drop=F])[2],i]+C
##         }
    
##     return(beta)

## }



# Sparse Own/Other (VAR)
.SparseGroupLassoVAROO<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,dual=FALSE,C1,intercept) 
{

        ## if(length(alpha)>1){browser()}

    m <- 0
    k <- ncol(Y)
    if(MN)
        {

            C <- matrix(0,nrow=k,ncol=k*p)        
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }

    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))    
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    ZZ <- kronecker(t(Z),diag(k))  
    M1f <- list()
    M2f <- list()
    jj <- .lfunction3(p,k)
    eigs <- c()
    q1 <- list()
    # get step size from inverse of max eigenvalue via power method
    for (j in 1:length(jj)) {
        M1f[[j]] <- ZZ[,jj[[j]] ]
        M2f[[j]] <- crossprod(ZZ[,jj[[j]]])
        gg1 <- powermethod(M2f[[j]], q1a[[j]])
        eigs[j] <- gg1$lambda
        q1[[j]] <- gg1$q1
    }    
    jj <- .lfunction3cpp(p,k)
    jjfull <- jj
    jjcomp <- .lfunctioncomp(p, k)
    dims <- dim(beta)
    beta <- array(beta[,2:ncol(beta[,,1]),],dim=c(dims[1],dims[2]-1,dims[3]))
    ## if(nrow(gamm)>1){browser()}
    ## browser()
    if(!dual){

        BB <- GamLoopSGLOO(beta,INIactive,gamm,alpha,Y,ZZ,jj,jj,jjcomp,eps,YMean,ZMean,k,p*k,M2f,eigs,m)


    }else{

        ## browser()
        BB <- GamLoopSGLOODP(beta,INIactive,gamm,alpha,Y,ZZ,jj,jj,jjcomp,eps,YMean,ZMean,k,p*k,M2f,eigs,m)
        
        }

    BB$q1 <- q1

    if (MN)
        {
            for(i in 1:(dim(BB$beta)[3]))

                {
                
                    BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
                }
        }
    return(BB)
}



# Sparse Lag (VAR)
.SparseGroupLassoVAR<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,C1,intercept) 
{
    k <- ncol(Y)

    if (MN)
        {

            C <- matrix(0,nrow=k,ncol=k*p)        
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }

    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))    
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()
    jj <- .groupfun(p,k)


                                        # get step size from inverse of max eigenvalue via power method
    for (j in 1:length(jj)) {
        M1f[[j]] <- Z[jj[[j]], ]
        M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])
        gg1 <- powermethod(M2f[[j]], q1a[[j]])
        eigs[j] <- gg1$lambda
        q1[[j]] <- gg1$q1
    }
    
    jj <- .groupfuncpp(p,k)
    jjfull <- jj
    jjcomp <- .groupfuncomp(p, k)
    dims <- dim(beta)
    beta <- array(beta[,2:ncol(beta[,,1]),],dim=c(dims[1],dims[2]-1,dims[3]))
    ## browser()    
    ## beta <- beta[,2:ncol(beta[,,1]),]
    
    BB <- GamLoopSGL(beta,INIactive,gamm,alpha,Y,Z,jj,jjfull,jjcomp,eps,YMean,as.matrix(ZMean),k,p*k,M1f,M2f,eigs)
    BB$q1 <- q1

        if (MN)
        {
            for(i in 1:(dim(BB$beta)[3]))
                {
                
                    BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
                }

        }
    
    return(BB)
}

# Sparse Lag (VAR) Dual Search
.SparseGroupLassoVARDual<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,C1,intercept) 
{
k <- ncol(Y)
    if (MN)
        {

            C <- matrix(0,nrow=k,ncol=k*p)
            diag(C) <- C1
            ## diag(C) <- rep(1,k)
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }

    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))    
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()
    jj <- .groupfun(p,k)
    # get step size from inverse of max eigenvalue via power method
    for (j in 1:length(jj)) {
        M1f[[j]] <- Z[jj[[j]], ]
        M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])
        gg1 <- powermethod(M2f[[j]], q1a[[j]])
        eigs[j] <- gg1$lambda
        q1[[j]] <- gg1$q1
    }
    
    jj <- .groupfuncpp(p,k)
    jjfull <- jj
    jjcomp <- .groupfuncomp(p, k)
    ngp <- length(alpha)*length(gamm)
    dims <- dim(beta)
    beta <- array(beta[,2:ncol(beta[,,1]),],dim=c(dims[1],dims[2]-1,dims[3]))
    ## browser()
    BB <- GamLoopSGLDP(beta,INIactive,gamm,alpha,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k,p*k,M1f,M2f,eigs)
    BB$q1 <- q1

        if (MN)
        {
            for(i in 1:(dim(BB$beta)[3]))
                {
                
                    BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
                }

        }
    
    return(BB)
}


# Lag Group (VAR/VARX-L)
.GroupLassoVAR1 <- 
function (beta, groups,jjcomp, Y, Z, gamm, INIactive, eps,p,MN,k,k1,s,C1,intercept) 
{

    
    
    if(!"matrix"%in%class(Y))
        {
            Y <- matrix(Y,ncol=1)
        }

    m <- k-k1

    if(MN)
        {
            C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
            ## diag(C) <- rep(1,k1)
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }
        
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))    
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }

    if(k>1){    

        Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))}

    else{Y <- Y-mean(Y)}    

    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    jj <- groups
    jjfull <- groups
    
    Eigsys <- Eigencomp(Z,jj,length(jj),k1)
    M2 <- Eigsys$M3
    eigvals <- Eigsys$eigval
    eigvecs <- Eigsys$eigvec
    
    beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,(k1)*p+s*m,length(gamm)))
    BB <- GamLoopGL2(beta,INIactive,gamm,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k1,p*k1+m*s,M2,eigvals,eigvecs)

    if(MN)
        {
            for(i in 1:(dim(BB$beta)[3]))
                {
                
                    BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
                }

        }
    return(BB)
}



# Group Lasso Own/Other (VARXL)
.GroupLassoOOX <- function (beta, groups, Y, Z, gamm, INIactive, eps,p,MN,k,k1,s,C1,intercept) 
{
    m <- k-k1

    if(MN)
        {
            C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
            ## diag(C) <- rep(1,k1)
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }
    

    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z

    if(intercept){
    if(k1>1){

        YMEAN <- apply(Y, 1, mean)

    }else{
        YMEAN <- mean(Y)
    }         

    ZMEAN <- apply(Z,1,mean)
    }else{

        YMEAN <- rep(0,nrow(Y))
        ZMEAN <- rep(0,nrow(Z))
        
        }
    if(intercept){
    if(k>1){    
        Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    }else{Y <- Y-mean(Y)}    

    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }
    kk=diaggroupfunVARX(p,k,k1,s)
    kkcomp <- diaggroupfunVARXcomp(p,k,k1,s)

    ZZ <- kronecker(t(Z), diag(k1))

    Eigsys <- EigencompOO(ZZ,kk,length(kk),k1)
    M2 <- Eigsys$M3
    eigvals <- Eigsys$eigval
    eigvecs <- Eigsys$eigvec

    kkfull <- kk

    beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,(k1)*p+m*s,length(gamm)))

    BB <- GamLoopGLOO(beta,INIactive,gamm,Y,ZZ,kk,kkfull,kkcomp,eps,YMEAN,ZMEAN,k1,p*(k1)+m*s,M2,eigvals,eigvecs,k1)

    if(MN)
         {
             for(i in 1:(dim(BB$beta)[3]))
                 {
                     BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
                 }
         }
    
    return(BB)
}


# Own/Other Group VAR-L
.GroupLassoOO <- function (beta, groups, Y, Z, gamm, INIactive, eps,p,MN,C1,intercept) 
{

    if(!"matrix"%in%class(Y))
        {
            Y <- matrix(Y,ncol=1)
        }

    k <- ncol(Y)

    if(MN)
        {
    C <- matrix(0,nrow=k,ncol=k*p)        
    ## diag(C) <- rep(1,k)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
        }

            
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z

    if(intercept){
    if(k>1){

        YMEAN <- apply(Y, 1, mean)

    }else{
        YMEAN <- mean(Y)
    }         

    ZMEAN <- apply(Z,1,mean)
    }else{

        YMEAN <- rep(0,nrow(Y))
        ZMEAN <- rep(0,nrow(Z))
        
        }
    if(intercept){
    if(k>1){    
        Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    }else{Y <- Y-mean(Y)}    

    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }
    kk <- groups
    kkfull <- groups
    kkcomp <- .lfunctioncomp(p, k)

    ZZ <- kronecker(t(Z), diag(k))

    Eigsys <- EigencompOO(ZZ,kk,length(kk),k)
    M2 <- Eigsys$M3
    eigvals <- Eigsys$eigval
    eigvecs <- Eigsys$eigvec

    beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k,k*p,length(gamm)))
            

    BB <- GamLoopGLOO(beta,INIactive,gamm,Y,ZZ,kk,kkfull,kkcomp,eps,YMEAN,ZMEAN,k,p*k,M2,eigvals,eigvecs,k)

        if(MN)
            {
                for(i in 1:(dim(BB$beta)[3]))
                    {
                        BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C
                    }

            }
            
    return(BB)
}



# Elementwise HVAR
.HVARElemAlg <- function (beta, Y, Z, lambda, eps,p,MN,C1,intercept,separate_lambdas=FALSE) 
{

    k <- ncol(Y)
    if(MN)
        {
            C <- matrix(0,nrow=k,ncol=k*p)        
            ## diag(C) <- rep(1,k)
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }
    
    k <- ncol(Y)
    betafin <- beta
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
    YMean <- c(apply(Y, 2, mean))
    ZMean <- c(apply(Z, 1, mean))

    Y <- Y - matrix(c(rep(1, nrow(Y))), ncol = 1) %*% matrix(c(apply(Y, 
        2, mean)), nrow = 1)
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,ncol(Y))
        ZMean <- rep(0,nrow(Z))
        }
    
    tk <- 1/max(Mod(eigen(Z %*% t(Z))$values))
    lambda <- as.matrix(lambda)
    betaini <- array(beta[, 2:ncol(beta[, , 1]), ],dim=c(k,k*p,nrow(lambda)))
    betafin <- gamloopElem(betaini,Y,Z,lambda,eps,YMean,ZMean,as.matrix(betaini[,,1]),k,p,separate_lambdas)

    if(MN)
        {
            for(i in 1:(dim(betafin)[3]))
                {
                    betafin[,2:ncol(betafin[,,i]),i] <- betafin[,2:ncol(betafin[,,i]),i]+C
                }

        }
   
    return(betafin)
}

.lassoVARFistX <- function (B, Z, Y, gamm, eps,p,MN,k,k1,s,m,C1,intercept,separate_lambdas=FALSE) 
{
    
    if(!"matrix"%in%class(Y))
   {
       Y <- matrix(Y,ncol=1)
   }
  

    if(MN)
        {
            C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }
    
    Y <- t(Y)
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
    Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    Y <- t(Y)

    tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))

    BFOO1 <- matrix(B[, 2:dim(B)[2], 1],ncol=k1)
    ## BFOO <- array(B[,2:ncol(as.matrix(B[,,1])),],dim=c(k1,k1*p+(k-k1)*s,length(gamm)))
    
    nc <- apply(B,3,ncol)[1]
    ## BFOO1 <-as.matrix(B[, 2:nc,1,drop=F])
    BFOO <- B[,2:nc,,drop=F]
    ## browser()
    beta <- gamloopFista(BFOO, Y, Z, as.matrix(gamm), eps, 
                         as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k1,s,separate_lambdas)

    if(MN)
        {
            for(i in 1:(dim(beta)[3]))             
                beta[,2:dim(beta[,,i,drop=F])[2],i] <- beta[,2:dim(beta[,,i,drop=F])[2],i]+C
        }
    
    return(beta)

}

.lassoVARFistXEN <- function (B, Z, Y, gamm,alpha, eps,p,MN,k,k1,s,m,C1,intercept,separate_lambdas=FALSE) 
{

    if(!"matrix"%in%class(Y))
   {
       Y <- matrix(Y,ncol=1)
   }
  

    if(MN)
        {
            C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }
    
    Y <- t(Y)
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
    Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    Y <- t(Y)

    tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))

    BFOO1 <- matrix(B[, 2:dim(B)[2], 1],ncol=k1)
    
    nc <- apply(B,3,ncol)[1]
    BFOO <- B[,2:nc,,drop=F]
    if(length(alpha)==1){alpha=rep(alpha,dim(B)[3])}
    
    beta <- gamloopFistaEN(BFOO, Y, Z, as.matrix(gamm),as.matrix(alpha), eps, 
                         as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k1,s,separate_lambdas)

    if(MN)
        {
            for(i in 1:(dim(beta)[3]))             
                beta[,2:dim(beta[,,i,drop=F])[2],i] <- beta[,2:dim(beta[,,i,drop=F])[2],i]+C
        }
    
    return(beta)

}


#Basic VAR Fista Implementation
.lassoVARFist <- function (B, Z, Y, gamm, eps,p,MN,C1,intercept,separate_lambdas=FALSE) 
{

    if(!"matrix"%in%class(Y))
        {
            Y <- matrix(Y,ncol=1)
        }
  
    k <- ncol(Y)

    if(MN)
        {
            C <- matrix(0,nrow=k,ncol=k*p)        
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }
    
    Y <- t(Y)
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
    Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    Y <- t(Y)
    
    tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))

    nc <- apply(B,3,ncol)[1]
    BFOO1 <- matrix(B[, 2:dim(B)[2], 1],nrow=k,ncol=k*p)
    BFOO <- B[,2:nc,,drop=F]
    ## browser()
    ## beta <- gamloopFistaEN(BFOO, Y, Z, as.matrix(gamm),alpha=alphas, eps, 
    ##     as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k,p,separate_lambdas)


    beta <- gamloopFista(BFOO, Y, Z, as.matrix(gamm), eps, 
        as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k,p,separate_lambdas)

    if(MN)
        {
            for(i in 1:(dim(beta)[3]))             
                beta[,2:dim(beta[,,i,drop=F])[2],i] <- beta[,2:dim(beta[,,i,drop=F])[2],i]+C
        }
    return(beta)
}


#Basic-elasticNET VAR Fista Implementation
.lassoVARFistEN <- function (B, Z, Y, gamm,alpha, eps,p,MN,C1,intercept,separate_lambdas=FALSE) 
{

    if(!"matrix"%in%class(Y))
        {
            Y <- matrix(Y,ncol=1)
        }
  
    k <- ncol(Y)

    if(MN)
        {
            C <- matrix(0,nrow=k,ncol=k*p)        
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }
    
    Y <- t(Y)
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
    Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    Y <- t(Y)
    
    tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))

    nc <- apply(B,3,ncol)[1]
    BFOO1 <- matrix(B[, 2:dim(B)[2], 1],nrow=k,ncol=k*p)
    BFOO <- B[,2:nc,,drop=F]
    if(length(alpha)==1){alpha=rep(alpha,dim(B)[3])}

    beta <- gamloopFistaEN(BFOO, Y, Z, as.matrix(gamm),alpha=alpha, eps, 
        as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k,p,separate_lambdas)

    if(MN)
        {
            for(i in 1:(dim(beta)[3]))             
                beta[,2:dim(beta[,,i,drop=F])[2],i] <- beta[,2:dim(beta[,,i,drop=F])[2],i]+C
        }
    return(beta)
}

# Componentwise HVAR
.HVARCAlg <- function(beta,Y,Z,lambda,eps,p,MN,C1,intercept,separate_lambdas=FALSE)
    {
    if(!"matrix"%in%class(Y))
        {
            Y <- matrix(Y,ncol=1)
        }

    k <- ncol(Y)

    if(MN)
        {
            C <- matrix(0,nrow=k,ncol=k*p)        
            ## diag(C) <- rep(1,k)
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }


    betafin <- beta
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
    YMean <- c(apply(Y, 2, mean))
    ZMean <- c(apply(Z, 1, mean))
    }else{
        YMean <- rep(0,ncol(Y))
        ZMean <- rep(0,nrow(Z))
    }
    if(intercept){
    if(k>1){
        Y <- Y -  matrix(c(rep(1, nrow(Y))),ncol=1)%*%matrix(c(apply(Y, 2, mean)),nrow=1)
    }else{
        Y <- Y-mean(Y)
    }

    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }
    tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))

    lambda <- as.matrix(lambda)
    if(separate_lambdas){
    betaini <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k,k*p,nrow(lambda)))
    }else{
    betaini <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k,k*p,length(lambda)))
    }
    betafin <- gamloopHVAR(betaini,Y,Z,lambda,eps,YMean,ZMean,as.matrix(betaini[,,1]),k,p,separate_lambdas)

    if(MN)
        {
            for(i in 1:(dim(betafin)[3]))
                {
                    betafin[,2:ncol(betafin[,,i]),i] <- betafin[,2:ncol(betafin[,,i]),i]+C
                }
        }    
    return(betafin)
}

# Endogenous First VARX-L
.EFVARX <- function(beta,Y,Z,lambda,eps,MN,k1,s,m,p,C1,intercept)
    {        

        

        if(MN)
        {


            C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
            ## diag(C) <- rep(1,k1)
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }
        
    betafin <- beta  
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    if(intercept){  
    YMEAN <- apply(Y, 1, mean)
    ZMEAN <- apply(Z, 1, mean)
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
        YMEAN <- rep(0,nrow(Y))
        ZMEAN <- rep(0,nrow(Z))
    }
    Y <- t(Y)
    tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))


        if (k1==1)
        {

            betaini <- array(beta[,2:(k1*p+m*s+1),],dim=c(1,k1*p+m*s,dim(beta)[3]))

        }else{

            betaini <- beta[,2:ncol(beta[,,1]),]

        }
   
        for(i in 1:length(lambda))
        {

            if(dim(beta)[3]>1){

                betaF <- fistaX(Y,Z,matrix(betaini[,,i],nrow=k1),p,k1,lambda[i],eps,tk,m,s)
            }

            else{ 

                betaF <- fistaX(Y,Z,matrix(betaini,nrow=k1),p,k1,lambda[i],eps,tk,m,s)

            }

            nu <- YMEAN - betaF %*% ZMEAN

            betafin[,,i] <- cbind(nu, betaF)

        }

        if(MN)
        {

            for(i in 1:(dim(betafin)[3]))
            {

                betafin[,2:ncol(betafin[,,i]),i] <- betafin[,2:ncol(betafin[,,i]),i]+C

            }

        }
     
        return(betafin)
        }

# HVAR Own/Other
.HVAROOAlg <- function(beta,Y,Z,lambda,eps,p,MN,C1,intercept,separate_lambdas=FALSE)
    {

        k <- ncol(Y)
        betafin <- beta
        YOLD <- Y
        ZOLD <- Z

        if(MN)
        {
            C <- matrix(0,nrow=k,ncol=k*p)        
            ## diag(C) <- rep(1,k)
            diag(C) <- C1
            Y <- t(Y)
            Y <- Y-C%*%Z
            Y <- t(Y)
        }

        weights <- sqrt(c(rep(c(1,k-1),length=2*p)))
        if(intercept){
        YMEAN <- c(apply(Y, 2, mean))
        ZMEAN <- c(apply(Z, 1, mean))
        Y <- Y -  matrix(c(rep(1, nrow(Y))),ncol=1)%*%matrix(c(apply(Y, 2, mean)),nrow=1)
        Z <- Z -  c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
        }else{
            YMEAN <- rep(0,ncol(Y))
            ZMEAN <- rep(0,nrow(Z))
            }
        groups <- list()
        for(i in 1:k)
        {

            groups[[i]] <- .vecoovarscpp(p,k,i)

        }

        lambda <- as.matrix(lambda)
        betaini <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k,k*p,nrow(lambda)))

        betafin <- gamloopOO(betaini,Y,Z,lambda,eps,YMEAN,ZMEAN,as.matrix(betaini[,,1]),k,p,weights,groups,separate_lambdas)


        if(MN)
        {
            for(i in 1:(dim(betafin)[3]))
            {

                betafin[,2:ncol(betafin[,,i]),i] <- betafin[,2:ncol(betafin[,,i]),i]+C

            }

        }


        return(betafin)

    }


# indexing for efx
vxsubs <- function(i,k,m,p,s)
{


        vv <- c(((i-1)*k+1):((i-1)*k+k),((i-1)*m+k*p+1):((i-1)*m+k*p+m))
                  
        vv

}

prox2 <- function(v,lambda,k,p,m,s)
{

    for(i in 1:p)
    {

        if(i<=s){

            vv <- vxsubs(i,k,m,p,s)
            F1 <- 0

        }

        if(i>s){
            vv <- ((i-1)*k+1):((i-1)*k+k)
            F1 <- 1
        }
        
        v2 <- proxvx2(v[vv],p,lambda,m,k,F1)
        v[vv] <- v2


    }
    
    v

}


fistaX <- function(Y,Z,beta,p,k1,lambda,eps,tk,m,s)
{

    for(i in 1:k1)
    {

        phiOLD <- beta[i,]
        phiOLDOLD <- beta[i,]
        j <- 1
        thresh <- 10*eps

        while(thresh>eps)
    {

        v <- matrix(phiOLD+((j-2)/(j+1))*(phiOLD-phiOLDOLD),nrow=1)
        phiR <- prox2(v+tk*as.vector((Y[,i]-v%*%Z)%*%t(Z)),tk*lambda,k1,p,m,s)
        thresh <- max(abs(phiR-v))
        phiOLDOLD <- phiOLD
        phiOLD <- phiR
        j <- j+1
    }

        beta[i,] <- phiR

    }

    beta
}

# Sparse Lag VARX-L
.SparseGroupLassoVARX<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,k,s,k1,C1,intercept) 
{

    m <- k-k1

    if(MN)
    {
        C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
        diag(C) <- C1
        Y <- t(Y)
        Y <- Y-C%*%Z
        Y <- t(Y)
    }

    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()
    jj <- groupfunVARXLG(p,k,k1,s)

    for (j in 1:length(jj)) {
        M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])

        if(j<=p){

            gg1 <- powermethod(M2f[[j]], q1a[[j]])
            eigs[j] <- gg1$lambda
            q1[[j]] <- gg1$q1
        }

        else{

            M2f[[j]] <- as.vector(Z[jj[[j]], ]) %*% as.vector(t(Z[jj[[j]], ]))
            eigs[j] <- M2f[[j]]

        }

    }


    jj <- groupfunVARX(p,k,k1,s)    
    jjfull <- jj
    jjcomp <- groupfunVARXcomp(p,k,k1,s)

    beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,(k1)*p+(s*m),length(gamm)))
    BB <- GamLoopSGLX(beta,INIactive,gamm,alpha,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k1,(k1)*p+(s*m),M2f,eigs,k1)

    BB$q1 <- q1
    
    if(MN)
    {

        for(i in 1:(dim(BB$beta)[3]))
        {

            BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C

        }

    }

    return(BB)
}



.SparseGroupLassoVARXDual<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,q1a,p,MN,k,s,k1,C1,intercept) 
{

    m <- k-k1

    if(MN)
    {
        C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
        ## diag(C) <- rep(1,k1)
        diag(C) <- C1
        Y <- t(Y)
        Y <- Y-C%*%Z
        Y <- t(Y)
    }

    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()
    jj <- groupfunVARXLG(p,k,k1,s)

    for (j in 1:length(jj)) {

        M2f[[j]] <- Z[jj[[j]], ] %*% t(Z[jj[[j]], ])

        if(j<=p){

            gg1 <- powermethod(M2f[[j]], q1a[[j]])
            eigs[j] <- gg1$lambda
            q1[[j]] <- gg1$q1
        }

        else{

            M2f[[j]] <- as.vector(Z[jj[[j]], ]) %*% as.vector(t(Z[jj[[j]], ]))
            eigs[j] <- M2f[[j]]

        }

    }


    jj <- groupfunVARX(p,k,k1,s)    
    jjfull <- jj
    jjcomp <- groupfunVARXcomp(p,k,k1,s)

    beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,(k1)*p+(s*m),nrow(gamm)*length(alpha)))

    BB <- GamLoopSGLXDP(beta,INIactive,gamm,alpha,Y,Z,jj,jjfull,jjcomp,eps,YMean,ZMean,k1,(k1)*p+(s*m),M2f,eigs,k1)

    BB$q1 <- q1

    if(MN)
    {

        for(i in 1:(dim(BB$beta)[3]))
        {

            BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C

        }

    }

    return(BB)
}


# Sparse Own/Other (VARX)
.SparseGroupLassoVAROOX<-function(beta,groups,Y,Z,gamm,alpha,INIactive,eps,p,MN,k1,s,k,dual=FALSE,C1,intercept) 
{

    m <- k-k1
    if(MN)
    {

        C <- matrix(0,nrow=k1,ncol=k1*p+s*m)        
        diag(C) <- C1
        Y <- t(Y)
        Y <- Y-C%*%Z
        Y <- t(Y)
    }

    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    }else{
        YMean <- rep(0,nrow(Y))
        ZMean <- rep(0,nrow(Z))
        }
    M1f <- list()
    M2f <- list()
    eigs <- c()
    q1 <- list()

    # function for R calculations
    jj<- diaggroupfunVARXLG(p,k,k1,s)
    # for c++ calculations
    kk <- diaggroupfunVARX(p,k,k1,s)
    jjcomp <- diaggroupfunVARXcomp(p,k,k1,s)

    ZZ <- kronecker(t(Z),diag(k1))  
    

    for (j in 1:length(jj)) {

        M2f[[j]] <- crossprod(ZZ[,jj[[j]]])

        eigs[j] <- max(Mod(eigen(M2f[[j]],only.values=TRUE)$values))
    }

    jjfull <- kk

    if(!dual){

        gran2 <- length(gamm)

    }else{

        gran2 <- nrow(gamm)*ncol(gamm)
        
        }
    
    beta <- array(beta[,2:ncol(as.matrix(beta[,,1])),],dim=c(k1,k1*p+m*s,gran2))
    ## browser()
    if(dual){
        BB <- GamLoopSGLOODP(beta,INIactive,gamm,alpha,Y,ZZ,kk,jjfull,jjcomp,eps,YMean,ZMean,k1,p*k1+m*s,M2f,eigs,m)
    }else{
        
        ## browser()}
        BB <- GamLoopSGLOO(beta,INIactive,gamm,alpha,Y,ZZ,kk,jjfull,jjcomp,eps,YMean,ZMean,k1,p*k1+m*s,M2f,eigs,m)
        ## if(length(gamm)>1){browser()}
    }
    if(MN)
    {
        for(i in 1:(dim(BB$beta)[3]))
        {

            BB$beta[,2:ncol(BB$beta[,,i]),i] <- BB$beta[,2:ncol(BB$beta[,,i]),i]+C

        }

    }
   

    return(BB)
}



# Lag weighted lasso: VAR only

.lassoVARTL <- function (B, Z, Y, gamm, eps,p,MN,alpha,C1,intercept) 
{

    if(!"matrix"%in%class(Y))
    {
        Y <- matrix(Y,ncol=1)
    }
  
    k <- ncol(Y)

   if(MN)
   {
       C <- matrix(0,nrow=k,ncol=k*p)        
       diag(C) <- C1
       Y <- t(Y)
       Y <- Y-C%*%Z
       Y <- t(Y)

   }
    
    Y <- t(Y)
if(intercept){
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
    Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
}else{
    YMean <- rep(0,nrow(Y))
    ZMean <- rep(0,nrow(Z))
    }
    Y <- t(Y)
    tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))
    BFOO1 <- as.matrix(B[, 2:ncol(B[, , 1]), 1])
    BFOO <- array(B[,2:ncol(as.matrix(B[,,1])),],dim=c(k,k*p,length(gamm)*length(alpha)))
           p2 <- 1:p

    gran2 <- length(gamm)

    for(i in 1:length(alpha))
    {

        W <- rep(p2^(-alpha[i]),each=k)
        ZADJ <- diag(W)%*%Z

        ## B[,,(1+(i-1)*gran2):(i*length(gamm))] <- gamloopFista(array(BFOO[,,(1+(i-1)*gran2):(i*length(gamm))],dim=c(k,k*p,length(gamm))), Y, ZADJ, gamm, eps,as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k,p)
        B[,,(1+(i-1)*gran2):(i*length(gamm))] <- gamloopFista(array(BFOO[,,(1+(i-1)*gran2):(i*length(gamm))],dim=c(k,k*p,length(gamm))), Y, ZADJ, as.matrix(gamm), eps,as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k,p)


        for(j in (1+(i-1)*gran2):(i*length(gamm)))
        {

            B[,2:(k*p+1),j] <- B[,2:(k*p+1),j]%*%diag(W)

        }               


    }
            if(MN)
        {
            for(i in 1:(dim(B)[3]))
             
                B[,2:ncol(B[,,i]),i] <- B[,2:ncol(B[,,i]),i]+C

        }

## browser()
    return(B)
}


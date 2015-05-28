### Support functions for BIGVAR Package
### These are mostly utility functions that will not be seen by the user


# mean benchmark
.evalMean <- function(Y,T1,T2)
{
  MSFE <- rep(0,T2-T1-1)
  k <- ncol(Y)
 for (u in (T1 + 1):T2) {
        trainY1 <- Y[1:(u - 1), ]
        if(k>1){
        ypred <- colMeans(trainY1)}
        else{ypred=mean(trainY1)}
        uhat <- matrix(Y[u, ] - ypred, 
            ncol = k)
        MSFE[u - T1 ] <- norm2(uhat)^2
    }
  return(list(Mean=mean(na.omit(MSFE)),SD=sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE)))))
  }

# random walk benchmark
.evalRW <- function(Y,T1,T2)
{
  MSFE <- rep(0,T2-T1-1)
  k <- ncol(Y)
 for (u in (T1 + 1):T2) {
        trainY1 <- Y[u - 1, ]
        ## ypred <- colMeans(trainY1)
        uhat <- matrix(Y[u, ] - trainY1, 
            ncol = k)
        MSFE[u - T1 ] <- norm2(uhat)^2
    }
  return(list(Mean=mean(na.omit(MSFE)),SD=sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE)))))
  }


# univariate AR benchmark
.evalAR <- function(Y,T1,T2,p)
{
  k <- ncol(Y)
  MSFE <- rep(0,T2-T1-1)
uhat <- matrix(ncol=k,nrow=1)
for(u in (T1+1):T2)

          {
            
        trainY1 <- matrix(Y[1:(u - 1), ],ncol=k)
          mod1 <- ar(trainY1,aic=T)
          uhat <- matrix(Y[u,]-as.numeric(predict(mod1)$pred),ncol=k)
          MSFE[u-T1+1] <- norm2(uhat)^2
          }

return(list(Mean=mean(na.omit(MSFE)),SD=sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE)))))
}


# Constructs the Grid of Lambda Values

.LambdaGridE<- function (gran1, gran2, jj = jj, Y, Z, group,p,k) 
{
    if (group == "Lag") {
        mat = list()
        for (i in 1:length(jj)) {
            if(k>1){
            mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y)
            }
            else{mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])}
        }
        gamstart <- max(unlist(mat))
    }
    if (group == "None"|group=="Tapered") {
        gamstart = max(t(Y) %*% t(Z))
    }
    if (group == "SparseLag") {
        mat = list()
        for (i in 1:length(jj)) {
           if(k>1){
            mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y)*(1/(k+1))
            }
            else{mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])}
         }
        gamstart <- max(unlist(mat))
    }
    if (group == "Diag") {
        mat = list()
        ZZ <- kronecker(t(Z), diag(k))
        for (i in 1:length(jj)) {
            mat[[i]] <- norm(as.vector(t(Y)) %*% ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), 
                "F")
        }
        gamstart <- max(unlist(mat))
    }
    if (group == "SparseDiag") {
        mat = list()
        ZZ <- kronecker(t(Z), diag(k))
        for (i in 1:length(jj)) {
            mat[[i]] <- norm(1/(k + 1) * as.vector(t(Y)) %*% 
                ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), "F")
        }
        gamstart <- max(unlist(mat))
    }
    if (group == "HVARC"|group=="HVAROO"|group=="HVARELEM") {
        gmax <- c()
        for (i in 1:k) {
            gmax[i] = norm2(Z %*% Y[, i])
        }
        gamstart <- max(gmax)
     }
    if(group=="Tapered"){
    beta=array(0,dim=c(k,k*p+1,10))
    }else{beta=array(0,dim=c(k,k*p+1,1))}
    gamstart <- LGSearch(gamstart,Y,Z,beta,group,k,p,jj)

    
   
    gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), 
        length = gran2))
    return(gamm)
}

.LambdaGridX<- function (gran1, gran2, jj = jj, Y, Z, group,p,k1,s,m,k) 
{
    ## k = ncol(Y)
    if (group == "Lag") {
        mat = list()
        for (i in 1:length(jj)) {
            if(k>1){
            mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y)
            }
            else{mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])}
        }
        gamstart <- max(unlist(mat))
    }
    if (group == "None") {
        gamstart = max(t(Y) %*% t(Z[(k1*p+1):nrow(Z),]))
    }
    if (group == "SparseLag") {
        mat = list()
        for (i in 1:length(jj)) {
           if(k>1){
            mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y*(1/(k1+1)))
            }
            else{mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])}
         }
        gamstart <- max(unlist(mat))
    }
    if (group == "Diag") {
        mat = list()
        ZZ <- kronecker(t(Z), diag(k1))
        for (i in 1:length(jj)) {
            mat[[i]] <- norm(as.vector(t(Y)) %*% ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), 
                "F")
        }
        gamstart <- max(unlist(mat))
    }
    if (group == "SparseDiag") {
        mat = list()
        ZZ <- kronecker(t(Z), diag(k1))
        for (i in 1:length(jj)) {
            mat[[i]] <- norm(1/(k1 + 1) * as.vector(t(Y)) %*% 
                ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), "F")
        }
        gamstart <- max(unlist(mat))
    }
    if (group=="EFX") {
        gmax <- c()
        for (i in 1:k1) {
            gmax[i] = norm2(Z %*% Y[, i])/sqrt(k*p)
        }
        
        gamstart <- max(gmax)

     }
    

    
        beta <- array(0,dim=c(k1,k1*p+s*m+1,1))

        gamstart <- LGSearchX(gamstart,Y,Z,beta,group,k1,p,s,m,jj,k)

        
       gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), 
        length = gran2))
    return(gamm)
    }





# Forecast evaluation (called in cv.bigvar)
.EvalLVARX <- function(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m)
  {
    gran2=1	
    gamm <- gamopt
    Y <- ZFull$Y
    MSFE <- rep(0,T-T2-1)
    alpha=1/(k1+1)
    beta=array(0,dim=c(k1,k1*p+(k-k1)*s+1,1))  
        if (group == "Lag") {
        jj <- groupfunVARX(p,k,k1,s)
        jjcomp <- groupfunVARXcomp(p,k,k1,s)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)
    }
    if (group == "SparseLag") {
        jj <- groupfunVARX(p, k,k1,s)
        q1a <- list()
        for (i in 1:(p+s)) {
            q1a[[i]] <- matrix(runif(k1, -1, 1), ncol = 1)
        }
            activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)

      }
    if (group == "Diag") {
        kk <- diaggroupfunVARX(p, k,k1,s)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
    }
         if (group == "SparseDiag") {
        kk <- diaggroupfunVARX(p, k,k1,s)
       activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
      }

    if(verbose==TRUE){
	print("Evaluation Stage")
    pb <- txtProgressBar(min = T2+1, max = T, style = 3)}
    for (v in (T2+1):T) {
        trainY <- ZFull$Y[1:(v-1), ]
       trainZ <- ZFull$Z[,1:(v-1)]
        if (group == "None") {
            beta <- .lassoVARFistX(beta, trainZ, trainY,gamm, 1e-05,p,MN,k,k1,s,m)
        }
        if (group == "Lag") {
            GG <- GroupLassoVAR(beta,trainY,trainZ,gamm,1e-04,k,p,activeset,jj,jjcomp,k1,s)
            beta <- GG$beta
            activeset <- GG$active
        }
        if (group == "SparseLag") {
            GG <- .SparseGroupLassoVARX(beta, jj, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN,k,s,k1)
            beta <- GG$beta
            activeset = GG$active
            q1a <- GG$q1
        }
        if (group == "Diag") {
            GG <- .GroupLassoOOX(beta, kk, trainY, trainZ, gamm, 
                activeset, 1e-04,p,MN,k,k1,s)
            beta <- GG$beta
            activeset <- GG$active
        }
          if (group == "SparseDiag") {
 
            GG <- .SparseGroupLassoVAROOX(beta, kk, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04,p,MN,k1,s,k)
            beta <- GG$beta
            activeset = GG$active
        }

            if(group=="EFX")
            {

              beta <- .EFVARX(beta,trainY,trainZ,gamm,1e-4,MN,k1,s,m,p)

                }

        
    betaEVAL <- matrix(beta[,,1],nrow=k1,ncol=(k1*p+(k-k1)*s+1))
            if (RVAR == TRUE) {
                         betaEVAL <- RelaxedLS(cbind(t(trainZ),trainY),betaEVAL,k,p,k1,s)
       
             }

                
        if(MN==TRUE){ eZ <- matrix(ZFull$Z[,v],ncol=1)}
else{
       eZ <- matrix(c(1,ZFull$Z[,v]),ncol=1)
   }
   

if(h==1){
                if(MN==TRUE){MSFE[v-T2] <- norm2(Y[v, ] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2}

               else{MSFE[v-T2] <- norm2(ZFull$Y[v,1:k1] - betaEVAL %*% eZ)^2

                    
                }
}

      else{
            if(v+h>T){break}
            eZ2 <- as.matrix(.Zmat2(Y[(v - p):v, ], p, k)$Z, ncol = 1)
 
            B <- VarptoVar1MC(betaEVAL[,2:ncol(betaEVAL)],p,k)

            MSFE[v - (T2 )] <- norm2(Y[v+h,]-(B%^%h)%*%eZ2)^2
            
            }
            
        




        if(verbose==TRUE){

	    setTxtProgressBar(pb, v)
            }
        }
eZ <- as.matrix(.Zmat2(Y[(nrow(Y) - p):nrow(Y), ], p, k)$Z, ncol = 1)

        if (group == "None") {
            betaPred <- .lassoVARFistX(beta, ZFull$Z, ZFull$Y,gamm, 1e-05,p,MN,k,k1,s,m)
        }
        if (group == "Lag") {
                 GG <- GroupLassoVAR(beta,ZFull$Y,ZFull$Z,gamm,1e-04,k,p,activeset,jj,jjcomp,k1,s)
      
            betaPred <- GG$beta
        }
        if (group == "SparseLag") {
            GG <- .SparseGroupLassoVARX(beta, jj, ZFull$Y, ZFull$Z, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN,k,s,k1)
            betaPred <- GG$beta
        }
        if (group == "Diag") {
            GG <- .GroupLassoOOX(beta, kk, ZFull$Y, ZFull$Z, gamm, 
                activeset, 1e-04,p,MN,k,k1,s)
            betaPred <- GG$beta
        }
          if (group == "SparseDiag") {
            GG <- .SparseGroupLassoVAROOX(beta, kk, ZFull$Y, ZFull$Z, 
                gamm, alpha, INIactive = activeset, 1e-04,p,MN,k1,s,k)
            betaPred <- GG$beta
          
        }
            if(group=="EFX")
            {

              betaPred <- .EFVARX(beta,ZFull$Y,ZFull$Z,gamm,1e-4,MN,k1,s,m,p)

                }
            
## 	if(group=="BVAR")
## {
## #C <- .1*diag(rep(1,k*p+1))
## #V0 <- diag(rep(1,k))
## betaPred <- .BVARMNPrior(Y,Z,.2,.5,p,k)
## }
	## if(group!="BVAR")
## {
    betaPred <- as.matrix(betaPred[,,1])
## }
            
return(list(MSFE=MSFE,betaPred=betaPred,zvals=eZ))
 
    }


# Forecast evaluation (called in cv.bigvar)
.EvalLVAR <- function(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T1,T2)
  {
      
    gran2=1	
    gamm <- gamopt
    Y <- ZFull$Y
    s=p;k1=k
    MSFE <- rep(0,T2-T1-1)
    alpha=1/(ncol(Y)+1)
    beta=array(0,dim=c(k,p*k+1,1))  
        if (group == "Lag") {
        jj <- .groupfuncpp(p,k)
        jjcomp <- .groupfuncomp(p,k)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)
    }
    if (group == "SparseLag") {
        jj <- .groupfun(p, k)
        q1a <- list()
        for (i in 1:(p)) {
            q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
        }
            activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)

      }
    if (group == "Diag") {
        kk <- .lfunction3cpp(p, k)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
    }
         if (group == "SparseDiag") {
        ## ZZ <- kronecker(t(ZFull$Z),diag(k))
       kk <- .lfunction3cpp(p, k)
       jjcomp=.lfunctioncomp(p,k)
       jj=.lfunction3(p,k)
       activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
          q1a <- list()
        for (i in 1:(2*p)) {
            q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)
        }
      }

    if(verbose==TRUE){
	print("Evaluation Stage")
    pb <- txtProgressBar(min = T1+1, max = T2, style = 3)}
    for (v in (T1+1):T2) {
        trainY <- ZFull$Y[1:(v-1), ]
       trainZ <- ZFull$Z[,1:(v-1)]
        if (group == "None") {
            beta <- .lassoVARFist(beta, trainZ, trainY,gamm, 1e-05,p,MN)
        }
        if (group == "Lag") {
            GG <- GroupLassoVAR(beta,trainY,trainZ,gamm,1e-04,k,p,activeset,jj,jjcomp,k,s)
            beta <- GG$beta
            activeset <- GG$active
        }
        if (group == "SparseLag") {
            GG <- .SparseGroupLassoVAR(beta, jj, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)
            beta <- GG$beta
            activeset = GG$active
            q1a <- GG$q1
        }
        if (group == "Diag") {
            GG <- .GroupLassoOO(beta, kk, trainY, trainZ, gamm, 
                activeset, 1e-04,p,MN)
            beta <- GG$beta
            activeset <- GG$active
        }
          if (group == "SparseDiag") {
 
            GG <- .SparseGroupLassoVAROO(beta, kk, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04,q1a,p,MN)
            beta <- GG$beta
            activeset = GG$active
            q1a <- GG$q1
        }
                 if (group=="HVARC")
           {

		beta <- .HVARCAlg(beta,trainY,trainZ,gamm,1e-05,p,MN)
		

             }

        if(group=="HVAROO")
          {
            beta <- .HVAROOAlg(beta,trainY,trainZ,gamm,1e-05,p,MN)

          }
	if(group=="HVARELEM")
	{
	    beta<-.HVARElemAlg(beta,trainY,trainZ,gamm,1e-05,p,MN) 	
	}
    if(group=="Tapered")
        {

        beta <- .lassoVARTL(beta, trainZ, trainY,gamm, 1e-4,p,MN,palpha)

            }
        
    betaEVAL <- matrix(beta[,,1],nrow=k,ncol=(k*p+1))
            if (RVAR == TRUE) {
                         betaEVAL <- RelaxedLS(cbind(t(trainZ),trainY),betaEVAL,k,p,k1,s)
       
             }

                
        if(MN==TRUE){ eZ <- matrix(ZFull$Z[,v],ncol=1)}
else{
       eZ <- matrix(c(1,ZFull$Z[,v]),ncol=1)
   }
   

if(h==1){
              # We don't consider an intercept for the MN lasso
                if(MN==TRUE){MSFE[v-T1+1] <- norm2(Y[v, ] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2}

               else{MSFE[v-T1] <- norm2(ZFull$Y[v,] - betaEVAL %*% eZ)^2

                    
                }
}

      else{
            if(v+h>T2){break}
            eZ2 <- as.matrix(.Zmat2(Y[(v - p):v, ], p, k)$Z, ncol = 1)
            B <- VarptoVar1MC(betaEVAL[,2:ncol(betaEVAL)],p,k)
            MSFE[v - (T1)] <- norm2(Y[v+h,]-(B%^%h)%*%eZ2)^2            
            }
                    
        if(verbose==TRUE){

	    setTxtProgressBar(pb, v)
            }
        }
        if (group == "None") {
            betaPred <- .lassoVARFist(beta, ZFull$Z, ZFull$Y,gamm, 1e-05,p,MN)
        }
        if (group == "Lag") {
                 GG <- GroupLassoVAR(beta,ZFull$Y,ZFull$Z,gamm,1e-04,k,p,activeset,jj,jjcomp,k,s)
      
            betaPred <- GG$beta
        }
        if (group == "SparseLag") {
            GG <- .SparseGroupLassoVAR(beta, jj, ZFull$Y, ZFull$Z, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)
            betaPred <- GG$beta
        }
        if (group == "Diag") {
            GG <- .GroupLassoOO(beta, kk, ZFull$Y, ZFull$Z, gamm, 
                activeset, 1e-04,p,MN)
            betaPred <- GG$beta
        }
          if (group == "SparseDiag") {
            GG <- .SparseGroupLassoVAROO(beta, kk, ZFull$Y, ZFull$Z, 
                gamm, alpha, INIactive = activeset, 1e-04,q1a,p,MN)
            betaPred <- GG$beta
          
        }
                 if (group=="HVARC")
           {
		betaPred <- .HVARCAlg(beta,ZFull$Y,ZFull$Z,gamm,1e-05,p,MN)

             }
        if(group=="HVAROO")
          {
            betaPred <- .HVAROOAlg(beta,ZFull$Y,ZFull$Z,gamm,1e-05,p,MN)

          }
	if(group=="HVARELEM")
	{
	    betaPred<-.HVARElemAlg(beta,ZFull$Y,ZFull$Z,gamm,1e-05,p,MN) 	
	}
    if(group=="Tapered")
        {

        betaPred <- .lassoVARTL(beta,ZFull$Z,ZFull$Y,gamm,1e-4,p,MN,palpha)
            
        }    
            
## 	if(group=="BVAR")
## {
## #C <- .1*diag(rep(1,k*p+1))
## #V0 <- diag(rep(1,k))
## betaPred <- .BVARMNPrior(Y,Z,.2,.5,p,k)
## }
	## if(group!="BVAR")
## {
    betaPred <- as.matrix(betaPred[,,1])
## }
            
return(list(MSFE=MSFE,betaPred=betaPred,zvals=eZ))
 
    }



#' Converts a VAR coefficient matrix of order p to multiple companion form
#' 
#' @param B a \eqn{k \times kp} coefficient matrix
#' @param p Lag order
#' @param k Number of Series
#' @return Returns a \eqn{kp \times kp} coefficient matrix representing all coefficient matrices contained in Ai as a VAR(1).
#' @references See page 15 of Lutkepohl, "A New Introduction to Multiple Time Series Analysis"
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
VarptoVar1MC <- function(B,p,k)
  {
Fp=matrix(0,nrow=k*p,ncol=k*p)      
Fp[1:k,] = B
Fp[-(1:k),1:(k*(p-1))] = diag(k*(p-1))
if(max(Mod(eigen(Fp)$values))>1){warning("Coefficient Matrix is not stationary")}
return(Fp)
}
    
# Simulate a VAR with p=1 
# k=number of time series
# A1 = Coefficient Matrix
# Sigma = Residual Matrix
# n=number of simulations

#' Simulate a VAR
#' 
#' @param k Number of Series
#' @param A1 Either a \eqn{k \times k} coefficient matrix or a \eqn{kp \times kp} matrix created using \code{\link{VarptoVar1MC}}. 
#' @param p Maximum Lag Order
#' @param Sigma Residual Covariance Matrix of dimension \eqn{k\times k}
#' @param n Number of simulations
#' @return Returns a \eqn{n \times k} of realizations from a VAR.
#' @references Lutkepohl, "A New Introduction to Multiple Time Series Analysis"
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
 MultVarSim <- function (k, A1, p, Sigma, n) 
{
    if(max(Mod(eigen(A1)$values))>1){stop("Error: Generator Matrix is not stationary")}

    Y <- matrix(0, nrow = n+500+p , ncol = k)
    YY <- as.vector(Y)
    for (i in seq(from = (k * p + 1), to = (nrow(Y) * k - 1), 
        by = k)) {
        u <- as.vector(c(mvrnorm(1, rep(0, k), Sigma), rep(0, 
            k * p - k)))
        YY[(i + k):(i - k * p + 1 + k)] <- A1 %*% YY[(i):(i - 
            k * p + 1)] + as.matrix(u, ncol = 1)
    }
    YY <- YY[YY!=0]
    Y <- matrix(YY, ncol = k, byrow = TRUE)
    Y <- Y[,c(ncol(Y):1)]
    Y <- Y[501:nrow(Y), ]
    return(Y)
}
# create Z matrix without intercept, used for modeling  These functions are in C++ now.
.Zmat2 <- function(Y,p2,k)
  {
if(class(Y)!="matrix")
    {
Y <- matrix(Y,ncol=1)
    }
T <- nrow(Y)
Zcol <- function(q)
  {
if(ncol(Y)>1)
    {
YPS <- as.vector(t(apply(as.matrix(Y[(1+q):(p2+q),]),2,rev)))
}
else{YPS <- as.vector(t(rev(as.matrix(Y[(1+q):(p2+q),]))))
}
return(YPS)
}
q=0:(T-1-(2*p2-p2))
YPS2 <- sapply(X=q,FUN=Zcol)


return(list(Z=YPS2,Y=as.matrix(Y[(p2+1):nrow(Y),])))
}

.ZScalar <- function(Y,p)
  {
Z <- c()

for(i  in 1:p)
  {
    Z1=Y[(p-i+1):(length(Y)-i)]
    Z <- cbind(Z,Z1)

    }
return(list(Y= as.matrix(Y[(p+1):length(Y)]),Z=as.matrix(Z)))
}

.ZScalar2 <- function(Y,p)
  {
Z <- c()

for(i  in 1:p)
  {
    Z1=c(1,Y[(p-i+1):(length(Y)-i)])
    Z <- cbind(Z,Z1)

    }
return(list(Y= matrix(Y[(p+1):length(Y)],ncol=1),Z=as.matrix(Z)))
}


# create z matrix with an intecept, used in forecast evaluation
.Zmat <- function(Y,p2,k)
{
if(class(Y)!="matrix")
    {Y <- matrix(Y,ncol=1)}
Z <- matrix(nrow=k*p2+1,ncol=(nrow(Y)-p2))
if(k>1)
    {
YPS <- c(1,as.vector(t(apply(as.matrix(Y[1:p2,]),2,rev))))
}
else{
YPS <- c(1,as.vector(t(rev(Y[1:(p2),]))))
    }
q=1
while(q<(nrow(Y)-(2*p2-p2)))
{
  if(k>1){  
  Y2 <- c(1,as.vector(t(apply(as.matrix(Y[(1+q):(p2+q),]),2,rev))))
}
  else{
   Y2 <- c(1,as.vector(t(rev(Y[(1+q):(p2+q),]))))   
      }
    YPS <- cbind(YPS,Y2)
  q=q+1
    }
return(list(Z=YPS,Y=as.matrix(Y[(p2+1):nrow(Y),])))

}

# function to create subsets for group lasso-will need to create a new one for each group structure- this is for the block group structure
.groupfun <- function(p,k=nrow(beta))
  {
    jjj <- list()
    jjj[[1]] <- 1:k
    if(p>1){
    for(i in 2:p){
      jjj[[i]] <- jjj[[i-1]]+k
}
}
return(jjj)
    }

.groupfuncpp <- function(p,k=nrow(beta))
  {
    jjj <- list()
    jjj[[1]] <- 0:(k-1)
    if(p>1)
    {
    for(i in 2:p){
      jjj[[i]] <- jjj[[i-1]]+k
}
}
return(jjj)
    }

# subsetting groups in rcpp
.groupfuncomp <- function(p,k)
  {

ownoth <- .groupfuncpp(p,k)

kk2=list()
pmax <- max(unlist(ownoth))
to <- 0:(pmax)
for(i in 1:length(ownoth))
  {
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]

    }


return(kk2)
}    


.lfunction2 <- function(p,k)
  {
kk=list()
kk[[1]] <- 1:(k^2)
if(p>1)
{
for(i in 2:p)
  {
    kk[[i]] <- 1:(k^2)+tail(kk[[i-1]],1)

    }
}    
return(kk)

}

.lfunction2cpp <- function(p,k)
  {
kk=list()
kk[[1]] <- 0:(k^2-1)
if(p>1)
{
for(i in 2:p)
  {
    kk[[i]] <- 0:(k^2-1)+tail(kk[[i-1]],1)+1

    }
}    
return(kk)

}

.lfunction3 <- function(p,k)
  {
kk <- .lfunction2(p,k)
oo <- list()
pp <- list()
for(i in 1:length(kk))
{
  j=0
  oo[[i]] <- kk[[i]][(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
  pp[[i]] <- kk[[i]][-(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
  j=j+1
}

ownoth <- c(oo,pp)

return(ownoth)
}    

.lfunction3cpp <- function(p,k)
  {
kk <- .lfunction2cpp(p,k)
oo <- list()
pp <- list()
for(i in 1:length(kk))
{
  j=0
  oo[[i]] <- kk[[i]][(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
  pp[[i]] <- kk[[i]][-(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
  j=j+1
}

ownoth <- c(oo,pp)

return(ownoth)
}    


.lfunctioncomp <- function(p,k)
  {

ownoth <- .lfunction3cpp(p,k)

kk2=list()
pmax <- max(unlist(ownoth))
to <- 0:(pmax)
for(i in 1:length(ownoth))
  {
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]

    }


return(kk2)
}    


# This function should work for arbitrary groups
.lfunction <- function(groups,p)
  {
     H <- as.vector(do.call('cbind',groups))
     kk <- list()
     kk[[1]] <- H
     if(p>1)
     {
    for(i in 2:p){
      kk[[i]] <-as.vector(do.call('cbind',groups))+tail(kk[[i-1]],1) }
      }
   return(kk)

    }
# Indexing for HVAR
.vsubs <- function(p,k)
    {
vi <- list()
for(i in p:1)
{
g <- max(k*i-k,0)
vi[[i]] <- g:(k*(p)-1)
    }
 return(vi)   
}

#indexing for HVAR OO
#lfunction3(p,k)
#oofun(p,k)
# new group function in hierarchical order
.oofun <- function(p,k)
 {
    kk <- .lfunction2(p, k)
    oo <- list()
    pp <- list()
    for (i in 1:length(kk)) {
        j = 0
        oo[[i]] <- kk[[i]][(seq(1, length(kk[[1]]), k + 1) + 
            (j * k^2))]
        pp[[i]] <- kk[[i]][-(seq(1, length(kk[[1]]), k + 1) + 
            (j * k^2))]
        j = j + 1
    }
    ownoth <- list()
   jj <- .lfunction3(p,k)
   for(i in 1:length(jj))
     {
       if(i==1)
         {
           ownoth[[i]] <- jj[[i]]
           oo[[1]] <- NULL
           }
       if(i==length(jj))
         {
           ownoth[[i]] <- tail(jj,1)
           pp[[1]] <- NULL
         }
       if(i!=1& i%%2!=0)
         {
           ownoth[[i]] <- head(oo,1)
           oo[[1]] <- NULL
         }
      if(i!=length(jj)&i%%2==0)
        {
          ownoth[[i]] <- head(pp,1)
           pp[[1]] <- NULL
        }

     }
    return(rev(ownoth))
   }
.oocumfun <-function(p,k)
{
kk <- rev(.oofun(p,k))
oogroups=list()
oogroups[[1]] <- unlist(kk)
for(i in 2:length(kk))
  {
    oogroups[[i]] <- unlist(kk[-(1:(i-1))])

  }
return(oogroups)
}


# indexing function for variable lag own/other
.vecoovars<-function(p,k,k1)
{
vv=list()
vv[[1]]=1:(p*k)
vv[[2]]=vv[[1]][-k1]
q1=1
for(i in 3:(2*p))
{
if(i%%2!=0)
{
vv[[i]]=(q1*k+1):(k*p)
q1=q1+1
}
else{
vv[[i]]=vv[[i-1]][-k1]
}
}
return(vv)
}

# indexing to start at zero for use within rcpp
.vecoovarscpp<-function(p,k,k1)
{
vv=list()
vv[[1]]=0:(p*k-1)
vv[[2]]=vv[[1]][-(k1)]
q1=1
for(i in 3:(2*p))
{
if(i%%2!=0)
{
vv[[i]]=(q1*k):(k*p-1)
q1=q1+1
}
else{
vv[[i]]=vv[[i-1]][-(k1)]
}
}
return(vv)
}



# Generates prior covariance matrix for Bayesian VAR

## .SigmaGen <- function(phi0,phi1,j)
##     {
## s <- apply(Y,2,sd)

## S1 <- c()
## s2 <- 1:(k*p)
## s3 <- seq(j,k*p,k)
## s4 <- (1:k)[-j]
## ppp=1
## for(i in s3 )
## {
## S1[i] <- phi0/ppp
## ppp=ppp+1
## }
## bb=1
## for( q in s4)
##     {
## pp=1
## for(i in s2[-s3][bb:(bb+p)])
## {
## S1[i] <- phi0*phi1*(s[q]/(s[j]))^2*1/pp^2
## pp=pp+1
##     }
## bb=bb+1
## }
## return(S1)
## }

# Sparsity Plot for VAR models
SparsityPlot <- function(B,p,k,title=NULL)
    {
text <- c()
for(i in 1:p)
    {
        text1 <- as.expression(bquote(bold(B)^.(i)))
        text <- append(text,text1)
          
    }
f <- function(m) t(m)[,nrow(m):1]

rgb.palette <- colorRampPalette(c("white", "blue" ),space = "Lab")


at <- seq(k/2+.5,p*(k)+.5,by=k)

se2=seq(1.75,by=k,length=k)

L2 <- levelplot(f(abs(B)),col.regions=rgb.palette,colorkey=NULL,xlab=NULL,ylab=NULL,main=list(label=title,cex=1),panel=function(...){ panel.levelplot(...);panel.abline(a=NULL,b=1,h=seq(1.5,p*k+.5,by=1),v=seq(1.5,by=1,length=p*k));panel.abline(a=NULL,b=1,v=seq(k+.5,p*k+.5,by=k),lwd=3)},scales=list(x=list(alternating=1,labels=text,cex=1,at=at,tck=c(0,0)),y=list(alternating=0,tck=c(0,0))))

return(L2)

        }

# grouping functions for varX

groupfunVARX <- function(p,k,k1,s)
    {
         jj=list()
        m=k-k1
        jj <- .groupfuncpp(p, k1)
        kp <- k1*p+m*s-1
        jj2 <- list()
        startjj=max(unlist(jj))+1
        for(i in seq(startjj,kp,by=1))
            {
                jj[[i]] <- i
                
            }
        jj[sapply(jj, is.null)] <- NULL

          return(jj)
      }

groupfunVARXcomp <- function(p,k,k1,s)
    {
ownoth <- groupfunVARX(p,k,k1,s)

kk2=list()
pmax <- max(unlist(ownoth))
to <- 0:(pmax)
for(i in 1:length(ownoth))
  {
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]

    }

return(kk2)
}
        
groupfunVARXLG <- function(p,k,k1,s)
    {
         jj=list()
        jj <- .groupfun(p, k1)
        m=k-k1
        kp <- k1*p+m*s
        jj2 <- list()
        startjj=max(unlist(jj))+1
        for(i in seq(startjj,kp,by=1))
            {
                jj[[i]] <- i
                
            }
        jj[sapply(jj, is.null)] <- NULL

          return(jj)
      }



diaggroupfunVARX <- function(p,k,k1,s)
    {
        m=k-k1
         jj=list()
        jj <- .lfunction3cpp(p, k1)
        kp <-k1*(p*k1+s*m)-1
        jj2 <- list()
        startjj=max(unlist(jj))+1
        for(i in seq(startjj,kp,by=k1))
            {
                jj[[i]] <- i:(i+k1-1)
                
            }
        jj[sapply(jj, is.null)] <- NULL

          return(jj)
      }

diaggroupfunVARXcomp <- function(p,k,k1,s)
    {
ownoth <- diaggroupfunVARX(p,k,k1,s)

kk2=list()
pmax <- max(unlist(ownoth))
to <- 0:(pmax)
for(i in 1:length(ownoth))
  {
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]

    }


return(kk2)
}



        
        
diaggroupfunVARXLG <- function(p,k,k1,s)
    {
        m=k-k1
         jj=list()
        jj <- .lfunction3(p, k1)
        kp <- k1*(p*k1+s*m)
        jj2 <- list()
        startjj=max(unlist(jj))+1
        for(i in seq(startjj,kp,by=k1))
            {
                jj[[i]] <- i:(i+(k1-1))
                
            }
        jj[sapply(jj, is.null)] <- NULL

          return(jj)
      }



diaggroupfunVARXLGL <- function(p,k,k1)
    {
         jj=list()
        jj <- .lfunction3(p, k1)
        kp <- k1*p*k
        jj2 <- list()
        startjj=max(unlist(jj))+1
        for(i in seq(startjj,kp,by=1))
            {
                jj[[i]] <- i
                
            }
        jj[sapply(jj, is.null)] <- NULL

          return(jj)
      }


diaggroupfunVARXL <- function(p,k,k1)
    {
         jj=list()
        jj <- .lfunction3cpp(p, k1)
        kp <- k1*p*k-1
        jj2 <- list()
        startjj=max(unlist(jj))+1
        for(i in seq(startjj,kp,by=1))
            {
                jj[[i]] <- i
                
            }
        jj[sapply(jj, is.null)] <- NULL

          return(jj)
      }

diaggroupfunVARXcompL <- function(p,k,k1)
    {
ownoth <- diaggroupfunVARXL(p,k,k1)

kk2=list()
pmax <- max(unlist(ownoth))
to <- 0:(pmax)
for(i in 1:length(ownoth))
  {
    kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]

    }


return(kk2)
}
# iterative procedure to find a less coarse bound for lambda starting value
LGSearchX <- function(gstart,Y,Z,BOLD,group,k1,p,s,m,gs,k)
    {
     tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
     lambdah <- gstart
     lambdal <- 0
    activeset <- list(rep(rep(list(0), length(gs))))

        while(max(lambdah-lambdal)>.00001)
            {
                lambda <- (lambdah+lambdal)/2
                if(group=="EFX"){
                BOLD <- .EFVARX(BOLD,Y,Z,lambda,1e-4,FALSE,k1,s,m,p)
                param=BOLD[,2:(k1*p+m*s+1),1]
            }
            if(group=="None"){
                param <- .lassoVARFistX(BOLD,Z,Y[,1:k1],lambda,1e-04,p,FALSE,k1+m,k1,s,m)[,2:(k1*p+m*s+1),]
            }
            if(group=="Lag"){
                  jj <- groupfunVARX(p,k,k1,s)
                  jjcomp <- groupfunVARXcomp(p,k,k1,s)
                  # need to use old Group Lasso function ; 
                  BB <- .GroupLassoVAR1(BOLD,jj,jjcomp,Y[,1:k1],Z,lambda,activeset,1e-4,p,FALSE,k,k1,s)
                  BOLD <- BB$beta
                  param=BB$beta[,2:(k1*p+m*s+1),]
                  activeset <- BB$active
                }
              if(group=="Diag")
                  {
                  kk <- diaggroupfunVARX(p, k,k1,s)
                  BB <- .GroupLassoOOX(BOLD, kk, Y, Z, lambda,activeset, 1e-04,p,FALSE,k,k1,s)
                  param=BB$beta[,2:(k1*p+m*s+1),]
                  BOLD=BB$beta
                  activeset=BB$active
                      }
         if(group=="SparseDiag")
                  {
                  kk <- diaggroupfunVARX(p, k,k1,s)
                  BB <- .SparseGroupLassoVAROOX(BOLD, kk, Y[,1:k1], Z, lambda,1/(k1+1),activeset, 1e-04,p,FALSE,k1,s,k)
                  param=BB$beta[,2:(k1*p+m*s+1),]
                  BOLD=BB$beta
                  activeset=BB$active
                      }
     

          if(group=="SparseLag"){
                  jj <- groupfunVARX(p, k,k1,s)
                  jjcomp <- groupfunVARXcomp(p,k,k1,s)
                  q1a=list()
                  for (i in 1:(p+s)) {
            q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
                   }
      BB <- .SparseGroupLassoVARX(BOLD,jj,Y[,1:k1],Z,lambda,1/(k1+1),activeset,1e-4,q1a,p,FALSE,k,s,k1)
                  param=BB$beta[,2:(k1*p+m*s+1),]
                  BOLD <- BB$beta
                  activeset=BB$active
                }
  

                if(max(abs(param))==0)
                   {
                     lambdah <- lambda}
                 else{lambdal <- lambda
                  
                    }

            }


lambdah
        }

# Grid search for VAR model starting values
LGSearch <- function(gstart,Y,Z,BOLD,group,k,p,gs)
    {
     tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
     lambdah <- gstart
     lambdal <- 0
activeset <- list(rep(rep(list(0), length(gs))))

        if(group=="SparseDiag")
                  {
                  kk <- .lfunction3cpp(p, k)
                  activeset <- rep(list(rep(rep(list(0), length(kk)))),1)
                                    q1a=list()
                  for (i in 1:(2*p)) {
            q1a[[i]] <- matrix(runif(length(gs[[i]]), -1, 1), ncol = 1)
                   }
}

        while(max(abs(lambdah-lambdal))>.00001)
            {
                lambda <- (lambdah+lambdal)/2
            if(group=="None"){
                param <- .lassoVARFist(BOLD,Z,Y,lambda,1e-04,p,FALSE)[,2:(k*p+1),]
            }
           if(group=="Tapered"){
                param <- .lassoVARTL(BOLD,Z,Y,lambda,1e-04,p,FALSE,rev(seq(0,1,length=10)))[,2:(k*p+1),]
            }
 
            if(group=="Lag"){
                  jj <- .groupfuncpp(p, k)
                  jjcomp <- .groupfuncomp(p,k)
                  # need to use old Group Lasso function ALSO need to specify groups outside of function 
                  BB <- .GroupLassoVAR1(BOLD,jj,jjcomp,Y,Z,lambda,activeset,1e-4,p,FALSE,k,k,p)
                  BOLD <- BB$beta
                  param=BB$beta[,2:(k*p+1),]
                  activeset <- BB$active
                }
            if(group=="SparseLag"){
                  jj <- .groupfuncpp(p, k)
                  jjcomp <- .groupfuncomp(p,k)
                  q1a=list()
                  for (i in 1:p) {
            q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
                   }
    BB <- .SparseGroupLassoVAR(BOLD,jj,Y,Z,lambda,1/(k+1),activeset,1e-4,q1a,p,FALSE)
                  param=BB$beta[,2:(k*p+1),]
                  BOLD <- BB$beta
                  activeset=BB$active
                }
              if(group=="Diag")
                  {
                  kk <- .lfunction3cpp(p, k)
                  BB <- .GroupLassoOO(BOLD, kk, Y, Z, lambda,activeset, 1e-04,p,FALSE)
                  param=BB$beta[,2:(k*p+1),]
                  BOLD=BB$beta
                  activeset=BB$active
                      }
                  if(group=="SparseDiag")
                  {

                  BB <- .SparseGroupLassoVAROO(BOLD, kk, Y, Z, lambda,1/(k+1),activeset, 1e-04,q1a,p,FALSE)
                  param=BB$beta[,2:(k*p+1),]
                  BOLD=BB$beta
                  activeset=BB$active
                
                       }
          if(group=="HVARC")
              {
                  BOLD <- .HVARCAlg(BOLD,Y,Z,lambda,1e-4,p,FALSE)
                  param=BOLD[,2:(k*p+1),]
                  }
          if(group=="HVAROO")
              {
                  BOLD <- .HVAROOAlg(BOLD,Y,Z,lambda,1e-4,p,FALSE)
                  param=BOLD[,2:(k*p+1),]
                  }
          if(group=="HVARELEM")
              {
                  BOLD <- .HVARElemAlg(BOLD,Y,Z,lambda,1e-4,p,FALSE)
                  param=BOLD[,2:(k*p+1),]
                  }
              if(max(abs(param))==0)
                   {
                     lambdah <- lambda}
                 else{lambdal <- lambda
                  
                    }

            }


lambdah
        }





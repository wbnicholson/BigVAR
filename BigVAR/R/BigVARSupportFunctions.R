### Support functions for BIGVAR Package
### These are mostly utility functions that will not be seen by the user
### Updated June 19, 2014

### For mean, random walk, and AIC benchmarks:

.evalMean <- function(Y,T1,T2)
{
  MSFE <- c()
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


.evalRW <- function(Y,T1,T2)
{
  MSFE <- c()
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

.evalBIC <- function(Y,T1,T2,verbose)
    {
   MSFE <- c()
   k <- ncol(Y)
  uhat <- matrix(ncol = k, nrow = 1)
    order <- c()
    if(verbose==TRUE){
    print("EvalBIC")}
   if(is.null(colnames(Y))==TRUE){colnames(Y)=1:k}
      pb <- txtProgressBar(min = T1+1, max = T2, style = 3)
    for (u in (T1 + 1):T2) {
        trainY1 <- matrix(Y[1:(u - 1), ],ncol=k)
	if(is.null(colnames(trainY1))==TRUE){colnames(trainY1)=1:k}

        if(dim(trainY1)[1]>(dim(trainY1)[2]+10))
        {
        mod1 <- (VAR(trainY1, ic="FPE",lag.max=4,p=4))   
        p1 <- predict(na.omit(mod1),n.ahead=1)$fcst
        ypred <- matrix(sapply(p1, "[[", 1))
        uhat <- matrix(Y[u, ] - ypred, 
            ncol = k)
        MSFE[u - T1+1] <- norm2(uhat)^2
     #     if(MSFE[u-T1+1]>100){break}
       }
        else{MSFE[u-T1+1]=NA}
        if(verbose==TRUE){
      	setTxtProgressBar(pb, u)}
    }
 #   print(mean(na.omit(MSFE)))
  return(list(Mean=mean(na.omit(MSFE)),SD=sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE)))))
  
}


#' @importFrom vars VAR
.evalAIC <- function(Y,T1,T2,verbose)
    {
   MSFE <- c()
   k <- ncol(Y)
  uhat <- matrix(ncol = k, nrow = 1)
    order <- c()
    if(verbose==TRUE){
    print("EvalAIC")}
   if(is.null(colnames(Y))==TRUE){colnames(Y)=1:k}
      pb <- txtProgressBar(min = T1+1, max = T2, style = 3)
    for (u in (T1 + 1):T2) {
        trainY1 <- matrix(Y[1:(u - 1), ],ncol=k)
	if(is.null(colnames(trainY1))==TRUE){colnames(trainY1)=1:k}

        if(dim(trainY1)[1]>(dim(trainY1)[2]+10))
        {
        mod1 <- (VAR(trainY1, ic="AIC"))   
        p1 <- predict(mod1,n.ahead=1)$fcst
        ypred <- matrix(sapply(p1, "[[", 1))
        uhat <- matrix(Y[u, ] - ypred, 
            ncol = k)
        MSFE[u - T1+1] <- norm2(uhat)^2
     #     if(MSFE[u-T1+1]>100){break}
       }
        else{MSFE[u-T1+1]=NA}
        if(verbose==TRUE){
      	setTxtProgressBar(pb, u)}
    }
  return(list(Mean=mean(na.omit(MSFE)),SD=sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE)))))

  }

.evalAR <- function(Y,T1,T2,p)
{
  k <- ncol(Y)
  MSFE <- c()
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

.LambdaGrid<- function (gran1, gran2, jj = jj, Y, Z = Z, group,p) 
{
    k = ncol(Y)
    if (group == "Group") {
        mat = list()
        for (i in 1:length(jj)) {
            if(k>1){
            mat[[i]] <- norm2(t(Y) %*% t(Z[jj[[i]], ]))
            }
            else{mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])}
        }
        gamstart <- max(unlist(mat))
    }
    if (group == "None") {
        gamstart = max(t(Y) %*% t(Z))
    }
    if (group == "Sparse") {
        mat = list()
        for (i in 1:length(jj)) {
           if(k>1){
            mat[[i]] <- norm2(t(Y) %*% t(Z[jj[[i]], ]))*(1/(k+1))
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

   
      
    gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), 
        length = gran2))
    return(gamm)
}



# Forecast evaluation (called in cv.bigvar)
.EvalLVAR <- function(Y,Z,gamopt,k,p,group,h,MN,verbose,RVAR)
  {
    gran2=1	
    gamm <- gamopt
    T1 <- floor(2 * nrow(Y)/3)
    T2 <- nrow(Y)
    MSFE <- c()
    alpha=1/(ncol(Y)+1)	
    beta=array(0,dim=c(k,k*p+1,1))
    ZFull <-.Zmat2(Y, p, k)
    if (group == "Group") {
        jj <- .groupfuncpp(p, k)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)
    }
    if (group == "Sparse") {
        jj <- .groupfuncpp(p, k)
        q1a <- list()
        for (i in 1:p) {
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
        kk <- .lfunction3cpp(p, k)
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
        trainY <- ZFull$Y[1:(v-p-1), ]
       trainZ <- ZFull$Z[,1:(v-p-1)]
        if (group == "None") {
            ## if(parallel==TRUE){
            beta <- .lassoVARFist(beta, trainZ, trainY,gamm, 1e-05,p,MN)
            ## }
            ## else{
            ## beta <- lassoVAROMP(beta, trainZ, trainY, gamm, 1e-04)
            ## }
        }
        if (group == "Group") {
            GG <- .GroupLassoVAROpt(beta, jj, trainY, trainZ, gamm, 
                activeset, 1e-04,p,MN)
            beta <- GG$beta
            activeset <- GG$active
        }
        if (group == "Sparse") {
            GG <- .SparseGroupLassoVAROpt(beta, jj, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)
            beta <- GG$beta
            activeset = GG$active
            q1a <- GG$q1
        }
        if (group == "Diag") {
            GG <- .GroupLassoOONew(beta, kk, trainY, trainZ, gamm, 
                activeset, 1e-04,p,MN)
            beta <- GG$beta
            activeset <- GG$active
        }
          if (group == "SparseDiag") {
            GG <- .SparseGroupLassoVAROptOO(beta, kk, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)
            beta <- GG$beta
            activeset = GG$active
            q1a <- GG$q1
        }
     
    betaEVAL <- matrix(beta[,,1],nrow=k,ncol=(k*p+1))
            if (RVAR == TRUE) {
                ## betaEVAL <- .relaxedLeastSquares(trainY, betaEVAL, 
                ##   trainZ,k,p)

                betaEVAL <- .relaxedVAR(trainY, betaEVAL, 
                  trainZ,k,p)
            
          
            }
        if(MN==TRUE){ eZ <- as.matrix(.Zmat2(Y[(v - p):v, ], p, k)$Z, ncol = 1)}
else{
        eZ <- as.matrix(.Zmat(Y[(v - p):v, ], p, k)$Z, ncol = 1)}
 

if(h==1){

                if(MN==TRUE){MSFE[v-T1+1] <- norm2(Y[v, ] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2}

               else{MSFE[v-T1+1] <- norm2(Y[v, ] - betaEVAL %*% eZ)^2}
}

      else{
            if(v+h>T1){break}
            eZ2 <- as.matrix(.Zmat2(Y[(v - p):v, ], p, k)$Z, ncol = 1)
            betaP <- betaEVAL
            nu=betaP[,1]
            betaP <- betaP[,2:ncol(betaP)]
            Bi=list()
            jj <- .groupfun(p,k)
            for(i in 1:p)
    {
        Bi[[i]] <- betaPred[,jj[[i]]]


    }
            B <- VarptoVar1(Bi,p,k)

            MSFE[v - (T1 - 1)] <- norm2(Y[v+h,]-(B%^%h)%*%eZ2)^2
            
            }
                

        if(verbose==TRUE){

	    setTxtProgressBar(pb, v)
            }
        }
eZ <- as.matrix(.Zmat2(Y[(nrow(Y) - p):nrow(Y), ], p, k)$Z, ncol = 1)

        if (group == "None") {

            betaPred <- .lassoVARFist(beta, Z, Y,gamm, 1e-05,p,MN)
        }
        if (group == "Group") {
            GG <- .GroupLassoVAROpt(beta, jj, Y, Z, gamm, 
                activeset, 1e-04,p,MN)
            betaPred <- GG$beta
        }
        if (group == "Sparse") {
            GG <- .SparseGroupLassoVAROpt(beta, jj, Y, Z, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)
            betaPred <- GG$beta
        }
        if (group == "Diag") {
            GG <- .GroupLassoOONew(beta, kk, Y, Z, gamm, 
                activeset, 1e-04,p,MN)
            betaPred <- GG$beta
        }
          if (group == "SparseDiag") {
            GG <- .SparseGroupLassoVAROptOO(beta, kk, Y, Z, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)
            betaPred <- GG$beta
          
        }



    betaPred <- as.matrix(betaPred[,,1])

return(list(MSFE=MSFE,betaPred=betaPred,zvals=eZ))
  }

#' Converts a VAR of order p to a VAR of order 1
#' 
#' @param Ai A list containing p k x k coefficient matrices
#' @param p Lag order
#' @param k Number of Series
#' @return Returns a \eqn{kp \times kp} coefficient matrix representing all matrices contained in Ai as a VAR(1).
#' @references See page 15 of Lutkepohl, "A New Introduction to Multiple Time Series Analysis"
#' @seealso \code{\link{MultVarSim}}
#' @examples
#' library(MASS)
#' A1 <- matrix(c(.4,-.02,.01,-.02,.3,.02,.01,.04,.3),ncol=3,nrow=3)
#' A2 <- matrix(c(.2,0,0,0,.3,0,0,0,.13),ncol=3,nrow=3)
#' Ai=list()
#' Ai[[1]]=A1
#' Ai[[4]]=A2
#' A <- VarptoVar1(Ai,6,3)
#' @export
VarptoVar1 <- function(Ai,p,k)
  {
Y <- matrix(0,nrow=k*p,ncol=1)
Y[1,] <- 0
A <- matrix(0,nrow=k,ncol=k*p)
A[1:nrow(Ai[[1]]),1:ncol(Ai[[1]])] <- Ai[[1]]


for(i in 2:length(Ai))
  {
    if(is.null(Ai[[i]])==TRUE)
    {Ai[[i]] <- matrix(0,nrow=k,ncol=k)}
    else{
    A[1:nrow(Ai[[1]]),1:ncol(Ai[[1]])+(i-1)*ncol(Ai[[1]])] <- Ai[[i]]  
    }
    }
d <- diag((k*p-k))
d1 <- matrix(0,nrow=nrow(d),ncol=k)
d <- cbind(d,d1)
A <- rbind(A,d)
return(A)
}

# Simulate a VAR with p=1 
# k=number of time series
# A1 = Coefficient Matrix
# Sigma = Residual Matrix
# n=number of simulations

#' Simulate a VAR
#' 
#' @param k Number of Series
#' @param A1 Either a \eqn{k \times k} coefficient matrix or a \eqn{kp \times kp} matrix created using \code{VarptoVar1}. 
#' @param p Maximum Lag Order
#' @param Sigma Residual Coariance Matrix of dimension \eqn{k\times k}
#' @param n Number of simulations
#' @return Returns a \eqn{n \times k} of realizations from a VAR.
#' @references Lutkepohl, "A New Introduction to Multiple Time Series Analysis"
#' @seealso \code{\link{VarptoVar1}}
#' @examples
#' A1 <- matrix(c(.4,-.02,.01,-.02,.3,.02,.01,.04,.3),ncol=3,nrow=3)
#' A2 <- matrix(c(.2,0,0,0,.3,0,0,0,.13),ncol=3,nrow=3)
#' Ai=list()
#' Ai[[1]]=A1
#' Ai[[4]]=A2
#' k=6;p=3
#' A <- VarptoVar1(Ai,k,p)
#' Y <-MultVarSim(k,A,p,.1*diag(k),100)
#' @export
#' @importFrom MASS mvrnorm
 MultVarSim <- function (k, A1, p, Sigma, n) 
{
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
    Y <- Y[500:nrow(Y), ]
    return(Y)
}
# create Z matrix without intercept, used for modeling
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


.relaxedLeastSquares <- function(Y,B2,Z,k,p)
{
sds<- apply(Y,2,sd)
sds2 <- rep(sds,p)

if(all(as.numeric(B2[,2:ncol(B2)])==0)){return(B2)}
    kp <- ncol(B2)-1
    A <- matrix(0,nrow=k,ncol=kp)

 nu <- B2[,1]
 B2 <- B2[,2:ncol(B2)]

for(i in 1:k)
 {  
B2a <- B2[i,]

R1 <- which(abs(B2a)>1e-8)

if(length(R1)<2){A[i,]=A[i,]}
else{

    R <- matrix(0,nrow=kp,ncol=length(R1))
Z1 <- as.matrix(Z)
Y1 <- as.matrix(Y)
# Create restriction matrix
jj=1
for( ii in 1:nrow(R)){
if(ii %in% R1)
  {
  R[ii,jj] <- 1
  jj=jj+1
}

}

RLS2a <- QRFact(Y1,Z1,R,i,kp,k,p)


    
A[i,] <- matrix(c(R%*%RLS2a),ncol=kp,nrow=1)

}
}
betaR <- cbind(nu,A)
return(betaR)
}


QRFact <- function(Y,Z,R,i,kp,k,p)
{

Ra <- matrix(0,nrow=kp+k,ncol=ncol(R)+1)
Ra[1:kp,1:ncol(R)] <- R
for(j in (kp+1):(kp+k))
    {
        if(j-(kp)==i)
            {
                Ra[j,ncol(Ra)] <- 1
            }

        }
Ra2<- Ra
 K <- cbind(t(Z),Y)
K2<-K
   K2 <- K%*%Ra
   q = ncol(K2)
   delta = (q^2 + q + 1) * (sqrt(.Machine$double.eps))
   scale = sqrt(delta) * sqrt(apply(K2^2, 2, sum))
   R1 = qr.R(qr((rbind(K2, diag(scale)))), complete = TRUE)
   R11 <- R1[1:ncol(R),1:ncol(R)]
   R12 <- as.matrix(R1[1:ncol(R),ncol(Ra)],ncol=1)
   return(solve(R11)%*%R12)

}


.relaxedVAR <- function(Y1,B2,Z1,k,p)
  {

if(all(as.numeric(B2[,2:ncol(B2)])==0)){return(B2)}
else{
nu <- B2[,1]
B2 <- B2[,2:ncol(B2)]
  
R1 <- which(abs(B2)>1e-8)
R <- matrix(0,nrow=length(B2),ncol=length(R1))
Z1 <- as.matrix(Z1)
Y1 <- as.matrix(Y1)
# Create restriction matrix
jj=1
for( ii in 1:nrow(R)){
if(ii %in% R1)
  {
  R[ii,jj] <- 1
  jj=jj+1
}

}


q <- ncol(crossprod(t(Z1)))
delta = (q^2 + q + 1) * (sqrt(.Machine$double.eps))
scale = sqrt(delta) * sqrt(apply(crossprod(t(Z1))^2, 2, sum))
Z1a <- crossprod(t(Z1))+diag(scale)


A <- .relaxedLeastSquares(Y1,cbind(nu,B2),Z1,k,p)[,2:(ncol(B2)+1)]

sigmaU <- 1/(nrow(Y1))*(t(Y1)-A%*%Z1)%*%t(t(Y1)-A%*%Z1)}

EGLS <- solve(t(R)%*%kronecker(Z1a,solve(sigmaU))%*%R)%*%t(R)%*%(kronecker(Z1,solve(sigmaU)))%*%as.vector(t(Y1))
RVAR <- matrix(R%*%EGLS,ncol=nrow(Z1),nrow=ncol(Y1))
betaR <- cbind(nu,RVAR)
return(betaR)
}


#' Sparsity Plot of a Coefficient Matrix 
#'
#' @param B \eqn{k \times kp} coefficient matrix
#' @param k Number of series
#' @param p Maximal Lag order
#' @param title (optional) Plot title
#' @return NA, side effect is graph
#' @details Similar to SparsityPlot.BigVAR, but the input object does not need to be of class \link{BigVAR.results}
#' @name SparsityPlot
#' @import lattice
#' @rdname SparsityPlot
#' @seealso \code{\link{SparsityPlot.BigVAR.results}},\code{\link{BigVAR.results}} 
#' @examples
#' data(Generator)
#' k=3;p=4
#' SparsityPlot(A[1:k,],p,k,title="Sparsity Plot of Example Generator Matrix")
#' @export
SparsityPlot <- function(B,p,k,title=NULL)
    {
text <- c()
for(i in 1:p)
    {

        text1 <- as.expression(bquote(bold(B)[.(i)]))
        text <- append(text,text1)
          
    }
f <- function(m) t(m)[,nrow(m):1]

rgb.palette <- colorRampPalette(c("white", "blue" ),space = "Lab")


at <- seq(k/2+.5,p*(k)+.5,by=k)

se2=seq(1.75,by=k,length=k)

L2 <- levelplot(f(abs(B)),col.regions=rgb.palette,colorkey=NULL,xlab=NULL,ylab=NULL,main=list(label=title,cex=1),panel=function(...){ panel.levelplot(...);panel.abline(a=NULL,b=1,h=seq(1.5,p*k+.5,by=1),v=seq(1.5,by=1,length=p*k));panel.abline(a=NULL,b=1,v=seq(k+.5,p*k+.5,by=k),lwd=3)},scales=list(x=list(alternating=1,labels=text,cex=1,at=at,tck=c(0,0)),y=list(alternating=0,tck=c(0,0))))

return(L2)

        }


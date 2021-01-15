### Support functions for BigVAR Package
### These are mostly utility functions that will not be seen by the user




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
#' \item{"Z"}{\eqn{kp+ms+1\times T-max(p,s)} VARX lag matrix}
#' \item{"Y"}{adjusted \eqn{k\times T-max(p,s)} endogenous series}
#' }
#' @details This function is not required unless you which to design your own cross validation routine. 
#' @references
#' See page 15 of Lutkepohl, "A New Introduction to Multiple Time Series Analysis
#' @seealso \code{\link{MultVarSim}}
#' @examples
#' data(Y)
#' # construct VAR lag matrix with p=4
#' ZZ<-VARXLagCons(Y,X=NULL,p=4,s=0)
#' @export
VARXLagCons <- function(Y,X=NULL,p,s=0,oos=FALSE,contemp=FALSE){
    if(is.null(X)){X <- matrix(0,nrow=nrow(Y))}
    if(nrow(Y)!=nrow(X)){
        stop("Y and X must have same dimensions")
    }
    if(s==0 & !contemp){
        m=0
    }else{
        m=ncol(X)
    }

    if(p<0|m<0){stop("lag orders must be positive")}
    k <- ifelse(is.null(Y),ncol(Y),0)    
    XX <- VARXCons(Y,X,k,p,m,s,oos,contemp)
    Y <- t(Y[(max(c(p,s))+1):nrow(Y),])
    return(list(Z=XX,Y=Y))
}

.huber_loss <- function(r,delta){
    l <- ifelse(abs(r)<delta,1/2*abs(r)^2,delta*(abs(r)-1/2*delta))
    return(l)
}

.calc.loss <-function(x,univ=FALSE,loss,delta){
    if(loss=="L1"){
        l <- sum(abs(x))
    }else if(loss=="L2"){
        if(univ){
            l=x^2
        }else{
            l=norm2(x)^2
        }

    }else if(loss=="Huber"){
        l=.huber_loss(x,delta)
        
    }
    return(sum(l))
}
                                        # mean benchmark
.evalMean <- function(Y,T1,T2,h=1,loss="L2",delta=2.5)
{

    ypredF <- NULL
    if(!"matrix"%in%class(Y)){
        Y <- matrix(Y,ncol=1)
    }

    MSFE <- c()

    k <- ncol(Y)

    for (u in (T1-h+2):T2) {

        if(h+u-1>T2){break}
        
        trainY1 <- Y[1:(u-1), ]

        if(k>1){

            ypred <- colMeans(trainY1)
            ypredF <- rbind(ypredF,ypred)
        }else{
            ypred <- mean(trainY1)
            ypredF <- c(ypredF,ypred)
        }

        uhat <- matrix(Y[u+h-1, ] - ypred, 
                       ncol = k)

        MSFE <- c(MSFE,.calc.loss(uhat,univ=FALSE,loss,delta))
        ## MSFE <- c(MSFE,norm2(uhat)^2)

    }
    return(list(Mean=mean(na.omit(MSFE)),SD=sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE))),preds=as.matrix(ypredF)))
}

                                        # random walk benchmark
.evalRW <- function(Y,T1,T2,h=1,loss="L2",delta=2.5)
{

    if(!"matrix"%in%class(Y)){
        Y <- matrix(Y,ncol=1)
    }
    ypredF <- NULL
    MSFE <- c()
    
    k <- ncol(Y)
    
    for (u in (T1-h+2):T2) {
        
        if(h+u-1>T2){break}
        
        trainY1 <- Y[u-1, ]
        ypredF <- rbind(ypredF,trainY1)
        uhat <- matrix(Y[u+h-1, ] - trainY1, 
                       ncol = k)
        
        ## MSFE <- c(MSFE,norm2(uhat)^2)
        MSFE <- c(MSFE,.calc.loss(uhat,univ=FALSE,loss,delta))
        
    }
    return(list(Mean=mean(na.omit(MSFE)),SD=sd(na.omit(MSFE))/sqrt(length(na.omit(MSFE))),preds=as.matrix(ypredF)))
}

                                        # Constructs the Grid of Lambda Values: VAR
.LambdaGridE<- function (gran1, gran2, jj = jj, Y, Z, group,p,k,MN,alpha,C,intercept,tol,separate_lambdas=FALSE,verbose=FALSE,linear=FALSE) 
{
    ## browser()
    if (group == "Lag") {
        mat <- list()
        for (i in 1:length(jj)) {
            if(k>1){
                mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y)
            }
            else{
                mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])
            }
        }
        gamstart <- max(unlist(mat))
    }
    if (group == "Basic"|group=="Tapered"|group=="BasicEN") {
        if(!separate_lambdas ){
            ## browser()
            if(group == "Basic"|group=="Tapered"){
                gamstart <- max(abs(t(Y) %*% t(Z)))
                }else{
                    gamstart <- max(abs(t(Y) %*% t(Z)))/alpha
                    }
            ## print(gamstart)
        }else{
            gamstart <- c()
            for(i in 1:k){

                if(group == "Basic"|group=="Tapered"){
                    gamstart[i] <- max(abs(t(Y[,i,drop=F]) %*% t(Z)))
                }else{
                    gamstart[i] <- max(abs(t(Y[,i,drop=F]) %*% t(Z)))/alpha
                }  
            }
        }
    }

    ## if (group == "MCP"|group=="SCAD") {
    ##     ## browser()
    ##     gamstart <- max(abs(crossprod(t(Z),Y)))/nrow(Y)
    ## }


    if (group == "SparseLag") {

        mat <- list()

        if(alpha>0){
            for (i in 1:length(jj)) {

                if(k>1){

                    mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y)*(1/(alpha))

                }

                else{mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])

                }

            }

            gamstart <- max(unlist(mat))

        }else{
            gamstart <- max(t(Y) %*% t(Z))

        }
    }

    if (group == "OwnOther") {

        mat <- list()

        ZZ <- kronecker(t(Z), diag(k))

        for (i in 1:length(jj)) {

            mat[[i]] <- norm(as.vector(t(Y)) %*% ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), 
                             "F")

        }

        gamstart <- max(unlist(mat))

    }
    if (group == "SparseOO") {

        mat <- list()

        ZZ <- kronecker(t(Z), diag(k))

        if(alpha>0){
            for (i in 1:length(jj)) {

                mat[[i]] <- norm(1/(k + 1) * as.vector(t(Y)) %*% 
                                 ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), "F")

            }

            gamstart <- max(unlist(mat))

        }else{
            gamstart <- max(t(Y) %*% t(Z))

        }
    }

    if (group == "HVARC"|group=="HVAROO"|group=="HVARELEM") {

        gmax <- c()

        for (i in 1:k) {

            gmax[i]  <- norm2(Z %*% Y[, i])

        }

        if(!separate_lambdas){
            gamstart <- max(gmax)
        }else{
            gamstart <- gmax
        }
    }

    if(group=="Tapered"){

        beta <- array(0,dim=c(k,k*p+1,1))

    }else{
        beta <- array(0,dim=c(k,k*p+1,1))
    }

    if(!separate_lambdas){
        gamstart <- LGSearch(gamstart,Y,Z,beta,group,k,p,jj,MN,alpha,C,intercept,tol)
        if(!linear){
            gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), 
                            length = gran2))
        }else{
            gamm <- exp(seq(from = gamstart, to = gamstart/gran1, 
                            length = gran2))
            
        }
    }else{

        gamm <- matrix(NA,nrow=c(gran2),ncol=k)

        
        for(i in 1:length(gamstart)){
            gamstart[i] <- LGSearch(gamstart[i],Y,Z,beta,group,k,p,jj,MN,alpha,C,intercept,tol)
            if(verbose & i %%20==0){
                print( sprintf('determined lambda grid for series %s',i))
            }
            
            gamm[,i] <- exp(seq(from = log(gamstart[i]), to = log(gamstart[i]/gran1), 
                                length = gran2))

        }
    }

    return(gamm)

}
                                        # Constructs penalty grid for each value of alpha in case of dual cv
.LambdaGridEDual <- function(gran1,gran2,jj,Y,Z,group,p,k,MN,alpha,C,intercept,tol){

    Lambda <- matrix(0,nrow=gran2,ncol=length(alpha))
    
    for(i in 1:length(alpha)){

        Lambda[,i] <- .LambdaGridE(gran1,gran2,jj,Y,Z,group,p,k,MN,alpha[i],C,intercept,tol)

    }
    return(Lambda)

    


}

                                        # Construct Lambda Grid: VARX
.LambdaGridX<- function (gran1, gran2, jj = jj, Y, Z, group,p,k1,s,m,k,MN,alpha,C,intercept,tol,separate_lambdas=FALSE,verbose=FALSE) 
{

    if (group == "Lag") {

        mat  <- list()

        for (i in 1:length(jj)) {

            if(k>1){

                mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y)

            }

            else{
                mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])
            }

        }

        gamstart <- max(unlist(mat))

    }

    if (group == "Basic"|group=="BasicEN") {
        if(!separate_lambdas ){
            ## browser()
            if(group == "Basic"){
                gamstart <- max(abs(t(Y) %*% t(Z)))
                }else{
                    gamstart <- max(abs(t(Y) %*% t(Z)))/alpha
                    }
            ## print(gamstart)
        }else{
            gamstart <- c()
            for(i in 1:k){

                if(group == "Basic"){
                    gamstart[i] <- max(abs(t(Y[,i,drop=F]) %*% t(Z)))
                }else{
                    gamstart[i] <- max(abs(t(Y[,i,drop=F]) %*% t(Z)))/alpha
                }  
            }
        }
    }


    ## if (group == "MCP"|group=="SCAD") {
    ##     ## browser()
    ##     gamstart <- max(abs(crossprod(t(Z),Y)))/nrow(Y)

    ## }


    if (group == "SparseLag") {

        mat  <- list()

        if(alpha>0){
            
            for (i in 1:length(jj)) {

                if(k>1){

                    mat[[i]] <- norm2(Z[jj[[i]], ]%*%Y*(1/(alpha)))

                }

                else{

                    mat[[i]] <- norm2(t(Y)%*%Z[jj[[i]],])

                }

            }

            gamstart <- max(unlist(mat))

        }else{

            gamstart <- max(t(Y)%*%t(Z))




        }


    }

    if (group == "OwnOther") {

        mat <- list()

        ZZ <- kronecker(t(Z), diag(k1))

        for (i in 1:length(jj)) {

            mat[[i]] <- norm(as.vector(t(Y)) %*% ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), 
                             "F")

        }

        gamstart <- max(unlist(mat))

    }

    if (group == "SparseOO") {

        mat <- list()

        ZZ <- kronecker(t(Z), diag(k1))

        if(alpha>0){
            for (i in 1:length(jj)) {

                mat[[i]] <- norm(1/(alpha) * as.vector(t(Y)) %*% 
                                 ZZ[, jj[[i]]]/sqrt(length(jj[[i]])), "F")
            }

            gamstart <- max(unlist(mat))

        }else{

            gamstart <- max(t(Y)%*%t(Z))

        }

    }

    if (group=="EFX") {

        gmax <- c()

        for (i in 1:k1) {

            gmax[i] <- norm2(Z %*% Y[, i])/sqrt(k*p)

        }
        
        gamstart <- max(gmax)

    }
    

    beta <- array(0,dim=c(k1,k1*p+s*m+1,1))
    ## gamstart <- LGSearchX(gamstart,Y,Z,beta,group,k1,p,s,m,jj,k,MN,alpha,C,intercept,tol)
    ## gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), 
    ##                 length = gran2))

    if(!separate_lambdas){
        gamstart <- LGSearchX(gamstart,Y,Z,beta,group,k1,p,s,m,jj,k,MN,alpha,C,intercept,tol)

        gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), 
                        length = gran2))

    }else{

        gamm <- matrix(NA,nrow=c(gran2),ncol=k1)

        ## browser()
        for(i in 1:length(gamstart)){
            gamstart[i] <- LGSearchX(gamstart[i],Y,Z,beta,group,k1,p,s,m,jj,k,MN,alpha,C,intercept,tol)
            gamstart[i] <- ifelse(gamstart[i]==0,1e-4,gamstart[i])
            if(verbose & i %%20==0){
                print( sprintf('determined lambda grid for series %s',i))
            }
            gamm[,i] <- exp(seq(from = log(gamstart[i]), to = log(gamstart[i]/gran1), 
                                length = gran2))

        }
        
    }

    
    return(gamm)

}

.LambdaGridXDual <- function(gran1,gran2,jj,Y,Z,group,p,k1,s,m,k,MN,alpha,C,intercept,tol){

    Lambda <- matrix(0,nrow=gran2,ncol=length(alpha))
    for(i in 1:length(alpha)){

        Lambda[,i] <- .LambdaGridX(gran1,gran2,jj,Y,Z,group,p,k1,s,m,k,MN,alpha[i],C,intercept,tol)

    }
    return(Lambda)

    

}

                                        # Forecast evaluation: VARX (called in cv.bigvar)
.BigVAREVALX <- function(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m,contemp,alpha,C,intercept,tol,window.size,separate_lambdas,loss="L2",delta=2.5)
{
    VARX <- TRUE
    if(contemp){
        s1 <- 1

    }else{

        s1 <- 0

    }

    preds <- matrix(NA,nrow=length((T2+1):T),ncol=k1)

    gran2 <- 1	

    gamm <- gamopt
    if(separate_lambdas){
        gamm <- matrix(gamm,nrow=1,ncol=k1 )
        
    }

    Y <- ZFull$Y

    ## MSFE <- rep(NA,length((T2+1):T))

    beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,1))  

    kk <- NULL
    jj <- NULL
    jjcomp <- NULL
    activeset <- NULL
    q1a <- NULL

    if (group == "Lag") {

        jj <- groupfunVARX(p,k,k1,s+s1)

        jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

    }else if (group == "SparseLag") {

        jj <- groupfunVARX(p, k,k1,s+s1)

        q1a <- list()

        for (i in 1:(p+s+s1)) {

            q1a[[i]] <- matrix(runif(k1, -1, 1), ncol = 1)

        }

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)


    }else if (group == "OwnOther") {

        kk <- diaggroupfunVARX(p, k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
    }else if (group == "SparseOO") {

        kk <- diaggroupfunVARX(p, k,k1,s+s1)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
    }else{

        jj <- NULL
        kk <- NULL
        activeset <- NULL
    }


    if(verbose){

        print("Evaluation Stage")

        pb <- txtProgressBar(min = T2-h+2, max = T, style = 3)
    }
    ## if(separate_lambdas){
    ##     MSFE <- matrix(NA,nrow=length((T1+1):T2),ncol=k1)
    ##     betaArray <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,nrow(MSFE)))

    ## }else{
    MSFE <- rep(NA,length((T2+1):T))
    betaArray <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,length(MSFE)))
    ## }

    for (v in (T2-h+2):T) {

        if(v+h-1>T){
            break
        }
        if(window.size!=0){
            ws1 <- max(c(v-window.size-h,1))
            trainY <- ZFull$Y[(ws1+h):(v-1), ]
            trainZ <- ZFull$Z[, (ws1+h):(v-h)]         
        }else{
            
            
            trainY <- ZFull$Y[h:(v-1), ]

            trainZ <- ZFull$Z[,1:(v-h),drop=F]

        }
        
        ## needed.objs <- c('group','beta','trainZ','trainY','gamm','tol','p','m','k1','k','s','s1','m','MN','C','intercept','separate_lambdas','activeset','alpha','jj','jjcomp','kk','palpha','q1a')
        dual <- FALSE
        if(!group%in%c("SparseLag","SparseOO")){
            q1a <- NULL
        }
        ## objs <- sapply(needed.objs,exists)
        ## objs <- setdiff(needed.objs,ls())
        ## if(length(objs)>0){
        
        ##     for(i in 1:length(objs)){               
        ##         assign(objs[i],NULL)
        ##                }
        ## }
        temp <- .BigVAR.fit(group,beta,trainZ,trainY,gamm,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
        beta <- temp$beta
        activeset <- temp$activeset
        q1a <- temp$q1a


        
        betaEVAL <- matrix(beta[,,1],nrow=k1,ncol=(k1*p+(k-k1)*(s+s1)+1))

        if (RVAR) {

            betaEVAL <- RelaxedLS(cbind(t(trainZ),trainY),betaEVAL,k,p,k1,s+s1)
            
        }
        

        if(MN ){
            
            eZ <- matrix(ZFull$Z[,v],ncol=1)

        }else{

            eZ <- matrix(c(1,ZFull$Z[,v]),ncol=1)

        }

        if(MN){
            preds[v-T2+h-1,] <-  betaEVAL[,2:ncol(betaEVAL)] %*% eZ

            ## if(separate_lambdas){

            ##         for(uu in 1:ncol(ZFull$Y)){
            ##             ## browser()
            ##             MSFE[v-T2+h-1,uu] <- norm2(ZFull$Y[v+h-1,uu] - preds[v-T2+h-1,uu])^2
            ##         }
            ## }else{
            ## MSFE[v-T2+h-1] <- norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
            MSFE[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1,1:k1] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ,univ=FALSE,loss,delta)

            
            ## preds[v-T2+h-1,] <-  betaEVAL[,2:ncol(betaEVAL)] %*% eZ

            diag(beta[,2:(k1+1),1]) <- diag(beta[,2:(k1+1),1])-C # subtract one for warm start purposes 
            ## }
        }else{
            ## browser()
            preds[v-T2+h-1,] <- betaEVAL %*% eZ
            ## if(separate_lambdas){

            ##         for(uu in 1:ncol(ZFull$Y)){
            ##             ## browser()
            ##             MSFE[v-T2+h-1,uu] <- norm2(ZFull$Y[v+h-1,uu] - preds[v-T2+h-1,uu])^2
            ##         }
            ## }else{

            ## MSFE[v-T2+h-1] <- norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL %*% eZ)^2

            MSFE[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1,1:k1] - betaEVAL %*% eZ,univ=FALSE,loss,delta)

            ## }
            
        }

        if(verbose){

            setTxtProgressBar(pb, v)

        }


        betaArray[,,v-T2+h-1] <- beta
    }

    temp <- .BigVAR.fit(group,beta,ZFull$Z,ZFull$Y,gamm,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
    betaPred <- temp$beta        
    ## browser()
    betaPred <- as.matrix(betaPred[,,1])

    
    return(list(MSFE=MSFE,betaPred=betaPred,predictions=preds,betaArray=betaArray))
    

}



                                        # BigVAR with incremental re-estimation 
.BigVAREVAL_rolling <- function(ZFull,MSFE,gamm,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,alpha,recursive,C,intercept,tol,window.size,separate_lambdas,loss,delta,ONESE)
{
    ## VARX <- TRUE
    ## if(contemp){
    ##     s1 <- 1

    ## }else{

    ##     s1 <- 0

    ## }
    s=0;k1=0;s1=0;m=0
    preds <- matrix(NA,nrow=length((T2+1):T),ncol=k)
    gran2 <- nrow(gamm)
    ## if(separate_lambdas){
    ##     gamm <- matrix(gamm,nrow=1,ncol=k )
    
    ## }


    if(verbose)
    {

        print("Evaluation Stage")

        pb <- txtProgressBar(min = T2-h+2, max = T, style = 3)

    }

    Y <- ZFull$Y


    betaArray <-  array(0,dim=c(k,p*k+1,dim(MSFE)[1]))

    beta <-  array(0,dim=c(k,p*k+1,nrow(gamm)))

    
    kk <- NULL
    jj <- NULL
    jjcomp <- NULL
    activeset <- NULL
    q1a <- NULL
    if (group == "Lag"){

        jj <- .groupfuncpp(p,k)

        jjcomp <- .groupfuncomp(p,k)

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

    }else if (group == "SparseLag"){

        jj <- .groupfun(p, k)

        q1a <- list()

        for (i in 1:(p)) {

            q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)

        }

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

        

    }else if (group == "OwnOther"){


        kk <- .lfunction3cpp(p, k)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)

    }else if (group == "SparseOO") {

        kk <- .lfunction3cpp(p, k)

        jjcomp <- .lfunctioncomp(p,k)

        jj=.lfunction3(p,k)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)

        q1a <- list()

        for (i in 1:(2*p))
        {

            q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)

        }

    }else{

        kk=NULL
        jj=NULL
        activeset=NULL
        jjcomp=NULL
        
    }
    if(!separate_lambdas){
        MSFE_oos <- c()
    }else{
        MSFE_oos <- matrix()
    }
    gamm_evolve <- c()
    MSFE_old <- MSFE
    ## browser()
    for(v in (T2-h+2):T){

        if(h>1 & !recursive){

            if(window.size!=0){
                ws1 <- max(c(v-window.size-h,1))
                trainY <- ZFull$Y[(ws1+h):(v-1), ]
                trainZ <- ZFull$Z[, (ws1+h):(v-h)]         
            }else{

                trainY <- ZFull$Y[(h):(v-1), ]
                
                trainZ <- ZFull$Z[, 1:(v-h)]
            }
            
        }else{
            if(window.size!=0){
                ws1 <- max(c(v-window.size,1))
                trainY <- ZFull$Y[(ws1):(v-1), ]
                trainZ <- ZFull$Z[, (ws1):(v-1)]         
            }else{
                trainY <- ZFull$Y[(1):(v-1), ]                       
                trainZ <- ZFull$Z[, (1):(v-1)]
            }
        }
        

        if(v+h-1>T){
            break
        }

        dual <- FALSE
        if(!group%in%c("SparseLag","SparseOO")){
            q1a <- NULL
        }
        ## browser()
        VARX <- FALSE
        temp <- .BigVAR.fit(group,beta,trainZ,trainY,gamm,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
        beta <- temp$beta
        activeset <- temp$activeset
        q1a <- temp$q1a
        

        if(MN ){
            
            eZ <- matrix(ZFull$Z[,v],ncol=1)

        }else{

            eZ <- matrix(c(1,ZFull$Z[,v]),ncol=1)

        }
        ## browser()
        ## abind::adrop(beta[,,1,drop=F],3)
        if(!ONESE & !separate_lambdas){
            optind <- max(which(colMeans(na.omit(MSFE))==min(colMeans(na.omit(MSFE)))))

            gamm_opt <- apply(MSFE,1,which.min)

            gamopt <- gamm[optind]
        }else if(ONESE & !separate_lambdas){

            MSFE2 <- MSFE 
            G2 <- colMeans(na.omit(MSFE2))
            G3 <- sd(na.omit(MSFE2))/sqrt(nrow(na.omit(MSFE2)))
            optind <- min(which(G2<(min(G2)+G3)))
            gamopt <- gamm[optind]

            
        }else if (separate_lambdas){
            ## browser()
            if(ONESE){
                MSFES <- t(apply(MSFE,3,colMeans))
                sds <- t(apply(MSFE,3,function(x)sd(na.omit(x))/sqrt(nrow(na.omit(x)))))
                ## for()
                ## optinds <- apply(MSFES,2,function(x)min(which(x<x+sds)))
                gamopt <- c()
                optinds <- c()
                for(i in 1:nrow(MSFES)){
                    optinds[i] <- min(which(MSFES[i,]<sds[i]+min(MSFES[i,])))
                    gamopt[i] <- gamm[optinds[i],i,drop=F]
                }
                optind=optinds
            }else{
                ## MSFES <- t(apply(MSFE,3,colMeans))
                ## browser()
                ## G3 <- sd(na.omit(MSFE2))/sqrt(nrow(na.omit(MSFE2)))
                ## optind <- min(which(G2<(min(G2)+G3)))
                ## gamopt <- gamm[optind]


                MSFES <- t(apply(MSFE,3,colMeans))
                optinds <- apply(MSFES,1,which.min)
                ## optinds <- sapply(MSFES,which.min)
                ## browser()
                gamopt <- c()
                for(i in 1:nrow(MSFES)){
                    gamopt[i] <- gamm[optinds[i],i]
                }
                optind=optinds


            }


        }
        
        
        gamm_evolve <- rbind(gamm_evolve,gamopt)
        if(!separate_lambdas){
            MSFE_temp <- matrix(0,nrow=1,ncol=length(gamm))
        }else{
            MSFE_temp <- array(0,dim=c(1,nrow(gamm),ncol(ZFull$Y)))
            temp_loss <- matrix(0,nrow=1,ncol=k)
        }
        ## browser()
        for(ii in 1:nrow(gamm)){
            betaEVAL <-  abind::adrop(beta[,,ii,drop=F],3)

            
            if (RVAR) {
                
                betaEVAL <- RelaxedLS(cbind(t(trainZ),trainY),betaEVAL,k,p,0,0)      
                
            }

            
            if(MN){
                if(h>1 & recursive){

                    ptemp <- betaEVAL[,2:ncol(betaEVAL)] %*% eZ

                    pred <- matrix(ptemp,nrow=1)

                    
                    if(ii%in%optind){
                        preds[v-T2+h-1,] <- predictMS(pred,trainY,h-1,betaEVAL[,2:ncol(betaEVAL)],p,TRUE)
                        if(separate_lambdas){
                            inds <- which(optind==ii)
                            for(j in 1:ncol(ZFull$Y)){
                                if(j%in%inds){
                                    ## MSFE_oos[v-T2+h-1,j] <- .calc.loss(ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j],univ=FALSE,loss,delta)
                                    temp_loss[,j] <- ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j]

                                    betaArray[j,,v-T2+h-1] <- betaEVAL[j,]    
                                    preds[v-T2+h-1,j] <- preds[v-T2+h-1,j]   ## betaEVAL[j,2:ncol(betaEVAL),drop=F] %*% eZ                

                                }
                            }


                        }else{
                            MSFE_oos[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1] - preds[v-T2+h-1,],univ=FALSE,loss,delta)
                        }
                    }
                    ptemp2 <- predictMS(pred,trainY,h-1,betaEVAL[,2:ncol(betaEVAL)],p,TRUE)
                    if(separate_lambdas){
                        for(j in 1:ncol(ZFull$Y)){

                            MSFE_temp[1,ii,j] <- .calc.loss(ZFull$Y[v+h-1,j] - ptemp[j,1],univ=FALSE,loss,delta)
                        }
                    }else{
                        MSFE_temp[1,ii] <- .calc.loss(ZFull$Y[v+h-1,] - ptemp,univ=FALSE,loss,delta)## norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                    }
                }else{     
                    
                    if(ii%in%optind){

                        preds[v-T2+h-1,] <-  betaEVAL[,2:ncol(betaEVAL)] %*% eZ
                        ## MSFE_oos[v-T2+h-1] <- norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                        if(separate_lambdas){
                            inds <- which(optind==ii)
                            for(j in 1:ncol(ZFull$Y)){
                                if(j%in%inds){
                                    temp_loss[,j] <- ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j]

                                    betaArray[j,,v-T2+h-1] <- betaEVAL[j,]    
                                    preds[v-T2+h-1,j] <- betaEVAL[j,2:ncol(betaEVAL),drop=F] %*% eZ                

                                }
                            }
                            ## MSFE_oos[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j],univ=FALSE,loss,delta)
                        }else{
                            MSFE_oos[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1,] - preds[v-T2+h-1,],univ=FALSE,loss,delta)
                        }
                    }
                    if(separate_lambdas){
                        for(j in 1:ncol(ZFull$Y)){

                            MSFE_temp[1,ii,j] <- .calc.loss(ZFull$Y[v+h-1,j] - betaEVAL[j,2:ncol(betaEVAL)] %*% eZ,univ=FALSE,loss,delta)
                        }
                    }else{
                        MSFE_temp[1,ii] <- .calc.loss(ZFull$Y[v+h-1,] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ,univ=FALSE,loss,delta)## norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                    }
                    diag(beta[,2:(dim(beta)[2]),ii]) <- diag(beta[,2:(dim(beta)[2]),ii])-C # subtract one for warm start purposes
                }
            }else{

                if(h>1&recursive){
                    ptemp <- betaEVAL %*% eZ

                    pred <- matrix(ptemp,nrow=1)

                    
                    if(ii%in%optind){
                        preds[v-T2+h-1,] <- predictMS(pred,trainY,h-1,betaEVAL,p,FALSE)

                        if(separate_lambdas){

                            inds <- which(optind==ii)
                            for(j in 1:ncol(ZFull$Y)){
                                if(j%in%inds){
                                    temp_loss <- ZFull$Y[v+h-1,j] - betaEVAL[j,] %*% eZ
                                    betaArray[j,,v-T2+h-1] <- betaEVAL[j,]    
                                    preds[v-T2+h-1,j] <- betaEVAL[j,] %*% eZ                

                                }
                            }
                        }else{
                            MSFE_oos[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1] - preds[v-T2+h-1,],univ=FALSE,loss,delta)
                        }
                    }
                    ptemp2 <- predictMS(pred,trainY,h-1,betaEVAL,p,FALSE)
                    if(separate_lambdas){


                        for(j in 1:ncol(ZFull$Y)){

                            MSFE_temp[1,ii,j] <- .calc.loss(ZFull$Y[v+h-1,j] - ptemp2[j],univ=FALSE,loss,delta)
                        }

                    }else{
                        
                        
                        MSFE_temp[1,ii] <- .calc.loss(ZFull$Y[v+h-1,] - ptemp2,univ=FALSE,loss,delta)## norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                    }
                }else{
                    
                    if(ii%in%optind){
                        ## browser()

                        if(separate_lambdas){

                            inds <- which(optind==ii)
                            for(j in 1:ncol(ZFull$Y)){
                                if(j%in%inds){
                                    temp_loss[,j] <- ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j] %*% eZ

                                    betaArray[j,,v-T2+h-1] <- betaEVAL[j,]    
                                    preds[v-T2+h-1,j] <- betaEVAL %*% eZ                
                                    ## MSFE_oos[v-T2+h-1] <-.calc.loss(ZFull$Y[v+h-1,] - preds[v-T2+h-1,],univ=FALSE,loss,delta) ## norm2(ZFull$Y[v+h-1,1:k1] - betaEVA                      
                                }
                            }
                        }else{

                            
                            if(!separate_lambdas){

                                betaArray[,,v-T2+h-1] <- betaEVAL    

                                preds[v-T2+h-1,] <- betaEVAL %*% eZ
                                
                                MSFE_oos[v-T2+h-1] <-.calc.loss(ZFull$Y[v+h-1,] - preds[v-T2+h-1,],univ=FALSE,loss,delta) ## norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL %*% eZ)^2
                            }

                            if(separate_lambdas){
                                for(j in 1:ncol(ZFull$Y)){
                                    MSFE_temp[1,ii,j] <- ZFull$Y[v+h-1,j] - betaEVAL[j,] %*% eZ
                                    
                                }
                            }else{          
                                MSFE_temp[1,ii] <- .calc.loss(ZFull$Y[v+h-1,] - betaEVAL %*% eZ,univ=FALSE,loss,delta)
                            }
                        }
                        
                        ## if(v>T2-h+2){browser()}
                        ## MSFE_temp[1,ii] <- norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL %*% eZ)^2
                    }

                }
            }
        }
        if(!separate_lambdas){    
            MSFE <- rbind(MSFE,MSFE_temp)
            MSFE <- MSFE[2:nrow(MSFE),]
        }else{
            MSFE <- abind::abind(MSFE,MSFE_temp,along=1)
            MSFE <- MSFE[2:dim(MSFE)[1],,]
        }
        if(verbose){

            setTxtProgressBar(pb, v)

        }
        
        if(separate_lambdas){
            MSFE_oos[v-T2+h-1] <- .calc.loss(temp_loss,univ=FALSE,loss,delta)
        }
    }
    ## browser()
    betaPred <- betaEVAL
    
    return(list(MSFE=MSFE_oos,betaPred=betaPred,predictions=preds,betaArray=betaArray,gamm_evolve=gamm_evolve))
    

}



.BigVAREVALX_rolling <- function(ZFull,MSFE,gamm,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m,contempalpha,recursive,C,intercept,tol,window.size,separate_lambdas,loss,delta,ONESE)
{
    VARX <- TRUE
    if(contemp){
        s1 <- 1

    }else{

        s1 <- 0

    }
    ## s=0;k1=0;s1=0;m=0
    preds <- matrix(NA,nrow=length((T2+1):T),ncol=k1)
    gran2 <- nrow(gamm)
    ## if(separate_lambdas){
    ##     gamm <- matrix(gamm,nrow=1,ncol=k )
    
    ## }


    if(verbose)
    {

        print("Evaluation Stage")

        pb <- txtProgressBar(min = T2-h+2, max = T-h, style = 3)

    }

    Y <- ZFull$Y


    beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2))  


    ## beta <-  array(0,dim=c(k,p*k+1,nrow(gamm)))

    
    kk <- NULL
    jj <- NULL
    jjcomp <- NULL
    activeset <- NULL
    q1a <- NULL
    if (group == "Lag") {

        jj <- groupfunVARX(p,k,k1,s+s1)

        jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

    }else if (group == "SparseLag") {

        jj <- groupfunVARX(p, k,k1,s+s1)

        q1a <- list()

        for (i in 1:(p+s+s1)) {

            q1a[[i]] <- matrix(runif(k1, -1, 1), ncol = 1)

        }

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)


    }else if (group == "OwnOther") {

        kk <- diaggroupfunVARX(p, k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
    }else if (group == "SparseOO") {

        kk <- diaggroupfunVARX(p, k,k1,s+s1)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
    }else{

        jj <- NULL
        kk <- NULL
        activeset <- NULL
    }
    ## if(!separate_lambdas){
        MSFE_oos <- c()
    ## }else{
    ##     MSFE_oos <- matrix()
    ## }
    gamm_evolve <- c()
    MSFE_old <- MSFE
    ## browser()
    for(v in (T2-h+2):T){

        if(h>1 & !recursive){

            if(window.size!=0){
                ws1 <- max(c(v-window.size-h,1))
                trainY <- ZFull$Y[(ws1+h):(v-1), ]
                trainZ <- ZFull$Z[, (ws1+h):(v-h)]         
            }else{

                trainY <- ZFull$Y[(h):(v-1), ]
                
                trainZ <- ZFull$Z[, 1:(v-h)]
            }
            
        }else{
            if(window.size!=0){
                ws1 <- max(c(v-window.size,1))
                trainY <- ZFull$Y[(ws1):(v-1), ]
                trainZ <- ZFull$Z[, (ws1):(v-1)]         
            }else{
                trainY <- ZFull$Y[(1):(v-1), ]                       
                trainZ <- ZFull$Z[, (1):(v-1)]
            }
        }
        

        if(v+h-1>T){
            break
        }

        dual <- FALSE
        if(!group%in%c("SparseLag","SparseOO")){
            q1a <- NULL
        }
        ## browser()
        temp <- .BigVAR.fit(group,beta,trainZ,trainY,gamm,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
        beta <- temp$beta
        activeset <- temp$activeset
        q1a <- temp$q1a
        

        if(MN ){
            
            eZ <- matrix(ZFull$Z[,v],ncol=1)

        }else{

            eZ <- matrix(c(1,ZFull$Z[,v]),ncol=1)

        }
        ## browser()
        ## abind::adrop(beta[,,1,drop=F],3)
        if(!ONESE & !separate_lambdas){
            optind <- max(which(colMeans(na.omit(MSFE))==min(colMeans(na.omit(MSFE)))))

            gamm_opt <- apply(MSFE,1,which.min)

            gamopt <- gamm[optind]
        }else if(ONESE & !separate_lambdas){

            MSFE2 <- MSFE 
            G2 <- colMeans(na.omit(MSFE2))
            G3 <- sd(na.omit(MSFE2))/sqrt(nrow(na.omit(MSFE2)))
            optind <- min(which(G2<(min(G2)+G3)))
            gamopt <- gamm[optind]

            
        }else if (separate_lambdas){
            ## browser()
            if(ONESE){
                MSFES <- t(apply(MSFE,3,colMeans))
                sds <- t(apply(MSFE,3,function(x)sd(na.omit(x))/sqrt(nrow(na.omit(x)))))
                ## for()
                ## optinds <- apply(MSFES,2,function(x)min(which(x<x+sds)))
                gamopt <- c()
                optinds <- c()
                for(i in 1:nrow(MSFES)){
                    optinds[i] <- min(which(MSFES[i,]<sds[i]+min(MSFES[i,])))
                    gamopt[i] <- gamm[optinds[i],i,drop=F]
                }
                optind=optinds
            }else{
                ## MSFES <- t(apply(MSFE,3,colMeans))
                ## browser()
                ## G3 <- sd(na.omit(MSFE2))/sqrt(nrow(na.omit(MSFE2)))
                ## optind <- min(which(G2<(min(G2)+G3)))
                ## gamopt <- gamm[optind]


                MSFES <- t(apply(MSFE,3,colMeans))
                optinds <- apply(MSFES,1,which.min)
                ## optinds <- sapply(MSFES,which.min)
                ## browser()
                gamopt <- c()
                for(i in 1:nrow(MSFES)){
                    gamopt[i] <- gamm[optinds[i],i]
                }
                optind=optinds


            }


        }
        
        
        gamm_evolve <- rbind(gamm_evolve,gamopt)
        if(!separate_lambdas){
            MSFE_temp <- matrix(0,nrow=1,ncol=length(gamm))
        }else{
            MSFE_temp <- array(0,dim=c(1,nrow(gamm),ncol(ZFull$Y)))
            temp_loss <- matrix(0,nrow=1,ncol=k1)
        }
        ## browser()
        for(ii in 1:nrow(gamm)){
            betaEVAL <-  abind::adrop(beta[,,ii,drop=F],3)

            
            if (RVAR) {
                
            betaEVAL <- RelaxedLS(cbind(t(trainZ),trainY),betaEVAL,k,p,k1,s)      
                
            }

            
            if(MN){
                if(h>1 & recursive){

                    ptemp <- betaEVAL[,2:ncol(betaEVAL)] %*% eZ

                    pred <- matrix(ptemp,nrow=1)

                    
                    if(ii%in%optind){
                        preds[v-T2+h-1,] <- predictMS(pred,trainY,h-1,betaEVAL[,2:ncol(betaEVAL)],p,TRUE)
                        if(separate_lambdas){
                            inds <- which(optind==ii)
                            for(j in 1:ncol(ZFull$Y)){
                                if(j%in%inds){
                                    ## MSFE_oos[v-T2+h-1,j] <- .calc.loss(ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j],univ=FALSE,loss,delta)
                                    temp_loss[,j] <- ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j]

                                    betaArray[j,,v-T2+h-1] <- betaEVAL[j,]    
                                    preds[v-T2+h-1,j] <- preds[v-T2+h-1,j]   ## betaEVAL[j,2:ncol(betaEVAL),drop=F] %*% eZ                

                                }
                            }


                        }else{
                            MSFE_oos[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1] - preds[v-T2+h-1,],univ=FALSE,loss,delta)
                        }
                    }
                    ptemp2 <- predictMS(pred,trainY,h-1,betaEVAL[,2:ncol(betaEVAL)],p,TRUE)
                    if(separate_lambdas){
                        for(j in 1:ncol(ZFull$Y)){

                            MSFE_temp[1,ii,j] <- .calc.loss(ZFull$Y[v+h-1,j] - ptemp[j,1],univ=FALSE,loss,delta)
                        }
                    }else{
                        MSFE_temp[1,ii] <- .calc.loss(ZFull$Y[v+h-1,] - ptemp,univ=FALSE,loss,delta)## norm2(ZFull$Y[v+h-1,1:k11] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                    }
                }else{     
                    
                    if(ii%in%optind){

                        preds[v-T2+h-1,] <-  betaEVAL[,2:ncol(betaEVAL)] %*% eZ
                        ## MSFE_oos[v-T2+h-1] <- norm2(ZFull$Y[v+h-1,1:k11] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                        if(separate_lambdas){
                            inds <- which(optind==ii)
                            for(j in 1:ncol(ZFull$Y)){
                                if(j%in%inds){
                                    temp_loss[,j] <- ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j]

                                    betaArray[j,,v-T2+h-1] <- betaEVAL[j,]    
                                    preds[v-T2+h-1,j] <- betaEVAL[j,2:ncol(betaEVAL),drop=F] %*% eZ                

                                }
                            }
                            ## MSFE_oos[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j],univ=FALSE,loss,delta)
                        }else{
                            MSFE_oos[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1,] - preds[v-T2+h-1,],univ=FALSE,loss,delta)
                        }
                    }
                    if(separate_lambdas){
                        for(j in 1:ncol(ZFull$Y)){

                            MSFE_temp[1,ii,j] <- .calc.loss(ZFull$Y[v+h-1,j] - betaEVAL[j,2:ncol(betaEVAL)] %*% eZ,univ=FALSE,loss,delta)
                        }
                    }else{
                        MSFE_temp[1,ii] <- .calc.loss(ZFull$Y[v+h-1,] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ,univ=FALSE,loss,delta)## norm2(ZFull$Y[v+h-1,1:k11] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                    }
                    diag(beta[,2:(dim(beta)[2]),ii]) <- diag(beta[,2:(dim(beta)[2]),ii])-C # subtract one for warm start purposes
                }
            }else{

                if(h>1&recursive){
                    ptemp <- betaEVAL %*% eZ

                    pred <- matrix(ptemp,nrow=1)

                    
                    if(ii%in%optind){
                        preds[v-T2+h-1,] <- predictMS(pred,trainY,h-1,betaEVAL,p,FALSE)

                        if(separate_lambdas){

                            inds <- which(optind==ii)
                            for(j in 1:ncol(ZFull$Y)){
                                if(j%in%inds){
                                    temp_loss <- ZFull$Y[v+h-1,j] - betaEVAL[j,] %*% eZ
                                    betaArray[j,,v-T2+h-1] <- betaEVAL[j,]    
                                    preds[v-T2+h-1,j] <- betaEVAL[j,] %*% eZ                

                                }
                            }
                        }else{
                            MSFE_oos[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1] - preds[v-T2+h-1,],univ=FALSE,loss,delta)
                        }
                    }
                    ptemp2 <- predictMS(pred,trainY,h-1,betaEVAL,p,FALSE)
                    if(separate_lambdas){


                        for(j in 1:ncol(ZFull$Y)){

                            MSFE_temp[1,ii,j] <- .calc.loss(ZFull$Y[v+h-1,j] - ptemp2[j],univ=FALSE,loss,delta)
                        }

                    }else{
                        
                        
                        MSFE_temp[1,ii] <- .calc.loss(ZFull$Y[v+h-1,] - ptemp2,univ=FALSE,loss,delta)## norm2(ZFull$Y[v+h-1,1:k11] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                    }
                }else{
                    
                    if(ii%in%optind){
                        ## browser()

                        if(separate_lambdas){

                            inds <- which(optind==ii)
                            for(j in 1:ncol(ZFull$Y)){
                                if(j%in%inds){
                                    temp_loss[,j] <- ZFull$Y[v+h-1,j] - preds[v-T2+h-1,j] %*% eZ

                                    betaArray[j,,v-T2+h-1] <- betaEVAL[j,]    
                                    preds[v-T2+h-1,j] <- betaEVAL %*% eZ                
                                    ## MSFE_oos[v-T2+h-1] <-.calc.loss(ZFull$Y[v+h-1,] - preds[v-T2+h-1,],univ=FALSE,loss,delta) ## norm2(ZFull$Y[v+h-1,1:k11] - betaEVA                      
                                }
                            }
                        }else{

                            
                            if(!separate_lambdas){

                                betaArray[,,v-T2+h-1] <- betaEVAL    

                                preds[v-T2+h-1,] <- betaEVAL %*% eZ
                                
                                MSFE_oos[v-T2+h-1] <-.calc.loss(ZFull$Y[v+h-1,] - preds[v-T2+h-1,],univ=FALSE,loss,delta) ## norm2(ZFull$Y[v+h-1,1:k11] - betaEVAL %*% eZ)^2
                            }

                            if(separate_lambdas){
                                for(j in 1:ncol(ZFull$Y)){
                                    MSFE_temp[1,ii,j] <- ZFull$Y[v+h-1,j] - betaEVAL[j,] %*% eZ
                                    
                                }
                            }else{          
                                MSFE_temp[1,ii] <- .calc.loss(ZFull$Y[v+h-1,] - betaEVAL %*% eZ,univ=FALSE,loss,delta)
                            }
                        }
                        
                        ## if(v>T2-h+2){browser()}
                        ## MSFE_temp[1,ii] <- norm2(ZFull$Y[v+h-1,1:k11] - betaEVAL %*% eZ)^2
                    }

                }
            }
        }
        if(!separate_lambdas){    
            MSFE <- rbind(MSFE,MSFE_temp)
            MSFE <- MSFE[2:nrow(MSFE),]
        }else{
            MSFE <- abind::abind(MSFE,MSFE_temp,along=1)
            MSFE <- MSFE[2:dim(MSFE)[1],,]
        }
        if(verbose){

            setTxtProgressBar(pb, v)

        }
        
        if(separate_lambdas){
            MSFE_oos[v-T2+h-1] <- .calc.loss(temp_loss,univ=FALSE,loss,delta)
        }
    }
    ## browser()
    betaPred <- betaEVAL
    
    return(list(MSFE=MSFE_oos,betaPred=betaPred,predictions=preds,betaArray=betaArray,gamm_evolve=gamm_evolve))
    

}




                                        # BigVAR with incremental re-estimation 
.BigVAREVALX_rolling <- function(ZFull,MSFE,gamm,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m,contemp,alpha,C,intercept,tol,window.size,separate_lambdas,loss,delta)
{
    VARX <- TRUE
    if(contemp){
        s1 <- 1

    }else{

        s1 <- 0

    }

    if(verbose)
    {

        print("Evaluation Stage")

        pb <- txtProgressBar(min = T2-h+2, max = T, style = 3)

    }
    
    preds <- matrix(NA,nrow=length((T2+1):T),ncol=k1)
    gran2 <- length(gamm)
    if(separate_lambdas){
        gamm <- matrix(gamm,nrow=1,ncol=k1 )
        
    }

    Y <- ZFull$Y

    beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2))  

    kk <- NULL
    jj <- NULL
    jjcomp <- NULL
    activeset <- NULL
    q1a <- NULL
    if (group == "Lag") {

        jj <- groupfunVARX(p,k,k1,s+s1)

        jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

    }else if (group == "SparseLag") {

        jj <- groupfunVARX(p, k,k1,s+s1)

        q1a <- list()

        for (i in 1:(p+s+s1)) {

            q1a[[i]] <- matrix(runif(k1, -1, 1), ncol = 1)

        }

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)


    }else if (group == "OwnOther") {

        kk <- diaggroupfunVARX(p, k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
    }else if (group == "SparseOO") {

        kk <- diaggroupfunVARX(p, k,k1,s+s1)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
    }else{

        jj <- NULL
        kk <- NULL
        activeset <- NULL
    }


    if(verbose){

        print("Evaluation Stage")

        pb <- txtProgressBar(min = T2-h+2, max = T, style = 3)
    }
    ## if(separate_lambdas){
    ##     MSFE <- matrix(NA,nrow=length((T1+1):T2),ncol=k1)
    ##     betaArray <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,nrow(MSFE)))

    ## }else{
    ## MSFE <- rep(NA,length((T2+1):T))
    betaArray <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,length((T2-h+2):T)))
    ## betaArray <- list()
    ## }
    MSFE_oos <- c()
    gamm_evolve <- c()
    MSFE_old <- MSFE
    ## browser()
    for(v in (T2-h+2):T){

        if(v+h-1>T){
            break
        }
        if(window.size!=0){
            ws1 <- max(c(v-window.size-h,1))
            trainY <- ZFull$Y[(ws1+h):(v-1), ]
            trainZ <- ZFull$Z[, (ws1+h):(v-h)]         
        }else{
            
            
            trainY <- ZFull$Y[h:(v-1), ]

            trainZ <- ZFull$Z[,1:(v-h)]

        }

        ## needed.objs <- c('group','beta','trainZ','trainY','gamm','tol','p','m','k1','k','s','s1','m','MN','C','intercept','separate_lambdas','activeset','alpha','jj','jjcomp','kk','palpha','q1a')
        dual <- FALSE
        if(!group%in%c("SparseLag","SparseOO")){
            q1a <- NULL
        }
        ## objs <- sapply(needed.objs,exists)
        ## objs <- setdiff(needed.objs,ls())
        ## if(length(objs)>0){
        
        ##     for(i in 1:length(objs)){               
        ##         assign(objs[i],NULL)
        ##                }
        ## }
        ## browser()
        temp <- .BigVAR.fit(group,beta,trainZ,trainY,gamm,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
        beta <- temp$beta
        activeset <- temp$activeset
        q1a <- temp$q1a
        
        
        
        ## betaEVAL <- matrix(beta[,,1],nrow=k1,ncol=(k1*p+(k-k1)*(s+s1)+1))

        if (RVAR) {

            betaEVAL <- RelaxedLS(cbind(t(trainZ),trainY),betaEVAL,k,p,k1,s+s1)
            
        }
        

        if(MN ){
            
            eZ <- matrix(ZFull$Z[,v],ncol=1)

        }else{

            eZ <- matrix(c(1,ZFull$Z[,v]),ncol=1)

        }
        gamm_opt <- apply(MSFE,1,which.min)
        ## browser()
        abind::adrop(beta[,,1,drop=F],3)
        optind <- max(which(colMeans(na.omit(MSFE))==min(colMeans(na.omit(MSFE)))))

        gamopt <- gamm[optind]

        gamm_evolve <- c(gamm_evolve,gamopt)

        MSFE_temp <- matrix(0,nrow=1,ncol=length(gamm))
        ## browser()
        for(ii in 1:length(gamm)){
            betaEVAL <-  abind::adrop(beta[,,ii,drop=F],3)
            
            if(MN){

                if(ii==optind){

                    preds[v-T2+h-1,] <-  betaEVAL[,2:ncol(betaEVAL)] %*% eZ

                    ## MSFE_oos[v-T2+h-1] <- norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                    MSFE_oos[v-T2+h-1] <- .calc.loss(ZFull$Y[v+h-1,1:k1] - betaEVAL,univ=FALSE,loss,delta)
                }
                ## MSFE_temp[v-T2+h-1,ii] <- norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL[,2:ncol(betaEVAL)] %*% eZ)^2
                MSFE_temp[v-T2+h-1,ii] <-  .calc.loss(ZFull$Y[v+h-1,1:k1] - preds[v-T2+h-1,],univ=FALSE,loss,delta)
                ## preds[v-T2+h-1,] <-  betaEVAL[,2:ncol(betaEVAL)] %*% eZ

                diag(beta[,2:(k1+1),ii]) <- diag(beta[,2:(k1+1),ii])-C # subtract one for warm start purposes 
                ## }
            }else{
                
                ## if(separate_lambdas){

                ##         for(uu in 1:ncol(ZFull$Y)){
                ##             ## browser()
                ##             MSFE[v-T2+h-1,uu] <- norm2(ZFull$Y[v+h-1,uu] - preds[v-T2+h-1,uu])^2
                ##         }
                ## }else{
                if(ii==optind){
                    ## browser()
                    betaArray[,,v-T2+h-1] <- betaEVAL    
                    preds[v-T2+h-1,] <- betaEVAL %*% eZ                
                    ## MSFE_oos[v-T2+h-1] <- norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL %*% eZ)^2
                    MSFE_oos[v-T2+h-1,ii] <-  .calc.loss(ZFull$Y[v+h-1,1:k1] - preds[v-T2+h-1,],univ=FALSE,loss,delta)
                }
                ## if(v>T2-h+2){browser()}
                MSFE_temp[1,ii] <- .calc.loss(ZFull$Y[v+h-1,1:k1] - betaEVAL%*%eZ,univ=FALSE,loss,delta)## norm2(ZFull$Y[v+h-1,1:k1] - betaEVAL %*% eZ)^2
                
            }

        }
        MSFE <- rbind(MSFE,MSFE_temp)
        MSFE <- MSFE[2:nrow(MSFE),]
        if(verbose){

            setTxtProgressBar(pb, v)

        }
    }
    
    betaPred <- betaEVAL
    
    return(list(MSFE=MSFE_oos,betaPred=betaPred,predictions=preds,betaArray=betaArray,gamm_evolve=gamm_evolve))
    

}



                                        # Forecast evaluation: VAR (called in cv.bigvar)
.BigVAREVAL <- function(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T1,T2,alpha,recursive,C,intercept,tol,window.size,separate_lambdas,loss="L2",delta=2.5)
{
    

    gran2 <- 1

    preds <- matrix(NA,nrow=length((T1+1):T2),ncol=k)


    gamm <- gamopt
    if(separate_lambdas){
        gamm <- matrix(gamm,nrow=1,ncol=k)
    }
    Y <- ZFull$Y

    s <- p
    k1 <- k
    m=0;
    s1=0;
    dual <- FALSE
    ## if(separate_lambdas){
    ##     MSFE <- matrix(NA,nrow=length((T1+1):T2),ncol=k)
    ##     betaArray <-  array(0,dim=c(k,p*k+1,nrow(MSFE)))
    ## }else{
    MSFE <- rep(NA,length((T1+1):T2))
    betaArray <-  array(0,dim=c(k,p*k+1,length(MSFE)))
    ## }
    
    beta <- array(0,dim=c(k,p*k+1,1))

    kk <- NULL
    jj <- NULL
    jjcomp <- NULL
    activeset <- NULL
    q1a <- NULL

    if (group == "Lag"){

        jj <- .groupfuncpp(p,k)

        jjcomp <- .groupfuncomp(p,k)

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

    }else if (group == "SparseLag"){

        jj <- .groupfun(p, k)

        q1a <- list()

        for (i in 1:(p)) {

            q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)

        }

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

        

    }else if (group == "OwnOther"){


        kk <- .lfunction3cpp(p, k)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)

    }else if (group == "SparseOO") {

        kk <- .lfunction3cpp(p, k)

        jjcomp <- .lfunctioncomp(p,k)

        jj=.lfunction3(p,k)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)

        q1a <- list()

        for (i in 1:(2*p))
        {

            q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)

        }

    }else{

        kk=NULL
        jj=NULL
        activeset=NULL
        jjcomp=NULL
        
    }


    if(verbose)
    {

        print("Evaluation Stage")

        pb <- txtProgressBar(min = T1-h+2, max = T2, style = 3)

    }
    
    for (v in (T1-h+2):T2)
    {


        
        if(h>1 & !recursive){

            if(window.size!=0){
                ws1 <- max(c(v-window.size-h,1))
                trainY <- ZFull$Y[(ws1+h):(v-1), ]
                trainZ <- ZFull$Z[, (ws1+h):(v-h)]         
            }else{

                trainY <- ZFull$Y[(h):(v-1), ]
                
                trainZ <- ZFull$Z[, 1:(v-h)]
            }
            
        }else{
            if(window.size!=0){
                ws1 <- max(c(v-window.size,1))
                trainY <- ZFull$Y[(ws1):(v-1), ]
                trainZ <- ZFull$Z[, (ws1):(v-1)]         
            }else{
                trainY <- ZFull$Y[(1):(v-1), ]                       
                trainZ <- ZFull$Z[, (1):(v-1)]
            }
        }
        

        if(v+h-1>T2){
            break
        }

        ## needed.objs <- c('group','beta','trainZ','trainY','gamm','tol','p','m','k1','s','m','MN','c','intercept','separate_lambdas','dual','activeset','q1a')
        VARX <- FALSE

        dual <- FALSE
        ## needed.objs <- c('group','beta','trainZ','trainY','gamm','tol','p','MN','C','intercept','separate_lambdas','activeset','alpha','jj','jjcomp','kk','palpha')
        if(!group%in%c("SparseLag","SparseOO")){
            q1a <- NULL
        }
        ## ## objs <- sapply(needed.objs,exists)
        ## objs <- setdiff(needed.objs,ls())
        ## if(length(objs)>0){
        
        ##     for(i in 1:length(objs)){               
        ##         assign(objs[i],NULL)
        ##                }
        ## }
        
        ## browser()
        temp <- .BigVAR.fit(group,beta,trainZ,trainY,gamm,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
        beta <- temp$beta
        activeset <- temp$activeset
        q1a <- temp$q1a

        
        ## betaEVAL <- matrix(beta[,,1],nrow=k,ncol=(k*p+1))
        ## if(group!="BGR"){
        betaEVAL <- matrix(beta[,,1],nrow=k,ncol=(k*p+1))
        
        if (RVAR) {

            betaEVAL <- RelaxedLS(cbind(t(trainZ),trainY),betaEVAL,k,p,k1,s)      
        }
        
        if(MN){
            eZ <- matrix(ZFull$Z[,v],ncol=1)

        }else{


            eZ <- matrix(c(1,ZFull$Z[,v]),ncol=1)


        }

                                        # We don't consider an intercept for the MN lasso
        if(MN){

            ptemp <- betaEVAL[,2:ncol(betaEVAL)] %*% eZ
            if(h>1 & recursive){

                pred <- matrix(ptemp,nrow=1)
                
                preds[v-T1+h-1,] <- predictMS(pred,trainY,h-1,betaEVAL[,2:ncol(betaEVAL)],p,TRUE)

                MSFE[v-T1+h-1] <- .calc.loss(ZFull$Y[v+h-1] - preds[v-T1+h-1,],univ=FALSE,loss,delta)
                ## }
            }else{
                preds[v-T1+h-1,] <- ptemp 
                MSFE[v-T1+h-1] <- .calc.loss(ZFull$Y[v+h-1] - preds[v-T1-h+1,],univ=FALSE,loss,delta)
                diag(beta[,2:(k+1),1]) <- diag(beta[,2:(k+1),1])-C # subtract one for warm start purposes 

            }
            
        }else{
            preds[v-T1+h-1,] <- betaEVAL %*% eZ                       
            if(h>1 & recursive){

                pred <- matrix(preds[v-T1+h-1,],nrow=1)

                preds[v-T1+h-1,] <- predictMS(pred,trainY,h-1,betaEVAL,p,FALSE)

            }
            
            MSFE[v-T1+h-1] <- .calc.loss(ZFull$Y[v+h-1,] - preds[v-T1+h-1,],univ=FALSE,loss,delta)

            
        }
        if(verbose){


            setTxtProgressBar(pb, v)

        }
        betaArray[,,v-T1+h-1] <- beta
    }
    ## browser()
    temp <- .BigVAR.fit(group,beta,ZFull$Z,ZFull$Y,gamm,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
    betaPred <- temp$beta
    
    
    ## if(group!="BGR"){
    betaPred <- as.matrix(betaPred[,,1],drop=F)

    
    ## betaPred <- matrix(betaPred,nrow=k)
    ## if(group=="BGR"){
    ##     betaPred <- t(betaPred)
    ##     }
    ## }else{
    ##     betaPred <- matrix(0,nrow=1,ncol=1)
    ##     }

    ## betaPred <- as.matrix(betaPred[,,1])

    ## betaPred <- matrix(betaPred,nrow=k)

    

    return(list(MSFE=MSFE,betaPred=betaPred,predictions=preds,betaArray=betaArray))
    

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
                                        # We require that the coefficient matrix generates a stationary VAR
    if(max(Mod(eigen(Fp)$values))>1){warning("Coefficient Matrix is not stationary")}

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
MultVarSim <- function (k, A1, p, Sigma, T) 
{

    if(max(Mod(eigen(A1)$values))>1){stop("Error: Generator Matrix is not stationary")}
    
                                        # add 500 observations for initialization purposes

    Y <- matrix(0, nrow = T+500+p , ncol = k)

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

                                        # function to create subsets for lag group VARX-L
.groupfun <- function(p,k)
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

                                        #C++ groupings to account for 0 indexing
.groupfuncpp <- function(p,k)
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

                                        # subsetting complement of groups in rcpp
.groupfuncomp <- function(p,k)
{


    ownoth <- .groupfuncpp(p,k)

    kk2 <- list()

    pmax <- max(unlist(ownoth))

    to <- 0:(pmax)

    for(i in 1:length(ownoth))
    {

        kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]
    }

    return(kk2)

}    

                                        # Group indexing for own/other VARX-L
.lfunction2 <- function(p,k)
{

    kk <- list()

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

    kk <- list()

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
        j <- 0
        oo[[i]] <- kk[[i]][(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
        pp[[i]] <- kk[[i]][-(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
        j <- j+1

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
        j <- 0
        oo[[i]] <- kk[[i]][(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
        pp[[i]] <- kk[[i]][-(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
        j <- j+1
    }

    ownoth <- c(oo,pp)

    return(ownoth)

}    


.lfunctioncomp <- function(p,k)
{

    ownoth <- .lfunction3cpp(p,k)
    kk2 <- list()
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
.oofun <- function(p,k)
{

    kk <- .lfunction2(p, k)
    oo <- list()
    pp <- list()

    for (i in 1:length(kk)) {
        j <- 0
        oo[[i]] <- kk[[i]][(seq(1, length(kk[[1]]), k + 1) + 
                            (j * k^2))]
        pp[[i]] <- kk[[i]][-(seq(1, length(kk[[1]]), k + 1) + 
                             (j * k^2))]
        j <- j + 1
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
    oogroups <- list()
    oogroups[[1]] <- unlist(kk)
    for(i in 2:length(kk))
    {
        oogroups[[i]] <- unlist(kk[-(1:(i-1))])

    }
    
    return(oogroups)
}


                                        # indexing function for Own/Other HVAR
.vecoovars<-function(p,k,k1)
{

    vv <- list()

    vv[[1]] <- 1:(p*k)

    vv[[2]] <- vv[[1]][-k1]

    q1 <- 1

    if(p>1){
        for(i in 3:(2*p))
        {
            if(i%%2!=0)
            {
                vv[[i]] <- (q1*k+1):(k*p)
                q1 <- q1+1
            }else{
                vv[[i]] <- vv[[i-1]][-k1]
            }
        }
    }
    return(vv)

}

                                        # indexing to start at zero for use within rcpp
.vecoovarscpp<-function(p,k,k1)
{

    vv <- list()
    vv[[1]] <- 0:(p*k-1)
    vv[[2]] <- vv[[1]][-(k1)]
    q1 <- 1

    if(p>1){
        for(i in 3:(2*p))
        {
            if(i%%2!=0)
            {
                vv[[i]] <- (q1*k):(k*p-1)
                q1 <- q1+1
            }else{
                
                vv[[i]] <- vv[[i-1]][-(k1)]

            }

        }

    }

    return(vv)

}


                                        #VARX Lag Group function
groupfunVARX <- function(p,k,k1,s)
{

    jj <- list()
    m <- k-k1
    jj <- .groupfuncpp(p, k1)
    kp <- k1*p+m*s-1
    jj2 <- list()
    startjj <- max(unlist(jj))+1
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
    kk2 <- list()
    pmax <- max(unlist(ownoth))
    to <- 0:(pmax)
    for(i in 1:length(ownoth))
    {
        kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]
    }
    return(kk2)

}


                                        # indexing starts at 1 for use in R        
groupfunVARXLG <- function(p,k,k1,s)
{

    jj <- list()
    jj <- .groupfun(p, k1)
    m <- k-k1
    kp <- k1*p+m*s
    jj2 <- list()
    startjj <- max(unlist(jj))+1
    for(i in seq(startjj,kp,by=1))
    {
        jj[[i]] <- i
        
    }
    jj[sapply(jj, is.null)] <- NULL

    return(jj)
}



diaggroupfunVARX <- function(p,k,k1,s)
{
    m <- k-k1
    jj <- list()
    jj <- .lfunction3cpp(p, k1)
    kp <-k1*(p*k1+s*m)-1
    jj2 <- list()
    startjj <- max(unlist(jj))+1

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
    kk2 <- list()
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

    m <- k-k1
    jj <- list()
    jj <- .lfunction3(p, k1)
    kp <- k1*(p*k1+s*m)
    jj2 <- list()
    startjj <- max(unlist(jj))+1
    for(i in seq(startjj,kp,by=k1))
    {
        jj[[i]] <- i:(i+(k1-1))                
    }
    
    jj[sapply(jj, is.null)] <- NULL

    return(jj)

}

diaggroupfunVARXLGL <- function(p,k,k1)
{

    jj <- list()
    jj <- .lfunction3(p, k1)
    kp <- k1*p*k
    jj2 <- list()
    startjj <- max(unlist(jj))+1

    for(i in seq(startjj,kp,by=1))
    {
        jj[[i]] <- i                
    }
    
    jj[sapply(jj, is.null)] <- NULL

    return(jj)
}


diaggroupfunVARXL <- function(p,k,k1)
{
    jj <- list()
    jj <- .lfunction3cpp(p, k1)
    kp <- k1*p*k-1
    jj2 <- list()
    startjj <- max(unlist(jj))+1

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
    kk2 <- list()
    pmax <- max(unlist(ownoth))

    to <- 0:(pmax)

    for(i in 1:length(ownoth))

    {

        kk2[[i]] <- to[is.na(pmatch(to,ownoth[[i]]))]

    }

    return(kk2)

}


                                        # iterative procedure to find a less coarse bound for lambda starting value via binary search
LGSearchX <- function(gstart,Y,Z,BOLD,group,k1,p,s,m,gs,k,MN,alpha,C,intercept,tol)
{
    s1=0;palpha=NULL
    tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
    lambdah <- gstart
    lambdal <- 0
    activeset <- list(rep(rep(list(0), length(gs))))

    ## if(group=="SparseLag"){
    

    ##     jj <- groupfunVARX(p, k,k1,s)

    ##     jjcomp <- groupfunVARXcomp(p,k,k1,s)

    ##     q1a=list()

    ##             for (i in 1:(p+s)) {

    ##                 q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)

    ##             }
    ## }

    gran2=1
    kk <- NULL
    jj <- NULL
    jjcomp <- NULL
    ## activeset <- NULL
    q1a <- NULL

    if (group == "Lag") {

        jj <- groupfunVARX(p,k,k1,s+s1)

        jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

    }else if (group == "SparseLag") {

        jj <- groupfunVARX(p, k,k1,s+s1)

        q1a <- list()

        for (i in 1:(p+s+s1)) {

            q1a[[i]] <- matrix(runif(k1, -1, 1), ncol = 1)

        }

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)


    }else if (group == "OwnOther") {

        kk <- diaggroupfunVARX(p, k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
    }else if (group == "SparseOO") {

        kk <- diaggroupfunVARX(p, k,k1,s+s1)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)
    }else{

        jj <- NULL
        kk <- NULL
        activeset <- NULL
    }
    
    
    while(max(lambdah-lambdal)>10*tol)
    {

        lambda <- (lambdah+lambdal)/2
        ## if(group=="EFX"){
        ##     BOLD <- .EFVARX(BOLD,Y,Z,lambda,tol,MN,k1,s,m,p,C,intercept)
        ##     param <- BOLD[,2:(k1*p+m*s+1),1]
        ## }

        ## if(group=="Basic"){
        ##     param <- .lassoVARFistX(BOLD,Z,Y[,1:k1],lambda,tol,p,MN,k1+m,k1,s,m,C,intercept)[,2:(k1*p+m*s+1),]
        ## }

        ## if(group=="Lag"){

        ##     jj <- groupfunVARX(p,k,k1,s)
        ##     jjcomp <- groupfunVARXcomp(p,k,k1,s)
        ##     BB <- .GroupLassoVAR1(BOLD,jj,jjcomp,Y[,1:k1],Z,lambda,activeset,tol,p,MN,k,k1,s,C,intercept)
        ##     BOLD <- BB$beta
        ##     param <- BB$beta[,2:(k1*p+m*s+1),]
        ##     activeset <- BB$active
        ## }
        

        ## if(group=="OwnOther")
        ##     {

        ##         kk <- diaggroupfunVARX(p, k,k1,s)
        ##         BB <- .GroupLassoOOX(BOLD, kk, Y, Z, lambda,activeset, tol,p,MN,k,k1,s,C,intercept)
        ##         param <- BB$beta[,2:(k1*p+m*s+1),]
        ##         BOLD <- BB$beta
        ##         activeset <- BB$active
        ##     }

        ## if(group=="SparseOO")
        ##     {

        ##         kk <- diaggroupfunVARX(p, k,k1,s)
        ##         BB <- .SparseGroupLassoVAROOX(BOLD, kk, Y[,1:k1], Z, lambda,alpha,activeset, tol,p,MN,k1,s,k,FALSE,C,intercept)
        ##         param <- BB$beta[,2:(k1*p+m*s+1),]
        ##         BOLD <- BB$beta
        ##         activeset <- BB$active

        ##     }
        

        ## if(group=="SparseLag"){

        ##     jj <- groupfunVARX(p, k,k1,s)
        ##     jjcomp <- groupfunVARXcomp(p,k,k1,s)
        ##     q1a=list()

        ##     for (i in 1:(p+s)) {

        ##         q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)

        ##     }
        
        ##     BB <- .SparseGroupLassoVARX(BOLD,jj,Y[,1:k1],Z,lambda,alpha,activeset,tol,q1a,p,MN,k,s,k1,C,intercept)
        ##     param <- BB$beta[,2:(k1*p+m*s+1),]
        ##     BOLD <- BB$beta
        ##     activeset <- BB$active
        ## }


        ## needed.objs <- c('group','BOLD','tol','p','m','k','k1','s','s1','m','MN','C','intercept','separate_lambdas','activeset','palpha')

        jj <- groupfunVARX(p, k,k1,s)

        kk <- diaggroupfunVARX(p, k,k1,s)

        jjcomp <- groupfunVARXcomp(p,k,k1,s)

        VARX <- TRUE
        dual <- FALSE
        separate_lambdas <- FALSE
        if(!group%in%c("SparseLag")){
            q1a <- NULL
        }## else{
        ##     q1a <- list()
        ##     for (i in 1:(p+s)) {

        ##         q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)

        ##     }
        ## }
        ## browser()
        ## ae=function(x){any(exists(x))}
        ## objs <- sapply(needed.objs,ae)
        ## if(!all(objs)){
        ##     oss <- names(objs)[!objs]
        ##     for(i in 1:length(oss)){               
        ##         assign(oss[i],NULL)
        ##     }
        ## }
        ## browser()
        temp <- .BigVAR.fit(group,BOLD,Z,Y,lambda,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
        BOLD <- temp$beta
        param <- BOLD[,2:(k1*p+m*s+1),]
        activeset <- temp$activeset
        q1a <- temp$q1a

        if(MN){
            ## diag(param[1:k1,1:k1]) <- 0
            ## diag(BOLD[1:(k1),2:(k1+1),1]) <- 0

            
            diag(param[1:k1,1:k1]) <- ifelse(C==0,diag(param[1:k1,1:k1]),diag(param[1:k1,1:k1])-C)
            diag(BOLD[,2:(k1*p+1+m*s),]) <- ifelse(C==0,diag(BOLD[,2:(k1*p+1+m*s),]),diag(BOLD[,2:(k1*p+1+m*s),])-C)

        }
        ## browser()
        if(max(abs(param))<sqrt(.Machine$double.eps))
        {
            lambdah <- lambda
        }else{

            lambdal <- lambda
            
        }

    }

    lambdah
}

                                        # Same as above, but for the VAR
                                        # no special handling for separate lambdas since it is called separately for each one
LGSearch <- function(gstart,Y,Z,BOLD,group,k,p,gs,MN,alpha,C,intercept,tol)
{

    m=0;s=0;s1=0
    tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
    lambdah <- gstart
    lambdal <- 0
    activeset <- list(rep(rep(list(0), length(gs))))

    ## if(group=="SparseOO")
    ##     {
    ##         kk <- .lfunction3cpp(p, k)
    ##         activeset <- rep(list(rep(rep(list(0), length(kk)))),1)
    ##         q1a=list()
    ##         for (i in 1:(2*p)) {
    ##             q1a[[i]] <- matrix(runif(length(gs[[i]]), -1, 1), ncol = 1)
    ##         }

    ##     }

    ## if(group=="SparseLag"){
    ##     q1a <- list()
    ##     jj <- .groupfuncpp(p, k)

    ##             for (i in 1:(p)) {

    ##                 q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)

    ##             }
    ## }
    gran2=1
    kk=NULL
    jj=NULL
    ## activeset=NULL
    jjcomp=NULL
    q1a <- NULL
    if (group == "Lag"){

        
        jj <- .groupfuncpp(p,k)

        jjcomp <- .groupfuncomp(p,k)

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

    }else if (group == "SparseLag"){

        jj <- .groupfun(p, k)

        q1a <- list()

        for (i in 1:(p)) {

            q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)

        }

        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                         gran2)

        

    }else if (group == "OwnOther"){


        kk <- .lfunction3cpp(p, k)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)

    }else if (group == "SparseOO") {

        kk <- .lfunction3cpp(p, k)

        jjcomp <- .lfunctioncomp(p,k)

        jj=.lfunction3(p,k)

        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                         gran2)

        q1a <- list()

        for (i in 1:(2*p))
        {

            q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)

        }

    }else{

        kk=NULL
        jj=NULL
        activeset=NULL
        jjcomp=NULL
        
    }

    
    while(max(abs(lambdah-lambdal))>10*tol)
    {

        lambda <- (lambdah+lambdal)/2
        ## if(group=="Basic"){
        ##     param <- .lassoVARFist(BOLD,Z,Y,lambda,tol,p,MN,C,intercept)[,2:(k*p+1),]
        ## }

        ## if(group=="Tapered"){

        ##     param <- .lassoVARTL(BOLD,Z,Y,lambda,tol,p,MN,rev(seq(0,1,length=10)),C,intercept)[,2:(k*p+1),]

        ##     param <- param[,,1]                    
        ## }
        
        ## if(group=="Lag"){

        ##     jj <- .groupfuncpp(p, k)
        ##     jjcomp <- .groupfuncomp(p,k)
        ##     BB <- .GroupLassoVAR1(BOLD,jj,jjcomp,Y,Z,lambda,activeset,tol,p,MN,k,k,p,C,intercept)
        ##     BOLD <- BB$beta
        ##     param <- BB$beta[,2:(k*p+1),]
        ##     activeset <- BB$active                  
        ## }

        ## if(group=="SparseLag"){

        ##     jj <- .groupfuncpp(p, k)
        ##     jjcomp <- .groupfuncomp(p,k)
        ##     q1a=list()
        ##     for (i in 1:p) {
        ##         q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
        ##     }
        ##     BB <- .SparseGroupLassoVAR(BOLD,jj,Y,Z,lambda,alpha,activeset,tol,q1a,p,MN,C,intercept)
        ##     param <- BB$beta[,2:(k*p+1),]
        ##     BOLD <- BB$beta
        ##     activeset <- BB$active
        ## }
        

        ## if(group=="OwnOther")
        ##     {
        ##         kk <- .lfunction3cpp(p, k)
        ##         BB <- .GroupLassoOO(BOLD, kk, Y, Z, lambda,activeset, tol,p,MN,C,intercept)
        ##         param <- BB$beta[,2:(k*p+1),]
        ##         BOLD <- BB$beta
        ##         activeset <- BB$active

        ##     }


        ## if(group=="SparseOO")
        ##     {
        ##         BB <- .SparseGroupLassoVAROO(BOLD, kk, Y, Z, lambda,alpha,activeset, tol,q1a,p,MN,FALSE,C,intercept)
        ##         param <- BB$beta[,2:(k*p+1),]
        ##         BOLD <- BB$beta
        ##         activeset <- BB$active                  
        ##     }

        ## if(group=="HVARC")
        ##     {
        ##         BOLD <- .HVARCAlg(BOLD,Y,Z,lambda,tol,p,MN,C,intercept)                  
        ##         param <- BOLD[,2:(k*p+1),]
        ##     }
        ## if(group=="HVAROO")
        ##     {
        ##         BOLD <- .HVAROOAlg(BOLD,Y,Z,lambda,tol,p,MN,C,intercept)
        ##         param <- BOLD[,2:(k*p+1),]
        ##     }

        ## if(group=="HVARELEM")
        ##     {
        ##         BOLD <- .HVARElemAlg(BOLD,Y,Z,lambda,tol,p,MN,C,intercept)
        ##         param <- BOLD[,2:(k*p+1),]
        ##     }

                                        # define just in case
        ## kk <- .lfunction3cpp(p, k)

        ## jj <- .groupfuncpp(p, k)

        ## jjcomp <- .groupfuncomp(p,k)
        ## needed.objs <- c('group','BOLD','lambda','tol','p','m','k1','k','s','s1','m','MN','C','intercept','activeset','palpha')

                                        # arguments used in generic BigVAR.fit not applicable here (dual, separate_lambdas,VARX) set to default values
        if(!group%in%c("SparseLag","SparseOO")){
            q1a <- NULL
        }else{


        }
        VARX <- FALSE
        k1=k
        separate_lambdas <- FALSE
        dual <- FALSE
        ## ae=function(x){any(exists(x))}
        ## objs <- sapply(needed.objs,ae)
        ## if(!all(objs)){
        ##     oss <- names(objs)[!objs]
        ##     for(i in 1:length(oss)){               
        ##         assign(oss[i],NULL)
        ##     }
        ## }
        ## browser()
        ## palpha <- rev(seq(0,1,length=10))
        palpha=0
        temp <- .BigVAR.fit(group,BOLD,Z,Y,lambda,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
        if(group=="Tapered"){
            param <- BOLD[,2:(k*p+1),1]
        }else{

            BOLD <- temp$beta
            param <- BOLD[,2:(k*p+1),]
        }
        activeset <- temp$activeset
        q1a <- temp$q1a

        
        if(MN){

            diag(param[1:k,1:k]) <- ifelse(C==0,diag(param[1:k,1:k]),0)
            if(group!="Tapered"){
                diag(BOLD[,2:(k*p+1),]) <- ifelse(C==0,diag(BOLD[,2:(k*p+1),]),0)
            }
        }
        ## browser()
        if(max(abs(param))< tol)
        {
            lambdah <- lambda

        }else{
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
#' @param IC specifies whether to select lag order according to "AIC" or "BIC"
#' @param h desired forecast horizon
#' @param loss loss function (default "L2", one of "L1","L2","Huber")
#' @param delta delta for Huber loss function (default 2.5)
#' @param iterated indicator as to whether to use iterated or direct multistep forecasts (if applicable, VAR context only)
#' @return Returns the one-step ahead MSFE as well as the forecasts over the evaluation period and lag order selected.
#' @details This function evaluates the one-step ahead forecasts of a VAR or VARX fit by least squares over an evaluation period.  At every point in time, lag orders for the endogenous and exogenous series are selected according to AIC or BIC.  This function is run automatically when \code{\link{cv.BigVAR}} is called unless \code{ic} is set to \code{FALSE} in \code{\link{constructModel}}.      
#' @references Neumaier, Arnold, and Tapio Schneider. "Estimation of parameters and eigenmodes of multivariate autoregressive models." ACM Transactions on Mathematical Software (TOMS) 27.1 (2001): 27-57.
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
#' BICMSFE <- VARXForecastEval(Y,X,p,0,T1,T2,"BIC",1)
#' 
#' @export
VARXForecastEval <- function(Y,X,p,s,T1,T2,IC,h,iterated=FALSE,loss="L2",delta=2.5)
{


    if(T1>nrow(Y) | T2>nrow(Y) |T2<T1){stop("Training dates exceed series length")}

    if(!IC%in%c("AIC","BIC") )
    {

        stop("IC must either be AIC or BIC")

    }

    MSFE <- c()
    predF <- NULL
    pvec <- NULL
    svec <- NULL
    k <- ncol(Y)
    m <- ifelse(s!=0,ncol(X),0)
    for(i in (T1-h+2):T2){

        if(h+i-1>T2){break}

        testY <- as.matrix(Y[1:(i-1),])
        testX <- as.matrix(X[1:(i-1),])

        if(!iterated){
            hd=h
            
        }else{

            hd=1
            
        }
        if(IC=="BIC"){
            popt <- ICX(testY,testX,k,p,s,m,"BIC",h=hd)
        }
        if(IC=="AIC"){

            popt <- ICX(testY,testX,k,p,s,m,"AIC",h=hd) 
        }
        B1 <- popt$B
        
        if(popt$p==0 &popt$s==0){
            
            eZ <- matrix(rep(1,1),ncol=1)

            pred <- B1%*%eZ
            
        }else{


            C <- max(popt$p,popt$s)

            ## print(C)
            ## if(C==1){

            ## # possibly memory leak in VARX lag matrix construction in Eigen if maxlag is 1.
            ## # to be on the safe side, we will perform it in R
            ## if(popt$s==0){
            ## eZ <- c(1,Y[i-1,])
            ## eZ <- matrix(eZ,ncol=1)
            ## }else{
            ##     eZ <- c(1,Y[i-1,],X[i-1,])               
            ##     }
            ## }else{
            eZ <- VARXCons(as.matrix(Y[(i-C):(i),]),as.matrix(X[(i-C):(i),]),k,popt$p,m,popt$s)

            ## }
            
            pred <- B1%*%eZ
            
                                        # iterated multistep forecasts (if VAR and horizon greater than 1)
            if(h>1 & s==0 & iterated ){
                
                pred <- predictMS(matrix(pred,nrow=1),Y,h-1,B1,C,FALSE)

            }

        }
        ## browser()
        predF <- rbind(predF,t(pred))
        MSFEi <-  .calc.loss(Y[i+h-1,] - pred,univ=FALSE,loss,delta)

        ## MSFEi <- norm2(Y[i+h-1,]-pred)^2

        MSFE <- c(MSFE,MSFEi)
        svec <- c(svec,popt$s)
        pvec <- c(pvec,popt$p)
    }

    return(list(MSFE=MSFE,pred=as.matrix(predF),p=pvec,s=svec))


}

#' Fit a VAR or VARX model by least squares
#' 
#' @param Y a \eqn{t \times k} multivariate time series
#' @param p maximum lag order
#' @param IC Information criterion indicator, if set to \code{NULL}, it will fit a least squares VAR(X) of orders p and s.  Otherwise, if set to "AIC" or "BIC" it return the model with lag orders that minimize the given IC. 
#' @param VARX a list of VARX specifications (as in \code{\link{constructModel}} (or NULL )
#' @return Returns a list with four entries:
#' \itemize{
#' \item{"Bhat"}{Estimated \eqn{k\times kp+ms} coefficient matrix}
#' \item{"SigmaU}{Estimated \eqn{k\times k} residual covariance matrix}
#' \item{"phat"}{Selected lag order for VAR component}
#' \item{"shat"}{Selected lag order for VARX component}
#' \item{"Y"}{multivariate time series retained for prediction purposes}
#' \item{"Y"}{number of endogenous (modeled) time series}
#' }
#' @details This function uses a modified form of the least squares technique proposed by Neumaier and Schneider (2001).  It fits a least squares VAR or VARX via a QR decomposition that does not require explicit matrix inversion.  This results in improved computational performance as well as numerical stability over the conventional least squares approach. 
#' @references Neumaier, Arnold, and Tapio Schneider. "Estimation of parameters and eigenmodes of multivariate autoregressive models." ACM Transactions on Mathematical Software (TOMS) 27.1 (2001): 27-57.
#' @seealso \code{\link{constructModel}}, \code{\link{cv.BigVAR}},\code{\link{BigVAR.fit}}
#' @examples
#' data(Y)
#' # fit a VAR_3(3)
#' mod <- VARXFit(Y,3,NULL,NULL)
#' # fit a VAR_3 with p= 6 and lag selected according to AIC
#' modAIC <- VARXFit(Y,6,"AIC",NULL)
#' # Fit a VARX_{2,1} with p=6, s=4 and lags selected by BIC
#' modXBIC <- VARXFit(Y,6,"BIC",list(k=1,s=4))
#' 
#' @export
VARXFit <- function(Y,p,IC,VARX=NULL)
{

    if(!is.null(VARX)){

        if(is.list(VARX) & !(exists('k',where=VARX) & exists('s',where=VARX)))
        {

            stop("VARX Specifications entered incorrectly")

        }

    }
    if(is.list(VARX) & (length(VARX)!=0)){

        k1 <- VARX$k
        s <- VARX$s
        Y1 <- matrix(Y[,1:k1],ncol=k1)
        m <- ncol(Y)-k1
        X <- matrix(Y[,(k1+1):ncol(Y)],ncol=m)
        Z <- VARXCons(Y1,X,k1,p,m,s)
        offset <- max(p,s)+1
        YT <- matrix(Y1[offset:nrow(Y),],ncol=k1)
        X <- matrix(X[offset:nrow(X),],ncol=m)
        
    }else{

        k <- ncol(Y)
        k1 <- k
        s=0;m=0
        offset <- p+1
        X <- matrix(0,nrow=nrow(Y))
        
        Z <- VARXCons(Y,X,k,p,m,s)
        YT <- matrix(Y[(offset):nrow(Y),],ncol=ncol(Y))
        
    }
    if(is.null(IC)){

        Res <- ARFitVARXR(cbind(t(Z),YT),k1,p,m,s)

        shat <- s
        phat <- p
        
    }else{
        if(!IC%in%c("AIC","BIC") )
        {

            stop("IC must either be AIC,BIC, or set to NULL")

        }

        Res <- ICX(YT,X,k1,p,s,m,IC)

        shat <- Res$s
        phat <- Res$p
        
    }
    if(is.null(VARX)){
        k=ncol(Y)
    }else{
        k=VARX$k
    }
    list(Bhat=Res$B,SigmaU=Res$SigmaU,phat=phat,shat=shat,Y=Y,k=k)
    
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
PredictVARX <- function(VARXRes){

    B <- VARXRes$Bhat
    Y <- VARXRes$Y
    k <- VARXRes$k
    m <- ncol(Y)-k
    ## browser() 
    
    if(k<ncol(Y)){
        Z <- VARXCons(Y[,1:k,drop=FALSE],Y[,(k+1):ncol(Y),drop=FALSE],k,VARXRes$phat,m,VARXRes$shat,oos=TRUE)
    }else{
        Z <- VARXCons(Y[,1:k,drop=FALSE],matrix(0,nrow=nrow(Y)),k,VARXRes$phat,m,VARXRes$shat,oos=TRUE)
    }
    
    return(as.numeric(tail(t(B%*%Z),1)))
    
    
}

                                        # Recursive multi-step predictions


predictMS <- function(pred,Y,n.ahead,B,p,MN=FALSE){

                                        # Augment Y with predictions, create lag matrix (no intercept if MN)
    Y <- rbind(Y,pred)

                                        # Can't call this function direclty from R due to assert errors
    ## Z <- ZmatF(Y,p,ncol(Y),oos=TRUE,intercept=!MN)
    
    Z <- VARXCons(Y,matrix(0,nrow=nrow(Y),ncol=1),ncol(Y),p,0,0,oos=TRUE)
    if(MN){
        Z <- Z[2:nrow(Z),,drop=F]
    }
    Z <- Z[,ncol(Z),drop=F]

    pred <- matrix(B%*%Z,ncol=ncol(Y),nrow=1)

    if(n.ahead==1){return(pred)}
    
    predictMS(pred,Y,n.ahead-1,B,p,MN)

}

                                        # Multi-step VARX with new data.
predictMSX <- function(pred,Y,n.ahead,B,p,newxreg,X,m,s,cumulative,MN,contemp=FALSE){

    Y <- rbind(Y,pred)
    X <- rbind(X,matrix(newxreg[cumulative,],ncol=m))
    
    if(nrow(Y)!=nrow(X)){stop("error, dimension issue")}
    if(!contemp){
        Z <- VARXCons(as.matrix(Y),X,ncol(Y),p,m,s,oos=TRUE)
    }else{
        ## ## browser()
        Z <- VARXCons(as.matrix(Y),as.matrix(X),ncol(Y),p,m,s,oos=FALSE,contemp=TRUE)
    }
    Z <- Z[,ncol(Z),drop=F]
    if(MN){

        Z <- as.matrix(Z[2:nrow(Z),drop=F])
        pred <- matrix(B[,2:ncol(B),drop=F]%*%Z,ncol=ncol(Y),nrow=1)

    }else{
        
        pred <- matrix(B%*%Z,ncol=ncol(Y),nrow=1)
    }
    if(n.ahead==1){return(pred)}

    predictMSX(pred,Y,n.ahead-1,B,p,newxreg,X,m,s,cumulative+1,MN)    

}




                                        # Find optimal values in 2-d gridsearch
findind <- function(opt,lambda1,lambda2)
{
    if(opt<length(lambda2)){
        lambda1ind <- 1
    }else{
        lambda1ind <- ceiling(opt/length(lambda2))
    }
    if(lambda1ind==1)
    {
        jind <- opt
    }else{
        jind <- opt-(length(lambda2))*(lambda1ind-1)
    }
    return(c(lambda1ind,jind))
}




BVARLitterman <- function(Y,Z,p,tau,mu,H,iRW)
{
    T <- nrow(Y); k <- ncol(Y)

    ## browser()
                                        # prior covariance based on univariate AR models
    sigmas <- c()
    for(i in 1:k){
        Z1 <- VARXCons(Y[,i,drop=F],matrix(0,nrow=nrow(Y),ncol=1),1,p,0,0)
        ## Y1 <- matrix(Y[(p+1):nrow(Y),i],ncol=1)
                                        # get the prior cov
        K <- cbind(t(Z1),Y[(p+1):nrow(Y),i])
        sigmas[i] <- sqrt(ARFitVARXR(K,1,p,0,0)$SigmaU)    
        ## print(sigmas[i])
    }
    ## browser()
    MMO <- colMeans(Y)



    ## Z <- VARX(Y,matrix(0,nrow=Y,ncol=1),k,p,0,0)
    
    ## Y1 <- matrix(Y[(p+1):nrow(Y),],ncol=k)

                                        # create prior random walk dummy
    Yrw1 <- diag(sigmas*iRW)
    Yrw2 <- matrix(0,nrow=k*(p-1),ncol=k)
    Yrw<- tau*(rbind(Yrw1,Yrw2))
    Zrw <- tau*cbind(kronecker(diag(1:p),diag(sigmas)),matrix(0,nrow=k*p,ncol=1))


                                        # create dummy for intercept
    epsilon=1e-5
    Ycs <- 1e-5*matrix(0,nrow=1,ncol=k)
    Zcs <- epsilon*cbind(matrix(0,ncol=k*p,nrow=1),1)

                                        # dummy on the sums of coefficients
    Ylr <- mu*diag(MMO*iRW)
    Zlr1 <- kronecker(matrix(1,nrow=1,ncol=p),diag(MMO)*iRW)
    Zlr <- mu*(cbind(Zlr1,matrix(0,nrow=k,ncol=1)))


                                        # Dummy for residual covariance matrix
    Ycv <- diag(sigmas)
    Zcv <- matrix(0,nrow=k,ncol=k*p+1)

    Yprior <- rbind(Yrw,Ylr,Ycv,Ycs)
    Zprior <- rbind(Zrw,Zlr,Zcv,Zcs)

    ## Zp2 <- rev(Zprior[(nrow(Zprior)),])
    ## Zp3 <- Zprior[1:(nrow(Zprior)-1),]
    ## Zp3 <- rbind(Zp2,Zp3)
    ## ZZI2 <- solve(t(Zp3)%*%Zp3+t(Z)%*%Z)
    ## dim(ZZI2)
    Tstar <- nrow(Yprior)
    ## dim(Zprior)
    ## # posterior
    ## dim(Z)
    ## browser()
    Z <- t(Z)
    Z <- cbind(Z[,2:ncol(Z)],1)
    ## dim(Yprior)
    ## browser()        

    ## kappa(crossprod(Z))
    ZZinv <- solve(t(Zprior)%*%Zprior+t(Z)%*%Z)
    ZY <- t(Zprior)%*%Yprior+t(Z)%*%Y
    beta <- ZZinv%*%ZY
    ## browser()
    ## preds <- Z[nrow(Z),]%*%beta

    return(t(beta))

}


BGRGridSearch <- function(Y,Z,p,grid,RWIND)
{
    preds <- list()
    for(i in 1:length(grid))
    {
        pi <- grid[i]
        mu=pi*.1 # used in BGR paper
        preds[[i]] <- BVARLitterman(Y,Z,p,pi,mu,-1,RWIND)

    }
    preds <- array(unlist(preds), dim = c(nrow(preds[[1]]), ncol(preds[[1]]), length(preds)))

    return(preds)

}

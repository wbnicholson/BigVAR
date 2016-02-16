#' BigVAR Object Class
#'
#' An object class to be used with cv.BigVAR
#' 
#' @slot Data a \eqn{T x k} multivariate time Series
#' @slot lagmax Maximal lag order for modeled series
#' @slot Structure Penalty Structure
#' @slot Relaxed Indicator for relaxed VAR
#' @slot Granularity Granularity of Penalty Grid
#' @slot horizon Desired Forecast Horizon
#' @slot crossval Cross-Validation Procedure
#' @slot Minnesota Minnesota Prior Indicator
#' @slot verbose Indicator for Verbose output
#' @slot ic Indicator for including AIC and BIC benchmarks
#' @slot VARX VARX Model Specifications
#' @slot T1 Index of time series in which to start cross validation
#' @slot T2  Index of times series in which to start forecast evaluation
#' @slot ONESE Indicator for "One Standard Error Heuristic"
#' @slot ownlambdas Indicator for user-supplied lambdas

#' @details To construct an object of class BigVAR, use the function "ConstructModel"
#' @seealso \code{\link{constructModel}}
#' @export
 setClass(
Class="BigVAR",
  representation(
    Data="matrix",
    lagmax="numeric",
    Structure="character",
    Relaxed="logical",
    Granularity="numeric",
    Minnesota="logical",  
    horizon="numeric",
    verbose="logical",  
    crossval="character",
    ic="logical",
    VARX="list",
    T1="numeric",
    T2="numeric",
    ONESE="logical",
    ownlambdas="logical",
    tf="logical"  
  )
  )


#' Construct an object of class BigVAR
#' @param Y \eqn{T x k} multivariate time series or Y \eqn{T x k+m} endogenous and exogenous series, respectively 
#' @param p Predetermined maximal lag order (for modeled series)
#' @param struct The choice of penalty structure (see details).
#' @param gran vector containing how deep to construct the penalty grid (parameter 1) and how many gridpoints to use (parameter 2)  If ownlambas is set to TRUE, gran denotes the user-supplied penalty parameters.
#' @param RVAR True or False: whether to refit using the Relaxed-VAR procedure
#' @param h Desired forecast horizon
#' @param cv Cross-validation approach, either "Rolling" for rolling cross-validation or "LOO" for leave-one-out cross-validation.
#' @param MN Minnesota Prior Indicator
#' @param verbose Verbose output while estimating
#' @param IC True or False: whether to include AIC and BIC benchmarks
#' @param VARX List containing VARX model specifications. 
#' @param T1 Index of time series in which to start cross validation
#' @param T2  Index of times series in which to start forecast evaluation
#' @param ONESE True or False: whether to use the "One Standard Error Heuristic"
#' @param ownlambdas True or False: Indicator for user-supplied penalty parameters"


#' @details The choices for "struct" are as follows
#' \itemize{
#' \item{  "None" (Lasso Penalty)}
#' \item{  "Lag" (Lag Group Lasso)} 
#' \item{  "SparseLag" (Lag Sparse Group Lasso)} 
#' \item{  "Diag" (Own/Other Group Lasso) }
#' \item{  "SparseDiag" (Own/Other Sparse Group Lasso) }
#' \item{  "EFX" (Endogenous First VARX)}
#' \item{  "HVARC" (Componentwise Hierarchical Group Lasso) }
#' \item{  "HVAROO" (Own/Other Hierarchical Group Lasso) }
#' \item{  "HVARELEM" (Elementwise Hierarchical Group Lasso)}
#' \item{  "Tapered" (Lag weighted Lasso)}
#' }
#'
#' VARX specifications consist of a list with entry k1 denoting the series that are to be modeled and entry s to denote the maximal lag order for exogenous series.

#' @note The specifications "None", "Lag," "SparseLag," "SparseDiag," and "Diag" can accomodate both VAR and VARX models.  EFX only applies to VARX models.  "HVARC," "HVAROO," "HVARELEM," and "Tapered" can only be used with VAR models.
#'
#' @references William B Nicholson, Jacob Bien, and David S Matteson. "Hierarchical vector autoregression." arXiv preprint 1412.5250, 2014.
#' 
#' William B Nicholson, David S. Matteson, and Jacob Bien (2015), "VARX-L Structured regularization for large vector
#' autoregressions with exogenous variables," \url{http://www.wbnicholson.com/Nicholsonetal2015.pdf}.
#' 
#' @examples
#' library(BigVAR)
#' # VARX Example
#' VARX=list()
#' VARX$k=2 # indicates that the first two series are modeled
#' VARX$s=2 # sets 2 as the maximal lag order for exogenous series
#' data(Y)
#' T1=floor(nrow(Y)/3)
#' T2=floor(2*nrow(Y)/3)
#' Model1=constructModel(Y,p=4,struct="None",gran=c(50,10),verbose=FALSE,VARX=VARX,T1=T1,T2=T2)
#' @export
constructModel <- function(Y,p,struct,gran,RVAR=FALSE,h=1,cv="Rolling",MN=FALSE,verbose=TRUE,IC=TRUE,VARX=list(),T1=floor(nrow(Y)/3),T2=floor(2*nrow(Y)/3),ONESE=FALSE,ownlambdas=FALSE)
  {
if(any(is.na(Y))){stop("Remove NA values before running ConstructModel")}      
if(dim(Y)[2]>dim(Y)[1] & length(VARX)==0){warning("k is greater than T, is Y formatted correctly (k x T)?")}      
if(p<0){stop("Maximal lag order must be at least 0")}
if(p==0& struct!="None"){stop("Only Basic VARX-L supports a transfer function")}
structures=c("None","Lag","SparseLag","Diag","SparseDiag","HVARC","HVAROO","HVARELEM","Tapered","EFX")
cond1=struct %in% structures
if(cond1==FALSE){stop(cat("struct must be one of",structures))}
if(h<1){stop("Forecast Horizon must be at least 1")}
if(cv!="Rolling" & cv!="LOO"){stop("Cross-Validation type must be one of Rolling or LOO")}
if(length(gran)!=2&ownlambdas==FALSE){stop("Granularity must have two parameters")}
structure2 <- c("None","Lag","HVARC")
cond2=struct %in% structure2
k1=0
if(length(VARX)!=0){k1 <- VARX$k}else{k=0}
nseries <- ncol(Y)-k1
if(p==0){tf=TRUE
     }else{
         tf=FALSE
         }
## if(!is.null(VARX$tf)){tf=TRUE}else{
## tf=FALSE
##     }
if(nseries==1 & cond2==FALSE ){stop("Univariate support is only available for Lasso, Lag Group, and Componentwise HVAR")}
if(length(VARX)==0 & struct=="EFX"){stop("EFX is only supported in the VARX framework")}
if(struct=="EFX" & !is.null(VARX$contemp)){stop("EFX does not support contemporaneuous dependence")}
structs=c("HVARC","HVAROO","HVARELEM")
if(length(VARX)!=0& struct %in% structs){stop("EFX is the only nested model supported in the VARX framework")}

(BV1 <- new(
      "BigVAR",
      Data=Y,
      lagmax=p,
      Structure=struct,
      Relaxed=RVAR,
      Granularity=gran,
      Minnesota=MN,
      verbose=verbose,
      horizon=h,
      crossval=cv,
      ic=IC,
      VARX=VARX,
      T1=T1,
      T2=T2,
      ONESE=ONESE,
      ownlambdas=ownlambdas,
      tf=tf    
       ))

  return(BV1)

  }
# show-default method to show an object when its name is printed in the console.
#' Default show method for an object of class BigVAR
#'
#' @param object \code{BigVAR} object created from \code{ConstructModel}
#' @return Displays the following information about the BigVAR object:
#' \itemize{
#' \item{Prints the first 10 rows of \code{Y}}
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
setMethod("show","BigVAR",
          function(object)
          {
            nrowShow <- min(10,nrow(object@Data))
            cat("*** BIGVAR MODEL *** \n")
            cat("Data (First 10 Observations):\n")
            print(formatC(object@Data[1:nrowShow,]),quote=FALSE)
            cat("Structure\n") ;print(object@Structure)
            ## cat("Forecast Horizon \n") ;print(object@horizon)
            cat("Relaxed VAR \n") ;print(object@Relaxed)
            cat("Minnesota Prior \n") ;print(object@Minnesota)
            cat("Maximum Lag Order \n") ;print(object@lagmax)
            if(length(object@VARX)!=0){
            cat("VARX Specs \n") ;print(object@VARX)}
            cat("Start of Cross Validation Period \n") ;print(object@T1-object@lagmax)
            cat("End of Cross Validation Period \n") ;print(object@T2-object@lagmax)
            
            }
)

#' Plot a BigVAR object
#'
#' @param x BigVAR object created from \code{ConstructModel}
#' @param y NULL
#' @param ... additional plot arguments
#' @return NA, side effect is graph
#' @details Uses plot.zoo to plot each individual series of \code{Y} on a single plot
#' @name plot.BigVAR
#' @import methods
#' @seealso \code{\link{constructModel}}
#' @aliases plot,BigVAR-method
#' @aliases plot-methods
#' @docType methods
#' @method plot
#' @rdname plot.BigVAR-methods
#' @export
#' @importFrom zoo plot.zoo
#' @importFrom zoo as.zoo
setMethod(f="plot",signature="BigVAR",
     definition= function(x,y=NULL,...)
          {
              g=ncol(x@Data)
              plot.zoo(as.zoo(x@Data),plot.type="single",col=1:g)

            }
)

#' Cross Validation for BigVAR
#' Performs rolling cross-validation on a BigVAR object
#' @usage cv.BigVAR(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details Will perform cross validation to select penalty parameters over a training sample, then evaluate them over a test set.  Compares against sample mean, random walk, AIC, and BIC benchmarks.  The resulting object is of class \code{BigVAR.results}
#' @return An object of class \code{BigVAR.results}.
#' @seealso \code{\link{constructModel}}, \code{\link{BigVAR.results}} 
#' @name cv.BigVAR
#' @aliases cv.BigVAR,BigVAR-method
#' @docType methods
#' @rdname cv.BigVAR-methods
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' # construct a Lasso VAR 
#' Model1=constructModel(Y,p=4,struct="None",gran=c(50,10))
#' results=cv.BigVAR(Model1)
#' @export
setGeneric(

  name="cv.BigVAR",
  def=function(object)
  {
    standardGeneric("cv.BigVAR")

    }
   

  )
# Cross-validation and evaluation function
setMethod(

  f="cv.BigVAR",
  signature="BigVAR",
  definition=function(object){
     p=object@lagmax
	s1=0
    Y <- object@Data
     k <- ncol(Y)
     alpha <- 1/(k+1)
    RVAR=object@Relaxed
    group=object@Structure
    cvtype=object@crossval
    MN <- object@Minnesota
     jj=0
     if(class(Y)!="matrix"){Y=matrix(Y,ncol=1)}
     if(object@crossval=="Rolling"){
     T1=object@T1
     }else{
     T1=p+2    
         }
     T2=object@T2
     s <- ifelse(length(object@VARX)!=0,object@VARX$s,0)
     ## T <- nrow(Y)-max(p,s)
if(object@ownlambdas==TRUE){gamm=object@Granularity
gran2 <- length(gamm)
}
   ONESE=object@ONESE
     VARX=object@VARX
     if(object@ownlambdas==FALSE){
  gran2 = object@Granularity[2]
    gran1=object@Granularity[1]
  }


     if(length(VARX)!=0){
     VARX=TRUE
     k1=object@VARX$k
     s=object@VARX$s
     if(!is.null(object@VARX$contemp)){
         contemp=TRUE
         s1=1
         }else{
             contemp=FALSE
             s1=0}
     m=k-k1
 ##  if(k1>1){
 ##     Z1 <- Zmat1(Y[,1:k1],p,k1)
 ##     Z1 <- Z1[2:nrow(Z1),]}else{
 ##         Z1<-matrix(t(.ZScalar(Y[,k1],p)$Z),nrow=p)}
 ##     if(m>1){
 ##     Z2 <- Zmat1(Y[,(k1+1):k],s,k-k1)
 ##     Z2 <- Z2[2:nrow(Z2),]
 ## }else{Z2<-matrix(t(.ZScalar(Y[,m],s)$Z),nrow=p)}

 ##     trainY <- matrix(Y[(max(c(p,s))+1):nrow(Y),1:k1],ncol=k1)
 ##     if(p<s){Z1 <- Z1[,(s-p+1):ncol(Z1)]
 ##         }
 ##     if(p>s){Z2 <- Z2[,(p-s+1):ncol(Z2)]}
 ##     trainZ <- rbind(Z1,Z2)

 ##     if(k1>1){
 ##     Z3 <- Zmat1(Y[,1:k1],p,k1)
 ##     } else{Z3=t(.ZScalar(Y[,k1],p)$Z)}
 ##     if(m>1){
 ##     Z4 <- Zmat1(Y[,(k1+1):k],s,k-k1)
 ##     Z4 <- Z4[2:nrow(Z4),]
 ## }else{Z4<-matrix(t(.ZScalar(Y[,m],s)$Z),nrow=s)}
 ##     if(p<s){Z3 <- Z3[,(s-p+1):ncol(Z3)]
 ##         }
 ##     if(p>s){Z4 <- Z4[,(p-s+1):ncol(Z4)]}
 ##     ZF2 <- rbind(Z3,Z4)

     Y1 <-matrix(Y[,1:k1],ncol=k1)
     ## if(VARX==TRUE){

     X <- matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m)
     ## }else{X=matrix(0,nrow=nrow(Y))}
     if(!object@tf){
     trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)
     }else{
       ## browser()
         ## print("TEST")
         trainZ <- VARXCons(matrix(0,ncol=1,nrow=nrow(X)),matrix(X,ncol=m),k=0,p=0,m=m,s=s,contemp=contemp,oos=FALSE)
         }
     ## browser()
     trainZ <- trainZ[2:nrow(trainZ),]
     trainY <- matrix(Y[(max(c(p,s))+1):nrow(Y),1:k1],ncol=k1)

    if(group=="Lag"|group=="SparseLag"){
        alpha=1/(k1+1)
        jj=groupfunVARXLG(p,k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)

        }
        if(group=="Diag"|group=="SparseDiag"){
            jj=diaggroupfunVARXLG(p,k,k1,s+s1)
            activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)

            }
     if(object@ownlambdas==FALSE){
     gamm <- .LambdaGridX(gran1, gran2, jj, as.matrix(trainY[1:T2,]), trainZ[,1:T2],group,p,k1,s+s1,m,k,MN)
     }

     
     beta=array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2))


        if (group == "Lag") {
          jj=groupfunVARX(p,k,k1,s+s1)
          jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)
    }
    if (group == "SparseLag") {
        jj <- groupfunVARX(p, k,k1,s+s1)
        q1a <- list()
        alpha=1/(k1+1)
        for (i in 1:(p+s+s1)) {
            q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
        }
             activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)
   
    }
    if (group == "Diag") {
        kk <- diaggroupfunVARX(p,k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
    }
         if (group == "SparseDiag") {
             alpha=1/(k1+1)
        kk <- diaggroupfunVARX(p,k,k1,s+s1)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
          q1a <- list()

        for (i in 1:(2*p)) {
            q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
        }
        
      
     }
     }else{
         ## s=p
       if(group=="Lag"|group=="SparseLag")
      {
      jj=.groupfun(p,k)
     }else{
      jj <- .lfunction3(p,k)
      }

        Z1=VARXCons(Y,matrix(0,nrow=nrow(Y)),k,p,0,0) 
           ## Z1=Zmat1(Y,p,k)
        trainZ=Z1[2:nrow(Z1),]   
        trainY=matrix(Y[(p+1):nrow(Y),],ncol=k)          
      GY=matrix(trainY[1:T2,],ncol=k)
      GZ=trainZ[,1:T2]
      if(object@ownlambdas==FALSE){   
      gamm <- .LambdaGridE(gran1, gran2, jj, GY, GZ,group,p,k,MN)
      }
        VARX=FALSE
        k1=k
        s=0   
     if (group == "Lag") {
        jj <- .groupfuncpp(p, k)
        jjcomp <- .groupfuncomp(p,k)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)
        k1=k
    }
    if (group == "SparseLag") {
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
        jjcomp=.lfunctioncomp(p,k)
        jj=.lfunction3(p,k)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
        q1a <- list()
        for (i in 1:(2*p)) {
            q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)
        }

        
     }
    beta=array(0,dim=c(k,k*p+1,gran2))
    }           
    h <- object@horizon
    verbose <- object@verbose         
    activeL <- rep(list(list()), gran2)
    activeset = activeL
     ZFull <- list()
     if(!is.matrix(trainZ)){trainZ <- matrix(trainZ,ncol=1)}
     if(!is.matrix(trainY)){trainY <- matrix(trainY,ncol=1)}
     
     ZFull$Z=trainZ
     ZFull$Y <- trainY
     T <- nrow(trainY)
## browser()
  if(object@ownlambdas==TRUE){gamm <- object@Granularity}    
     if(group=="Tapered")
         {
palpha=seq(0,1,length=10)
palpha <- rev(palpha)

gran2=length(gamm)*length(palpha)
beta=array(0,dim=c(k,k*p+1,gran2))
}

     if(class(ZFull$Y)!="matrix" ){ZFull$Y=matrix(ZFull$Y,ncol=1)}
  MSFE <- matrix(0, nrow = T2 - T1+1, ncol = gran2)
  if(verbose==TRUE){
    pb <- txtProgressBar(min = T1, max = T2, style = 3)
    cat("Cross-Validation Stage:",group)}
    YT <- Y[1:T2,]
    for (v in T1:T2) {
      if(cvtype=="Rolling")
        {
        trainY <- ZFull$Y[1:(v-1), ]
        trainZ <- ZFull$Z[, 1:(v-1)]
        }else{
            if(VARX==TRUE)
                {
                   YT2 <- YT[-v,]
                    Y1 <-matrix(YT2[,1:k1],ncol=k1)
                    X <- matrix(YT2[,(ncol(YT2)-m+1):ncol(YT2)],ncol=m)
                    trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)
                    trainZ <- trainZ[2:nrow(trainZ),]
                    trainY <- matrix(YT2[(max(c(p,s))+1):nrow(YT2),1:k1],ncol=k1)
                    }else{

                        YT2 <- YT[-v,]
                        Z1=VARXCons(YT2,matrix(0,nrow=nrow(YT2)),k,p,0,0) 
                        trainZ=Z1[2:nrow(Z1),]        
                        trainY=matrix(YT2[(p+1):nrow(YT2),],ncol=k)                                  
                        }
            }

        if (group == "None") {
            if(VARX==TRUE){
           beta <- .lassoVARFistX(beta, trainZ, trainY,gamm, 1e-04,p,MN,k,k1,s+s1,m)
       }else{beta <- .lassoVARFist(beta, trainZ, trainY,gamm, 1e-04,p,MN)}
           }
        if (group == "Lag") {
            # Should be the same function for both cases
            ## GG <- GroupLassoVAR(beta,trainY,trainZ,gamm,1e-04,k,p,activeset,jj,jjcomp,k1,s,MN)
             GG <- .GroupLassoVAR1(beta,jj,jjcomp,trainY,trainZ,gamm,activeset,1e-04,p,MN,k,k1,s+s1)
            beta <- GG$beta
            activeset <- GG$active
        }
        if (group == "SparseLag") {
            if(VARX==TRUE){
            GG <- .SparseGroupLassoVARX(beta, jj, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN,k,s+s1,k1)
            }else{
               GG <- .SparseGroupLassoVAR(beta, jj, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)
            }
            beta <- GG$beta
            activeset = GG$active
            q1a <- GG$q1
        }
        if (group == "Diag") {
            if(VARX==TRUE){
            GG <- .GroupLassoOOX(beta, kk, trainY, trainZ, gamm, 
                activeset, 1e-04,p,MN,k,k1,s+s1)
            }else{
            GG <- .GroupLassoOO(beta, kk, trainY, trainZ, gamm, 
                activeset, 1e-04,p,MN)
                }
            beta <- GG$beta
            activeset <- GG$active
        }
          if (group == "SparseDiag") {
              if(VARX==TRUE){
            GG <- .SparseGroupLassoVAROOX(beta, kk, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04,p,MN,k1,s+s1,k)
                  }else{
            GG <- .SparseGroupLassoVAROO(beta, kk, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04,q1a,p,MN)
            q1a <- GG$q1
                      }
            beta <- GG$beta
            activeset = GG$active
        }
      if(group=="Tapered")
          {
              beta <- .lassoVARTL(beta,trainZ,trainY,gamm,1e-4,p,MN,palpha)
              
              }
      if(group=="EFX")
          {
              beta <- .EFVARX(beta,trainY,trainZ,gamm,1e-5,MN,k1,s,m,p)
              }
      if(group=="HVARC")
          {
              beta <- .HVARCAlg(beta,trainY,trainZ,gamm,1e-5,p,MN)

              }
      if(group=="HVAROO")
          {
              beta <- .HVAROOAlg(beta,trainY,trainZ,gamm,1e-5,p,MN)


              }
      if(group=="HVARELEM")
          {
              beta <- .HVARElemAlg(beta,trainY,trainZ,gamm,1e-5,p,MN)


              }

      
      eZ <- c(1,ZFull$Z[,v])


        ## eZh <- as.matrix(.Zmat2(Y[(v - p):v, ], p, k)$Z, ncol = 1)
        for (ii in 1:gran2) {
            if (RVAR == TRUE) {
                # Relaxed Least Squares (intercept ignored)
                beta[,,ii] <- RelaxedLS(cbind(t(trainZ),trainY),beta[,,ii],k,p,k1,s+s1)
            }
                        if(MN==TRUE){
            MSFE[v - (T1 - 1), ii] <- norm2(ZFull$Y[v,1:k1] - beta[,2:dim(beta)[2],ii] %*% eZ[2:length(eZ)])^2

                }else{
                    if(object@crossval=="Rolling"){
            MSFE[v - (T1 - 1), ii] <- norm2(ZFull$Y[v,1:k1] - beta[,,ii] %*% eZ)^2
                    }else{

                        ## browser()
                        ## sp <- max(p,s)
              if(VARX==TRUE){          
              eZ<- VARXCons(matrix(Y[(v-p):(v),1:k1],ncol=k1),Y[(v-p):(v),(ncol(Y)-k1):(ncol(Y))],k1,p,m,s,contemp=contemp)
              }else{
              eZ<- VARXCons(Y[(v-p):(v),1:k1],matrix(0,nrow=length((v-p):v)),k1,p,0,0)
                  }
        ## eZ <- as.matrix(.Zmat2(Y[(v - p-1):(v-1), ], p, k)$Z, ncol = 1)
## head(eZ)
##                    head(Y[v-1,])
              ## eZ <- c(1,ZFull$Z[,v-1])
             MSFE[v - (T1 - 1), ii] <- norm2(Y[v,1:k1] - beta[,,ii] %*% eZ)^2
      


                    }
                    
            }
            }
        if(verbose==TRUE){
         setTxtProgressBar(pb, v)}

      }

if(group=="Tapered")
    {
indopt <- which.min(colMeans(MSFE))

if(indopt<length(gamm))
    {
alphaind <- 1
alphaopt <- palpha[1]
        }
else{
alphaopt <- palpha[floor(indopt/length(gamm))]
alphaind <- floor(indopt/length(gamm))

}
if(alphaind==1)
{
    gamopt <- gamm[indopt]
}else if(indopt %% length(gamm)==0)
         {
             
           gamopt <- gamm[length(gamm)]

             }else{
    gamind <- indopt-length(gamm)*alphaind
    gamopt <- gamm[gamind]
    }

palpha<-alphaopt


}


     
     if(ONESE==TRUE){
     MSFE2 <-MSFE 
     G2 <- colMeans(na.omit(MSFE2))
     G3 <- sd(na.omit(MSFE2))/sqrt(nrow(na.omit(MSFE2)))
     gamopt <- gamm[min(which(G2<(min(G2)+G3)))]
     G2diag<-G2
     gamopt2<-gamopt
     G32<-G3
     }
     else{
if(group!="Tapered")
  {
            gamopt <- gamm[max(which(colMeans(na.omit(MSFE))==min(colMeans(na.omit(MSFE)))))]}

     }
     if(VARX==TRUE){
     OOSEval <- .EvalLVARX(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m,contemp)

     }else{
    OOSEval <- .EvalLVAR(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T)
     }
     MSFEOOSAgg <- na.omit(OOSEval$MSFE)
     betaPred <- OOSEval$betaPred
     preds <- OOSEval$predictions
Y <- object@Data

if(VARX==TRUE){

if(contemp){OOS=FALSE}else{OOS=TRUE}
## browser()
if(!object@tf){
    Zvals <- VARXCons(matrix(Y[,1:k1],ncol=k1),matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m),k1,p,m,s,oos=OOS,contemp=contemp)}else{
    Zvals <- VARXCons(matrix(0,ncol=1,nrow=nrow(Y)),matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m),0,0,m,s,oos=FALSE,contemp=contemp)
        }
}else{
m=0;s=0
Zvals <- VARXCons(matrix(Y[,1:k1],ncol=k1),matrix(0,nrow=nrow(Y)),k1,p,m,s,oos=TRUE)}

Zvals <- matrix(Zvals[,ncol(Zvals)],ncol=1)
## tail(Y)
     ## Zvals <- OOSEval$zvals
     ## resids <- OOSEval$resids
if(ncol(Y)==1| k1==1){betaPred <- matrix(betaPred,nrow=1)}

resids <- t(t(ZFull$Y)-betaPred%*%rbind(rep(1,ncol(ZFull$Z)),ZFull$Z))
                                                               
     MSFEOOS<-mean(na.omit(MSFEOOSAgg))
     seoos=sd(na.omit(MSFEOOSAgg))/sqrt(length(na.omit(MSFEOOSAgg)))
     ## if(class(Y)!="matrix"){Y=matrix(Y,ncol=1)}
if(VARX==FALSE){k1=k}
      meanbench=.evalMean(ZFull$Y[,1:k1],T2,T)
	RWbench=.evalRW(ZFull$Y[,1:k1],T2,T)
 
     if(object@ic==FALSE|object@tf){AICbench=list()
       AICbench$Mean=0
       AICbench$SD=0
       BICbench=list()
       BICbench$Mean=0
       BICbench$SD=0                          
                      }else{
    AICbench1 <- EvalIC(matrix(ZFull$Y[,1:k1],ncol=k1),T2,k1,p,"AIC")
    AICbench <- list()
    AICbench$Mean <- mean(AICbench1$MS)
    AICbench$SD <- sd(AICbench1$MS)/sqrt(length(AICbench1$MS))
    BICbench1 <- EvalIC(matrix(ZFull$Y[,1:k1],ncol=k1),T2,k1,p,"BIC")
    BICbench <- list()
    BICbench$Mean <- mean(BICbench1$MS)
    BICbench$SD <- sd(BICbench1$MS)/sqrt(length(BICbench1$MS))
 
}

    results <- new("BigVAR.results",InSampMSFE=colMeans(MSFE),InSampSD=apply(MSFE,2,sd)/sqrt(nrow(MSFE)),LambdaGrid=gamm,index=which.min(colMeans(na.omit(MSFE))),OptimalLambda=gamopt,OOSMSFE=MSFEOOSAgg,seoosmsfe=seoos,MeanMSFE=meanbench$Mean,AICMSFE=AICbench$Mean,RWMSFE=RWbench$Mean,MeanSD=meanbench$SD,AICSD=AICbench$SD,BICMSFE=BICbench$Mean,BICSD=BICbench$SD,RWSD=RWbench$SD,Data=object@Data,lagmax=object@lagmax,Structure=object@Structure,Minnesota=object@Minnesota,Relaxed=object@Relaxed,Granularity=object@Granularity,horizon=object@horizon,betaPred=betaPred,Zvals=Zvals,resids=resids,VARXI=VARX,VARX=list(k=k1,s=s),preds=preds)
     
    return(results)
     }
 )


#' BigVAR Estimation
#' @description
#' Fit a BigVAR object with a structured penalty
#' @usage BigVAR.est(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details Fits a using a structured penalty on a BigVAR object.  Does not perform cross-validation.
#' @return An array of \eqn{k \times kp} or \eqn{k\times kp+ms} coefficient matrices; one for each value of lambda.
#' @seealso \code{\link{constructModel}}, \code{\link{BigVAR.results}} 
#' @name BigVAR.est
#' @aliases BigVAR.est,BigVAR-method
#' @docType methods
#' @rdname BigVAR.est-methods
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' # construct a Lasso VAR 
#' Model1=constructModel(Y,p=4,struct="None",gran=c(50,10))
#' BigVAR.est(Model1)
#' @export
setGeneric(

  name="BigVAR.est",
  def=function(object)
  {
    standardGeneric("BigVAR.est")

    }
   

  )
setMethod(

  f="BigVAR.est",
  signature="BigVAR",
  definition=function(object){
     p=object@lagmax
if(object@ownlambdas==TRUE){gamm=object@Granularity
gran2 <- length(gamm)
}      
    Y <- object@Data
     k <- ncol(Y)
     alpha <- 1/(k+1)
    RVAR=object@Relaxed
    group=object@Structure
    MN <- object@Minnesota
     jj=0
     if(class(Y)!="matrix"){Y=matrix(Y,ncol=1)}
     s <- ifelse(length(object@VARX)!=0,object@VARX$s,0)
     T <- nrow(Y)-max(p,s)
     VARX=object@VARX
if(object@ownlambdas==FALSE){
     gran2 = object@Granularity[2]
    gran1=object@Granularity[1]}
     if(length(VARX)!=0){
     VARX=TRUE
     k1=object@VARX$k
     s=object@VARX$s
     m=k-k1
  if(k1>1){
     Z1 <- Zmat1(Y[,1:k1],p,k1)
     Z1 <- Z1[2:nrow(Z1),]}else{
         Z1<-matrix(t(.ZScalar(Y[,k1],p)$Z),nrow=p)}

     Z2 <- Zmat1(Y[,(k1+1):k],s,k-k1)
     Z2 <- Z2[2:nrow(Z2),]

     trainY <- matrix(Y[(max(c(p,s))+1):nrow(Y),1:k1],ncol=k1)
     if(p<s){Z1 <- Z1[,(s-p+1):ncol(Z1)]
         }
     if(p>s){Z2 <- Z2[,(p-s+1):ncol(Z2)]}
     trainZ <- rbind(Z1,Z2)

     if(k1>1){
     Z3 <- Zmat1(Y[,1:k1],p,k1)
     }
     else{Z3=t(.ZScalar(Y[,k1],p)$Z)}
     Z4 <- Zmat1(Y[,(k1+1):k],s,k-k1)
     Z4 <- Z4[2:nrow(Z4),]
     if(p<s){Z3 <- Z3[,(s-p+1):ncol(Z3)]
         }
     if(p>s){Z4 <- Z4[,(p-s+1):ncol(Z4)]}
     ZF2 <- rbind(Z3,Z4)
    if(group=="Lag"|group=="SparseLag"){
        alpha=1/(k1+1)
        jj=groupfunVARXLG(p,k,k1,s)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)

        }
        if(group=="Diag"|group=="SparseDiag"){
            jj=diaggroupfunVARXLG(p,k,k1,s)
            activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)

            }
if(object@ownlambdas==FALSE){
     gamm <- .LambdaGridX(gran1, gran2, jj, trainY, trainZ,group,p,k1,s,m,k)
}

     beta=array(0,dim=c(k1,k1*p+(k-k1)*s+1,gran2))

        if (group == "Lag") {
          jj=groupfunVARX(p,k,k1,s)
          jjcomp <- groupfunVARXcomp(p,k,k1,s)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)
    }
    if (group == "SparseLag") {
        jj <- groupfunVARX(p, k,k1,s)
        q1a <- list()
        alpha=1/(k1+1)
        for (i in 1:(p+s)) {
            q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
        }
             activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)
   
    }
    if (group == "Diag") {
        kk <- diaggroupfunVARX(p,k,k1,s)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
    }
         if (group == "SparseDiag") {
             alpha=1/(k1+1)
        kk <- diaggroupfunVARX(p,k,k1,s)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
          q1a <- list()

        for (i in 1:(2*p)) {
            q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
        }
        
      
     }
     }else{
         s=p
       if(group=="Lag"|group=="SparseLag")
      {
      jj=.groupfun(p,k)
     }else{
      jj <- .lfunction3(p,k)
      }
           Z1=Zmat1(Y,p,k)
        trainZ=Z1[2:nrow(Z1),]   
        trainY=Y[(p+1):nrow(Y),]          
if(object@ownlambdas==FALSE){
      gamm <- .LambdaGridE(gran1, gran2, jj, trainY, trainZ,group,p,k)
}
        VARX=FALSE
        k1=k
        s=0   
     if (group == "Lag") {
        jj <- .groupfuncpp(p, k)
        jjcomp <- .groupfuncomp(p,k)
        activeset <- rep(list(rep(rep(list(0), length(jj)))), 
            gran2)
        k1=k
    }
    if (group == "SparseLag") {
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
        jjcomp=.lfunctioncomp(p,k)
        jj=.lfunction3(p,k)
        activeset <- rep(list(rep(rep(list(0), length(kk)))), 
            gran2)
        q1a <- list()
        for (i in 1:(2*p)) {
            q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)
        }

        
     }
    beta=array(0,dim=c(k,k*p+1,gran2))
    }           
     if(group=="Tapered")
         {
palpha=seq(0,1,length=10)
palpha <- rev(palpha)
gran2=length(gamm)*length(palpha)
beta=array(0,dim=c(k,k*p+1,gran2))
}

     ## if(class(ZFull$Y)!="matrix" ){ZFull$Y=matrix(ZFull$Y,ncol=1)}
        if (group == "None") {
            if(VARX==TRUE){
           beta <- .lassoVARFistX(beta, trainZ, trainY,gamm, 1e-04,p,MN,k,k1,s,m)}
            else{beta <- .lassoVARFist(beta, trainZ, trainY,gamm, 1e-04,p,MN)}
           }
        if (group == "Lag") {
            GG <- .GroupLassoVAR1(beta,jj,jjcomp,trainY,trainZ,gamm,activeset,1e-04,p,MN,k,k1,s)
            beta <- GG$beta
            activeset <- GG$active
        }
        if (group == "SparseLag") {
            if(VARX==TRUE){
            GG <- .SparseGroupLassoVARX(beta, jj, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN,k,s,k1)
            }else{
               GG <- .SparseGroupLassoVAR(beta, jj, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)
            }
            beta <- GG$beta
            activeset = GG$active
            q1a <- GG$q1
        }
        if (group == "Diag") {
            if(VARX==TRUE){
            GG <- .GroupLassoOOX(beta, kk, trainY, trainZ, gamm, 
                activeset, 1e-04,p,MN,k,k1,s)
            }else{
            GG <- .GroupLassoOO(beta, kk, trainY, trainZ, gamm, 
                activeset, 1e-04,p,MN)
                }
            beta <- GG$beta
            activeset <- GG$active
        }
          if (group == "SparseDiag") {
              if(VARX==TRUE){
            GG <- .SparseGroupLassoVAROOX(beta, kk, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04,p,MN,k1,s,k)
                  }else{
            GG <- .SparseGroupLassoVAROO(beta, kk, trainY, trainZ, 
                gamm, alpha, INIactive = activeset, 1e-04,q1a,p,MN)
            q1a <- GG$q1
                      }
            beta <- GG$beta
            activeset = GG$active
        }
      if(group=="Tapered")
          {
              beta <- .lassoVARTL(beta,trainZ,trainY,gamm,1e-4,p,MN,palpha)
              
              }
      if(group=="EFX")
          {
              beta <- .EFVARX(beta,trainY,trainZ,gamm,1e-5,MN,k1,s,m,p)
              }
      if(group=="HVARC")
          {
              beta <- .HVARCAlg(beta,trainY,trainZ,gamm,1e-5,p,MN)

              }
      if(group=="HVAROO")
          {
              beta <- .HVAROOAlg(beta,trainY,trainZ,gamm,1e-5,p,MN)


              }
      if(group=="HVARELEM")
          {
              beta <- .HVARElemAlg(beta,trainY,trainZ,gamm,1e-5,p,MN)


              }      
     
    return(list(B=beta,lambdas=gamm))
     }
 )


## New object class: bigVAR results, inherits class bigVAR, prints results from cv.bigVAR

# Prints the results, but stores a lot more information if needed
# compares everything against AIC, Mean, Random Walk

#' BigVAR.results
#'
#' This class contains the results from cv.BigVAR.
#'
#' It inherits the class BigVAR, but contains substantially more information. 
#' 
#'    @field InSampMSFE In-sample MSFE from optimal value of lambda
#'    @field LambdaGrid Grid of candidate lambda values
#'    @field index Rank of optimal lambda value 
#'    @field OptimalLambda Value of lambda which minimizes MSFE
#'    @field OOSMSFE Average Out of sample MSFE of BigVAR model with Optimal Lambda
#'   @field seoosfmsfe Standard Error of Out of sample MSFE of BigVAR model with Optimal Lambda
#'   @field MeanMSFE Average Out of sample MSFE of Unconditional Mean Forecast
#'   @field MeanSD Standard Error of out of sample MSFE of Unconditional Mean Forecast
#' @field RWMSFE Average Out of sample MSFE of Random Walk Forecast
#' @field RWSD Standard Error of out of sample MSFE of Random Walk Forecast
#' @field AICMSFE Average Out of sample MSFE of AIC Forecast
#' @field AICSD Standard Error of out of sample MSFE of AIC Forecast
#' @field BICMSFE Average Out of sample MSFE of BIC Forecast
#' @field BICSD Standard Error of out of sample MSFE of BIC Forecast
#' @field betaPred The final out of sample coefficient matrix of \code{B}, to be used for prediction
#' @field Zvals The final lagged values of \code{Y}, to be used for prediction
#' @field resids residuals obtained from betaPred
#' @field Data a \eqn{T x k} multivariate time Series
#' @field lagmax Maximal lag order
#' @field Structure Penalty Structure
#' @field Relaxed Indicator for relaxed VAR
#' @field Granularity Granularity of Penalty Grid
#' @field horizon Desired Forecast Horizon
#' @field crossval Cross-Validation Procedure
#' @field alpha penalty for Sparse Group Lasso
#' @field VARXI VARX Indicator 
#' @field Minnesota Minnesota Prior Indicator
#' @field verbose  verbose indicator
#' 
#' @note One can also access any object of class BigVAR from BigVAR.results
#' @name BigVAR.results 
#' @rdname BigVAR.results
#' @aliases BigVAR.results-class
#' @exportClass BigVAR.results
#' @author Will Nicholson
#' @export
setClass("BigVAR.results",
         representation(InSampMSFE="numeric",InSampSD="numeric",LambdaGrid="numeric",index="numeric",OptimalLambda="numeric",OOSMSFE="numeric",seoosmsfe="numeric",MeanMSFE="numeric",AICMSFE="numeric",BICMSFE="numeric",BICSD="numeric",RWMSFE="numeric",MeanSD="numeric",AICSD="numeric",RWSD="numeric",betaPred="matrix",Zvals="matrix",VARXI="logical",resids="matrix",preds="matrix"),
         contains="BigVAR"
         )


#' Plot an object of class BigVAR.results
#' 
#' @param x BigVAR.results object created from \code{cv.BigVAR}
#' @param y NULL
#' @param ... additional arguments
#' @details Plots the in sample MSFE of all values of lambda
#' @name plot
#' @import methods
#' @aliases plot,BigVAR.results-method
#' @aliases plot-methods
#' @docType methods
#' @method plot
#' @rdname BigVAR.results-plot-methods
#' @export
setMethod(f="plot",signature="BigVAR.results",
     definition= function(x,y=NULL,...)
          {
            plot(x@LambdaGrid,x@InSampMSFE,type="l",xlab="Value of Lambda",ylab="MSFE",log="x")
            abline(v=x@OptimalLambda,col="green")
            }
)

#' Default show method for an object of class BigVAR.results
#'
#' @param object BigVAR.results object created from \code{cv.BigVAR}
#' @details prints forecast results and additional diagonostic information as well as comparisons with mean, random walk, and AIC, and BIC benchmarks
#' @seealso \code{\link{cv.BigVAR}} 
#' @name show
#' @aliases show,BigVAR.results-method
#' @docType methods
#' @method show
#' @rdname show-methods-BigVAR.results
#' @export
setMethod("show","BigVAR.results",
          function(object)
          {
            nrowShow <- min(10,nrow(object@Data))
            ## ncolShow <- min(10,ncol(object@data))

            cat("*** BIGVAR MODEL Results *** \n")
            cat("Structure\n") ;print(object@Structure)
            ## cat("Forecast Horizon \n") ;print(object@horizon)
            if(object@Relaxed==TRUE){
            cat("Relaxed VAR \n") ;print(object@Relaxed)}
            if(object@VARXI){
            cat("VARX Specs \n") ;print(object@VARX)}
            cat("Maximum Lag Order \n") ;print(object@lagmax)
            cat("Optimal Lambda \n"); print(round(object@OptimalLambda,digits=4))
            cat("Grid Depth \n") ;print(object@Granularity[1])
            cat("Index of Optimal Lambda \n");print(object@index)
 
            cat("In-Sample MSFE\n");print(round(min(object@InSampMSFE),digits=3))
            cat("BigVAR Out of Sample MSFE\n");print(round(mean(object@OOSMSFE),digits=3))
            cat("*** Benchmark Results *** \n")
	    ## cat("BigVAR Out of Sample Standard Error\n");print(min(object@seoosmsfe))
            cat("Conditional Mean Out of Sample MSFE\n");print(round(object@MeanMSFE,digits=3))
	    ## cat("Unconditional Mean Standard Error\n");print(min(object@MeanSD))
            cat("AIC Out of Sample MSFE\n");print(round(object@AICMSFE,digits=3))
            ## cat("AIC Out of Sample Standard Error\n");print(min(object@AICSD))
            cat("BIC Out of Sample MSFE\n");print(round(object@BICMSFE,digits=3))
            ## cat("BIC Out of Sample Standard Error\n");print(min(object@BICSD))
	    cat("RW Out of Sample MSFE\n");print(round(object@RWMSFE,digits=3))
	    ## cat("RW Out of Sample Standard Error\n");print(min(object@RWSD))
            }
)
#' Forecast using a BigVAR.results object
#' @usage predict(object,...)
#' @param object BigVAR.results object from \code{cv.BigVAR}
#' @param ... additional arguments affecting the predictions produced (e.g. \code{n.ahead})
#' @details Provides \code{n.ahead} step forecasts using the model produced by cv.BigVAR. 
#' @seealso \code{\link{cv.BigVAR}} 
#' @name predict
#' @aliases predict,BigVAR.results-method
#' @docType methods
#' @method predict
#' @rdname predict-methods-BigVAR.results
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' Model1=constructModel(Y,p=4,struct="None",gran=c(50,10),verbose=FALSE)
#' results=cv.BigVAR(Model1)
#' predict(results,n.ahead=1)
#' @export
setMethod("predict","BigVAR.results",
          function(object,n.ahead,...)
{
    if(object@Minnesota==TRUE){
eZ <- object@Zvals
}else{
eZ <- object@Zvals[2:nrow(object@Zvals),]

    }

## object@VARX
betaPred <- object@betaPred
## str(object)    
k <-object@VARX$k
m <- ncol(object@Data)-k
p <- object@lagmax
s <- object@VARX$s
contemp=object@VARX$contemp
s1 <- 0
    if(!is.null(contemp)){
if(contemp){
    s1=1}
}else{s1=0}

nu=betaPred[,1]
betaPred <- betaPred[,2:ncol(betaPred)]
if(n.ahead==1)
    {
fcst <- betaPred%*%eZ+nu
        }
else{
if(object@VARXI){stop("multi step ahead predictions not supported for VARX models")}    
B <- VarptoVAR1MC(betaPred,p,k,m,s+s1)
    
fcst <- ((B%^%n.ahead)%*%eZ)[1:k,]+nu
    }
return(fcst)

    }
)

#' Sparsity Plot of a BigVAR.results object 
#'
#' @param object BigVAR.results object
#' @return NA, side effect is graph
#' @details Uses \code{levelplot} from the \code{lattice} package to plot the magnitude of each coefficient
#' @name SparsityPlot.BigVAR.results
#' @aliases SparsityPlot.BigVAR.results,BigVAR.results-method
#' @docType methods
#' @rdname SparsityPlot.BigVAR.results-methods
#' @export
#' @importFrom lattice levelplot
#' @importFrom lattice panel.abline
#' @importFrom lattice panel.levelplot

setGeneric(

  name="SparsityPlot.BigVAR.results",
  def=function(object)
  {
    standardGeneric("SparsityPlot.BigVAR.results")

    }
   )
setMethod(
f="SparsityPlot.BigVAR.results",
  signature="BigVAR.results",
  definition=function(object){

      B <- object@betaPred
      if(nrow(B)==1){
      B <- matrix(B[,2:ncol(B)],nrow=1)
         }else{
      B <- B[,2:ncol(B)]}

      k <- nrow(B)
      p <- object@lagmax
      s1=0
      if(length(object@VARX!=0)){
          s=object@VARX$s
          m=ncol(object@Data)-object@VARX$k
          contemp=object@VARX$contempt
        if(!is.null(contemp)){
	if(contemp){
	s1=1}
	}else{s1=0}

      }
      else{m=0;s=0}
	s=s+s1
    text <- c()
    for (i in 1:p) {
        text1 <- as.expression(bquote(bold(Phi)^(.(i))))
        text <- append(text, text1)
    }
    ## text <- c()
    if(m>0){
     for (i in (p+1):(p+s+1)) {
        text1 <- as.expression(bquote(bold(beta)^(.(i-p-s1))))
        text <- append(text, text1)
    }
     }
    f <- function(m) t(m)[, nrow(m):1]

    rgb.palette <- colorRampPalette(c("white", "blue" ),space = "Lab")
    ## rgb.palette <- colorRampPalette(c("white", "blue"), space = "Lab")
    at <- seq(k/2 + 0.5, p * (k)+ 0.5, by = k)
    if(m>0){
    at2 <- seq(p*k+m/2+.5,p*k+s*m+.5,by=m)}else{at2=c()}
    at <- c(at,at2)
    se2 = seq(1.75, by = k, length = k)

    L2 <- levelplot(as.matrix(f(abs(B))), col.regions = rgb.palette, colorkey = NULL, 
        xlab = NULL, ylab = NULL, main = list(label = "Sparsity Pattern Generated by BigVAR", 
            cex = 1), panel = function(...) {
            panel.levelplot(...)
            panel.abline(a = NULL, b = 1, h = seq(1.5, m*s+p* k + 
                0.5, by = 1), v = seq(1.5, by = 1, length = p * 
                k+m*s),lwd=.5)
            bl1 <- seq(k + 0.5, p * 
                k + 0.5, by = k)
            ## if(m>0){

            ## bl2 <- seq(p*k + 0.5, p * 
            ##     k + 0.5+s*m, by = 1)}else(bl2=c())
            b23 <- seq(p*k + 0.5, p * 
                k + 0.5+s*m, by = m)
    
            b1 <- c(bl1,b23)
            panel.abline(a = NULL, b = 1, v = p*k+.5, lwd = 3)
      
            panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
            ## panel.abline(a = NULL, b = 1, v = b23, lwd = 1.5)

        }, scales = list(x = list(alternating = 1, labels = text, 
            cex = 1, at = at, tck = c(0, 0)), y = list(alternating = 0, 
            tck = c(0, 0))))
    return(L2)
}
    )
      
      

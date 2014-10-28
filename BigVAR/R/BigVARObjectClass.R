
#' BigVAR Object Class
#'
#' An object class to be used with cv.BigVAR
#' 
#' @name BigVAR-class
#' @rdname BigVAR-class
#' @slot Data a \eqn{T x k} multivariate time Series
#' @slot lagmax Maximal lag order
#' @slot Structure Penalty Structure
#' @slot Relaxed Indicator for relaxed VAR
#' @slot Granularity Granularity of Penalty Grid
#' @slot horizon Desired Forecast Horizon
#' @slot crossval Cross-Validation Procedure
#' @slot alpha penalty for Sparse Group Lasso
#' @slot nseries Number of Series
#' @slot Minnesota Minnesota Prior Indicator
#' @slot verbose Indicator for Verbose output
#' @details Construct an object of class BigVAR via the function "ConstructModel"
#' @seealso \code{\link{constructModel}}
#' @exportClass BigVAR
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
    alpha="numeric",
    nseries="numeric"  
  )
  )


#' Construct an object of class BigVAR
#' @param Y T x K multivariate time series
#' @param p Predetermined maximal lag order
#' @param struct The choice of penalty structure (see details).
#' @param gran vector containing how deep to construct the penalty grid (parameter 1) and how many gridpoints to use (parameter 2)
#' @param RVAR True or False: whether to refit using the Relaxed-VAR procedure
#' @param h Desired forecast horizon
#' @param cv Cross-validation approach, either "Rolling" for rolling cross-validation or "LOO" for leave-one-out cross-validation.
#' @param MN Minnesota Prior Indicator
#' @param verbose, Verbose output while estimating
#' @details The choices for "struct" are as follows
#' \itemize{
#' \item{  "None" (Lasso Penalty)}
#' \item{  "Group" (Block Group Lasso)} 
#' \item{  "Sparse" (Block Sparse Group Lasso)} 
#' \item{  "Diag" (Own/Other Group Lasso) }
#' \item{  "SparseDiag" (Own/Other Sparse Group Lasso) }
#' }
#' @examples
#' library(BigVAR)
#' data(Y)
#` Y=Y[1:100,]
#' m1=constructModel(Y,p=4,struct="None",gran=c(50,10),
#' RVAR=FALSE,MN=FALSE,verbose=FALSE,h=1,cv="Rolling")
#' @export
constructModel <- function(Y,p,struct,gran,RVAR,h,cv,MN,verbose)
  {
if(dim(Y)[2]>dim(Y)[1]){warning("k is greater than T, is Y formatted correctly (k x T)?")}      
if(p<1){stop("Maximal lag order must be at least 1")}
structures=c("None","Group","Sparse","Diag","SparseDiag")
cond1=struct %in% structures
if(cond1==FALSE){stop(cat("struct must be one of",structures))}
if(h<1){stop("Forecast Horizon must be at least 1")}
if(cv!="Rolling" & cv!="LOO"){stop("Cross-Validation type must be one of Rolling or LOO")}
if(length(gran)!=2){stop("Granularity must have two parameters")}
structure2 <- c("None","Group")
cond2=struct %in% structure2
if(ncol(Y)==1 & cond2==FALSE ){stop("Univariate support is only available for Lasso")}


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
      alpha=1/(ncol(Y)+1),
      nseries=ncol(Y)

       ))

  return(BV1)

  }
# show-default method to show an object when its name is printed in the console.
#' Default show method for an object of class BigVAR
#'
#' @param object BigVAR object created from \code{ConstructModel}
#' @return Data prints the first 10 rows of \code{Y}
#' @return Structure prints desired structure'
#' @return Forecast Horizon
#' @return Relaxed Indicator
#' @return Maximum lag order p
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
            ## ncolShow <- min(10,ncol(object@data))

            cat("*** BIGVAR MODEL *** \n")
            cat("Data (First 10 Observations):\n")
            print(formatC(object@Data[1:nrowShow,]),quote=FALSE)
            cat("Structure\n") ;print(object@Structure)
            cat("Forecast Horizon \n") ;print(object@horizon)
            cat("Relaxed VAR \n") ;print(object@Relaxed)
            cat("Minnesota Prior \n") ;print(object@Minnesota)
            cat("Maximum Lag Order \n") ;print(object@lagmax)



            }
)

#' Plot a BigVAR object
#' 
#' @param x BigVAR object created from \code{ConstructModel}
#' @param y needed to mantain compatibility with generic, otherwise ignored
#' @param ... additional arguments
#' @return NA, side effect is graph
#' @details Uses plot.zoo to plot each individual series of \code{Y} on a single plot
#' @name plot.BigVAR
#' @import methods
#' @import graphics
#' @seealso \code{\link{constructModel}}
#' @aliases plot,BigVAR,ANY-method
#' @docType methods
#' @rdname plot-methods
#' @export
setMethod(f="plot",signature="BigVAR",
     definition= function(x,y=NULL,...)
          {
              g=ncol(x@Data)
              plot.zoo(as.zoo(x@Data),plot.type="single",col=1:g)

            }
)

#' Cross Validation for BigVAR
#' @description
#' Performs rolling or leave-one-out cross-validation on a BigVAR object
#' @usage cv.BigVAR(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details Will perform cross validation to select penalty parameters over a training sample, then evaluate them over a test set.  Compares against sample mean, random walk, and AIC benchmarks.  The resulting object is of class \code{BigVAR.results}
#' @return An object of class \code{BigVAR.results}.
#' @seealso \code{\link{constructModel}}, \code{\link{BigVAR.results}} 
#' @name cv.BigVAR
#' @aliases cv.BigVAR,BigVAR-method
#' @docType methods
#' @rdname cv.BigVAR-methods
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' m1=constructModel(Y,p=4,struct="None",gran=c(50,10),
#' RVAR=FALSE,MN=FALSE,h=1,cv="Rolling",verbose=FALSE)
#' results=cv.BigVAR(m1)
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
    k=object@nseries
    YYY <- .Zmat2(object@Data,p,ncol(object@Data))

    Y <- YYY$Y
    Z <- YYY$Z
  
    alpha <- object@alpha
    RVAR=object@Relaxed
    group=object@Structure
    cvtype=object@crossval
    MN <- object@Minnesota

    if(group=="Group"|group=="Sparse")
      {
      jj=.groupfun(p,k)
     }
    else{
      jj <- .lfunction3(p,k)
      }
    T0 <- floor(nrow(Y)/3)
    T1 <- floor(2 * nrow(Y)/3)
    T2 <- nrow(Y)

    ZZ1 <- .Zmat2(Y[1:(T1 - 1), ], p, k)
    trainZ <- ZZ1$Z
    trainY1 <- ZZ1$Y
    gran2 = object@Granularity[2]
    gran1=object@Granularity[1]	
    gamm <- .LambdaGrid(gran1, gran2, jj, trainY1, trainZ,group,p)
    h <- object@horizon
    verbose <- object@verbose 
        
    activeL <- rep(list(list()), gran2)
    beta=array(0,dim=c(k,k*p+1,gran2))

    ZFull <- .Zmat2(Y, p, k)
    activeset = activeL
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
            q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
        }

      }
     if(cvtype=="LOO")
       {
         T0=p+k
         }
         MSFE <- matrix(0, nrow = T1 - T0 + 1, ncol = gran2)
  if(verbose==TRUE){
    pb <- txtProgressBar(min = T0, max = T1, style = 3)
    cat("Cross-Validation Stage:",group)}
    for (v in T0:T1) {
      if(cvtype=="Rolling")
        {
        trainY <- ZFull$Y[1:(v-p-1), ]
        trainZ <- ZFull$Z[, 1:(v-p-1)]
        }
      if(cvtype=="LOO")
        {
          YTrainF <- trainY1[-v,]
          ZZ <- .Zmat2(YTrainF,p,k)
          trainY <- ZZ$Y
          trainZ <- ZZ$Z
          
          
        }  
        if (group == "None") {
         
           beta <- .lassoVARFist(beta, trainZ, trainY,gamm, 1e-05,p,MN)
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

        eZ <- as.matrix(.Zmat(Y[(v - p):v, ], p, k)$Z, ncol = 1)
        eZh <- as.matrix(.Zmat2(Y[(v - p):v, ], p, k)$Z, ncol = 1)
  
        for (ii in 1:gran2) {
            if (RVAR == TRUE) {
                  beta[,,ii] <- .relaxedVAR(trainY, beta[,,ii], 
                  trainZ,k,p)
                }
            
        if(h==1)
            {
            MSFE[v - (T0 - 1), ii] <- norm2(Y[v, ] - beta[,,ii] %*% eZ)^2


            if(MN==TRUE){ MSFE[v - (T0 - 1), ii] <- norm2(Y[v, ] - beta[,2:ncol(beta[,,ii]),ii] %*% eZh)^2
}
            }
        else{
            if(v+h>T1){break}
            betaP <- beta[,,ii]
            nu=betaP[,1]
            betaP <- betaP[,2:ncol(betaP)]
            Bi=list()
            jj <- .groupfun(p,k)
            for(i in 1:p)
    {
        Bi[[i]] <- betaP[,jj[[i]]]


    }
            B <- VarptoVar1(Bi,p,k)

            MSFE[v - (T0 - 1), ii] <- norm2(Y[v+h,]-(B%^%h)%*%eZh)^2
            }
            
        }

        if(verbose==TRUE){
         setTxtProgressBar(pb, v)}

      }
     gamopt <- gamm[which.min(colMeans(na.omit(MSFE)))]
     OOSEval <- .EvalLVAR(Y,Z,gamopt,k,p,group,h,MN,verbose,RVAR)
     MSFEOOSAgg <- OOSEval$MSFE
     betaPred <- OOSEval$betaPred
     Zvals <- OOSEval$zvals
     MSFEOOS<-mean(na.omit(MSFEOOSAgg))
     seoos=sd(na.omit(MSFEOOSAgg))/sqrt(length(na.omit(MSFEOOSAgg)))

	meanbench=.evalMean(Y,T1,T2)
     if(k>1){
	AICbench=.evalAIC(Y,T1,T2,verbose)}
     else{AICbench=.evalAR(Y,T1,T2,p)}
	RWbench=.evalRW(Y,T1,T2)
    results <- new("BigVAR.results",InSampMSFE=colMeans(MSFE),LambdaGrid=gamm,index=which.min(colMeans(na.omit(MSFE))),OptimalLambda=gamopt,OOSMSFE=MSFEOOS,seoosmsfe=seoos,MeanMSFE=meanbench$Mean,AICMSFE=AICbench$Mean,RWMSFE=RWbench$Mean,MeanSD=meanbench$SD,AICSD=AICbench$SD,RWSD=RWbench$SD,Data=object@Data,lagmax=object@lagmax,Structure=object@Structure,Relaxed=object@Relaxed,Granularity=object@Granularity,horizon=object@horizon,betaPred=betaPred,Zvals=Zvals)

    return(results)
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
#' @field betaPred The final out of sample coefficient matrix of \code{B}, to be used for prediction
#' @field Zvals The final lagged values of \code{Y}, to be used for prediction  
#' @field Data a \eqn{T x k} multivariate time Series
#' @field lagmax Maximal lag order
#' @field Structure Penalty Structure
#' @field Relaxed Indicator for relaxed VAR
#' @field Granularity Granularity of Penalty Grid
#' @field horizon Desired Forecast Horizon
#' @field crossval Cross-Validation Procedure
#' @field alpha penalty for Sparse Group Lasso
#' @field nseries Number of Series
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
         representation(InSampMSFE="numeric",LambdaGrid="numeric",index="numeric",OptimalLambda="numeric",OOSMSFE="numeric",seoosmsfe="numeric",MeanMSFE="numeric",AICMSFE="numeric",RWMSFE="numeric",MeanSD="numeric",AICSD="numeric",RWSD="numeric",betaPred="matrix",Zvals="matrix"),
         contains="BigVAR"
         )


#' Plot an object of class BigVAR.results
#' 
#' @param x BigVAR.results object created from \code{cv.BigVAR}
#' @param y Needed for compatibility with generic, otherwise ignored
#' @param ... additional arguments passed to generic
#' @details Plots optimal lambda value
#' @name plot.BigVAR.results
#' @import methods
#' @aliases plot,BigVAR.results,ANY-method
#' @docType methods
#' @rdname plot-methods-BigVAR.results
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
#' @param object BigVAR.results object created from \code{ConstructModel} 
#' @details prints forecast information and comparisons with mean, random walk, and AIC benchmarks
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
            cat("Forecast Horizon \n") ;print(object@horizon)
            cat("Relaxed VAR \n") ;print(object@Relaxed)
            cat("Maximum Lag Order \n") ;print(object@lagmax)
            cat("Optimal Lambda \n"); print(object@OptimalLambda)
            cat("In-Sample MSFE\n");print(min(object@InSampMSFE))
            cat("BigVAR Out of Sample MSFE\n");print(min(object@OOSMSFE))
	    cat("BigVAR Out of Sample Standard Error\n");print(min(object@seoosmsfe))

            cat("Unconditional Mean Out of Sample MSFE\n");print(min(object@MeanMSFE))
	    cat("Unconditional Mean Standard Error\n");print(min(object@MeanSD))
            cat("AIC Out of Sample MSFE\n");print(min(object@AICMSFE))
            cat("AIC Out of Sample Standard Error\n");print(min(object@AICSD))
	    cat("RW Out of Sample MSFE\n");print(min(object@RWMSFE))
	    cat("RW Out of Sample Standard Error\n");print(min(object@RWSD)) 	
            cat("Index of Optimal Lambda \n");print(object@index)
            }
)
#' Forecast using a BigVAR.results object
#' @usage predict(object,...)
#' @param object BigVAR.results object from \code{cv.BigVAR}
#' @param ... additional arguments affecting the predictions produced 
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
#' m1=constructModel(Y,p=4,struct="None",gran=c(50,10),MN=FALSE,
#' RVAR=FALSE,h=1,cv="Rolling",verbose=FALSE)
#' results=cv.BigVAR(m1)
#' predict(results,n.ahead=1)
#' @export
setMethod("predict","BigVAR.results",
          function(object,n.ahead,...)
{
eZ <- object@Zvals
betaPred <- object@betaPred
k <- nrow(betaPred)
p <- object@lagmax
nu=betaPred[,1]
betaPred <- betaPred[,2:ncol(betaPred)]
Bi=list()
jj <- .groupfun(p,k)
for(i in 1:p)
    {
Bi[[i]] <- betaPred[,jj[[i]]]


        }
B <- VarptoVar1(Bi,p,k)
if(n.ahead==1)
    {
fcst <- betaPred%*%eZ+nu
        }
else{
fcst <- ((B%^%n.ahead)%*%eZ)[1:k,]+nu
    }
return(fcst)

    }
)

#' Sparsity Plot of a BigVAR.results object 
#'
#' @param object BigVAR.results object
#' @return NA, side effect is graph
#' @details Uses levelplot from the lattice package to plot the magnitude of each coefficient
#' @name SparsityPlot.BigVAR.results
#' @import lattice
#' @aliases SparsityPlot.BigVAR.results,BigVAR.results-method
#' @docType methods
#' @seealso \code{\link{SparsityPlot}}
#' @rdname SparsityPlot.BigVAR.results-methods
#' @export
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
      B <- B[,2:ncol(B)]
      k <- nrow(B)
      p <- object@lagmax
text <- c()
for(i in 1:p)
    {
        ## text1 <- (paste("expression(bold(B)[",i,"])",sep="",collapse="")))
        text1 <- as.expression(bquote(bold(B)[.(i)]))
        text <- append(text,text1)
          
    }
f <- function(m) t(m)[,nrow(m):1]

rgb.palette <- colorRampPalette(c("white", "blue" ),space = "Lab")


at <- seq(k/2+.5,p*(k)+.5,by=k)

se2=seq(1.75,by=k,length=k)


L2 <- levelplot(f(abs(B)),col.regions=rgb.palette,colorkey=NULL,xlab=NULL,ylab=NULL,main="Sparsity Pattern Generated by BigVAR",panel=function(...){ panel.levelplot(...);panel.abline(a=NULL,b=1,h=seq(1.5,p*k+.5,by=1),v=seq(1.5,by=1,length=p*k));panel.abline(a=NULL,b=1,v=seq(k+.5,p*k+.5,by=k),lwd=3)},scales=list(x=list(alternating=1,labels=text,at=at,tck=c(0,0)),y=list(alternating=0,tck=c(0,0))))

return(print(L2))
}
    )
      

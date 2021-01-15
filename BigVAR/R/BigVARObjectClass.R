                                        # Ensures that the created BigVAR object is valid
check.BigVAR <- function(object){

    errors <- character()

    VARX <- object@VARX
    Y <- object@Data
    
    if(any(is.na(Y))){msg <- c("Remove NA values before running ConstructModel")

        errors <- c(errors,msg)
    }      
    if(dim(Y)[2]>dim(Y)[1] & length(VARX)==0){msg <- paste("k is",ncol(Y),"which is greater than T, is Y formatted correctly (k x T)?")}      
    if(object@lagmax<0){msg <- c("Maximal lag order must be at least 0")
        errors <- c(errors,msg)
    }
    if(object@lagmax==0& !object@Structure%in%c("Basic","BasicEN")){
        msg <- c("Only Basic VARX-L supports a transfer function")
        errors <- c(errors,msg)
    }
    structures=c("Basic","Lag","SparseLag","OwnOther","SparseOO","HVARC","HVAROO","HVARELEM","Tapered","EFX","BGR","BasicEN","MCP","SCAD")
    cond1=object@Structure%in% structures
    if(cond1==FALSE){
        msg <- paste("struct must be one of",structures)
        errors <- c(errors,msg)      
    }
    if(object@horizon<1){msg <- paste("Forecast Horizon is ",object@horizon, " must be at least 1")

    }

    
    if(object@crossval!="Rolling" & object@crossval!="LOO"){msg <- c("Cross-Validation type must be one of Rolling or LOO")
        errors <- c(errors,msg)
    }
    if(length(object@Granularity)!=2&object@ownlambdas==FALSE){msg("Granularity must have two parameters")
        errors <- c(errors,msg)
        
    }
    if(any(object@Granularity<=0)){
        msg <- c("Granularity parameters must be positive")
        errors <- c(errors,msg)
    }
    structure2 <- c("Basic","Lag","HVARC","BasicEN","MCP","SCAD")
    cond2=object@Structure %in% structure2
    k1=0
    if(length(VARX)!=0){

        k1 <- VARX$k
        if(k1>ncol(Y)){msg <- c("k is greater than the number of columns in Y")

            errors <- c(errors,msg)
        }        
    }else{k=0}
    m <- ncol(Y)-k1

    if(object@tf & object@lagmax!=0){
        msg <- c("p must be 0 if fitting a transfer function")
        errors <- c(errors,msg)
    }
    nseries <- ncol(Y)-ifelse(m<ncol(Y),m,0)    
    if(nseries==1 & cond2==FALSE ){
        msg <- c("Univariate support is only available for Basic VARX-L, Lag Group VARX-L, and Componentwise HVAR")

        errors <- c(errors,msg)
        
    }
    if(length(VARX)==0 & object@Structure=="EFX"){

        msg <- c("EFX is only supported in the VARX framework")

        errors <- c(errors,msg)          

    }

    if(is.list(VARX) & length(VARX)>0 & !(exists('k',where=VARX) & exists('s',where=VARX)))
    {

        msg <- c("VARX Specifications entered incorrectly")

        errors <- c(errors,msg)
    }

    
    if(object@Structure=="EFX" & !is.null(VARX$contemp)){
        if(VARX$contemp){
            msg <- c("EFX does not support contemporaneous dependence")
            errors <- c(errors,msg)
        }

    }
    structs=c("HVARC","HVAROO","HVARELEM")
    if(length(VARX)!=0& object@Structure %in% structs){
        msg <- c("EFX is the only nested model supported in the VARX framework")

        errors <- c(errors,msg)

    }
    if(object@T1>nrow(Y) | object@T2>nrow(Y) |object@T2<object@T1){
        msg <- c("Training dates exceed series length")
        errors <- c(errors,msg)

    }

    if(any(object@alpha<0) || any(object@alpha>1)){
        msg <- c("alpha must be between zero and 1")
        errors <- c(errors,msg)
    }
    if(object@recursive & length(VARX)>0){
        msg <- c("recursive forecasts can only be used with VAR models")
        errors <- c(errors,msg)
    }    
    
    if(length(errors)==0) TRUE else errors
    
}

#' BigVAR Object Class
#'
#' An object class to be used with cv.BigVAR
#' 
#' @slot Data a \eqn{T \times k} multivariate time Series
#' @slot lagmax Maximal lag order for modeled series
#' @slot intercept Indicator as to whether an intercept should be included 
#' @slot Structure Penalty Structure
#' @slot Relaxed Indicator for relaxed VAR
#' @slot Granularity Granularity of Penalty Grid
#' @slot horizon Desired Forecast Horizon
#' @slot crossval Cross-Validation Procedure
#' @slot Minnesota Minnesota Prior Indicator
#' @slot verbose Indicator for Verbose output
#' @slot dates dates extracted from an xts object 
#' @slot ic Indicator for including AIC and BIC benchmarks
#' @slot VARX VARX Model Specifications
#' @slot T1 Index of time series in which to start cross validation
#' @slot T2  Index of times series in which to start forecast evaluation
#' @slot ONESE Indicator for "One Standard Error Heuristic"
#' @slot ownlambdas Indicator for user-supplied lambdas
#' @slot tf Indicator for transfer function
#' @slot alpha Grid of candidate alpha values (applies only to Sparse VARX-L and Elastic Net models)
#' @slot recursive Indicator as to whether recursive multi-step forecasts are used (applies only to multiple horizon VAR models)
#' @slot constvec vector indicating variables to shrink toward a random walk instead of toward zero (valid only if Minnesota is \code{TRUE})
#' @slot tol optimization tolerance
#' @slot window.size size of rolling window.  If set to NULL an expanding window will be used.
#' @slot separate_lambdas indicator to use separate penalty parameter for each time series (default \code{FALSE})
#' @slot loss Loss function to select penalty parameter (one of "L1","L2","Huber").
#' @slot delta delta parameter for Huber loss (default 2.5)
#' @slot rolling_oos True or False: indicator to update the penalty parameter over the evaluation period (default \code{False})  

#' @details To construct an object of class BigVAR, use the function \code{\link{constructModel}}
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
        intercept="logical",
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
        tf="logical",
        alpha="numeric",
        recursive="logical",
        dates="character",
        constvec="numeric",
        tol="numeric",
        window.size="numeric",
        separate_lambdas="logical",
        loss="character",
        delta='numeric',
        rolling_oos="logical"
    ),validity=check.BigVAR
)


#' Construct an object of class BigVAR
#' @param Y \eqn{T \times k} multivariate time series or Y \eqn{T \times (k+m)} endogenous and exogenous series, respectively. 
#' @param p Predetermined maximal lag order (for modeled series).
#' @param struct The choice of penalty structure (see details).
#' @param gran vector of penalty parameter specifications.
#' @param intercept True or False: option to fit an intercept.
#' @param RVAR True or False: option to refit based upon the support selected using the Relaxed-VAR procedure.
#' @param h Desired forecast horizon.
#' @param cv Cross-validation approach, either "Rolling" for rolling cross-validation or "LOO" for leave-one-out cross-validation.
#' @param MN Minnesota Prior Indicator.
#' @param verbose Verbose output while estimating.
#' @param IC True or False: whether to include AIC and BIC benchmarks.
#' @param VARX List containing VARX model specifications. 
#' @param T1 Index of time series in which to start cross validation.
#' @param T2  Index of times series in which to start forecast evaluation.
#' @param ONESE True or False: whether to use the "One Standard Error Heuristic."
#' @param ownlambdas True or False: Indicator for user-supplied penalty parameters.
#' @param alpha grid of candidate parameters for the alpha in the Sparse Lag and Sparse Own/Other VARX-L.
#' @param recursive True or False: Indicator as to whether iterative multi-step predictions are desired in the VAR context if the forecast horizon is greater than 1.
#' @param C vector of coefficients to shrink toward a random walk (if \code{MN} is \code{TRUE}).
#' @param tol optimization tolerance (default 1e-4).
#' @param dates optional vector of dates corresponding to \eqn{Y}.
#' @param separate_lambdas indicator for separate penalty parameters for each time series (default \code{FALSE}).
#' @param window.size size of rolling window.  If set to 0 an expanding window will be used. 
#' @param linear indicator for linearly decrementing penalty grid (FALSE is log-linear).
#' @param loss Loss function to select penalty parameter (one of "L1","L2","Huber").
#' @param delta delta parameter for Huber loss (default 2.5)
#' @param rolling_oos True or False: indicator to update the penalty parameter over the evaluation period (default \code{False})  
#' @details The choices for "struct" are as follows
#' \itemize{
#' \item{  "Basic" (Basic VARX-L)}
#' \item{  "BasicEN" (Elastic Net VARX-L)}
#' \item{  "Lag" (Lag Group VARX-L)} 
#' \item{  "SparseLag" (Lag Sparse Group VARX-L)} 
#' \item{  "OwnOther" (Own/Other Group VARX-L) }
#' \item{  "SparseOO" (Own/Other Sparse Group VARX-L) }
#' \item{  "EFX" (Endogenous First VARX-L)}
#' \item{  "HVARC" (Componentwise HVAR) }
#' \item{  "HVAROO" (Own/Other HVAR) }
#' \item{  "HVARELEM" (Elementwise HVAR)}
#' \item{  "Tapered" (Lag weighted Lasso VAR)}
#' \item{  "BGR" (Bayesian Ridge Regression (cf. Banbura et al))}
#' }
#'
#' The first number in the vector "gran" specifies how deep to construct the penalty grid and the second specifies how many penalty parameters to use  If ownlambas is set to TRUE, gran should contain the user-supplied penalty parameters.
#' 
#' VARX specifications consist of a list with entry k denoting the series that are to be modeled and entry s to denote the maximal lag order for exogenous series.
#'
#' The argument alpha is ignored unless the structure choice is "SparseLag" or "Lag."  By default "alpha" is set to \code{NULL} and will be initialized as 1/(k+1) in \code{cv.BigVAR} and \code{BigVAR.est}.  Any user supplied values must be between 0 and 1.  

#' @note The specifications "Basic","BasicEN", "Lag," "SparseLag," "SparseOO," and "OwnOther" can accommodate both VAR and VARX models.  EFX only applies to VARX models.  "HVARC," "HVAROO," "HVARELEM," and "Tapered" can only be used with VAR models.
#'
#' @seealso \code{\link{cv.BigVAR}},\code{\link{BigVAR.est}}
#' 
#' @references  Nicholson, William, I. Wilms, J. Bien, and D. S. Matteson. High dimensional forecasting via interpretable vector autoregression. Journal of Machine Learning Research, 21(166):1–52, 2020.
#' William B. Nicholson, David S. Matteson, Jacob Bien,VARX-L: Structured regularization for large vector autoregressions with exogenous variables, International Journal of Forecasting, Volume 33, Issue 3, 2017, Pages 627-651,
#' William B Nicholson, David S. Matteson, and Jacob Bien (2016), "BigVAR: Tools for Modeling Sparse High-Dimensional Multivariate Time Series" arxiv:1702.07094
#'
#' Banbura, Marta, Domenico Giannone, and Lucrezia Reichlin. "Large Bayesian vector auto regressions." Journal of Applied Econometrics 25.1 (2010): 71-92.
#' @examples
#' # VARX Example
#' # Create a Basic VARX-L with k=2, m=1, s=2, p=4
#' VARX=list()
#' VARX$k=2 # indicates that the first two series are modeled
#' VARX$s=2 # sets 2 as the maximal lag order for exogenous series
#' data(Y)
#' T1=floor(nrow(Y)/3)
#' T2=floor(2*nrow(Y)/3)
#' Model1=constructModel(Y,p=4,struct="Basic",gran=c(50,10),verbose=FALSE,VARX=VARX,T1=T1,T2=T2)
#' @export
constructModel <- function(Y,p,struct,gran,RVAR=FALSE,h=1,cv="Rolling",MN=FALSE,verbose=TRUE,IC=TRUE,VARX=list(),T1=floor(nrow(Y)/3),T2=floor(2*nrow(Y)/3),ONESE=FALSE,ownlambdas=FALSE,alpha=as.double(NULL),recursive=FALSE,C=as.double(NULL),dates=as.character(NULL),intercept=TRUE,tol=1e-4,window.size=0,separate_lambdas=FALSE,linear=TRUE,loss="L2",delta=2.5,rolling_oos=FALSE)
{
    if(any(is.na(Y))){stop("Remove NA values before running constructModel")}      
    if(dim(Y)[2]>dim(Y)[1] & length(VARX)==0){warning("k is greater than T, is Y formatted correctly (k x T)?")}      
    if(p<0){stop("Maximal lag order must be at least 0")}
    if(p==0& !struct%in%c("Basic","BasicEN")){stop("Only Basic VARX-L supports a transfer function")}
    oldnames <- c("None","Diag","SparseDiag")
    if(struct%in%oldnames) stop("Naming Convention for these structures has changed. Use Basic, OwnOther, and SparseOO.")
    structures=c("Basic","Lag","SparseLag","OwnOther","SparseOO","HVARC","HVAROO","HVARELEM","Tapered","EFX","BGR","BasicEN","MCP","SCAD")
    if(struct=="BasicEN"&length(alpha)>1&separate_lambdas){stop("Cannot perform separate lambdas per series and range of alphas simultaneously")}
    cond1=struct %in% structures
    if(!loss%in%c("L1","L2","Huber")){stop("loss must be one of L1,L2,Huber")}
    if(!cond1){stop(cat("struct must be one of",structures))}
    if(h<1){stop("Forecast Horizon must be at least 1")}
    if(cv!="Rolling" & cv!="LOO"){stop("Cross-Validation type must be one of Rolling or LOO")}
    if(length(gran)!=2&ownlambdas==FALSE){stop("Granularity must have two parameters")}
    if(any(gran<=0)){stop("Granularity parameters must be positive")}
    if(tol<0 | tol>1e-1){stop("Tolerance must be positive")}
    if(window.size>nrow(Y) | window.size<0){stop("window size must be shorter than the series length")}
    if(delta<0){stop("huber delta must be positive")}
    bss <- c("Basic","HVARC","HVAROO","HVARELEM","Tapered")
                                        # check
    ## if(rolling_oos & separate_lambdas){stop("Cannot estimate rolling out of sample and separate lambdas jointly.")}

    if(separate_lambdas & !struct%in%c("Basic","HVARC","HVAROO","HVARELEM","Tapered","BasicEN")){stop(print(cat("separate lambda estimation only available for ",bss)))}

    start_ind <- (T1-p-h+1)
    ## print(start_ind)
    if(cv=="Rolling" & start_ind<5){
        stop("too few values for rolling validation, try running BigVAR.fit")
    }
    
    if( MN  & intercept){ intercept=FALSE }
    
    structure2 <- c("Basic","Lag","HVARC","BasicEN","MCP","SCAD")
    cond2=struct %in% structure2
    ## k <- ncol(Y)
    
    if(length(VARX)!=0){

        k <- VARX$k
        if(k>ncol(Y)){stop("k is greater than the number of columns in Y")}
    }else{k=ncol(Y)}
    m <- ncol(Y)-k
    nseries <- ncol(Y)-ifelse(m<ncol(Y),m,0)
    if(p==0){tf=TRUE
    }else{
        tf=FALSE
    }
    if(nseries==1 & cond2==FALSE ){stop("Univariate support is only available for Basic, Elastic Net, Lag Group, and Componentwise HVAR")}
    if(length(VARX)==0 & struct=="EFX"){stop("EFX is only supported in the VARX framework")}
    if(struct=="EFX" & !is.null(VARX$contemp)){
        if(VARX$contemp){
            stop("EFX does not support contemporaneous dependence")}
    }
    structs=c("HVARC","HVAROO","HVARELEM")
    if(length(VARX)!=0& struct %in% structs){stop("EFX is the only nested model supported in the VARX framework")}
    if(length(VARX)!=0& struct =="BGR"){stop("BGR is only available in the VAR framework")}
    if(length(VARX)!=0& struct =="Tapered"){stop("Lag-Weighted Lasso only available in VAR framework")}
    ## if(gran[2]==1 & separate_lambdas){
    ##     stop("Must have more than one lambda if fitting separate lambdas")
    ## }
    
    if(T1>nrow(Y) | T2>nrow(Y) |T2<T1){stop("Training dates exceed series length")}

    if(is.list(VARX) & length(VARX)>0 & !(exists('k',where=VARX) & exists('s',where=VARX)))
    {

        stop("VARX Specifications entered incorrectly")

    }

    if(!is.null(alpha)){
        if(any(alpha<0) || any(alpha>1)){stop("alpha must be between 0 and 1")}
    }
    if(length(C)!=0){

        if(length(C)!=k){stop("C must have length k")}
        if(!all(C%in%c(0,1))){stop("Values of C must be either 0 or 1")}
        
    }else{

        C <- rep(1,k)
        
    }
    ## if("xts"%in%class(Y)){
    ##     ind <- as.character(index(Y))
    ##     Y <- as.matrix(Y)

    if(length(dates)!=0){

        ind <- dates

    }else{
        ind <- as.character(NULL)
    }
    
                                        # Can't have a class named C
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
         tf=tf,
         alpha=alpha,
         recursive=recursive,
         dates=ind,
         constvec=C,
         intercept=intercept,
         tol=tol,
         window.size=window.size,
         separate_lambdas=separate_lambdas,
         loss=loss,
         delta=delta,
         rolling_oos=rolling_oos
     ))

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
setMethod("show","BigVAR",
          function(object)
          {
              
              
              T1P <- ifelse(length(object@dates)!=0,object@dates[object@T1],object@T1)

              T2P <- ifelse(length(object@dates)!=0,object@dates[object@T2],object@T2)

              nrowShow <- min(5,nrow(object@Data))
              cat("*** BIGVAR MODEL *** \n")
              cat("Data (First 5 Observations):\n")
              print(formatC(object@Data[1:nrowShow,]),quote=FALSE)
              cat("Structure\n") ;print(object@Structure)
              cat("Forecast Horizon \n") ;print(object@horizon)
              cat("Relaxed VAR \n") ;print(object@Relaxed)
              cat("Minnesota Prior \n") ;print(object@Minnesota)
              cat("Maximum Lag Order \n") ;print(object@lagmax)
              if(length(object@VARX)!=0){
                  cat("VARX Specs \n") ;print(object@VARX)}
              cat("Start of Cross Validation Period \n") ;print(T1P)
              cat("End of Cross Validation Period \n") ;print(T2P)
              
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
setMethod(f="plot",signature="BigVAR",
          definition= function(x,y=NULL,...)
          {

              T1P <- ifelse(length(x@dates)!=0,x@dates[x@T1],x@T1)

              T2P <- ifelse(length(x@dates)!=0,x@dates[x@T2],x@T2)

              g=ncol(x@Data)
              names <- ifelse(rep(!is.null(colnames(x@Data)),ncol(x@Data)),colnames(x@Data),as.character(1:g))
              if(length(x@dates)!=0){

                  dates <- as.yearqtr(x@dates)
              }else{

                  dates <- 1:nrow(x@Data)
              }

              Yzoo <- zoo(as.matrix(x@Data),order.by=dates)
              plot.zoo(Yzoo,plot.type="single",col=1:g)
              legend('topright',names,lty=1,col=1:g)

              abline(v=index(Yzoo[as.yearqtr(T1P)]))
              abline(v=index(Yzoo[as.yearqtr(T2P)]))
              
          }
          )

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
#' Model1=constructModel(Y,p=4,struct="Basic",gran=c(50,10))
#' results=cv.BigVAR(Model1)
#' @importFrom abind adrop
#' @importFrom abind abind
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
        p <- object@lagmax
        s1 <- 0
        Y <- object@Data
        k <- ncol(Y)
        alpha <- object@alpha
        RVAR <- object@Relaxed
        group <- object@Structure
        cvtype <- object@crossval
        intercept=object@intercept
        recursive <- object@recursive
        VARX <- object@VARX
        tol=object@tol
        window.size=object@window.size
        verbose <- object@verbose
        loss <- object@loss
        delta <- object@delta
        rolling_oos=object@rolling_oos
        if(length(alpha)==0){

            if(length(VARX)>0){    
                alpha <- 1/(VARX$k+1)
                
            }else{

                alpha <- 1/(k+1)
                
            }
        }

        C <- object@constvec
        
        if(length(alpha)>1 & group%in%c("SparseLag","SparseOO","BasicEN"))
        {
            ## browser()
            dual <- TRUE

        }else{

            dual <- FALSE
        }

        MN <- object@Minnesota
        h <- object@horizon
        jj <- 0
        separate_lambdas <- object@separate_lambdas
        if(!"matrix"%in%class(Y)){Y=matrix(Y,ncol=1)}

        if(object@crossval=="Rolling"){
            T1 <- object@T1
                                        # ensure it's an integer
            T1 <- floor(T1)

        }else{

            T1 <- p+2    

        }
        T2 <- object@T2
        T2 <- floor(T2)
        s <- ifelse(length(object@VARX)!=0,object@VARX$s,0)
        ONESE <- object@ONESE
        if(object@ownlambdas){
            gamm <- object@Granularity
            gran2 <- length(gamm)
            if(gran2==1){
                ONESE <-FALSE
            }
            
        }     
        if(object@Granularity[2]==1){
            stop("only one penalty parameter; run BigVAR.est instead of cv.BigVAR")
        }
        
        


                                        # Adjust T1, T2 by maximum lag order to account for initialization
        if(object@crossval=="Rolling"){
            T1 <- T1-max(p,s)
            T2 <- T2-max(p,s)
        }
        if(!object@ownlambdas){

            gran2 <- object@Granularity[2]
            gran1 <- object@Granularity[1]

        }
        


                                        # Constructing lag matrix in VARX setting
        if(length(VARX)!=0){

            VARX <- TRUE
            k1 <- object@VARX$k
            s <- object@VARX$s

            if(!is.null(object@VARX$contemp)){

                contemp <- TRUE
                s1 <- 1

            }else{

                contemp <- FALSE
                s1 <- 0
            }
            if(window.size==0){
                window.size=0
            }else{
                window.size <- window.size+1
            }

            m <- k-k1
            Y1 <-matrix(Y[,1:k1],ncol=k1)
            X <- matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m)

            if(!object@tf){
                trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)

            }else{

                trainZ <- VARXCons(matrix(0,ncol=1,nrow=nrow(X)),matrix(X,ncol=m),k=0,p=0,m=m,s=s,contemp=contemp,oos=FALSE)

            }
            
            trainZ <- trainZ[2:nrow(trainZ),,drop=F]

            trainY <- matrix(Y[(max(c(p,s))+1):nrow(Y),1:k1],ncol=k1)

            if(group=="Lag"|group=="SparseLag"){

                                        # Lag based groupings
                jj <- groupfunVARXLG(p,k,k1,s+s1)

                                        # initialize activeset as empty
                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2*length(alpha))

            }
            if(group=="OwnOther"|group=="SparseOO"){
                                        # Own other based groupings
                jj <- diaggroupfunVARXLG(p,k,k1,s+s1)
                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2)

            }
            if(object@ownlambdas==FALSE){

                
                if(dual){
                                        # Constructs penalty grid if both alpha and lambda are selected
                    gamm <- .LambdaGridXDual(gran1, gran2, jj, trainY, trainZ,group,p,k1,s,m,k,MN,alpha,C,intercept,tol)

                }else{

                                        # Penalty parameter grid for just lambda
                    gamm <- .LambdaGridX(gran1, gran2, jj, as.matrix(trainY[1:T2,]), trainZ[,1:T2],group,p,k1,s+s1,m,k,MN,alpha,C,intercept,tol,separate_lambdas,verbose)
                }
            }

                                        # Coefficient matrices
            ## if(!dual){

            ##     beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2))

            ## }else{

            beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2*length(alpha)))
            
            ## }

            q1a <- NULL
            kk <- NULL
            jj <- NULL
            jjcomp <- NULL
            activeset <- NULL
            q1a <- NULL

                                        # Groupings in accordance with C++ indexing standards
            if (group == "Lag") {
                jj <- groupfunVARX(p,k,k1,s+s1)
                jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)
                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2)


            }else if (group == "SparseLag") {
                
                jj <- groupfunVARX(p, k,k1,s+s1)
                q1a <- list()
                                        # Initializing warm start vectors for the power method
                for (i in 1:(p+s+s1)) {
                    q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
                }

                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2*length(alpha))
            }else if (group == "OwnOther") {
                kk <- diaggroupfunVARX(p,k,k1,s+s1)
                activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                                 gran2)
            }else if (group == "SparseOO") {
                kk <- diaggroupfunVARX(p,k,k1,s+s1)
                activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                                 gran2*length(alpha))
                q1a <- list()

                for (i in 1:(2*p+s+s1)) {
                    q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)
                }
                
            }else{
                kk <- NULL
                jj <- NULL
                jjcomp <- NULL
                activeset <- NULL
                q1a <- NULL


            }
        }else{
                                        # VAR estimation
            contemp <- FALSE
            if(group=="Lag"|group=="SparseLag")
            {
                jj <- .groupfun(p,k)
            }else{
                jj <- .lfunction3(p,k)
            }
            Z1 <- VARXCons(Y,matrix(0,nrow=nrow(Y)),k,p,0,0) 

            trainZ <- Z1[2:nrow(Z1),]   

            trainY <- matrix(Y[(p+1):nrow(Y),],ncol=k)          

            GY <- matrix(trainY[1:T2,],ncol=k)

                                        # We only use training period data to construct the penalty grid
            GZ <- trainZ[,1:T2]

            if(object@ownlambdas==FALSE){


                                        # Find starting values for penalty grid

                if(dual){

                    gamm <- .LambdaGridEDual(gran1, gran2, jj, GY, GZ,group,p,k,MN,alpha,C,intercept,tol)

                }else{

                    if(group!="BGR"){
                        
                        gamm <- .LambdaGridE(gran1, gran2, jj, GY, GZ,group,p,k,MN,alpha,C,intercept,tol,separate_lambdas =separate_lambdas ,verbose)

                    }else{

                                        # BGR operates on a fixed grid
                        gamm <- seq(1,5,by=.025)
                        gamm <- gamm*sqrt(k*p)
                        
                    }

                }
                
            }
            VARX <- FALSE
            k1 <- k
            s <- 0
            kk <- NULL
            jj <- NULL
            jjcomp <- NULL
            activeset <- NULL
            q1a <- NULL

            if (group == "Lag") {
                jj <- .groupfuncpp(p, k)
                jjcomp <- .groupfuncomp(p,k)
                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2)
                
            }else if (group == "SparseLag") {
                jj <- .groupfuncpp(p, k)
                q1a <- list()
                for (i in 1:p) {
                    q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
                }
                
                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2*length(alpha))
                
            }else if (group == "OwnOther") {
                kk <- .lfunction3cpp(p, k)
                activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                                 gran2)
            }else if (group == "SparseOO") {
                kk <- .lfunction3cpp(p, k)
                jjcomp <- .lfunctioncomp(p,k)
                jj <- .lfunction3(p,k)
                activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                                 gran2*length(alpha))
                q1a <- list()
                for (i in 1:(2*p)) {
                    q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)
                }

                
            }else{
                kk <- NULL
                jj <- NULL
                jjcomp <- NULL
                activeset <- NULL
                q1a <- NULL

            }
            
            beta <- array(0,dim=c(k,k*p+1,gran2*length(alpha)))
        }           
        h <- object@horizon
        
        ZFull <- list()

        if(!is.matrix(trainZ)){trainZ <- matrix(trainZ,ncol=1)}

        if(!is.matrix(trainY)){trainY <- matrix(trainY,ncol=1)}
        
        ZFull$Z <- trainZ
        ZFull$Y <- trainY
        T <- nrow(trainY)

        if(object@ownlambdas){gamm <- object@Granularity}    

                                        # Constructs grid for lag weighted lasso
        if(group=="Tapered")
        {

            palpha <- seq(0,1,length=10)
            palpha <- rev(palpha)
            gran2 <- length(gamm)*length(palpha)
            beta <- array(0,dim=c(k,k*p+1,gran2))

        }else{
            palpha <- NULL
        }

        if(!"matrix"%in%class(ZFull$Y)){
            ZFull$Y <- matrix(ZFull$Y,ncol=1)
        }

        if(!dual){
            if(separate_lambdas  ){
                if(!VARX){
                    MSFE <- array(0,dim=c(T2-T1+1,gran2,k))
                }else{
                    MSFE <- array(0,dim=c(T2-T1+1,gran2,k1))
                }
            }else{
                MSFE <- matrix(0, nrow = T2 - T1+1, ncol = gran2)
                gamm <- as.matrix(gamm)
            }

        }else{

            nalpha <- length(alpha)
            MSFE <- matrix(0, nrow = T2 - T1+1, ncol = gran2*nalpha)
            
        }
        if(verbose){
            pb <- txtProgressBar(min = T1, max = T2, style = 3)
            cat("Cross-Validation Stage:",group)}
        YT <- Y[1:T2,]
        
                                        # Start of penalty parameter selection     
        betaWS <- beta
        for (v in (T1-h+1):T2) {
            if(cvtype=="Rolling")

            {

                if(v+h-1>T2){
                    break
                }

                if(h>1 & !recursive){

                    if(window.size!=0){
                        ws1 <- max(c(v-window.size-h,1))
                        trainY <- ZFull$Y[(ws1+h):(v-1), ]
                        trainZ <- ZFull$Z[, (ws1+h):(v-h)]         
                    }else{

                        trainY <- ZFull$Y[(h):(v-1), ]
                        
                        trainZ <- ZFull$Z[, 1:(v-h),drop=F]
                    }
                    
                }else{
                    if(window.size!=0){
                        ws1 <- max(c(v-window.size,1))
                        trainY <- ZFull$Y[(ws1):(v-1), ]
                        trainZ <- ZFull$Z[, (ws1):(v-1)]         
                    }else{
                        trainY <- ZFull$Y[(1):(v-1), ]                       
                        trainZ <- ZFull$Z[, (1):(v-1),drop=F]
                    }
                }
            }else{

                if(VARX)

                {

                    YT2 <- YT[-v,]

                    Y1 <- matrix(YT2[,1:k1],ncol=k1)

                    X <- matrix(YT2[,(ncol(YT2)-m+1):ncol(YT2)],ncol=m)

                    trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)

                    trainZ <- trainZ[2:nrow(trainZ),,drop=FALSE]

                    trainY <- matrix(YT2[(max(c(p,s))+1):nrow(YT2),1:k1],ncol=k1)

                }else{


                    YT2 <- YT[-v,]

                    Z1 <- VARXCons(YT2,matrix(0,nrow=nrow(YT2)),k,p,0,0) 

                    trainZ <- Z1[2:nrow(Z1),]        

                    trainY <- matrix(YT2[(p+1):nrow(YT2),],ncol=k)                                  

                }

            }

            needed.objs <- c('group','beta','trainZ','trainY','gamm','tol','p','m','k1','s','s1','m','MN','C','intercept','separate_lambdas','dual','activeset','alpha','jj','jjcomp','kk','palpha')
            if(!group%in%c("SparseLag","SparseOO")){
                q1a <- NULL
            }
            ## objs <- sapply(needed.objs,exists)
            objs <- setdiff(needed.objs,ls())
            if(length(objs)>0){
                
                for(i in 1:length(objs)){               
                    assign(objs[i],NULL)
                }
            }
            temp <- .BigVAR.fit(group,betaWS,trainZ,trainY,gamm,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
            beta <- temp$beta
            betaWS <- temp$beta
            activeset <- temp$activeset
            q1a <- temp$q1a
            eZ <- c(1,ZFull$Z[,v])


            if(group!="BGR"){
                if(!dual)
                {

                    
                                        # Calculate h-step MSFE for each penalty parameter
                    for (ii in 1:gran2) {

                        if (RVAR)
                        {
                                        # Relaxed Least Squares (intercept ignored)
                            beta[,,ii] <- RelaxedLS(cbind(t(trainZ),trainY),beta[,,ii],k,p,k1,s+s1)
                        }
                        if (MN){
                                        # shrink toward vector random walk

                            pred <- beta[,2:dim(beta)[2],ii] %*% eZ[2:length(eZ)]

                            if(h>1 & recursive){
                                pred <- matrix(pred,nrow=1)
                                pred <- predictMS(pred,trainY,h-1,beta[,2:dim(beta)[2],ii],p,MN)
                            }
                            
                            if(separate_lambdas)
                            {
                                pred <- matrix(pred,ncol=1)  
                                for(uu in 1:ncol(ZFull$Y)){
                                    ## if(loss=="L2"){ 
                                    ##     MSFE[v - (T1 -h), ii,uu] <- (ZFull$Y[v+h-1,uu] - pred[uu,1])^2
                                    ## }else if(loss=="L1"){
                                    ##     MSFE[v - (T1 -h), ii,uu] <- abs(ZFull$Y[v+h-1,uu] - pred[uu,1])
                                    ## }else if(loss=="huber"){
                                    ##     MSFE[v - (T1 -h), ii,uu] <- .huber_loss(ZFull$Y[v+h-1,uu] - pred[uu,1],delta)
                                    ## }
                                    MSFE[v - (T1 -h), ii,uu] <- .calc.loss(ZFull$Y[v+h-1,uu] - pred[uu,1],univ=TRUE,loss,delta)
                                    ## }else{
                                }
                            }else{
                                MSFE[v - (T1 -h), ii] <- .calc.loss(ZFull$Y[v+h-1,1:k1] - pred,univ=FALSE,loss,delta)
                                ## if(loss=="L2"){ 
                                ##     MSFE[v - (T1 -h), ii] <- norm2(ZFull$Y[v+h-1,1:k1] - pred)^2
                                ## }else if(loss=="L1"){
                                ##     MSFE[v - (T1 -h), ii] <- sum(abs(ZFull$Y[v+h-1,1:k1] - pred)^2)
                                ## }else if(loss=="huber"){
                                ##     MSFE[v - (T1 -h), ii] <- sum(.huber_loss(ZFull$Y[v+h-1,1:k1] - pred,delta))
                                ## }

                            }
                                        # Subtract one from diagonal for warm start purposes
                            
                            diag(beta[,2:(k1+1),ii]) <- diag(beta[,2:(k1+1),ii])-C
                            
                        }else{

                            if(object@crossval=="Rolling"){

                                pred <- beta[,,ii] %*% eZ

                                if(h>1 & recursive){
                                    
                                    pred <- matrix(pred,nrow=1)
                                    pred <- predictMS(pred,trainY,h-1,beta[,,ii],p)
                                    ## pred <- 
                                }
                                if(separate_lambdas)
                                {
                                    pred <- matrix(pred,ncol=1)  
                                    for(uu in 1:ncol(ZFull$Y)){

                                        ## if(loss=="L2"){ 
                                        MSFE[v - (T1 -h), ii,uu] <- .calc.loss(ZFull$Y[v+h-1,uu]-pred[uu,1],univ=TRUE,loss,delta)
                                        ## }else if(loss=="L1"){
                                        ##     MSFE[v - (T1 -h), ii,uu] <- abs(ZFull$Y[v+h-1,uu] - pred[uu,1])
                                        ## }else if(loss=="huber"){
                                        ##     MSFE[v - (T1 -h), ii,uu] <- huber_loss(ZFull$Y[v+h-1,uu] - pred[uu,1],delta)
                                        ## }
                                    }
                                }else{
                                    ## ## browser()
                                    ## if(ii==3){browser()}
                                    MSFE[v - (T1 -h), ii] <- .calc.loss(ZFull$Y[v+h-1,1:k1] - pred,univ=FALSE,loss,delta)
                                    ## if(loss=="L2"){ 
                                    ##     MSFE[v - (T1 -h), ii] <- norm2(ZFull$Y[v+h-1,1:k1] - pred)^2
                                    ## }else if(loss=="L1"){
                                    ##     MSFE[v - (T1 -h), ii] <- sum(abs(ZFull$Y[v+h-1,1:k1] - pred))
                                    ## }else if(loss=="huber"){
                                    ##     MSFE[v - (T1 -h), ii] <- sum(huber_loss(ZFull$Y[v+h-1,1:k1] - pred,delta))
                                    ##     }
                                    
                                }


                            }else{
                                if(VARX){          

                                    eZ <- VARXCons(matrix(Y[(v-p):(v),1:k1],ncol=k1),Y[(v-p):(v),(ncol(Y)-m+1):(ncol(Y))],k1,p
                                                  ,m,s,contemp=contemp)
                                    
                                    pred <- beta[,,ii] %*% eZ


                                }else{

                                    eZ<- VARXCons(Y[(v-p):(v),1:k1],matrix(0,nrow=length((v-p):v)),k1,p,0,0)

                                    pred <- beta[,,ii] %*% eZ

                                }
                                MSFE[v - (T1 - h), ii] <- .calc.loss(ZFull$Y[v+h-1,1:k1] - pred,univ=FALSE,loss,delta)
                                ## MSFE[v - (T1 - h), ii] <- norm2(ZFull$Y[v+h-1,1:k1] - pred)^2     


                            }

                        }

                    }
                }else{
                                        # If alpha and lambda are jointly fit, calculate MSFE for each alpha, lambda combination
                    for (ii in 1:gran2) {
                        for(j in 1:length(alpha)){
                            if (RVAR) {

                                        # Relaxed Least Squares (intercept ignored)
                                beta[,,(ii-1)*nalpha+j] <- RelaxedLS(cbind(t(trainZ),trainY),beta[,,(ii-1)*nalpha+j],k,p,k1,s+s1)
                            }

                            if(MN){
                                if(!intercept){
                                    pred <- beta[,2:dim(beta)[2],(ii-1)*nalpha+j] %*% eZ[2:length(eZ)]
                                }else{
                                    pred <- beta[,,(ii-1)*nalpha+j] %*% eZ
                                }
                                if(h>1 & recursive){
                                    pred <- matrix(pred,nrow=1)
                                    pred <- predictMS(pred,trainY,h-1,beta[,2:dim(beta)[2],(ii-1)*nalpha+j],p,MN)
                                }

                                
                                MSFE[v - (T1 - h), (ii-1)*nalpha+j] <- .calc.loss(ZFull$Y[v+h-1,1:k1] - beta[,2:dim(beta)[2],(ii-1)*nalpha+j] %*% eZ[2:length(eZ)],univ=FALSE,loss,delta)
                                ## MSFE[v - (T1 - h), (ii-1)*nalpha+j] <- norm2(ZFull$Y[v+h-1,1:k1] - beta[,2:dim(beta)[2],(ii-1)*nalpha+j] %*% eZ[2:length(eZ)])^2
                                
                            }else{
                                pred <- beta[,,(ii-1)*nalpha+j] %*% eZ

                                if(h>1 & recursive){
                                    pred <- matrix(pred,nrow=1)
                                    
                                    pred <- predictMS(pred,trainY,h-1,beta[,,(ii-1)*nalpha+j],p)
                                }
                                
                                if(object@crossval=="Rolling"){
                                    ## pred <- beta[,,ii] %*% eZ

                                    ## if(h>1 & recursive){
                                    ##     pred <- matrix(pred,nrow=1)
                                    ##     pred <- predictMS(pred,trainY,h-1,beta[,,ii],p)
                                    ## }
                                    
                                    ## if(separate_lambdas)
                                    ## {
                                    
                                    ## for(uu in 1:k){
                                    ##     MSFE[v - (T1 -h), ii,uu] <- (ZFull$Y[v+h-1,uu] - pred[uu,1])^2
                                    ##     }
                                    ## }else{
                                    ## MSFE[v - (T1 -h), ii] <- norm2(ZFull$Y[v+h-1,1:k1] - pred)^2

                                    MSFE[v - (T1 - h), (ii-1)*nalpha+j] <- .calc.loss(ZFull$Y[v,1:k1] - pred,univ=FALSE,loss,delta)
                                    ## MSFE[v - (T1 - h), (ii-1)*nalpha+j] <- norm2(ZFull$Y[v,1:k1] - pred)^2

                                    ## browser()
                                    ## beta2[,,(ii-1)*nalpha+j] <<-
                                    ## if((ii-1)*nalpha+j==51){browser()}
                                    ## }

                                }else{
                                    if(VARX){          
                                        eZ<- VARXCons(matrix(Y[(v-p):(v),1:k1],ncol=k1),Y[(v-p):(v),(ncol(Y)-m+1):(ncol(Y))],k1,p,m,s,contemp=contemp)
                                    }else{
                                        eZ<- VARXCons(Y[(v-p):(v),1:k1],matrix(0,nrow=length((v-p):v)),k1,p,0,0)
                                    }
                                    MSFE[v - (T1 - h), (ii-1)*nalpha+j] <- .calc.loss(ZFull$Y[v,1:k1] - beta[,,(ii-1)*nalpha+j] %*% eZ,univ=FALSE,loss,delta)
                                    ## MSFE[v - (T1 - h), (ii-1)*nalpha+j] <- norm2(ZFull$Y[v,1:k1] - beta[,,(ii-1)*nalpha+j] %*% eZ)^2
                                    
                                    
                                    


                                }
                                
                            }


                        }
                    }
                }

            }else{
                for (ii in 1:ncol(MSFE)) {
                    MSFE[v - (T1 - h), ii] <- .calc.loss(Y[v+h-1,1:k1] - beta[,,ii],univ=FALSE,loss,delta)

                    ## MSFE[v - (T1 - h), ii] <- norm2(Y[v+h-1,1:k1] - beta[,,ii])^2
                }
            }


            if(verbose){
                setTxtProgressBar(pb, v)}
            ## browser()
            ## beta2[,,v-(T1-h)] <<- beta[,,13]
        }

                                        # Sort out indexing for 2-d gridsearch
        ## browser()
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
            optind <- indopt
        }

        ## browser()

                                        # one standard error correction     
        if(ONESE & !dual &!separate_lambdas){
            MSFE2 <- MSFE 
            G2 <- colMeans(na.omit(MSFE2))
            G3 <- sd(na.omit(MSFE2))/sqrt(nrow(na.omit(MSFE2)))
            optind <- min(which(G2<(min(G2)+G3)))
            gamopt <- gamm[optind]
        }else{

            if(group!="Tapered" & !dual){
                                        # in rare cases in which MSFE is equal, the smaller penalty parameter is chosen.
                                        # This prevents extremely sparse solutions

                if(separate_lambdas){
                    ## browser()
                    if(ONESE){
                        ## MSFEs_FULL <- MSFE
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
                }else{

                    optind <- max(which(colMeans(na.omit(MSFE))==min(colMeans(na.omit(MSFE)))))
                    gamopt <- gamm[optind]
                }
            }else if(dual){
                if(!ONESE){

                    optind <- max(which(colMeans(na.omit(MSFE))==min(colMeans(na.omit(MSFE)))))
                    ## browser()
                    inds <- findind(optind,gamm[,1],alpha)                                   
                }else{
                    ## browser()
                    G2 <- colMeans(na.omit(MSFE))

                    G3 <- sd(na.omit(MSFE))/sqrt(nrow(na.omit(MSFE)))

                    optind <- min(which(G2<(min(G2)+G3)))
                    inds <- findind(optind,gamm[,1],alpha)
                    
                }
                gamopt <- gamm[inds[1],inds[2]]
                gamm <- gamm[,inds[2]]
                alphaopt <- alpha[inds[2]]
                optind <- inds

            }
        }
        if(!dual){
            alphaopt <- alpha
        }

        
        if(VARX){

            if(rolling_oos){


                OOSEval <- .BigVAREVALX_rolling(ZFull,MSFE,gamm,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m,contemp,alphaopt,C,intercept,tol,window.size,separate_lambdas,loss=loss,delta=delta)


            }else{
                                        # Out of sample forecast evaluation: VARX
                OOSEval <- .BigVAREVALX(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m,contemp,alphaopt,C,intercept,tol,window.size,separate_lambdas,loss=loss,delta=delta)
            }
        }else{
            if(rolling_oos){
                ## browser()
                OOSEval <- .BigVAREVAL_rolling(ZFull,MSFE,gamm,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,alphaopt,recursive,C,intercept,tol,window.size,separate_lambdas,loss=loss,delta=delta,ONESE=ONESE)
            }else{                        # Out of sample evaluation for VAR    
                OOSEval <- .BigVAREVAL(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,alphaopt,recursive,C,intercept,tol,window.size,separate_lambdas,loss=loss,delta=delta)
            }
        }
        MSFEOOSAgg <- na.omit(OOSEval$MSFE)
        betaPred <- OOSEval$betaPred
        betaArray <- OOSEval$betaArray

        preds <- OOSEval$predictions
        Y <- object@Data # preserve Y for BigVAR.results object
        if(VARX){

                                        # Construct lag matrix for OOS predictions
                                        # need to be careful with OOS predictions if contemporaneous exogenous covariates are included
            if(contemp){OOS <- FALSE}else{OOS <- TRUE}
                                        # special case of transfer function                                      
            if(!object@tf){
                Zvals <- VARXCons(matrix(Y[,1:k1],ncol=k1),matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m),k1,p,m,s,oos=OOS,contemp=contemp)
            }else{

                Zvals <- VARXCons(matrix(0,ncol=1,nrow=nrow(Y)),matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m),0,0,m,s,oos=FALSE,contemp=contemp)
            }

        }else{
            m <- 0;s <- 0
            Zvals <- VARXCons(matrix(Y[,1:k1],ncol=k1),matrix(0,nrow=nrow(Y)),k1,p,m,s,oos=TRUE)
        }

        Zvals <- matrix(Zvals[,ncol(Zvals)],ncol=1)
        if(ncol(Y)==1| k1==1){betaPred <- matrix(betaPred,nrow=1)}
        lagmatrix <- rbind(rep(1,ncol(ZFull$Z)),ZFull$Z)
        
        fitted <- t(betaPred%*%lagmatrix)
                                        #Residuals
        resids <- ((ZFull$Y)-fitted)

        ## lagmatrix <- ZFull$Z
        
        MSFEOOS<-mean(na.omit(MSFEOOSAgg))

        seoos <- sd(na.omit(MSFEOOSAgg))/sqrt(length(na.omit(MSFEOOSAgg)))

        if(!VARX){k1 <- k}


                                        # naive benchmarks     
        meanbench <- .evalMean(ZFull$Y[,1:k1],T2,T,h=h,loss=loss,delta=delta)
        RWbench <- .evalRW(ZFull$Y[,1:k1],T2,T,h=h,loss=loss,delta=delta)


        if(object@ic==FALSE|object@tf){

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

                                        # Information Criterion Benchmarks    

        }else{

            if(!VARX){

                X <- matrix(0,nrow=nrow(Y),ncol=k)

                AICbench1 <- VARXForecastEval(matrix(ZFull$Y,ncol=k),X,p,0,T2,T,"AIC",h,loss=loss,delta=delta)

                AICbench <- list()

                AICbench$Mean <- mean(AICbench1$MSFE)

                AICbench$SD <- sd(AICbench1$MSFE)/sqrt(length(AICbench1$MSFE))
                AICbench$preds <- AICbench1$pred
                AICbench$pvec <- AICbench1$p
                AICbench$svec <- AICbench1$s

                BICbench1 <- VARXForecastEval(matrix(ZFull$Y,ncol=k),X,p,0,T2,T,"BIC",h,loss=loss,delta=delta)
                
                BICbench <- list()

                BICbench$Mean <- mean(BICbench1$MSFE)

                BICbench$SD <- sd(BICbench1$MSFE)/sqrt(length(BICbench1$MSFE))
                BICbench$preds <- BICbench1$pred
                BICbench$pvec <- BICbench1$p
                BICbench$svec <- BICbench1$s

            }else{

                offset <- max(c(p,s))

                X <- matrix(Y[(offset+1):nrow(Y),(k1+1):ncol(Y)],ncol=m)

                AICbench1 <- VARXForecastEval(matrix(ZFull$Y[,1:k1],ncol=k1),as.matrix(X),p,s,T2,T,"AIC",h=h,loss=loss,delta=delta)

                AICbench <- list()

                AICbench$Mean <- mean(AICbench1$MSFE)

                AICbench$SD <- sd(AICbench1$MSFE)/sqrt(length(AICbench1$MSFE))
                AICbench$preds <- AICbench1$pred
                AICbench$pvec <- AICbench1$p
                AICbench$svec <- AICbench1$s

                BICbench1 <- VARXForecastEval(matrix(ZFull$Y[,1:k1],ncol=k1),X,p,s,T2,T,"BIC",h=h,loss=loss,delta=delta)

                BICbench <- list()

                BICbench$Mean <- mean(BICbench1$MSFE)

                BICbench$SD <- sd(BICbench1$MSFE)/sqrt(length(BICbench1$MSFE))  
                BICbench$preds <- BICbench1$pred
                BICbench$pvec <- BICbench1$p
                BICbench$svec <- BICbench1$s

            }

        }

        if(!VARX){contemp=FALSE}

        if(VARX & contemp){
            VARXL <- list(k=k1,s=s,contemp=contemp)
        }else if(VARX & !contemp){

            VARXL <- list(k=k1,s=s,contemp=FALSE)
            
        }else{
            VARXL <- list()
        }

        if(separate_lambdas){
            tmean <- t(apply(MSFE,3,colMeans))
        }
        
        ## isMSFE <- ifelse(separate_lambdas,apply(tmean,1,mean),colMeans(MSFE))
        if(separate_lambdas){
            ## isMSFE <- tmean
            ## isMSFE <- as.matrix(MSFE[,,1])
            ## isMS
            isMSFE <- MSFE
        }else{
            
            isMSFE <- array(MSFE,dim=c(nrow(MSFE),ncol(MSFE),1))
        }

        ## ## browser()
        ## if(group%in%c("Basic","BasicEN","Lag","HVARC","HVAROO","HVARELEM")){
        sparse_count <- function(x){
            x_ss <- x[,2:ncol(x)]
            sc <- length(x_ss[x_ss!=0])/length(x)
            sc
            
        }

        sc <- mean(apply(betaArray,3,sparse_count))
        ## }else{
        ##     ## browser()
        ##     ## jj
        ##     sparse_count <- function(x){
        ##          x_ss <- x[,2:ncol(x)]
        ##          sc <- length(x_ss[x_ss!=0])/length(x)
        ##          sc
        
        ## }
        ## jj=groupfunVARXLG(p,k,k1,s)
        ## jj
        ## b1=betaArray[,2:ncol(betaArray),1]
        ## b1[,unlist(jj)+1]
        ##     ## betaArray[,unlist(jj)+1,]
        ##     sparse_count <- mean(apply(betaArray,3,sparse_count))
        
        ## }

        
        ## browser()
                                        # Create a new BigVAR.Results Object
        results <- new("BigVAR.results",InSampMSFE=isMSFE,InSampSD=apply(MSFE,2,sd)/sqrt(nrow(MSFE)),LambdaGrid=as.matrix(gamm),index=optind,OptimalLambda=gamopt,OOSMSFE=as.matrix(MSFEOOSAgg),seoosmsfe=seoos,MeanMSFE=meanbench$Mean,AICMSFE=AICbench$Mean,AICpvec=AICbench$pvec,AICsvec=AICbench$svec,AICPreds=AICbench$preds,BICpvec=BICbench$pvec,BICsvec=BICbench$svec,BICPreds=BICbench$preds,RWMSFE=RWbench$Mean,RWPreds=RWbench$preds,MeanSD=meanbench$SD,MeanPreds=meanbench$preds,AICSD=AICbench$SD,BICMSFE=BICbench$Mean,BICSD=BICbench$SD,RWSD=RWbench$SD,sparse_count=sc,betaPred=betaPred,Zvals=Zvals,resids=resids,VARXI=VARX,preds=preds,alpha=alphaopt,fitted=fitted,lagmatrix=lagmatrix,betaArray=betaArray,dual=dual,contemp=contemp,object)
        ## results <- new("BigVAR.results",InSampMSFE=isMSFE,InSampSD=apply(MSFE,2,sd)/sqrt(nrow(MSFE)),LambdaGrid=as.matrix(gamm),index=optind,OptimalLambda=gamopt,OOSMSFE=as.matrix(MSFEOOSAgg),seoosmsfe=seoos,MeanMSFE=meanbench$Mean,AICMSFE=AICbench$Mean,AICpvec=AICbench$pvec,AICsvec=AICbench$svec,AICPreds=AICbench$preds,BICpvec=BICbench$pvec,BICsvec=BICbench$svec,BICPreds=BICbench$preds,RWMSFE=RWbench$Mean,RWPreds=RWbench$preds,MeanSD=meanbench$SD,MeanPreds=meanbench$preds,AICSD=AICbench$SD,BICMSFE=BICbench$Mean,BICSD=BICbench$SD,RWSD=RWbench$SD,Data=object@Data,lagmax=object@lagmax,Structure=object@Structure,Minnesota=object@Minnesota,Relaxed=object@Relaxed,Granularity=object@Granularity,horizon=object@horizon,betaPred=betaPred,Zvals=Zvals,resids=resids,VARXI=VARX,VARX=VARXL,preds=preds,T1=T1,T2=T2,dual=dual,alpha=alphaopt,crossval=object@crossval,ownlambdas=object@ownlambdas,tf=object@tf,recursive=recursive,constvec=C,intercept=intercept,tol=tol,fitted=fitted,lagmatrix=lagmatrix,betaArray=betaArray,window.size=object@window.size,separate_lambdas=object@separate_lambdas,sparse_count=sparse_count,object) 
        
        return(results)
    }
)


#' BigVAR Estimation
#' @description
#' Fit a BigVAR object with a structured penalty (VARX-L or HVAR).
#' @usage BigVAR.est(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details Fits HVAR or VARX-L model on a BigVAR object.  Does not perform cross-validation.  This method allows the user to construct their own penalty parameter selection procedure.
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
#' Model1=constructModel(Y,p=4,struct="Basic",gran=c(50,10))
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
        tol=object@tol
        if(object@ownlambdas==TRUE){
            gamm=object@Granularity
            gran2 <- length(gamm)

        }      
        separate_lambdas <- object@separate_lambdas
        C <- object@constvec
        group <- object@Structure
        Y <- object@Data
        k <- ncol(Y)
        VARX <- object@VARX
        intercept <- object@intercept
        if(length(object@alpha)==0){

            if(length(VARX)>0){    
                alpha <- 1/(VARX$k+1)
                
            }else{

                alpha <- 1/(k+1)
                
            }
        }else{

            alpha <- object@alpha

        }

        
        if(length(alpha)>1 & group%in%c("SparseLag","SparseOO"))
        {
            dual <- TRUE

        }else{

            dual <- FALSE
        }

        RVAR <- object@Relaxed
        group <- object@Structure
        MN <- object@Minnesota
        jj <- 0
        if(!"matrix"%in%class(Y)){Y=matrix(Y,ncol=1)}
        s <- ifelse(length(object@VARX)!=0,object@VARX$s,0)
        T <- nrow(Y)-max(p,s)
        VARX <- object@VARX

        if(!object@ownlambdas){
            gran2 <- object@Granularity[2]
            gran1 <- object@Granularity[1]}

        if(length(VARX)!=0){

            VARX <- TRUE
            k1 <- object@VARX$k
            s <- object@VARX$s

            if(!is.null(object@VARX$contemp)){

                contemp <- TRUE
                s1 <- 1

            }else{

                contemp <- FALSE
                s1 <- 0
            }

            m <- k-k1
            Y1 <- matrix(Y[,1:k1],ncol=k1)
            X <- matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m)

            if(!object@tf){
                
                trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)

            }else{

                trainZ <- VARXCons(matrix(0,ncol=1,nrow=nrow(X)),matrix(X,ncol=m),k=0,p=0,m=m,s=s,contemp=contemp,oos=FALSE)

            }

            
            trainZ <- trainZ[2:nrow(trainZ),]

            trainY <- matrix(Y[(max(c(p,s))+1):nrow(Y),1:k1],ncol=k1)



            if(group=="Lag"|group=="SparseLag"){


                jj=groupfunVARXLG(p,k,k1,s)

                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2*length(alpha))


            }

            if(group=="OwnOther"|group=="SparseOO"){

                jj=diaggroupfunVARXLG(p,k,k1,s)
                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2)


            }


            if(group=="BGR"){
                Grid <- seq(1,5,by=.025)
                grid <- Grid*sqrt(k*p)
                MSFE <- matrix(0, nrow = 1, ncol = length(grid))
            }

            
            if(!object@ownlambdas){

                
                if(dual){
                                        # Constructs penalty grid if both alpha and lambda are selected
                    gamm <- .LambdaGridXDual(gran1, gran2, jj, trainY, trainZ,group,p,k1,s,m,k,MN,alpha,C,intercept,tol)

                }else{

                                        # Penalty parameter grid for just lambda
                    ## gamm <- .LambdaGridX(gran1, gran2, jj, as.matrix(trainY), trainZ,group,p,k1,s+s1,m,k,MN,alpha,C,intercept,tol)
                    if(group!="BGR"){
                        gamm <- .LambdaGridX(gran1, gran2, jj, trainY, trainZ,group,p,k1,s+s1,m,k,MN,alpha,C,intercept,tol)
                    }else{
                        gamm <- seq(1,5,by=.025)
                        gamm <- gamm*sqrt(k*p)
                        
                    }
                }
            }else{

                gamm <- object@Granularity

            }

                                        # Initial Coefficient Matrix         
            beta=array(0,dim=c(k1,k1*p+(k-k1)*s+1,gran2*length(alpha)))

                                        # Initialize groups, active sets, power method calculations, etc
            if (group == "Lag") {

                jj=groupfunVARX(p,k,k1,s+s1)

                jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)

                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2)

            }

            if (group == "SparseLag") {

                jj <- groupfunVARX(p, k,k1,s+s1)

                q1a <- list()


                for (i in 1:(p+s+s1)) {

                    q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)

                }
                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2*length(alpha))
                
            }
            if (group == "OwnOther") {

                kk <- diaggroupfunVARX(p,k,k1,s+s1)

                activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                                 gran2)

            }
            if (group == "SparseOO") {


                kk <- diaggroupfunVARX(p,k,k1,s+s1)

                activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                                 gran2*length(alpha))

                q1a <- list()

                for (i in 1:(2*p+s+s1)) {

                    q1a[[i]] <- matrix(runif(length(jj[[i]]), -1, 1), ncol = 1)

                }
                
            }
        }else{

            s=p

            if (group=="Lag"|group=="SparseLag")

            {

                jj=.groupfun(p,k)

            }else{

                jj <- .lfunction3(p,k)

            }


            Z1 <- VARXCons(Y,matrix(0,nrow=nrow(Y)),k,p,0,0) 


            trainZ <- Z1[2:nrow(Z1),]   


            trainY <- matrix(Y[(p+1):nrow(Y),],ncol=k)          


            if(!object@ownlambdas){

                if(dual){
                    
                    gamm <- .LambdaGridEDual(gran1, gran2, jj, trainY, trainZ,group,p,k,MN,alpha,C,intercept,tol)
                    

                }else{
                    verbose=object@verbose
                    if(group!="BGR"){
                        
                        
                        gamm <- .LambdaGridE(gran1, gran2, jj, trainY, trainZ,group,p,k,MN,alpha,C,intercept,tol,separate_lambdas = separate_lambdas,verbose)
                    }else{
                        gamm <- seq(1,5,by=.025)
                        gamm <- gamm*sqrt(k*p)
                        
                    }

                }
            }else{

                gamm <- object@Granularity
                
            }

            VARX <- FALSE
            k1 <- k
            s <- 0   

            if (group == "Lag") {

                jj <- .groupfuncpp(p, k)

                jjcomp <- .groupfuncomp(p,k)

                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2)

                k1 <- k

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

            if (group == "OwnOther") {

                kk <- .lfunction3cpp(p, k)

                activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                                 gran2)

            }

            if (group == "SparseOO") {

                kk <- .lfunction3cpp(p, k)
                jjcomp <- .lfunctioncomp(p,k)
                jj <- .lfunction3(p,k)
                activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                                 gran2)
                q1a <- list()

                for (i in 1:(2*p)) {

                    q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)
                }
                
            }
            beta <- array(0,dim=c(k,k*p+1,gran2*length(alpha)))

        }           


        
        if(group=="Tapered")

        {
            palpha <- seq(0,1,length=10)
            palpha <- rev(palpha)
            gran2 <- length(gamm)*length(palpha)
            beta <- array(0,dim=c(k,k*p+1,gran2))
        }


        ## if (group == "BGR") {
        
        ##     trainZ <- rbind(1,trainZ)
        ##     beta <- BGRGridSearch(trainY,trainZ,p,gamm,as.numeric(MN))
        ##  }
        
        needed.objs <- c('group','beta','trainZ','trainY','gamm','tol','p','m','k1','s',
                         's1','m','MN','C','intercept','separate_lambdas','dual','activeset','alpha','jj','jjcomp','kk','palpha')
        if(!group%in%c("SparseLag","SparseOO")){
            q1a <- NULL
        }
        ## objs <- sapply(needed.objs,exists)
        objs <- setdiff(needed.objs,ls())
        if(length(objs)>0){
            
            for(i in 1:length(objs)){               
                assign(objs[i],NULL)
            }
        }
        temp <- .BigVAR.fit(group,beta,trainZ,trainY,gamm,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk,palpha)
        beta <- temp$beta
        activeset <- temp$activeset
        q1a <- temp$q1a
        
        return(list(B=beta,lambdas=gamm))

    }

)


## New object class: bigVAR results, inherits class bigVAR, prints results from cv.bigVAR

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

#'
#' @note One can also access any object of class BigVAR from BigVAR.results
#' @name BigVAR.results 
#' @rdname BigVAR.results
#' @aliases BigVAR.results-class
#' @exportClass BigVAR.results
#' @author Will Nicholson
#' @export
setClass("BigVAR.results",
         representation(InSampMSFE="array",InSampSD="numeric",LambdaGrid="matrix",index="numeric",OptimalLambda="numeric",OOSMSFE="matrix",seoosmsfe="numeric",MeanMSFE="numeric",AICMSFE="numeric",AICPreds="matrix",BICMSFE="numeric",BICpvec="numeric",BICsvec="numeric",AICpvec="numeric",AICsvec="numeric",BICSD="numeric",BICPreds="matrix",RWMSFE="numeric",RWPreds="matrix",MeanSD="numeric",MeanPreds="matrix",AICSD="numeric",RWSD="numeric",betaPred="matrix",Zvals="matrix",VARXI="logical",resids="matrix",preds="matrix",dual="logical",contemp="logical",fitted="matrix",lagmatrix="matrix",betaArray="array",sparse_count="numeric"),
         contains="BigVAR"
         )


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
setMethod(f="plot",signature="BigVAR.results",
          definition= function(x,y=NULL,...)
          {

              if(!x@separate_lambdas){
                  plot(x@LambdaGrid,colMeans(x@InSampMSFE[,,1]),type="o",xlab="Value of Lambda",ylab="MSFE",log="x")
              }else{
                  k <- ncol(x@Data)
                  par(mfrow=c(k,1))
                  for(i in 1:k){
                      plot(x@LambdaGrid[,i],colMeans(x@InSampMSFE[,,i]),type="o",xlab="Value of Lambda",ylab="MSFE",log="x")   
                      abline(v=x@OptimalLambda[i],col="green")
                  }
              }
          }
          )

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
setMethod("show","BigVAR.results",
          function(object)
          {
              cat("*** BIGVAR MODEL Results *** \n")
              cat("Structure\n") ;print(object@Structure)
              if(object@Relaxed==TRUE){
                  cat("Relaxed VAR \n") ;print(object@Relaxed)}
              cat("Loss \n") ;print(object@loss)
              cat("Forecast Horizon \n") ;print(object@horizon)
              cat("Minnesota VAR\n") ;print(object@Minnesota)
              

              if(object@VARXI){
                  cat("VARX Specs \n") ;print(object@VARX)}
              cat("Maximum Lag Order \n") ;print(object@lagmax)
              cat("Optimal Lambda \n"); print(signif(object@OptimalLambda,digits=4))
              if(object@dual){

                  cat("Optimal Alpha \n"); print(signif(object@alpha,digits=2))
                  
              }        
              cat("Grid Depth \n") ;print(object@Granularity[1])
              cat("Index of Optimal Lambda \n");print(object@index)
              cat("Fraction of active coefficients \n");print(signif(object@sparse_count,digits=4))
              if(!object@separate_lambdas){
                  cat("In-Sample Loss\n");print(signif(mean(object@InSampMSFE[,object@index,]),digits=3))
              }else{
                  ## browser()

                  cat("In-Sample Loss\n");print(signif(apply(object@InSampMSFE[,object@index,],2,mean),digits=3))
              }
              cat("BigVAR Out of Sample Loss\n");print(signif(mean(object@OOSMSFE),digits=3))
              cat("*** Benchmark Results *** \n")
              cat("Conditional Mean Out of Sample Loss\n");print(signif(object@MeanMSFE,digits=3))
              cat("AIC Out of Sample Loss\n");print(signif(object@AICMSFE,digits=3))
              cat("BIC Out of Sample Loss\n");print(signif(object@BICMSFE,digits=3))
              cat("RW Out of Sample Loss\n");print(signif(object@RWMSFE,digits=3))
          }
          )



#' Forecast using a BigVAR.results object
#' 
#' @usage predict(object,...)
#' @param object BigVAR.results object from \code{cv.BigVAR}
#' @param ... additional arguments affecting the predictions produced (e.g. \code{n.ahead})
#' @details Provides \code{n.ahead} step forecasts using the model produced by cv.BigVAR. 
#' @seealso \code{\link{cv.BigVAR}} 
#' @name predict
#' @aliases predict,BigVAR.results-method
#' @docType methods
#' @method predict method
#' @rdname predict-methods-BigVAR.results
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' Model1=constructModel(Y,p=4,struct="Basic",gran=c(50,10),verbose=FALSE)
#' results=cv.BigVAR(Model1)
#' predict(results,n.ahead=1)
#' @export
setMethod("predict","BigVAR.results",
          function(object,n.ahead,newxreg=NULL,...)
          {
                                        # MN option removes intercept
              MN <- object@Minnesota
              eZ <- object@Zvals
              ## if(MN & !object@intercept){
              ##     eZ <- eZ[2:nrow(eZ),,drop=F]
              ##     }
              betaPred <- object@betaPred
              Y <- object@Data
              k <-object@VARX$k
              m <- ncol(object@Data)-k
              p <- object@lagmax
              s <- object@VARX$s
              VARX <- object@VARXI
              contemp <- object@contemp
              s1 <- 0
              ## if(VARX){
              fcst <- matrix(betaPred%*%eZ,ncol=1)

              if(n.ahead==1)
              {
                  return(fcst)
              }else{
                  if(!VARX){
                                        # iterative multistep forecasts
                      fcst <- matrix(predictMS(matrix(fcst,nrow=1),Y[(nrow(Y)-p+1):nrow(Y),],n.ahead-1,betaPred,p,MN),ncol=1)

                  }else{

                                        # experimental VARX multistep forecasts
                      if(is.null(newxreg))
                      {
                          stop("Need new data for multi-step VARX forecasts.  Re-run with new data in newxreg")

                      }else{
                          if(nrow(newxreg)<n.ahead-1){
                              stop(paste("Need at least ",n.ahead-1,"rows of new data"))
                          }
                          
                          C <- max(p,s)
                          if(contemp){C <- C+3}
                          ## browser()
                          fcst <- matrix(predictMSX(matrix(fcst,nrow=1),as.matrix(Y[(nrow(Y)-C+1):nrow(Y),1:(k)]),n.ahead-1,betaPred,p,newxreg,matrix(Y[(nrow(Y)-C+1):nrow(Y),(ncol(Y)-m+1):ncol(Y)],ncol=m),m,s,1,MN,contemp),ncol=1)
                      }
                  }

              }

              return(fcst)
          }
          )

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
#' Model1 <- constructModel(Y,p=4,struct="Basic",gran=c(50,10),verbose=FALSE)
#' SparsityPlot.BigVAR.results(cv.BigVAR(Model1))
#' @export
#' @importFrom lattice levelplot
#' @importFrom lattice panel.abline
#' @importFrom lattice panel.levelplot
#' @importFrom grDevices colorRampPalette

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
            B <- B[,2:ncol(B)]
        }
        k <- nrow(B)
        p <- object@lagmax
        s1 <- 0
        if(length(object@VARX!=0)){
            s <- object@VARX$s
            m <- ncol(object@Data)-object@VARX$k
            contemp <- object@VARX$contemp
            if(!is.null(contemp)){
                if(contemp){

                    s1 <- 1
                }            
            }else{
                s1 <- 0

            }
        }
        else{
            m <- 0;s <- 0
        }

        s <- s+s1

        text <- c()

        for (i in 1:p) {

            text1 <- as.expression(bquote(bold(Phi)^(.(i))))

            text <- append(text, text1)

        }

        if(m>0){

            for (i in (p+1):(p+s+1)) {

                text1 <- as.expression(bquote(bold(beta)^(.(i-p-s1))))

                text <- append(text, text1)

            }

        }

        f <- function(m) t(m)[, nrow(m):1]

        rgb.palette <- colorRampPalette(c("white", "blue" ),space = "Lab")

        at <- seq(k/2 + 0.5, p * (k)+ 0.5, by = k)

        if(m>0){

            at2 <- seq(p*k+m/2+.5,p*k+s*m+.5,by=m)

        }else{
            at2=c()

        }

        at <- c(at,at2)

        se2 = seq(1.75, by = k, length = k)

        L2 <- levelplot(as.matrix(f(abs(B))), col.regions = rgb.palette, colorkey = NULL, 
                        xlab = NULL, ylab = NULL, main = list(label = "Sparsity Pattern Generated by BigVAR", 
                                                              cex = 1), panel = function(...) {
                                                                  panel.levelplot(...)
                                                                  panel.abline(a = NULL, b = 1, h = seq(1.5, m*s+p* k + 
                                                                                                             0.5, by = 1), v = seq(1.5, by = 1, length = p*k+m*s),lwd=.5)
                                                                  bl1 <- seq(k + 0.5,p*k + 0.5, by = k)
                                                                  b23 <- seq(p*k + 0.5, p * 
                                                                                        k + 0.5+s*m, by = m)
                                                                  b1 <- c(bl1,b23)
                                                                  panel.abline(a = NULL, b = 1, v = p*k+.5, lwd = 3)
                                                                  
                                                                  panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
                                                                  

                                                              }, scales = list(x = list(alternating = 1, labels = text, 
                                                                                        cex = 1, at = at, tck = c(0, 0)), y = list(alternating = 0,tck
                                                                                                                                   = c(0, 0))))

        return(L2)

    }

)



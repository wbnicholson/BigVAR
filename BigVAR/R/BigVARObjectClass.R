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
    if(object@lagmax==0& object@Structure!="Basic"){
        msg <- c("Only Basic VARX-L supports a transfer function")
        errors <- c(errors,msg)
    }
    structures=c("Basic","Lag","SparseLag","OwnOther","SparseOO","HVARC","HVAROO","HVARELEM","Tapered","EFX")
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
    structure2 <- c("Basic","Lag","HVARC")
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
        msg <- c("EFX does not support contemporaneous dependence")
        errors <- c(errors,msg)

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
    
    if(length(errors)==0) TRUE else errors
    
}

#' BigVAR Object Class
#'
#' An object class to be used with cv.BigVAR
#' 
#' @slot Data a \eqn{T \times k} multivariate time Series
#' @slot lagmax Maximal lag order for modeled series
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
#' @slot alpha Grid of candidate alpha values (applies only to Sparse VARX-L models)
#' @slot recursive Indicator as to whether recursive multi-step forecasts are used (applies only to multiple horizon VAR models)

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
        tf="logical",
        alpha="numeric",
        recursive="logical",
        dates="character"
        ),validity=check.BigVAR
    )


#' Construct an object of class BigVAR
#' @param Y \eqn{T x k} multivariate time series or Y \eqn{T x (k+m)} endogenous and exogenous series, respectively 
#' @param p Predetermined maximal lag order (for modeled series)
#' @param struct The choice of penalty structure (see details).
#' @param gran vector containing how deep to construct the penalty grid (parameter 1) and how many gridpoints to use (parameter 2)  If ownlambas is set to TRUE, gran denotes the user-supplied penalty parameters.
#' @param RVAR True or False: whether to refit based upon the support selected using the Relaxed-VAR procedure
#' @param h Desired forecast horizon
#' @param cv Cross-validation approach, either "Rolling" for rolling cross-validation or "LOO" for leave-one-out cross-validation.
#' @param MN Minnesota Prior Indicator
#' @param verbose Verbose output while estimating
#' @param IC True or False: whether to include AIC and BIC benchmarks
#' @param VARX List containing VARX model specifications. 
#' @param T1 Index of time series in which to start cross validation
#' @param T2  Index of times series in which to start forecast evaluation
#' @param ONESE True or False: whether to use the "One Standard Error Heuristic"
#' @param ownlambdas True or False: Indicator for user-supplied penalty parameters
#' @param alpha grid of candidate parameters for the alpha in the Sparse Lag and Sparse Own/Other VARX-L 
#' @param recursive True or False: Indicator as to whether iterative multi-step predictions are desired in the VAR context if the forecast horizon is greater than 1
#' 


#' @details The choices for "struct" are as follows
#' \itemize{
#' \item{  "Basic" (Basic VARX-L)}
#' \item{  "Lag" (Lag Group VARX-L)} 
#' \item{  "SparseLag" (Lag Sparse Group VARX-L)} 
#' \item{  "OwnOther" (Own/Other Group VARX-L) }
#' \item{  "SparseOO" (Own/Other Sparse Group VARX-L) }
#' \item{  "EFX" (Endogenous First VARX-L)}
#' \item{  "HVARC" (Componentwise HVAR) }
#' \item{  "HVAROO" (Own/Other HVAR) }
#' \item{  "HVARELEM" (Elementwise HVAR)}
#' \item{  "Tapered" (Lag weighted Lasso VAR)}
#' }
#'
#' VARX specifications consist of a list with entry k denoting the series that are to be modeled and entry s to denote the maximal lag order for exogenous series.
#'
#' The argument alpha is ignored unless the structure choice is "SparseLag" or "Lag."  By default "alpha" is set to \code{NULL} and will be initialized as 1/(k+1) in \code{cv.BigVAR} and \code{BigVAR.est}.  Any user supplied values must be between 0 and 1.  

#' @note The specifications "Basic", "Lag," "SparseLag," "SparseOO," and "OwnOther" can accommodate both VAR and VARX models.  EFX only applies to VARX models.  "HVARC," "HVAROO," "HVARELEM," and "Tapered" can only be used with VAR models.
#'
#' @seealso \code{\link{cv.BigVAR}},\code{\link{BigVAR.est}}
#' 
#' @references William B Nicholson, Jacob Bien, and David S Matteson. "High Dimensional Forecasting via Interpretable Vector Autoregression." arXiv preprint 1412.5250, 2016.
#' 
#' William B Nicholson, David S. Matteson, and Jacob Bien (2015), "VARX-L Structured regularization for large vector
#' autoregressions with exogenous variables," \url{http://www.wbnicholson.com/Nicholsonetal2015.pdf}.
#' 
#' @examples
#' library(BigVAR)
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
constructModel <- function(Y,p,struct,gran,RVAR=FALSE,h=1,cv="Rolling",MN=FALSE,verbose=TRUE,IC=TRUE,VARX=list(),T1=floor(nrow(Y)/3),T2=floor(2*nrow(Y)/3),ONESE=FALSE,ownlambdas=FALSE,alpha=as.double(NULL),recursive=FALSE)
{
    if(any(is.na(Y))){stop("Remove NA values before running ConstructModel")}      
    if(dim(Y)[2]>dim(Y)[1] & length(VARX)==0){warning("k is greater than T, is Y formatted correctly (k x T)?")}      
    if(p<0){stop("Maximal lag order must be at least 0")}
    if(p==0& struct!="Basic"){stop("Only Basic VARX-L supports a transfer function")}
    oldnames <- c("None","Diag","SparseDiag")
    if(struct%in%oldnames) stop("Naming Convention for these structures has changed. Use Basic, OwnOther, and SparseOO.")
    structures=c("Basic","Lag","SparseLag","OwnOther","SparseOO","HVARC","HVAROO","HVARELEM","Tapered","EFX")
    cond1=struct %in% structures
    if(!cond1){stop(cat("struct must be one of",structures))}
    if(h<1){stop("Forecast Horizon must be at least 1")}
    if(cv!="Rolling" & cv!="LOO"){stop("Cross-Validation type must be one of Rolling or LOO")}
    if(length(gran)!=2&ownlambdas==FALSE){stop("Granularity must have two parameters")}
    if(any(gran<=0)){stop("Granularity parameters must be positive")}

    structure2 <- c("Basic","Lag","HVARC")
    cond2=struct %in% structure2
    k1=0
    if(length(VARX)!=0){

        k1 <- VARX$k
        if(k1>ncol(Y)){stop("k is greater than the number of columns in Y")}
    }else{k=0}
    m <- ncol(Y)-k1
    nseries <- ncol(Y)-ifelse(m<ncol(Y),m,0)
    if(p==0){tf=TRUE
         }else{
             tf=FALSE
         }
    if(nseries==1 & cond2==FALSE ){stop("Univariate support is only available for Lasso, Lag Group, and Componentwise HVAR")}
    if(length(VARX)==0 & struct=="EFX"){stop("EFX is only supported in the VARX framework")}
    if(struct=="EFX" & !is.null(VARX$contemp)){stop("EFX does not support contemporaneous dependence")}
    structs=c("HVARC","HVAROO","HVARELEM")
    if(length(VARX)!=0& struct %in% structs){stop("EFX is the only nested model supported in the VARX framework")}
    if(T1>nrow(Y) | T2>nrow(Y) |T2<T1){stop("Training dates exceed series length")}

    if(is.list(VARX) & length(VARX)>0 & !(exists('k',where=VARX) & exists('s',where=VARX)))
        {

            stop("VARX Specifications entered incorrectly")

        }

    if(!is.null(alpha)){
        if(any(alpha<0) || any(alpha>1)){stop("alpha must be between zero and 1")}
    }
    if("xts"%in%class(Y)){
        ind <- as.character(index(Y))
        Y <- as.matrix(Y)
    }else{
        ind <- as.character(NULL)
    }
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
        dates=ind
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
              
              T1P <- ifelse(!is.null(object@dates),object@dates[object@T1],object@T1)
              T2P <- ifelse(!is.null(object@dates),object@dates[object@T2],object@T2)
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
#' @aliases plot-methods
#' @docType methods
#' @method plot
#' @rdname plot.BigVAR-methods
#' @export
#' @importFrom zoo plot.zoo
#' @importFrom zoo as.zoo
#' @importFrom zoo zoo
#' @importFrom zoo as.yearqtr
#' @importFrom zoo index
setMethod(f="plot",signature="BigVAR",
          definition= function(x,y=NULL,...)
          {

              T1P <- ifelse(length(x@dates)!=0,x@dates[x@T1],x@T1)

              T2P <- ifelse(length(x@dates)!=0,x@dates[x@T2],x@T2)

              g=ncol(x@Data)
              names <- ifelse(rep(!is.null(colnames(x@Data)),ncol(x@Data)),colnames(x@Data),as.character(1:g))
              ## browser()
              ## dates <- ifelse(!is.null(x@dates),x@dates,1:nrow(Y))
              if(length(x@dates)!=0){

                  dates <- as.yearqtr(x@dates)
              }else{

                  dates <- 1:nrow(x@Data)
              }

              Yzoo <- zoo(as.matrix(x@Data),order.by=dates)
              plot.zoo(Yzoo,plot.type="single",col=1:g)
              legend('topright',names,lty=1,col=1:g)

              ## index()
              abline(v=index(Yzoo[as.yearqtr(T1P)]))
              abline(v=index(Yzoo[as.yearqtr(T2P)]))
              
          }
          )

#' Cross Validation for BigVAR
#' 
#' @usage cv.BigVAR(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details Will perform cross validation to select penalty parameters over a training sample, then evaluate them over a test set.  Compares against sample mean, random walk, AIC, and BIC benchmarks.  Creates an object of class \code{BigVAR.results}
#' @return An object of class \code{BigVAR.results}.
#' @seealso \code{\link{constructModel}}, \code{\link{BigVAR.results}} 
#' @name cv.BigVAR
#' @aliases cv.BigVAR,BigVAR-method
#' @docType methods
#' @rdname cv.BigVAR-methods
#' @examples
#' data(Y)
#' Y=Y[1:100,]
#' # Fit a Basic VARX-L with rolling cross validation 
#' Model1=constructModel(Y,p=4,struct="Basic",gran=c(50,10))
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
        p <- object@lagmax
        s1 <- 0
        Y <- object@Data
        k <- ncol(Y)
        alpha <- object@alpha
        RVAR <- object@Relaxed
        group <- object@Structure
        cvtype <- object@crossval

        recursive <- object@recursive
        VARX <- object@VARX

        if(length(alpha)==0){

            if(length(VARX)>0){    
                alpha <- 1/(VARX$k+1)
                
            }else{

                alpha <- 1/(k+1)
                
            }
        }

        
        if(length(alpha)>1 & group%in%c("SparseLag","SparseOO"))
            {
                dual <- TRUE

            }else{

                dual <- FALSE
            }

        MN <- object@Minnesota
        h <- object@horizon
        jj <- 0

        if(class(Y)!="matrix"){Y=matrix(Y,ncol=1)}

        if(object@crossval=="Rolling"){
            T1 <- object@T1

        }else{

            T1 <- p+2    

        }
        T2 <- object@T2
        s <- ifelse(length(object@VARX)!=0,object@VARX$s,0)

        if(object@ownlambdas){
            gamm <- object@Granularity
            gran2 <- length(gamm)
        }     

        ONESE <- object@ONESE



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

            m <- k-k1
            Y1 <-matrix(Y[,1:k1],ncol=k1)
            X <- matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m)

            if(!object@tf){
                trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)

            }else{

                trainZ <- VARXCons(matrix(0,ncol=1,nrow=nrow(X)),matrix(X,ncol=m),k=0,p=0,m=m,s=s,contemp=contemp,oos=FALSE)

            }

            trainZ <- trainZ[2:nrow(trainZ),]

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
                    gamm <- .LambdaGridXDual(gran1, gran2, jj, trainY, trainZ,group,p,k1,s,m,k,MN,alpha)

                }else{

                                        # Penalty parameter grid for just lambda
                    gamm <- .LambdaGridX(gran1, gran2, jj, as.matrix(trainY[1:T2,]), trainZ[,1:T2],group,p,k1,s+s1,m,k,MN,alpha)
                }
            }

                                        # Coefficient matrices
            ## if(!dual){

            ##     beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2))

            ## }else{

            beta <- array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2*length(alpha)))
            
            ## }
                                        # Groupings in accordance with C++ indexing standards
            if (group == "Lag") {
                jj <- groupfunVARX(p,k,k1,s+s1)
                jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)
                activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                                 gran2)
            }
            if (group == "SparseLag") {
                
                jj <- groupfunVARX(p, k,k1,s+s1)
                q1a <- list()
                                        # Initializing warm start vectors for the power method
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

                    gamm <- .LambdaGridEDual(gran1, gran2, jj, trainY, trainZ,group,p,k,MN,alpha)

                }else{

                    gamm <- .LambdaGridE(gran1, gran2, jj, GY, GZ,group,p,k,MN,alpha)
                }
                
            }
            VARX <- FALSE
            k1 <- k
            s <- 0   
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
                                 gran2*length(alpha))
                
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
                                 gran2*length(alpha))
                q1a <- list()
                for (i in 1:(2*p)) {
                    q1a[[i]] <- matrix(runif(length(kk[[i]]), -1, 1), ncol = 1)
                }

                
            }
            beta <- array(0,dim=c(k,k*p+1,gran2))
        }           
        h <- object@horizon
        verbose <- object@verbose         
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

            }

        if(class(ZFull$Y)!="matrix" ){
            ZFull$Y <- matrix(ZFull$Y,ncol=1)
        }

        if(!dual){
            MSFE <- matrix(0, nrow = T2 - T1+1, ncol = gran2)
        }else{

            nalpha <- length(alpha)
            MSFE <- matrix(0, nrow = T2 - T1+1, ncol = gran2*nalpha)
            
        }
        if(verbose){
            pb <- txtProgressBar(min = T1, max = T2, style = 3)
            cat("Cross-Validation Stage:",group)}
        YT <- Y[1:T2,]

                                        # Start of penalty parameter selection     

        for (v in (T1-h+1):T2) {

            if(cvtype=="Rolling")

                {

                    if(v+h-1>T2){
                        break
                    }

                    if(h>1 & !recursive){
                        
                        trainY <- ZFull$Y[(h):(v-1), ]
                        
                        trainZ <- ZFull$Z[, 1:(v-h)]
                        
                    }else{

                        trainY <- ZFull$Y[1:(v-1), ]

                        
                        trainZ <- ZFull$Z[, 1:(v-1)]         
                    }
                }else{

                    if(VARX)

                        {

                            YT2 <- YT[-v,]

                            Y1 <- matrix(YT2[,1:k1],ncol=k1)

                            X <- matrix(YT2[,(ncol(YT2)-m+1):ncol(YT2)],ncol=m)

                            trainZ <- VARXCons(Y1,X,k1,p,m,s,contemp=contemp)

                            trainZ <- trainZ[2:nrow(trainZ),]

                            trainY <- matrix(YT2[(max(c(p,s))+1):nrow(YT2),1:k1],ncol=k1)

                        }else{


                            YT2 <- YT[-v,]

                            Z1 <- VARXCons(YT2,matrix(0,nrow=nrow(YT2)),k,p,0,0) 

                            trainZ <- Z1[2:nrow(Z1),]        

                            trainY <- matrix(YT2[(p+1):nrow(YT2),],ncol=k)                                  

                        }

                }


            if (group == "Basic") {

                if(VARX){

                    beta <- .lassoVARFistX(beta, trainZ, trainY,gamm, 1e-04,p,MN,k,k1,s+s1,m)

                }else{

                    beta <- .lassoVARFist(beta, trainZ, trainY,gamm, 1e-04,p,MN)
                }

            }

            if (group == "Lag") {

                GG <- .GroupLassoVAR1(beta,jj,jjcomp,trainY,trainZ,gamm,activeset,1e-04,p,MN,k,k1,s+s1)

                beta <- GG$beta

                activeset <- GG$active

            }

            if (group == "SparseLag") {

                if(VARX){

                    if(!dual){

                        GG <- .SparseGroupLassoVARX(beta, jj, trainY, trainZ, 
                                                    gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN,k,s+s1,k1)

                    }else{

                        GG <- .SparseGroupLassoVARXDual(beta, jj, trainY, trainZ, 
                                                        gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN,k,s+s1,k1)

                    }
                }else{

                    if(!dual){
                        GG <- .SparseGroupLassoVAR(beta, jj, trainY, trainZ, 
                                                   gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)

                    }else{
                        GG <- .SparseGroupLassoVARDual(beta, jj, trainY, trainZ, 
                                                       gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)


                    }
                }

                beta <- GG$beta

                activeset <- GG$active

                q1a <- GG$q1

            }

            if (group == "OwnOther") {

                if(VARX){

                    GG <- .GroupLassoOOX(beta, kk, trainY, trainZ, gamm, 
                                         activeset, 1e-04,p,MN,k,k1,s+s1)

                }else{

                    GG <- .GroupLassoOO(beta, kk, trainY, trainZ, gamm, 
                                        activeset, 1e-04,p,MN)
                }

                beta <- GG$beta

                activeset <- GG$active

            }

            if (group == "SparseOO") {
                if(VARX){


                    GG <- .SparseGroupLassoVAROOX(beta, kk, trainY, trainZ, 
                                                  gamm, alpha, INIactive = activeset, 1e-04,p,MN,k1,s+s1,k,dual)

                }else{

                    GG <- .SparseGroupLassoVAROO(beta, kk, trainY, trainZ, 
                                                 gamm, alpha, INIactive = activeset, 1e-04,q1a,p,MN,dual)

                    q1a <- GG$q1

                }

                beta <- GG$beta
                
                activeset <- GG$active

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
                            
                            MSFE[v - (T1 - h), ii] <- norm2(ZFull$Y[v+h-1,1:k1] - pred)^2
                                        # Subtract one from diagonal for warm start purposes
                            diag(beta[,2:(k1+1),ii]) <- diag(beta[,2:(k1+1),ii])-1
                            
                        }else{

                            if(object@crossval=="Rolling"){

                                pred <- beta[,,ii] %*% eZ

                                if(h>1 & recursive){
                                    pred <- matrix(pred,nrow=1)
                                    pred <- predictMS(pred,trainY,h-1,beta[,,ii],p)
                                }
                                MSFE[v - (T1 -h), ii] <- norm2(ZFull$Y[v+h-1,1:k1] - pred)^2


                            }else{
                                if(VARX){          

                                    eZ <- VARXCons(matrix(Y[(v-p):(v),1:k1],ncol=k1),Y[(v-p):(v),(ncol(Y)-m+1):(ncol(Y))],k1,p
                                                   ,m,s,contemp=contemp)

                                    pred <- beta[,,ii] %*% eZ


                                }else{

                                    eZ<- VARXCons(Y[(v-p):(v),1:k1],matrix(0,nrow=length((v-p):v)),k1,p,0,0)

                                    pred <- beta[,,ii] %*% eZ

                                }
                                MSFE[v - (T1 - h), ii] <- norm2(Y[v+h-1,1:k1] - pred)^2     


                            }                    

                        }

                    }
                }else{

                                        # If alpha and lambda are jointly fit, calculate MSFE for each alpha, lambda combination
                    for (ii in 1:gran2) {
                        for(jj in 1:length(alpha)){
                            if (RVAR) {

                                        # Relaxed Least Squares (intercept ignored)
                                beta[,,(ii-1)*nalpha+jj] <- RelaxedLS(cbind(t(trainZ),trainY),beta[,,(ii-1)*nalpha+jj],k,p,k1,s+s1)
                            }

                            if(MN){

                                pred <- beta[,2:dim(beta)[2],(ii-1)*nalpha+jj] %*% eZ[2:length(eZ)]                            
                                if(h>1 & recursive){
                                    pred <- matrix(pred,nrow=1)
                                    pred <- predictMS(pred,trainY,h-1,beta[,2:dim(beta)[2],(ii-1)*nalpha+jj],p,MN)
                                }

                                
                                MSFE[v - (T1 - h), (ii-1)*nalpha+jj] <- norm2(ZFull$Y[v+h-1,1:k1] - beta[,2:dim(beta)[2],(ii-1)*nalpha+jj] %*% eZ[2:length(eZ)])^2
                                
                            }else{
                                pred <- beta[,,(ii-1)*nalpha+jj] %*% eZ

                                if(h>1 & recursive){
                                    pred <- matrix(pred,nrow=1)
                                    pred <- predictMS(pred,trainY,h-1,beta[,,(ii-1)*nalpha+jj],p)
                                }
                                
                                if(object@crossval=="Rolling"){
                                    MSFE[v - (T1 - h), (ii-1)*nalpha+jj] <- norm2(ZFull$Y[v+h-1,1:k1] - beta[,,(ii-1)*nalpha+jj] %*% eZ)^2
                                    err <- norm2(ZFull$Y[v,1:k1] - beta[,,ii] %*% eZ)^2          
                                }else{
                                    if(VARX){          
                                        eZ<- VARXCons(matrix(Y[(v-p):(v),1:k1],ncol=k1),Y[(v-p):(v),(ncol(Y)-m+1):(ncol(Y))],k1,p,m,s,contemp=contemp)
                                    }else{
                                        eZ<- VARXCons(Y[(v-p):(v),1:k1],matrix(0,nrow=length((v-p):v)),k1,p,0,0)
                                    }
                                    MSFE[v - (T1 - h), (ii-1)*nalpha+jj] <- norm2(Y[v,1:k1] - beta[,,(ii-1)*nalpha+jj] %*% eZ)^2
                                    


                                }
                            }


                        }
                    }
                }



            

            if(verbose){
                setTxtProgressBar(pb, v)}
        }

                                        # Sort out indexing for 2-d gridsearch

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



                                        # one standard error correction     
        if(ONESE & !dual){
            MSFE2 <- MSFE 
            G2 <- colMeans(na.omit(MSFE2))
            G3 <- sd(na.omit(MSFE2))/sqrt(nrow(na.omit(MSFE2)))
            gamopt <- gamm[min(which(G2<(min(G2)+G3)))]
        }else{

            if(group!="Tapered" & !dual){
                                        # in rare cases in which MSFE is equal, the smaller penalty parameter is chosen.
                                        # This prevents extremely sparse solutions
                gamopt <- gamm[max(which(colMeans(na.omit(MSFE))==min(colMeans(na.omit(MSFE)))))]
            }else if(dual){
                
                if(!ONESE){
                    inds <- findind(which.min(colMeans(MSFE)),gamm[,1],alpha)                                   
                }else{

                    G2 <- colMeans(na.omit(MSFE))

                    G3 <- sd(na.omit(MSFE))/sqrt(nrow(na.omit(MSFE)))

                    inds <- findind(min(which(G2<(min(G2)+G3))),gamm[,1],alpha)
                    
                }
                gamopt <- gamm[inds[1],inds[2]]
                gamm <- gamm[,inds[2]]
                alphaopt <- alpha[inds[2]]
                ## print(paste("alphaopt",alphaopt))
            }
        }
        if(!dual){
            alphaopt <- alpha
        }
        
        if(VARX){
                                        # Out of sample forecast evaluation: VARX
            OOSEval <- .BigVAREVALX(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,k1,s,m,contemp,alphaopt)

        }else{
                                        # Out of sample evaluation for VAR    
            OOSEval <- .BigVAREVAL(ZFull,gamopt,k,p,group,h,MN,verbose,RVAR,palpha,T2,T,alphaopt,recursive)
        }
        MSFEOOSAgg <- na.omit(OOSEval$MSFE)
        betaPred <- OOSEval$betaPred
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

                                        #Residuals
        resids <- t(t(ZFull$Y)-betaPred%*%rbind(rep(1,ncol(ZFull$Z)),ZFull$Z))
        
        MSFEOOS<-mean(na.omit(MSFEOOSAgg))

        seoos <- sd(na.omit(MSFEOOSAgg))/sqrt(length(na.omit(MSFEOOSAgg)))

        if(!VARX){k1 <- k}


                                        # naive benchmarks     
        meanbench <- .evalMean(ZFull$Y[,1:k1],T2,T,h=h)
        RWbench <- .evalRW(ZFull$Y[,1:k1],T2,T,h=h)


        if(object@ic==FALSE|object@tf){

            AICbench <- list()
            AICbench$Mean <- as.double(NA)
            AICbench$SD <- as.double(NA)
            BICbench <- list()
            BICbench$Mean <- as.double(NA)
            BICbench$SD <- as.double(NA)                          


                                        # Information Criterion Benchmarks    

        }else{

            if(!VARX){

                X <- matrix(0,nrow=nrow(Y),ncol=k)

                AICbench1 <- VARXForecastEval(matrix(ZFull$Y[,1:k1],ncol=k1),X,p,0,T2,T,"AIC",h)

                AICbench <- list()

                AICbench$Mean <- mean(AICbench1)

                AICbench$SD <- sd(AICbench1)/sqrt(length(AICbench1))

                BICbench1 <- VARXForecastEval(matrix(ZFull$Y[,1:k1],ncol=k1),X,p,0,T2,T,"BIC",h)

                BICbench <- list()

                BICbench$Mean <- mean(BICbench1)

                BICbench$SD <- sd(BICbench1)/sqrt(length(BICbench1))

            }else{

                offset <- max(c(p,s))

                X <- matrix(Y[(offset+1):nrow(Y),(k1+1):ncol(Y)],ncol=m)

                AICbench1 <- VARXForecastEval(matrix(ZFull$Y[,1:k1],ncol=k1),as.matrix(X),p,s,T2,T,"AIC",h=h)

                AICbench <- list()

                AICbench$Mean <- mean(AICbench1)

                AICbench$SD <- sd(AICbench1)/sqrt(length(AICbench1))

                BICbench1 <- VARXForecastEval(matrix(ZFull$Y[,1:k1],ncol=k1),X,p,s,T2,T,"BIC",h=h)

                BICbench <- list()

                BICbench$Mean <- mean(BICbench1)

                BICbench$SD <- sd(BICbench1)/sqrt(length(BICbench1))  

            }

        }

        if(!VARX){contemp=FALSE}
                                        # Create a new BigVAR.Results Object
        results <- new("BigVAR.results",InSampMSFE=colMeans(MSFE),InSampSD=apply(MSFE,2,sd)/sqrt(nrow(MSFE)),LambdaGrid=gamm,index=which.min(colMeans(na.omit(MSFE))),OptimalLambda=gamopt,OOSMSFE=MSFEOOSAgg,seoosmsfe=seoos,MeanMSFE=meanbench$Mean,AICMSFE=AICbench$Mean,RWMSFE=RWbench$Mean,MeanSD=meanbench$SD,AICSD=AICbench$SD,BICMSFE=BICbench$Mean,BICSD=BICbench$SD,RWSD=RWbench$SD,Data=object@Data,lagmax=object@lagmax,Structure=object@Structure,Minnesota=object@Minnesota,Relaxed=object@Relaxed,Granularity=object@Granularity,horizon=object@horizon,betaPred=betaPred,Zvals=Zvals,resids=resids,VARXI=VARX,VARX=list(k=k1,s=s,contemp=contemp),preds=preds,T1=T1,T2=T2,dual=dual,alpha=alphaopt,crossval=object@crossval,ownlambdas=object@ownlambdas,tf=object@tf)
        
        return(results)
    }
    )


#' BigVAR Estimation
#' @description
#' Fit a BigVAR object with a structured penalty
#' @usage BigVAR.est(object)
#' @param object BigVAR object created from \code{ConstructModel}
#' @details Fits HVAR or VARX-L model on a BigVAR object.  Does not perform cross-validation.
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

        if(object@ownlambdas==TRUE){
            gamm=object@Granularity
            gran2 <- length(gamm)

        }      

        group <- object@Structure
        Y <- object@Data
        k <- ncol(Y)
        VARX <- object@VARX
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
        if(class(Y)!="matrix"){Y=matrix(Y,ncol=1)}
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

            if(!object@ownlambdas){

            
            if(dual){
                                        # Constructs penalty grid if both alpha and lambda are selected
                gamm <- .LambdaGridXDual(gran1, gran2, jj, trainY, trainZ,group,p,k1,s,m,k,MN,alpha)

            }else{

                                        # Penalty parameter grid for just lambda
                gamm <- .LambdaGridX(gran1, gran2, jj, as.matrix(trainY), trainZ,group,p,k1,s+s1,m,k,MN,alpha)
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

            gamm <- .LambdaGridEDual(gran1, gran2, jj, trainY, trainZ,group,p,k,MN,alpha)
            

        }else{
            gamm <- .LambdaGridE(gran1, gran2, jj, trainY, trainZ,group,p,k,MN)
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


    if (group == "Basic") {

        if(VARX){
            beta <- .lassoVARFistX(beta, trainZ, trainY,gamm, 1e-04,p,MN,k,k1,s,m)}
        else{
            beta <- .lassoVARFist(beta, trainZ, trainY,gamm, 1e-04,p,MN)
        }
    }

    if (group == "Lag") {

        GG <- .GroupLassoVAR1(beta,jj,jjcomp,trainY,trainZ,gamm,activeset,1e-04,p,MN,k,k1,s)
        beta <- GG$beta
        activeset <- GG$active
    }

    if (group == "SparseLag") {

        if(VARX){

            if(!dual){

                GG <- .SparseGroupLassoVARX(beta, jj, trainY, trainZ, 
                                            gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN,k,s+s1,k1)

            }else{

                GG <- .SparseGroupLassoVARXDual(beta, jj, trainY, trainZ, 
                                                gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN,k,s+s1,k1)

            }
        }else{

            if(!dual){
                GG <- .SparseGroupLassoVAR(beta, jj, trainY, trainZ, 
                                           gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)

            }else{
                GG <- .SparseGroupLassoVARDual(beta, jj, trainY, trainZ, 
                                               gamm, alpha, INIactive = activeset, 1e-04, q1a,p,MN)


            }
        }

        beta <- GG$beta

        activeset <- GG$active

        q1a <- GG$q1

    }

    if (group == "OwnOther") {
        if(VARX){
            GG <- .GroupLassoOOX(beta, kk, trainY, trainZ, gamm, 
                                 activeset, 1e-04,p,MN,k,k1,s)
        }else{
            GG <- .GroupLassoOO(beta, kk, trainY, trainZ, gamm, 
                                activeset, 1e-04,p,MN)
        }
        beta <- GG$beta
        activeset <- GG$active
    }

    if (group == "SparseOO") {
        if(VARX){


            GG <- .SparseGroupLassoVAROOX(beta, kk, trainY, trainZ, 
                                          gamm, alpha, INIactive = activeset, 1e-04,p,MN,k1,s+s1,k,dual)

        }else{

            GG <- .SparseGroupLassoVAROO(beta, kk, trainY, trainZ, 
                                         gamm, alpha, INIactive = activeset, 1e-04,q1a,p,MN,dual)

            q1a <- GG$q1

        }

        beta <- GG$beta
        
        activeset <- GG$active

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
#' This class contains the results from cv.BigVAR.
#'
#' It inherits the class BigVAR, but contains substantially more information. 
#' 
#' @field InSampMSFE In-sample MSFE from optimal value of lambda
#' @field LambdaGrid Grid of candidate lambda values
#' @field index Rank of optimal lambda value 
#' @field OptimalLambda Value of lambda which minimizes MSFE
#' @field OOSMSFE Average Out of sample MSFE of BigVAR model with Optimal Lambda
#' @field seoosfmsfe Standard Error of Out of sample MSFE of BigVAR model with Optimal Lambda
#' @field MeanMSFE Average Out of sample MSFE of Unconditional Mean Forecast
#' @field MeanSD Standard Error of out of sample MSFE of Unconditional Mean Forecast
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
#' @field dual indicator as to whether dual cross validation was conducted
#' @field contemp indicator if contemporaneous exogenous predictors are used
#'
#' @note One can also access any object of class BigVAR from BigVAR.results
#' @name BigVAR.results 
#' @rdname BigVAR.results
#' @aliases BigVAR.results-class
#' @exportClass BigVAR.results
#' @author Will Nicholson
#' @export
setClass("BigVAR.results",
representation(InSampMSFE="numeric",InSampSD="numeric",LambdaGrid="numeric",index="numeric",OptimalLambda="numeric",OOSMSFE="numeric",seoosmsfe="numeric",MeanMSFE="numeric",AICMSFE="numeric",BICMSFE="numeric",BICSD="numeric",RWMSFE="numeric",MeanSD="numeric",AICSD="numeric",RWSD="numeric",betaPred="matrix",Zvals="matrix",VARXI="logical",resids="matrix",preds="matrix",dual="logical",contemp="logical"),
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
#' @method plot
#' @rdname BigVAR.results-plot-methods
#' @export
setMethod(f="plot",signature="BigVAR.results",
definition= function(x,y=NULL,...)
    {

        plot(x@LambdaGrid,x@InSampMSFE,type="o",xlab="Value of Lambda",ylab="MSFE",log="x")

        abline(v=x@OptimalLambda,col="green")

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
#' @method show
#' @rdname show-methods-BigVAR.results
#' @export
setMethod("show","BigVAR.results",
function(object)
    {
        cat("*** BIGVAR MODEL Results *** \n")
        cat("Structure\n") ;print(object@Structure)
        if(object@Relaxed==TRUE){
            cat("Relaxed VAR \n") ;print(object@Relaxed)}

        cat("Forecast Horizon \n") ;print(object@horizon)

        if(object@VARXI){
            cat("VARX Specs \n") ;print(object@VARX)}
        cat("Maximum Lag Order \n") ;print(object@lagmax)
        cat("Optimal Lambda \n"); print(round(object@OptimalLambda,digits=4))
        if(object@dual){

            cat("Optimal Alpha \n"); print(round(object@alpha,digits=2))
            
        }
        
        cat("Grid Depth \n") ;print(object@Granularity[1])
        cat("Index of Optimal Lambda \n");print(object@index)
        cat("In-Sample MSFE\n");print(round(min(object@InSampMSFE),digits=3))
        cat("BigVAR Out of Sample MSFE\n");print(round(mean(object@OOSMSFE),digits=3))
        cat("*** Benchmark Results *** \n")
        cat("Conditional Mean Out of Sample MSFE\n");print(round(object@MeanMSFE,digits=3))
        cat("AIC Out of Sample MSFE\n");print(round(object@AICMSFE,digits=3))
        cat("BIC Out of Sample MSFE\n");print(round(object@BICMSFE,digits=3))
        cat("RW Out of Sample MSFE\n");print(round(object@RWMSFE,digits=3))
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
#' @method predict
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
        betaPred <- object@betaPred
        Y <- object@Data
        k <-object@VARX$k
        m <- ncol(object@Data)-k
        p <- object@lagmax
        s <- object@VARX$s
        VARX <- object@VARXI
        contemp <- object@VARX$contemp
        s1 <- 0
        if(VARX){
            if(object@contemp){
                s1 <- 1
            }
        }else{
            s1=0
        }
        fcst <- betaPred%*%eZ

        if(n.ahead==1)
            {
                return(fcst)
            }else{
                if(!VARX){
                                        # iterative multistep forecasts
                    fcst <- predictMS(matrix(fcst,nrow=1),Y[(nrow(Y)-p+1):nrow(Y),],n.ahead-1,betaPred,p,MN)

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

                            fcst <- predictMSX(matrix(fcst,nrow=1),as.matrix(Y[(nrow(Y)-C+1):nrow(Y),1:(k)]),n.ahead-1,betaPred,p,newxreg,matrix(Y[(nrow(Y)-C+1):nrow(Y),(ncol(Y)-m+1):ncol(Y)],ncol=m),m,s,1,MN)
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
#' @details Uses \code{levelplot} from the \code{lattice} package to plot the magnitude of each coefficient
#' @name SparsityPlot.BigVAR.results
#' @aliases SparsityPlot.BigVAR.results,BigVAR.results-method
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



    
 # BigVAR fit used in rolling cv, out of sample forecast evaluation and BigVAR.est
.BigVAR.fit <- function(group,beta,trainZ,trainY,gamm,tol,p,m=0,k1,k,s=0,s1=0,MN=FALSE,C,intercept,separate_lambdas,dual,activeset=NULL,q1a=NULL,jj=NULL,jjcomp=NULL,VARX=list(),alpha=NULL,kk,palpha)
{
    if(is.null(s)){
        s=0
    }
    if(is.null(s1)){
        s1=0
    }
        if(is.null(m)){
        m=0
    }
    if (group == "Basic") {

        if(VARX){

            beta <- .lassoVARFistX(beta, trainZ, trainY,gamm, tol,p,MN,k,k1,s+s1,m,C,intercept,separate_lambdas)

        }else{

            beta <- .lassoVARFist(beta, trainZ, trainY,gamm, tol,p,MN,C,intercept,separate_lambdas)
        }

    }


    if (group == "BasicEN") {
        ## browser()
        if(VARX){

            beta <- .lassoVARFistXEN(beta, trainZ, trainY,gamm,alpha, tol,p,MN,k,k1,s+s1,m,C,intercept,separate_lambdas)

        }else{

            beta <- .lassoVARFistEN(beta, trainZ, trainY,gamm,alpha, tol,p,MN,C,intercept,separate_lambdas)
        }

    }

    
    if (group == "Lag") {

        GG <- .GroupLassoVAR1(beta,jj,jjcomp,trainY,trainZ,gamm,activeset,tol,p,MN,k,k1,s+s1,C,intercept)

        beta <- GG$beta

        activeset <- GG$active

    }

    if (group == "SparseLag") {

        if(VARX){

            if(!dual){

                GG <- .SparseGroupLassoVARX(beta, jj, trainY, trainZ, 
                                            gamm, alpha, INIactive = activeset, tol, q1a,p,MN,k,s+s1,k1,C,intercept)

            }else{

                GG <- .SparseGroupLassoVARXDual(beta, jj, trainY, trainZ, 
                                                gamm, alpha, INIactive = activeset, tol, q1a,p,MN,k,s+s1,k1,C,intercept)

            }
        }else{

            if(!dual){
                GG <- .SparseGroupLassoVAR(beta, jj, trainY, trainZ, 
                                           gamm, alpha, INIactive = activeset, tol, q1a,p,MN,C,intercept)

            }else{
                GG <- .SparseGroupLassoVARDual(beta, jj, trainY, trainZ, 
                                               gamm, alpha, INIactive = activeset, tol, q1a,p,MN,C,intercept)


            }
        }

        beta <- GG$beta

        activeset <- GG$active

        q1a <- GG$q1

    }

    if (group == "OwnOther") {

        if(VARX){

            GG <- .GroupLassoOOX(beta, kk, trainY, trainZ, gamm, 
                                 activeset, tol,p,MN,k,k1,s+s1,C,intercept)

        }else{

            
            GG <- .GroupLassoOO(beta, kk, trainY, trainZ, gamm, 
                                activeset, tol,p,MN,C,intercept)
        }

        beta <- GG$beta

        activeset <- GG$active

    }

    if (group == "SparseOO") {
        if(VARX){

## browser()
            GG <- .SparseGroupLassoVAROOX(beta, kk, trainY, trainZ, 
                                          gamm, alpha, INIactive = activeset, tol,p,MN,k1,s+s1,k,dual,C,intercept)

        }else{

            GG <- .SparseGroupLassoVAROO(beta, kk, trainY, trainZ, 
                                         gamm, alpha, INIactive = activeset, tol,q1a,p,MN,dual,C,intercept)

            q1a <- GG$q1

        }

        beta <- GG$beta
        
        activeset <- GG$active

    }

    if(group=="Tapered")
    {
        ## palpha <- seq(0,1,length=10)
        ## palpha <- rev(palpha)
        
        beta <- .lassoVARTL(beta,trainZ,trainY,gamm,tol,p,MN,palpha,C,intercept)    
    }

    if(group=="EFX")
    {

        beta <- .EFVARX(beta,trainY,trainZ,gamm,tol,MN,k1,s,m,p,C,intercept)
        
    }

    if(group=="HVARC")
    {
        beta <- .HVARCAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept,separate_lambdas)

    }

    if(group=="HVAROO")
    {
        beta <- .HVAROOAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept,separate_lambdas)
    }

    if(group=="HVARELEM")
    {

        beta <- .HVARElemAlg(beta,trainY,trainZ,gamm,tol,p,MN,C,intercept,separate_lambdas)

    }

    if(group=="BGR"){
        ## browser()
        trainZ <- rbind(1,trainZ)
        beta <- BGRGridSearch(trainY,trainZ,p,gamm,as.numeric(MN))
    }

    if(group=="MCP"|group=="SCAD")
    {

        beta <- .MCPFit(beta,trainZ,trainY,gamm,tol,p,MN,k,k1,s,m,C,intercept,group,gamma)

    }
    
    
    
    if(!exists('activeset')){
        activeset <- NULL
    }
    if(!exists('q1a')){
        q1a <- NULL
    }
    return(list(beta=beta,activeset=activeset,q1a=q1a))

}



#' Simple function to fit BigVAR model with fixed penalty parameter
#' @param Y \eqn{T \times k} multivariate time series or Y \eqn{T \times (k+m)} endogenous and exogenous series, respectively 
#' @param p Predetermined maximal lag order (for modeled series)
#' @param struct The choice of penalty structure (see details).
#' @param lambda vector or matrix of penalty parameters.
#' @param intercept True or False: option to fit an intercept
#' @param RVAR True or False: option to refit based upon the support selected using the Relaxed-VAR procedure
#' @param MN Minnesota Prior Indicator
#' @param VARX List containing VARX model specifications. 
#' @param alpha grid of candidate parameters for the alpha in the Sparse Lag and Sparse Own/Other VARX-L 
#' @param C vector of coefficients to shrink toward a random walk (if \code{MN} is \code{TRUE})
#' @param tf transfer function indicator (i.e. VARX in which p=0 & s>0) (default false)
#' @param tol optimization tolerance (default 1e-4)
#' @param separate_lambdas indicator for separate penalty parameters for each time series (default \code{FALSE})
#' @param beta optional \eqn{k\times (k\times p + m\times s +1)} coefficient matrix to use as a "warm start" (default \code{FALSE})
#' 
#'  @details The choices for "struct" are as follows
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
#' \item{  "BGR" (Bayesian Ridge Regression (cf. Banbura et al))}
#' }
#'
#' VARX specifications consist of a list with entry k denoting the series that are to be modeled and entry s to denote the maximal lag order for exogenous series.
#'
#' The argument alpha is ignored unless the structure choice is "SparseLag" or "Lag."  By default "alpha" is set to \code{NULL} and will be initialized as 1/(k+1) in \code{cv.BigVAR} and \code{BigVAR.est}.  Any user supplied values must be between 0 and 1.  

#' @note The specifications "Basic", "Lag," "SparseLag," "SparseOO," and "OwnOther" can accommodate both VAR and VARX models.  EFX only applies to VARX models.  "HVARC," "HVAROO," "HVARELEM," and "Tapered" can only be used with VAR models.
#'
#' @seealso \code{\link{cv.BigVAR}},\code{\link{BigVAR.est}},\code{\link{constructModel}}
#' 
#' @references  William B Nicholson, Jacob Bien, and David S Matteson. "High Dimensional Forecasting via Interpretable Vector Autoregression." arXiv preprint arXiv:1412.5250, 2016.
#' William B. Nicholson, David S. Matteson, Jacob Bien,VARX-L: Structured regularization for large vector autoregressions with exogenous variables, International Journal of Forecasting, Volume 33, Issue 3, 2017, Pages 627-651,
#' William B Nicholson, David S. Matteson, and Jacob Bien (2016), "BigVAR: Tools for Modeling Sparse High-Dimensional Multivariate Time Series" arxiv:1702.07094
#'
#' Banbura, Marta, Domenico Giannone, and Lucrezia Reichlin. "Large Bayesian vector auto regressions." Journal of Applied Econometrics 25.1 (2010): 71-92.
#' @examples
#' # VARX Example
#' # Fit a Basic VARX-L with k=2, m=1, s=2, p=4, lambda=1e-2
#' VARX=list()
#' VARX$k=2 # indicates that the first two series are modeled
#' VARX$s=2 # sets 2 as the maximal lag order for exogenous series
#' data(Y)
#' BigVAR.fit(Y,p=4,"Basic",lambda=1e-2,VARX=VARX)
#' @export
BigVAR.fit <- function(Y,p,struct,lambda,alpha=NULL,VARX=list(),separate_lambdas=F,MN=F,C=as.double(NULL),intercept=TRUE,tf=F,tol=1e-4,RVAR=F,beta=NULL){
    if(!is.matrix(Y)){stop("Y needs to be a matrix")}
    if(is.null(lambda)){stop("Must include penalty parameter")}
    if(is.null(alpha)){alpha=1/(ncol(Y)+1)}
    dual=FALSE
    k <- ncol(Y)
                                        # runs constructModel just to check for errors
    temp.bv <- constructModel(Y,p,struct=struct,gran=c(lambda),ownlambdas=TRUE,MN,intercept=intercept,C=C,VARX=VARX,cv="LOO")
    gran2 <- length(lambda)
    ## structures=c("Basic","Lag","SparseLag","OwnOther","SparseOO","HVARC","HVAROO","HVARELEM","Tapered","EFX","BGR")
    group=struct
    if(separate_lambdas){
        if(is.vector(lambda)){
            lambda <- matrix(lambda,nrow=1)
            }
        gran2 <- nrow(lambda)
    }
    
    jj <- 0
    if(!"matrix"%in%class(Y)){Y=matrix(Y,ncol=1)}
    s <- ifelse(VARX,VARX$s,0)
    T <- nrow(Y)-max(p,s)
    if(length(VARX)!=0){

        
        k1 <- VARX$k
        s <- VARX$s

        if(!is.null(VARX$contemp)){

            contemp <- TRUE
            s1 <- 1

        }else{

            contemp <- FALSE
            s1 <- 0
        }
        VARX <- TRUE
        m <- k-k1
        Y1 <- matrix(Y[,1:k1],ncol=k1)
        X <- matrix(Y[,(ncol(Y)-m+1):ncol(Y)],ncol=m)

        if(!tf){
            
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

        
                                        # Initial Coefficient Matrix         
        beta1=array(0,dim=c(k1,k1*p+(k-k1)*(s+s1)+1,gran2*length(alpha)))
        if(!is.null(beta)&all(dim(beta)==dim(beta1))){
            beta=beta
        }else{
            beta=beta1
        }
                                        # Initialize groups, active sets, power method calculations, etc
            kk <- NULL
            jj <- NULL
            jjcomp <- NULL
            activeset <- NULL
            q1a <- NULL

        if (group == "Lag") {

            jj=groupfunVARX(p,k,k1,s+s1)

            jjcomp <- groupfunVARXcomp(p,k,k1,s+s1)

            activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                             gran2)

        }else if (group == "SparseLag") {

            jj <- groupfunVARX(p, k,k1,s+s1)

            q1a <- list()


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
                                        # No VARX
        VARX <- FALSE
        s=p
        s1=0
        m=0
        if (group=="Lag"|group=="SparseLag")

        {

            jj=.groupfun(p,k)

        }else{

            jj <- .lfunction3(p,k)

        }


        Z1 <- VARXCons(Y,matrix(0,nrow=nrow(Y)),k,p,0,0) 

            trainZ <- Z1[2:nrow(Z1),,drop=F]   

        trainY <- matrix(Y[(p+1):nrow(Y),],ncol=k)          

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

            k1 <- k

        }else if (group == "SparseLag") {

            jj <- .groupfuncpp(p, k)

            q1a <- list()

            for (i in 1:p) {

                q1a[[i]] <- matrix(runif(k, -1, 1), ncol = 1)
            }


            activeset <- rep(list(rep(rep(list(0), length(jj)))), 
                             gran2)
            

        }else if (group == "OwnOther") {

            kk <- .lfunction3cpp(p, k)

            activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                             gran2)

        }else if (group == "SparseOO") {

            kk <- .lfunction3cpp(p, k)
            jjcomp <- .lfunctioncomp(p,k)
            jj <- .lfunction3(p,k)
            activeset <- rep(list(rep(rep(list(0), length(kk)))), 
                             gran2)
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
    
    if(group=="Tapered")

    {
        palpha <- seq(0,1,length=10)
        palpha <- rev(palpha)
        gran2 <- length(lambda)*length(palpha)
        beta <- array(0,dim=c(k,k*p+1,gran2))
    }


    if (group == "BGR") {
        
        trainZ <- rbind(1,trainZ)
        beta <- BGRGridSearch(trainY,trainZ,p,lambda,as.numeric(MN))
    }else{
    ## browser()
    
    ## needed.objs <- c('group','beta','trainZ','trainY','lambda','tol','p','m','k1','s','s1','m','MN','C','intercept','separate_lambdas','dual','activeset','q1a')
    ## ae=function(x){any(exists(x))}
    ## objs <- sapply(needed.objs,ae)
    ## if(!all(objs)){
    ##     oss <- names(objs)[!objs]
    ##     for(i in 1:length(oss)){               
    ##         assign(oss[i],NULL)
    ##     }
    ## }
    ## browser()
    temp <- .BigVAR.fit(group,beta,trainZ,trainY,lambda,tol,p,m,k1,k,s,s1,MN,C,intercept,separate_lambdas,dual,activeset,q1a,jj,jjcomp,VARX,alpha,kk)
    beta <- temp$beta
    }

    
                                        # refit if varx
    if(RVAR){
        beta <- RelaxedLS(cbind(t(trainZ),trainY),beta,k,p,k1,s+s1)
    }
    
    ## activeset <- temp$activeset
    ## q1a <- temp$q1a
    
    return(beta)
}



IFGLS <- function(Y,Z,B2,p,k,SigmaU,zerothresh,eps,oracle,WLS)
{

    pk <- k*p+1
    pkk <- k^2*p+k
    
    K <- cbind(t(Z),Y)
    K <- cbind(1,K)

                                        # add ridge penalty only if necessary
    if(nrow(Y)<(pk+k)){
        q = ncol(K)
        delta = (q^2 + q + 1) * (sqrt(.Machine$double.eps))
        scale = sqrt(delta) * sqrt(apply(K^2, 2, sum))
        K <- rbind(K,diag(scale))
        Y <- K[,(ncol(K)-k+1):ncol(K)]
    }

    Z <- (K[,1:(ncol(K)-k)])

    # get zero restrictions from VARX coefficient matrix 

    RB <- Matrix(RestMatFull(t(B2),pkk,zerothresh),sparse=TRUE)


    RR <- ncol(as.matrix(RB))


    T1 <- nrow(Y)

    # generate QR matrices for use later
    Q1Full <- QRConsM2(Z,T1,k,B2,RR,p,zerothresh)

    Q1AA <- Q1Full$Q1A
    Q1AR <- Q1Full$Q1R
    QQ <- Q1Full$QQ

    R22 <- qr.R(qr(K), complete = TRUE)
    R12 <- R22[1:(pk),(pk+1):(pk+k)]
    Y22 <- R22[(pk+1):(pk+k),(pk+1):(pk+k)]

    RZ <- qr.R(qr(Z,complete=F))

    Rres <- Matrix(kronecker(diag(k),RZ)%*%RB,sparse=T)

    Y22C <- crossprod(Y22)

    if(qr(SigmaU)$rank<k){
        print("Warning, low rank matrix")
        udv <- svd(SigmaU)
        C <- qr.R(qr(diag(sqrt(udv$d))%*%t(udv$u)))

    }else{
        C <- chol(SigmaU)
    }

    YstarFull <- (Q1AA)%*%as.vector(Y)
                                        # iterative portion
    thresh=10*eps
    j=0
    if(oracle|WLS)
    {
        test1 <- KronMatcppEigen(k,T1,Z,RR,B2,C,p,RR,QQ,zerothresh,as.matrix(YstarFull),as.matrix(Q1AR))

        BFGLS <- matrix(RB%*%test1$BFGLS,nrow=k,ncol=k*p+1,byrow=T)

        SigmaUOld <- SigmaU
        resid <- as.vector(R12)-Rres%*%test1$BFGLS

        resid <- matrix(resid,ncol=k,byrow=F)

        SigmaU <- (crossprod((resid))+Y22C)/(T-pk)

    }else{
        jmax=5
        # iterative portion
        while(thresh>eps)
        {

            test1 <- KronMatcppEigen(k=as.integer(k),T1=as.integer(T1),Z=Z,RR=as.integer(RR),BHAT=B2,C=C,p=p,R2i=as.integer(RR),QQ=QQ,zerothresh,as.matrix(YstarFull),as.matrix(Q1AR))

            BFGLS <- matrix(RB%*%test1$BFGLS,nrow=k,ncol=k*p+1,byrow=T)

            SigmaUOld <- SigmaU
            resid <- as.vector(R12)-Rres%*%test1$BFGLS

            resid <- matrix(resid,ncol=k,byrow=F)

            SigmaU <- (crossprod((resid))+Y22C)/T

            if(qr(SigmaU)$rank<k){
                udv <- svd(SigmaU)
                C <- qr.R(qr(diag(sqrt(udv$d))%*%t(udv$u)))
            }else{
                C <- chol(SigmaU)
            }


            thresh=onorm(SigmaUOld-SigmaU)

            j=j+1
            if(j>jmax){
                break
            }

        }

    }
    return(list(BFGLS=BFGLS,SigmaU=SigmaU))
}

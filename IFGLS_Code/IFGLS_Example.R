library(Rcpp)
# source IFGLS C++ file
sourceCpp('./IFGLS_Standalone.cpp')
# source VARX Lag matrix construction
sourceCpp('./BigVAR/BigVAR/src/DataCons.cpp')
library(Matrix)
library(BigVAR)
data(Y)
p=4;k=ncol(Y)
# run cv.BigVAR to get an example coefficient matrix
A <- constructModel(Y,p=p,struct="Basic",gran=c(50,10))

res <- cv.BigVAR(A)

Z <- VARXCons(Y,matrix(0,nrow=nrow(Y)),k=ncol(Y),p=p,m=0,s=0)

Y <- Y[(p+1):nrow(Y),]
# extract coefficient matrix
B2=res@betaPred
zerothresh=eps=1e-4
Z <- Z[2:nrow(Z),]
# initial value for SigmaU-sample covariance
SigmaU <- crossprod(Y)/nrow(Y)

source('/home/will/Dropbox/BigVAR_ec/IFGLS_function.R')
IFGLS_res <- IFGLS(Y,Z,B2,p,k,SigmaU=SigmaU,zerothresh=1e-4,eps=1e-4,oracle=FALSE,WLS=FALSE)


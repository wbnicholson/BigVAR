
#include <math.h>
#include <vector>
#include <limits>
#include <algorithm>

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace std;
using namespace Rcpp;


//[[Rcpp::export]]
MatrixXd ZmatF(MatrixXd Y, const int p, const int k,bool intercept=true,bool oos=false)
{
  int T=Y.rows();
    if(oos){
	  T+=1;
  }

  MatrixXd Y2=Y.rowwise().reverse();
  MatrixXd Y2a=Y2.topLeftCorner(p,k);
  Y2a.transposeInPlace();
  VectorXd Y1(Map<VectorXd>(Y2a.data(),Y2a.cols()*Y2a.rows()));
  MatrixXd Z(k*p,T-p);
  Z.setZero();
  Z.col(0)=Y1.reverse();
  VectorXd Y1c(Y1.size());
  for(int i=1;i<(T-p);++i){
    MatrixXd Y1M=Y2.block(i,0,p,k);
    Y1M.transposeInPlace();
    VectorXd Y1N(Map<VectorXd>(Y1M.data(),Y1M.cols()*Y1M.rows()));
    Z.col(i)=Y1N.reverse();
  }
  MatrixXd ones(1,T-p);
  ones.setOnes();
  if(p==0){return(ones);

  }
  if(intercept){
  MatrixXd ZF(k*p+1,T-p);
  ZF <<ones,Z;
    return(ZF);
  }else{
    return(Z);}
	      
}
//[[Rcpp::export]]
MatrixXd VARXCons(MatrixXd Y, MatrixXd X, const int k, const int p, const int m, const int s,bool oos=false)
{
  // int T=Y.cols();
  int T=Y.rows();
  MatrixXd Z1=ZmatF(Y,p,k,true,oos);
  MatrixXd Z2=ZmatF(X,s,m,false,oos);
  if(s==0){return(Z1);
      }
  if(p==0){
  MatrixXd ones(1,T-s);
  ones.setOnes();
  MatrixXd ZF(m*s+1,T-s);
  ZF <<ones,Z2;
return(ZF);}
  if(p!=0 & s!=0){
    if(p>s){
      int T1=Z2.cols();      
Z2=Z2.rightCols(T1-(p-s));}
    if(p<s){
      int T1=Z1.cols();
Z1=Z1.rightCols(T1-(s-p));}
    MatrixXd ZZ(Z1.rows()+Z2.rows(),Z1.cols());
    ZZ<< Z1,Z2;
    return(ZZ);
}
}

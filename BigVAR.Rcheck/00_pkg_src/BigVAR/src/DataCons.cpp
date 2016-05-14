
#include <math.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;
using namespace std;

//triangular solver 
//[[Rcpp::export]]
MatrixXd backsolve(const MatrixXd R2, const MatrixXd R){

	int k = R.rows();
	MatrixXd RI(k,k);
	RI=R.triangularView<Upper>().solve(R2);
	return RI;

}


//[[Rcpp::export]]
MatrixXd ZmatF(MatrixXd Y,  int p, const int k,bool intercept=true,bool oos=false,bool contemp=false)
{
	int T=Y.rows();
	// some issues with out of sample predictions and contemporaneous dependence
    if(oos & !contemp){
		T+=1;
	}

	MatrixXd Y2=Y.rowwise().reverse();
	if(contemp){p+=1;}
	MatrixXd Y2a=Y2.topLeftCorner(p,k);
	Y2a.transposeInPlace();
	VectorXd Y1(Map<VectorXd>(Y2a.data(),Y2a.cols()*Y2a.rows()));
	int M;
	if(contemp){M=T-p+1;
	}else{M=T-p;}

	MatrixXd Z(k*p,M);
	Z.setZero();
	Z.col(0)=Y1.reverse();
	VectorXd Y1c(Y1.size());
	for(int i=1;i<M;++i){
		MatrixXd Y1M=Y2.block(i,0,p,k);
		Y1M.transposeInPlace();
		VectorXd Y1N(Map<VectorXd>(Y1M.data(),Y1M.cols()*Y1M.rows()));
		Z.col(i)=Y1N.reverse();
	}
	MatrixXd ones(1,M);
	ones.setOnes();
	if(p==0){return(ones);

	}
	if(intercept){
		MatrixXd ZF(k*p+1,M);
		ZF <<ones,Z;
		return(ZF);
	}else{
		return(Z);}
	      
}
//[[Rcpp::export]]
MatrixXd VARXCons(MatrixXd Y, MatrixXd X, const int k, const int p, const int m,  int s,bool oos=false,bool contemp=false)
{
	int T=Y.rows();

	MatrixXd  Z1=ZmatF(Y,p,k,true,oos,false);
	MatrixXd Z2=ZmatF(X,s,m,false,oos,contemp);
	if(s==0 & !(contemp)){return(Z1);
	}
	if(p==0){
		MatrixXd ZF=ZmatF(X,s,m,true,oos,contemp);
	  
		return(ZF);
	}
	  
	//adjusting for contemp. dependence
	if(contemp & oos){Z1=Z1.rightCols(Z1.cols()-1);}
	if(p!=0 & (s!=0||contemp)){
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

//Fitting AIC/BIC using Neumaier's least square technique


//[[Rcpp::export]]
List ARFitVARX(const MatrixXd K2, const int k, const int p,int m, int s)
{
	// int RC=k*p+m*s+1;
	int RC=K2.cols()-k;
	int T=K2.rows();
	int q=K2.cols();
	double delta = (pow(q,2)+q+1)*sqrt(std::numeric_limits<double>::epsilon());
	VectorXd D1(K2.cols());
	for(int i=0;i<q;++i)
		{
			D1(i)=K2.col(i).norm();

		}
	D1=sqrt(delta)*D1;
	MatrixXd AA=D1.asDiagonal();
	MatrixXd K3(K2.rows()+AA.rows(),K2.cols());
	K3 << K2, AA;
  
	HouseholderQR<MatrixXd> QR1(K3);
	// HouseholderQR<MatrixXd> QR1(K2);
	MatrixXd R = QR1.matrixQR().triangularView<Upper>();
	MatrixXd R11=R.topLeftCorner(RC,RC);
	MatrixXd R22=R.block(RC,RC,k,k);
	MatrixXd R12=R.topRightCorner(RC,k);
	MatrixXd Test= backsolve(R12,R11);
	MatrixXd Sigma=(R22.transpose()*R22)/T;

	MatrixXd NewR = R.block(0,0,R.cols(),R.cols());

	List Results=List::create(Named("B")=Test.transpose(),Named("SigmaU")=Sigma,Named("R")=NewR);
	return(Results);			      
			   		      
}			    


//[[Rcpp::export]]
List ICX(MatrixXd Y, MatrixXd X, double k, int pmax,int smax,double m,std::string pen,int h=1)
{

	std::vector<double> crit;
	int T=Y.rows();
	double criti=0;
	for(int i =0;i<=pmax;++i)
		{
			for(int j=0;j<=smax;++j)    
				{
					if((int) (k*i+j*m)>T)
						{
							crit.push_back(1000000);
							// Rcout<<"Undetermined"<<std::endl;
							break;
						}else{
		  
						MatrixXd sigmaU;
						MatrixXd Z=VARXCons(Y,X,k,i,m,j);
						int c=max(i,j)+h-1;
						MatrixXd Y2=Y.bottomRows(T-c);
						Z=Z.leftCols(Z.cols()-h+1);
						MatrixXd K(Y2.rows(),Z.rows()+Y2.cols());
						K<<Z.transpose(),Y2;
						List Bres = ARFitVARX(K,k,i,m,j);
						sigmaU = as<MatrixXd>(Bres("SigmaU"));
						double T2=Y2.rows();

						double p= (double) i;
						double s= (double) j;
						if(pen=="BIC"){
							criti=log(sigmaU.determinant())+log(T2)*(pow(k,2)*p+k*m*s+k)/T2;
						}else if(pen=="AIC"){
							criti=log(sigmaU.determinant())+(2*(pow(k,2)*p+k*m*s+k))/T2;
						}
						// Rcout<< "p is " << i<<endl;
						// Rcout<< "s is " << j<<endl;
						// Rcout<<criti<<endl;

						crit.push_back(criti);
					}
				}
		}
	int min = std::distance(crit.begin(),std::min_element(crit.begin(),crit.end()));
	// Rcout<<min<<std::endl;
	int phat;
	int shat;
	if(smax==0){
		phat=min;
		shat=0;
	}else{
 
		if(min!=0){
			if(smax!=0){
				shat=min % (smax+1);
			}else{shat=0;}
			if(smax!=0){
				phat=min/(smax+1);
			}
		}else{
			phat=0;
			shat=0;
		}
	}
	MatrixXd BResF;
	MatrixXd ZF=VARXCons(Y,X,k,phat,m,shat);
	int cf=max(phat,shat);
	MatrixXd Y2F=Y.bottomRows(T-cf);
	MatrixXd KF(Y2F.rows(),ZF.rows()+Y2F.cols());
	KF<<ZF.transpose(),Y2F;
	List Bres= ARFitVARX(KF,k,phat,m,shat);
	BResF =Bres("B");

	MatrixXd Sigma=Bres("SigmaU");
	
	List Results=List::create(Named("B")=BResF,Named("p")=phat,Named("s")=shat,Named("SigmaU")=Sigma);

	return(Results);
  
}




//[[Rcpp::export]]
List ARFitVARXDD(MatrixXd R,int T,int k,int p,int pmax)
{
  int RC=R.cols()-k;	
  MatrixXd R11=R.topLeftCorner(RC-k*p,RC-k*p);
  // MatrixXd R22=R.bottomRightCorner(2*k*p,k);

  MatrixXd R12prime(p*k,k);
  if(p==pmax){
	  R12prime.setZero();
		  }else{
  R12prime=R.block(R.cols()-(p+1)*k,R.cols()-k,p*k,k);
  }
  MatrixXd R22=R.bottomRightCorner(k,k);

  MatrixXd R12=R.topRightCorner(RC-k*p,k);


  MatrixXd Test= backsolve(R12,R11);
  MatrixXd Sigma=(R22.transpose()*R22+R12prime.transpose()*R12prime)/T;
  // MatrixXd NewR = R.block(0,0,RC,RC);
  List Results=List::create(Named("B")=Test.transpose(),Named("SigmaU")=Sigma,Named("R12prime")=R12prime);
  return(Results);			      
			   		      
}			    




//[[Rcpp::export]]
List ICDD(MatrixXd Y, double k, int pmax,std::string pen)
{

	std::vector<double> crit;
	int T=Y.rows();
	double criti=0;

	MatrixXd sigmaU;
	MatrixXd Z=ZmatF(Y,pmax,k);
	MatrixXd Y2=Y.bottomRows(T-pmax);
	MatrixXd K(Y2.rows(),Z.rows()+Y2.cols());
	K<<Z.transpose(),Y2;
	List Bres = ARFitVARX(K,k,pmax,0,0);
	double T2=Y2.rows();
	sigmaU = as<MatrixXd>(Bres("SigmaU"));
	double p= (double) pmax;

	if(pen=="BIC"){
		criti=log(sigmaU.determinant())+log(T2)*(pow(k,2)*p+k)/T2;
	}else if(pen=="AIC"){
		criti=log(sigmaU.determinant())+(2*(pow(k,2)*p+k))/T2;
	}
	crit.push_back(criti);

	MatrixXd R = as<MatrixXd>(Bres("R"));
	int decrement=1;		
	for(int i =pmax-1;i>0;--i)
		{
			Bres=ARFitVARXDD(R,T2,k,decrement,pmax);
			sigmaU=as<MatrixXd>(Bres("SigmaU"));
			decrement+=1;
			p= (double) i;

			if(pen=="BIC"){
				criti=log(sigmaU.determinant())+log(T2)*(pow(k,2)*p+k)/T2;
			}else if(pen=="AIC"){
				criti=log(sigmaU.determinant())+(2*(pow(k,2)*p+k))/T2;
			}
			// Rcout<<criti<<std::endl;
			crit.push_back(criti);

		}


	MatrixXd ones(Y2.rows(),1);
	ones.setOnes();
	MatrixXd K1a(Y2.rows(),1+Y2.cols());
	K1a<<ones,Y2;
	Bres = ARFitVARX(K1a,k,0,0,0);
	sigmaU = as<MatrixXd>(Bres("SigmaU"));
	// Rcout<<sigmaU<<std::endl;

	if(pen=="BIC"){
		criti=log(sigmaU.determinant())+log(T2)*(k)/T2;
	}else if(pen=="AIC"){
		criti=log(sigmaU.determinant())+(2*(k))/T2;
	}
	// Rcout<<criti<<std::endl;
	crit.push_back(criti);

	
	int min = std::distance(crit.begin(),std::min_element(crit.begin(),crit.end()));
	int phat=pmax-min;
	MatrixXd BResF;

	Bres=ARFitVARXDD(R,T2,k,pmax-phat,pmax);


	// MatrixXd ZF=ZmatF(Y,phat,k);
	// MatrixXd Y2F=Y.bottomRows(T-phat);
	// MatrixXd KF(Y2F.rows(),ZF.rows()+Y2F.cols());
	// KF<<ZF.transpose(),Y2F;
	// Bres= ARFitVARX(KF,k,phat,0,0);
	BResF =Bres("B");

	sigmaU=Bres("SigmaU");
	
	List Results=List::create(Named("B")=BResF,Named("p")=phat,Named("SigmaU")=sigmaU);

	return(Results);
  
}
